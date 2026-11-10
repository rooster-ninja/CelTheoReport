"""Ciddor ray-trace refraction engine.

Vendored from rooster-ninja/AL_RefractionAnalysis (AstronomyLive_Refraction/refraction.py).
The physics and numerics are unchanged. What changed:

* The atmosphere is a per-event object (``CiddorAtmosphere``) built from a ``Sounding``
  instead of module-level globals, so many events can be processed in one run.
* The star (infinite-distance) trace can use either the Gauss radius (original behaviour)
  or the same normal-section radius used for the peak trace. Both are reported.
* Optional ray-path sampling for plotting.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from ambiance import Atmosphere
from geographiclib.geodesic import Geodesic


# --------------------------------------------------------------------------- sounding
@dataclass
class Sounding:
    """Vertical profile, ascending height. Units: m, hPa, °C, %."""
    height_m: List[float]
    pres_hpa: List[float]
    temp_c: List[float]
    rh_pct: List[float]
    station: str = ""
    wmo: str = ""
    valid_utc: str = ""
    source: str = ""
    distance_km: float = float("nan")

    @property
    def max_alt(self) -> float:
        return float(self.height_m[-1]) if self.height_m else -1e9

    def interp(self, alt_m: float, values: Sequence[float]) -> float:
        h = self.height_m
        if alt_m <= h[0]:
            return float(values[0])
        if alt_m >= h[-1]:
            return float(values[-1])
        lo, hi = 0, len(h) - 1
        while hi - lo > 1:  # binary search (original used a linear scan; same result)
            mid = (lo + hi) // 2
            if h[mid] <= alt_m:
                lo = mid
            else:
                hi = mid
        z0, z1 = h[lo], h[hi]
        t = (alt_m - z0) / (z1 - z0) if z1 != z0 else 0.0
        return float(values[lo] + t * (values[hi] - values[lo]))


EMPTY_SOUNDING = Sounding([], [], [], [], source="none")


# --------------------------------------------------------------------------- Ciddor
def ciddor_Z(T_k, p_pa, xw):
    t_c = T_k - 273.15
    a0, a1, a2 = 1.58123e-6, -2.9331e-8, 1.1043e-10
    b0, b1 = 5.707e-6, -2.051e-8
    c0, c1 = 1.9898e-4, -2.376e-6
    d, e = 1.83e-11, -0.765e-8
    return 1 - (p_pa / T_k) * (
        a0 + a1 * t_c + a2 * t_c**2 + (b0 + b1 * t_c) * xw + (c0 + c1 * t_c) * xw**2
    ) + (p_pa / T_k) ** 2 * (d + e * xw**2)


def ciddor(lambda_nm, t_c, p_pa, rh_percent, xc_ppm):
    """Ciddor (1996) refractive index of moist air."""
    lam_um = lambda_nm / 1000.0
    h = rh_percent / 100.0
    sigma = 1.0 / lam_um
    T = t_c + 273.15
    R = 8.314510
    k0, k1, k2, k3 = 238.0185, 5792105, 57.362, 167917
    w0, w1, w2, w3 = 295.235, 2.6422, -0.032380, 0.004028
    A, B, C, D = 1.2378847e-5, -1.9121316e-2, 33.93711047, -6.3431645e3
    alpha, beta, gamma = 1.00062, 3.14e-8, 5.6e-7
    svp = math.exp(A * T**2 + B * T + C + D / T) if t_c >= 0 else 10 ** (-2663.5 / T + 12.537)
    f = alpha + beta * p_pa + gamma * t_c**2
    xw = (f * h * svp / p_pa) if p_pa > 0 else 0.0
    nas = 1 + (k1 / (k0 - sigma**2) + k3 / (k2 - sigma**2)) * 1e-8
    naxs = 1 + (nas - 1) * (1 + 0.534e-6 * (xc_ppm - 450))
    nws = 1 + 1.022 * (w0 + w1 * sigma**2 + w2 * sigma**4 + w3 * sigma**6) * 1e-8
    Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc_ppm - 400))
    Mw = 0.018015
    Za = ciddor_Z(288.15, 101325, 0.0)
    Zw = ciddor_Z(293.15, 1333, 1.0)
    paxs = 101325 * Ma / (Za * R * 288.15)
    pws = 1333 * Mw / (Zw * R * 293.15)
    Z = ciddor_Z(T, p_pa, xw)
    pa = (p_pa * Ma / (Z * R * T) * (1 - xw)) if p_pa > 0 else 0.0
    pw = (p_pa * Mw / (Z * R * T) * xw) if p_pa > 0 else 0.0
    return 1.0 + (pa / paxs) * (naxs - 1.0) + (pw / pws) * (nws - 1.0)


# --------------------------------------------------------------------------- geometry
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
SPHERE_RADIUS_M = 6371000.0
AMB_MIN_ALT_M = -5_000.0
AMB_MAX_ALT_M = 80_000.0


def wgs84_radii_of_curvature(lat_deg):
    s = math.sin(math.radians(lat_deg))
    denom = math.sqrt(1.0 - WGS84_E2 * s * s)
    return WGS84_A * (1.0 - WGS84_E2) / denom**3, WGS84_A / denom  # M, N


def wgs84_normal_section_radius(lat_deg, azimuth_deg):
    M, N = wgs84_radii_of_curvature(lat_deg)
    A = math.radians(azimuth_deg)
    return 1.0 / ((math.cos(A) ** 2) / M + (math.sin(A) ** 2) / N)


def wgs84_gauss_radius(lat_deg):
    M, N = wgs84_radii_of_curvature(lat_deg)
    return math.sqrt(M * N)


def haversine_central_angle_rad(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dlat, dlon = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlon / 2) ** 2
    return 2.0 * math.asin(math.sqrt(min(1.0, max(0.0, a))))


@dataclass
class RefractionSettings:
    use_wgs84_ellipsoid: bool = True
    wavelength_nm: float = 550.0
    co2_ppm: float = 450.0
    segment_m: float = 50.0
    top_of_atmosphere_m: float = 80000.0
    max_total_path_m: float = 500_000.0


# --------------------------------------------------------------------------- atmosphere
class CiddorAtmosphere:
    """n(alt): sounding T/P/RH up to its ceiling, ambiance standard atmosphere above (RH 0),
    vacuum above 80 km. With an empty sounding, ambiance is used throughout."""

    def __init__(self, sounding: Sounding, settings: RefractionSettings):
        self.s = sounding
        self.settings = settings
        self.q = settings.segment_m
        self._cache = {}

    def TP_RH(self, alt_m):
        if not self.s.height_m or alt_m > self.s.max_alt:
            a = Atmosphere([min(max(alt_m, AMB_MIN_ALT_M), AMB_MAX_ALT_M - 1)])
            return float(a.temperature[0]), float(a.pressure[0]), 0.0
        T = self.s.interp(alt_m, self.s.temp_c) + 273.15
        P = self.s.interp(alt_m, self.s.pres_hpa) * 100
        RH = self.s.interp(alt_m, self.s.rh_pct)
        return T, P, RH

    def n_at(self, alt_m):
        z = float(alt_m)
        if z >= AMB_MAX_ALT_M:
            return 1.0
        if z < AMB_MIN_ALT_M:
            z = AMB_MIN_ALT_M
        k = int(round(z / self.q))
        if k in self._cache:
            return self._cache[k]
        T, P, RH = self.TP_RH(z)
        n = ciddor(self.settings.wavelength_nm, T - 273.15, P, RH, self.settings.co2_ppm)
        self._cache[k] = float(n)
        return float(n)


# --------------------------------------------------------------------------- tracing
def _trace_to_radius(atm, b, observer_alt_m, segment_m, r_target, top_m, s_max, R, path=None, path_every=1):
    r = R + observer_alt_m
    phi = s = 0.0
    ds = float(segment_m)
    i = 0
    while True:
        if path is not None and i % path_every == 0:
            path.append((phi, r))
        i += 1
        alt = r - R
        if alt >= top_m:
            return phi, r, s, False
        if s_max is not None and s >= s_max:
            return phi, r, s, False
        if r >= r_target:
            return phi, r, s, True
        n = atm.n_at(alt)
        arg = b / (n * r)
        if arg > 1.0:
            return phi, r, s, False
        alpha = math.acos(max(-1.0, min(1.0, arg)))
        dr = ds * math.sin(alpha)
        dphi = (ds * math.cos(alpha)) / r
        if r + dr >= r_target:
            frac = (r_target - r) / max(1e-12, dr)
            phi += dphi * frac
            r = r_target
            s += ds * frac
            if path is not None:
                path.append((phi, r))
            return phi, r, s, True
        r += dr
        phi += dphi
        s += ds
        if phi > 800.0 or r > (R + top_m) * 2:
            return phi, r, s, False


def geometric_from_apparent_infinite(atm, apparent_elev_deg, observer_alt_m, R, path=None):
    """True (vacuum) elevation of a star seen at ``apparent_elev_deg``."""
    st = atm.settings
    r0 = R + observer_alt_m
    b = atm.n_at(observer_alt_m) * r0 * math.cos(math.radians(apparent_elev_deg))
    phi_top, r_hit, _s, ok = _trace_to_radius(atm, b, observer_alt_m, st.segment_m, R + st.top_of_atmosphere_m,
                                             st.top_of_atmosphere_m, None, R, path, 20)
    if not ok:
        return float("nan")
    arg = b / r_hit
    if arg > 1.0:
        return float("nan")
    alpha_exit = math.acos(max(-1.0, min(1.0, arg)))
    rad = (math.cos(phi_top), math.sin(phi_top))
    tan = (-math.sin(phi_top), math.cos(phi_top))
    dx = math.sin(alpha_exit) * rad[0] + math.cos(alpha_exit) * tan[0]
    dy = math.sin(alpha_exit) * rad[1] + math.cos(alpha_exit) * tan[1]
    return math.degrees(math.atan2(dx, dy))


@dataclass
class TargetSolution:
    apparent_deg: float
    geometric_deg: float
    uplift_deg: float
    ground_m: float
    phi_deg: float
    azimuth_deg: float
    radius_m: float
    slant_m: float


def apparent_geometric_target(atm, obs_lat, obs_lon, obs_alt, tgt_lat, tgt_lon, tgt_alt) -> TargetSolution:
    """Apparent elevation at which a refracted ray from the observer lands on the target."""
    st = atm.settings
    if st.use_wgs84_ellipsoid:
        inv = Geodesic.WGS84.Inverse(obs_lat, obs_lon, tgt_lat, tgt_lon)
        ground_m, az = float(inv["s12"]), float(inv["azi1"])
        R = wgs84_normal_section_radius(obs_lat, az)
        phi_t = ground_m / R
    else:
        R = SPHERE_RADIUS_M
        phi_t = haversine_central_angle_rad(obs_lat, obs_lon, tgt_lat, tgt_lon)
        ground_m, az = R * phi_t, float("nan")
    r0, r_t = R + obs_alt, min(R + tgt_alt, R + st.top_of_atmosphere_m)
    vx, vy = r_t * math.cos(phi_t) - r0, r_t * math.sin(phi_t)
    alpha_geom = math.atan2(vx, vy)
    slant = math.hypot(vx, vy)
    n0 = atm.n_at(obs_alt)
    eps = math.radians(0.0001)

    def g(a):
        a = max(a, eps)
        ph, _r, _s, ok = _trace_to_radius(atm, n0 * r0 * math.cos(a), obs_alt, st.segment_m, r_t,
                                          st.top_of_atmosphere_m, st.max_total_path_m, R)
        return (ph - phi_t) if ok else float("nan")

    center, span, max_span = max(alpha_geom, eps), math.radians(0.5), math.radians(20.0)
    lo, hi = max(eps, center - span), min(math.radians(89.9), center + span)
    g_lo, g_hi = g(lo), g(hi)
    while (not math.isfinite(g_lo) or not math.isfinite(g_hi) or g_lo * g_hi > 0) and span < max_span:
        span *= 1.7
        lo, hi = max(eps, center - span), min(math.radians(89.9), center + span)
        g_lo, g_hi = g(lo), g(hi)
    if not math.isfinite(g_lo) or not math.isfinite(g_hi) or g_lo * g_hi > 0:
        samples = [max(eps, center - k * math.radians(0.5)) for k in range(40)]
        samples += [min(math.radians(89.9), center + k * math.radians(0.5)) for k in range(1, 100)]
        prev = None
        bracket = None
        for a in samples:
            ga = g(a)
            if not math.isfinite(ga):
                continue
            if prev is not None and prev[1] * ga <= 0:
                bracket = (prev[0], a, prev[1], ga)
                break
            prev = (a, ga)
        if bracket is None:
            gd = math.degrees(alpha_geom)
            return TargetSolution(float("nan"), gd, float("nan"), ground_m, math.degrees(phi_t), az, R, slant)
        lo, hi, g_lo, g_hi = bracket
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        g_mid = g(mid)
        if not math.isfinite(g_mid):
            lo = mid
            continue
        if abs(g_mid) < 1e-12:
            lo = hi = mid
            break
        if g_lo * g_mid <= 0:
            hi, g_hi = mid, g_mid
        else:
            lo, g_lo = mid, g_mid
    app = math.degrees(0.5 * (lo + hi))
    geom = math.degrees(alpha_geom)
    return TargetSolution(app, geom, app - geom, ground_m, math.degrees(phi_t), az, R, slant)


def ray_path_xy(atm, apparent_deg, obs_alt, R, until_m, every=10):
    """Sample a refracted ray (observer-local x along ground, y height) for plotting."""
    st = atm.settings
    r0 = R + obs_alt
    b = atm.n_at(obs_alt) * r0 * math.cos(math.radians(apparent_deg))
    path: List[Tuple[float, float]] = []
    _trace_to_radius(atm, b, obs_alt, st.segment_m, R + st.top_of_atmosphere_m, st.top_of_atmosphere_m,
                     until_m, R, path, every)
    return [(r * math.sin(p), r * math.cos(p) - R) for p, r in path]


def bennett_arcmin(apparent_deg, pressure_hpa=1010.0, temp_c=10.0):
    """Bennett (1982) refraction for an apparent altitude, with the usual P/T factor."""
    h = apparent_deg
    R = 1.0 / math.tan(math.radians(h + 7.31 / (h + 4.4)))
    R -= 0.06 * math.sin(math.radians(14.7 * R + 13))
    return R * (pressure_hpa / 1010.0) * (283.0 / (273.0 + temp_c))
