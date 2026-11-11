"""Per-event computation: astropy star positions + Ciddor ray trace + comparisons."""
from __future__ import annotations

import math
import warnings
from datetime import timedelta
from typing import Optional

from .data import Event
from .refraction import (CiddorAtmosphere, RefractionSettings, apparent_geometric_target, bennett_arcmin,
                         geometric_from_apparent_infinite, ray_path_xy, wgs84_gauss_radius)
from .soundings import sounding_for_event
from .stars import resolve
from .terrain import compute_terrain, to_peak_frame_xy

NaN = float("nan")


def _f(x, nd=None):
    """JSON-safe float (NaN -> None)."""
    if x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    return round(x, nd) if nd is not None else x


def pressure_hpa(alt_m):
    return 1013.25 * (1.0 - 2.25577e-5 * alt_m) ** 5.25588


def configure_astropy(offline: bool):
    from astropy.utils import iers
    iers.conf.iers_degraded_accuracy = "warn"
    if offline:
        iers.conf.auto_download = False
    warnings.filterwarnings("ignore", module="astropy")
    warnings.filterwarnings("ignore", message=".*ERFA.*")


def star_altitudes(e: Event, star: dict):
    """True (vacuum) altitude/azimuth at observer and peak, altitude rate, and astropy's own
    refracted altitude at the observer."""
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
    from astropy.time import Time

    t = Time(e.utc)
    obs = EarthLocation(lat=e.obs_lat * u.deg, lon=e.obs_lon * u.deg, height=e.obs_alt * u.m)
    peak = EarthLocation(lat=e.peak_lat * u.deg, lon=e.peak_lon * u.deg, height=e.peak_alt * u.m)

    def target(tt, loc):
        if star["kind"] == "body":
            return get_body(star["body"], tt, loc)
        c = SkyCoord(ra=star["ra_deg"] * u.deg, dec=star["dec_deg"] * u.deg,
                     pm_ra_cosdec=star.get("pmra_masyr", 0.0) * u.mas / u.yr,
                     pm_dec=star.get("pmdec_masyr", 0.0) * u.mas / u.yr,
                     distance=1e6 * u.pc, obstime=Time("J2000"), frame="icrs")
        return c.apply_space_motion(new_obstime=tt)

    def altaz(tt, loc, **kw):
        return target(tt, loc).transform_to(AltAz(obstime=tt, location=loc, **kw))

    a_obs = altaz(t, obs, pressure=0 * u.hPa)
    a_peak = altaz(t, peak, pressure=0 * u.hPa)
    dt = 30 * u.s
    rate = (altaz(t + dt, obs, pressure=0 * u.hPa).alt.deg - altaz(t - dt, obs, pressure=0 * u.hPa).alt.deg) / 60.0
    temp = e.csv.get("temp_c")
    temp = 10.0 if temp is None or not math.isfinite(temp) else temp
    a_app = altaz(t, obs, pressure=pressure_hpa(e.obs_alt) * u.hPa, temperature=temp * u.deg_C,
                  relative_humidity=0.5, obswl=0.55 * u.micron)
    return dict(true_alt_obs_deg=a_obs.alt.deg, az_obs_deg=a_obs.az.deg, true_alt_peak_deg=a_peak.alt.deg,
                alt_rate_deg_per_s=rate, erfa_apparent_alt_obs_deg=a_app.alt.deg, temp_c_used=temp)


CLOCK_CHECK_ARCSEC = 60.0


def implied_clock_offset(e: Event, star: dict, target_alt_deg: float, window_h=6.0, step_s=20.0):
    """Clock shift (minutes, + = later) at which astropy reproduces the CSV true altitude
    with the star within 3° of the peak azimuth. None if nothing fits."""
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time
    from geographiclib.geodesic import Geodesic
    paz = Geodesic.WGS84.Inverse(e.obs_lat, e.obs_lon, e.peak_lat, e.peak_lon)["azi1"] % 360
    loc = EarthLocation(lat=e.obs_lat * u.deg, lon=e.obs_lon * u.deg, height=e.obs_alt * u.m)
    t0 = Time(e.utc)
    c = SkyCoord(ra=star["ra_deg"] * u.deg, dec=star["dec_deg"] * u.deg, pm_ra_cosdec=star.get("pmra_masyr", 0) * u.mas / u.yr,
                 pm_dec=star.get("pmdec_masyr", 0) * u.mas / u.yr, distance=1e6 * u.pc, obstime=Time("J2000"), frame="icrs")
    c = c.apply_space_motion(new_obstime=t0)
    dts = np.arange(-window_h * 3600, window_h * 3600 + 1, step_s)
    aa = c.transform_to(AltAz(obstime=t0 + dts * u.s, location=loc, pressure=0 * u.hPa))
    alt, az = aa.alt.deg, aa.az.deg
    daz = np.abs(((az - paz + 540) % 360) - 180)
    cost = np.where(daz < 3.0, np.abs(alt - target_alt_deg), np.inf)
    if not np.isfinite(cost).any():
        return None
    i = int(np.argmin(cost))
    # refine to 1 s around the best coarse step
    fine = dts[i] + np.arange(-step_s, step_s + 1, 1.0)
    aa = c.transform_to(AltAz(obstime=t0 + fine * u.s, location=loc, pressure=0 * u.hPa))
    j = int(np.argmin(np.abs(aa.alt.deg - target_alt_deg)))
    return {"minutes": fine[j] / 60.0, "alt_residual_arcsec": (aa.alt.deg[j] - target_alt_deg) * 3600}


def compute_event(e: Event, settings: RefractionSettings, offline=False, log=print, profiles=True, refresh_stars=False,
                  terrain=True, terrain_settings=None) -> dict:
    warn = list(e.warnings)
    out = {
        "id": e.id, "row": e.row,
        "flags": {"target_type": e.target_type, "is_test": e.test_kind is not None, "test_kind": e.test_kind,
                  "warnings": warn},
        "peak": {"name": e.peak, "lat": _f(e.peak_lat), "lon": _f(e.peak_lon), "alt_m": _f(e.peak_alt)},
        "observer": {"name": e.observer, "lat": _f(e.obs_lat), "lon": _f(e.obs_lon), "alt_m": _f(e.obs_alt)},
        "time": {"local_date": e.local_date, "local_time": e.local_time, "utc_offset_h": _f(e.utc_offset_h),
                 "utc": e.utc.strftime("%Y-%m-%dT%H:%M:%SZ") if e.utc else None},
        "csv": {k: (v if isinstance(v, (str, dict)) else _f(v)) for k, v in e.csv.items()},
    }
    out["csv"]["star_name"] = e.star

    # ---- star
    star = resolve(e.star, offline=offline, log=log, refresh_no_pm=refresh_stars) if e.star else {"kind": "unresolved", "query": None, "source": "blank"}
    out["star"] = {"name": e.star, **{k: (_f(v) if isinstance(v, float) else v) for k, v in star.items()}}
    astro = None
    if star["kind"] != "unresolved" and e.utc:
        try:
            astro = star_altitudes(e, star)
        except Exception as ex:  # noqa: BLE001
            warn.append(f"astropy failed: {ex}")
    elif star["kind"] == "unresolved":
        warn.append(f"star not resolved ({star.get('source')})")
    out["astropy"] = {k: _f(v) for k, v in (astro or {}).items()}

    # ---- geometry + refraction
    planar = math.degrees(math.atan2(e.peak_alt - e.obs_alt, e.csv["distance"])) if math.isfinite(e.csv["distance"]) else NaN
    snd = sounding_for_event(e.id, e.obs_lat, e.obs_lon, e.utc, offline=offline, log=log) if e.utc else None
    out["refraction"] = {"sounding": None}
    sol = None
    if snd is not None and all(math.isfinite(v) for v in (e.obs_lat, e.peak_lat, e.obs_alt, e.peak_alt)):
        if not snd.height_m:
            warn.append("no radiosonde available; standard atmosphere used")
        atm = CiddorAtmosphere(snd, settings)
        sol = apparent_geometric_target(atm, e.obs_lat, e.obs_lon, e.obs_alt, e.peak_lat, e.peak_lon, e.peak_alt)
        star_geom_gauss = geometric_from_apparent_infinite(atm, sol.apparent_deg, e.obs_alt, wgs84_gauss_radius(e.obs_lat))
        star_geom_ns = geometric_from_apparent_infinite(atm, sol.apparent_deg, e.obs_alt, sol.radius_m)
        P = pressure_hpa(e.obs_alt)
        T = out["astropy"].get("temp_c_used") or (e.csv["temp_c"] if math.isfinite(e.csv["temp_c"]) else 10.0)
        out["refraction"] = {
            "sounding": {"station": snd.station, "wmo": snd.wmo, "valid_utc": snd.valid_utc, "source": snd.source,
                         "distance_km": _f(snd.distance_km, 1), "levels": len(snd.height_m),
                         "ceiling_m": _f(snd.max_alt) if snd.height_m else None},
            "peak_apparent_deg": _f(sol.apparent_deg), "peak_geometric_deg": _f(sol.geometric_deg),
            "peak_uplift_arcmin": _f(sol.uplift_deg * 60),
            "star_geometric_deg": _f(star_geom_ns), "star_geometric_gauss_deg": _f(star_geom_gauss),
            "star_uplift_arcmin": _f((sol.apparent_deg - star_geom_ns) * 60),
            "star_uplift_gauss_arcmin": _f((sol.apparent_deg - star_geom_gauss) * 60),
            "bennett_arcmin": _f(bennett_arcmin(sol.apparent_deg, P, T)),
            "pressure_hpa_obs": _f(P, 2), "temp_c_obs": _f(T),
        }
        if profiles:
            out["profiles"] = {
                "ray_to_peak": [[_f(x, 1), _f(y, 2)] for x, y in to_peak_frame_xy(
                    ray_path_xy(atm, sol.apparent_deg, e.obs_alt, sol.radius_m, sol.slant_m * 1.02, every=max(1, int(sol.ground_m / 50 / 120))),
                    math.radians(sol.phi_deg), sol.radius_m)],
                "atmosphere": [
                    {"alt_m": a, "n_minus_1_ppm": _f((atm.n_at(a) - 1) * 1e6, 3),
                     **dict(zip(("T_K", "P_Pa", "RH_pct"), (_f(v, 2) for v in atm.TP_RH(a))))}
                    for a in [e.obs_alt + k * 500 for k in range(int((40000 - e.obs_alt) / 500) + 1)]
                ],
            }
    out["geometry"] = {
        "csv_distance_m": _f(e.csv["distance"]),
        "geodesic_distance_m": _f(sol.ground_m) if sol else None,
        "slant_range_m": _f(sol.slant_m) if sol else None,
        "central_angle_deg": _f(sol.phi_deg) if sol else None,
        "azimuth_deg": _f(sol.azimuth_deg) if sol else None,
        "earth_radius_m": _f(sol.radius_m) if sol else None,
        "planar_angle_deg": _f(planar),
        "peak_height_above_observer_m": _f(e.peak_alt - e.obs_alt),
        "frame": "peak-anchored: origin at the peak foot, tangent plane at the peak; only the globe observer drops",
    }
    if sol:
        R, D = sol.radius_m, sol.ground_m
        gz = (R + e.obs_alt) * math.cos(D / R) - R
        out["geometry"].update({
            "observer_planar_xz_m": [_f(-D, 1), _f(e.obs_alt, 1)],
            "observer_globe_xz_m": [_f((R + e.obs_alt) * math.sin(-D / R), 1), _f(gz, 2)],
            "globe_observer_drop_m": _f(e.obs_alt - gz, 2),
            "peak_xz_m": [0.0, _f(e.peak_alt, 1)],
        })
    if sol and e.csv["distance"] and math.isfinite(e.csv["distance"]) and abs(sol.ground_m - e.csv["distance"]) > 200:
        warn.append(f"CSV distance differs from geodesic by {sol.ground_m - e.csv['distance']:.0f} m")

    # ---- comparisons (positive = star truly higher than the model predicts at occlusion)
    r = out["refraction"]
    a = out["astropy"]
    c = {}
    if a.get("true_alt_obs_deg") is not None:
        ta = a["true_alt_obs_deg"]
        rate = a.get("alt_rate_deg_per_s")
        sec = lambda d: _f(d / rate) if rate else None
        if r.get("star_geometric_deg") is not None:
            d = ta - r["star_geometric_deg"]
            c["globe_ciddor_arcsec"] = _f(d * 3600, 2)
            c["globe_ciddor_seconds"] = sec(d)
        if math.isfinite(planar):
            d = ta - planar
            c["planar_geometric_arcsec"] = _f(d * 3600, 2)
            c["planar_geometric_seconds"] = sec(d)
            if a.get("erfa_apparent_alt_obs_deg") is not None:
                d2 = a["erfa_apparent_alt_obs_deg"] - planar
                c["planar_refracted_arcsec"] = _f(d2 * 3600, 2)
                c["planar_refracted_seconds"] = sec(d2)
        if e.csv.get("true_obs") is not None and math.isfinite(e.csv["true_obs"]):
            c["csv_true_obs_minus_astropy_arcsec"] = _f((e.csv["true_obs"] - ta) * 3600, 2)
            if abs(e.csv["true_obs"] - ta) * 3600 > CLOCK_CHECK_ARCSEC and star["kind"] == "star":
                try:
                    k = implied_clock_offset(e, star, e.csv["true_obs"])
                    if k:
                        c["csv_implied_clock_offset_min"] = _f(k["minutes"], 1)
                        c["csv_implied_alt_residual_arcsec"] = _f(k["alt_residual_arcsec"], 1)
                        if abs(k["alt_residual_arcsec"]) < 30 and abs(round(k["minutes"] / 60) * 60 - k["minutes"]) < 2 and round(k["minutes"] / 60) != 0:
                            warn.append(f"CSV true altitude matches astropy only if the clock is shifted {k['minutes']:+.0f} min "
                                        f"(UTC offset {e.utc_offset_h - k['minutes'] / 60:+.0f} instead of {e.utc_offset_h:+.0f}); "
                                        "check the UTC column (config/time_corrections.json)")
                        else:
                            warn.append(f"CSV true altitude differs from astropy by {(e.csv['true_obs'] - ta) * 3600:+.0f}″ and no clean clock shift explains it "
                                        f"(best {k['minutes']:+.1f} min, residual {k['alt_residual_arcsec']:.0f}″): check the star identity, date and time")
                    else:
                        warn.append(f"CSV true altitude differs from astropy by {(e.csv['true_obs'] - ta) * 3600:+.0f}″; no time within ±6 h fits with the star behind the peak")
                except Exception as ex:  # noqa: BLE001
                    warn.append(f"clock check failed: {ex}")
        if e.csv.get("true_peak") is not None and math.isfinite(e.csv["true_peak"]) and a.get("true_alt_peak_deg") is not None:
            c["csv_true_peak_minus_astropy_arcsec"] = _f((e.csv["true_peak"] - a["true_alt_peak_deg"]) * 3600, 2)
    if r.get("star_uplift_arcmin") is not None and r.get("bennett_arcmin") is not None:
        c["ciddor_minus_bennett_arcmin"] = _f(r["star_uplift_arcmin"] - r["bennett_arcmin"], 4)
    out["comparison"] = c

    # ---- terrain (real surface along the line of sight)
    if terrain and sol is not None and e.target_type == "aircraft":
        out["terrain"] = {"skipped": "aircraft target: the occluder is the aircraft, not terrain"}
    elif terrain and sol is not None:
        try:
            doc, summ = compute_terrain(out, atm, terrain_settings, offline=offline, log=log)
            out["terrain"] = summ
            out["_terrain_doc"] = doc
            if summ.get("checks", {}).get("dtm_at_peak_m") is not None and abs(summ["checks"]["dtm_at_peak_m"] - e.peak_alt) > 30:
                warn.append(f"DEM height at the peak coordinate differs from CSV by {summ['checks']['dtm_at_peak_m'] - e.peak_alt:+.0f} m")
            for nm, cf in (summ.get("canopy") or {}).items():
                if cf.get("flag"):
                    warn.append(f"canopy check: DSM skyline is {cf['dsm_minus_dtm_skyline_arcsec']:+.0f}″ vs bare earth on the {nm.replace('_', ' ')} line (review)")
                    break
            for nm in ("star_azimuth", "to_peak"):
                ob = (summ.get(nm) or {}).get("obstruction") or {}
                gr = ((ob.get("dtm") or ob.get("dsm") or {}).get("globe_refracted")) or {}
                if gr and not gr.get("near_peak", True):
                    warn.append(f"terrain {gr['u_m'] / 1000:+.1f} km from the peak rises above the summit on the {nm.replace('_', ' ')} line")
                    break
        except Exception as ex:  # noqa: BLE001
            warn.append(f"terrain failed: {ex}")
    return out
