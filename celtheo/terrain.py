"""Terrain along the line of sight: DEM sampling, peak-anchored planar/globe frames, obstruction.

Frame (matches the Desmos graph): origin at the peak's foot, tangent plane at the peak.
Both models put the peak at the same place; only the globe observer (and terrain toward it) drops.

    u = s - D                                    along-track distance from the peak (negative toward observer)
    planar:  x = u,                    z = h
    globe:   x = (R+h) sin(u/R),       z = (R+h) cos(u/R) - R

Elevation angles are physical (frame independent) and measured from the observer's local horizon.
Two surfaces are sampled: DTM (bare earth, primary: most summits are bare rock above tree line) and
DSM (top of canopy/buildings). Where the DSM skyline differs from the DTM by more than
``canopy_flag_arcsec`` the event is flagged for case-by-case review.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from geographiclib.geodesic import Geodesic

ROOT = Path(__file__).resolve().parent.parent
SOURCES = json.loads((ROOT / "config" / "dem_sources.json").read_text())
STAC_CACHE = ROOT / "data" / "dem_stac_cache.json"
G = Geodesic.WGS84

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.TIF,.tiff,.vrt")
os.environ.setdefault("GDAL_HTTP_MULTIRANGE", "YES")
os.environ.setdefault("GDAL_HTTP_MERGE_CONSECUTIVE_RANGES", "YES")
os.environ.setdefault("VSI_CACHE", "TRUE")
os.environ.setdefault("GDAL_CACHEMAX", "512")


def _close_rasters():
    Raster.close_all()


@dataclass
class TerrainSettings:
    profile_step_m: float = 10.0
    beyond_peak_m: float = 20000.0      # profiles continue past the peak: anything on the line can hide the star
    corridor: bool = False              # 3D corridor grid (off: this analysis uses 2D profiles)
    corridor_half_width_m: float = 2500.0
    corridor_step_m: float = 100.0
    candidate_window_arcsec: float = 90.0
    max_candidates: int = 8
    canopy_flag_arcsec: float = 30.0   # DSM vs DTM skyline difference worth a manual look
    near_field_frac: float = 0.05      # samples closer than this fraction of D (clamped below) are the observer's
    near_field_min_m: float = 100.0    # foreground: reported separately, excluded from the obstruction search
    near_field_max_m: float = 1000.0


# --------------------------------------------------------------------------- tiles
def _tile_urls(src: dict, lat: float, lon: float) -> List[str]:
    """URL(s) of the tile covering a point for templated tile sources."""
    kind = src["tiling"]
    if kind == "usgs_1deg":      # named by NW corner: n39w106 covers 38..39N, 106..105W
        n, w = math.floor(lat) + 1, math.floor(-lon) + 1 if lon < 0 else None
        if w is None:
            return []
        return [t.format(n=n, w=w) for t in src["url_templates"]]
    if kind == "copernicus_1deg":  # named by SW corner
        la, lo = math.floor(lat), math.floor(lon)
        ns, ew = ("N" if la >= 0 else "S"), ("E" if lo >= 0 else "W")
        return [t.format(ns=ns, lat=abs(la), ew=ew, lon=abs(lo)) for t in src["url_templates"]]
    if kind == "single":
        return list(src["url_templates"])
    return []


def _stac_items(src: dict, bbox, log=print) -> List[str]:
    """COG hrefs (asset src['asset']) intersecting bbox from a STAC API, cached."""
    import requests
    key = f"{src['collection']}|{src['asset']}|{','.join(f'{b:.3f}' for b in bbox)}"
    cache = json.loads(STAC_CACHE.read_text()) if STAC_CACHE.exists() else {}
    if key in cache:
        return cache[key]
    try:
        r = requests.get(f"{src['stac_api'].rstrip('/')}/search", params={"collections": src["collection"],
                         "bbox": ",".join(str(b) for b in bbox), "limit": 100}, timeout=60)
        r.raise_for_status()
        hrefs = [f["assets"][src["asset"]]["href"] for f in r.json().get("features", []) if src["asset"] in f.get("assets", {})]
    except Exception as e:  # noqa: BLE001
        log(f"    STAC search failed for {src['name']}: {e}")
        return []
    cache[key] = hrefs
    STAC_CACHE.parent.mkdir(parents=True, exist_ok=True)
    STAC_CACHE.write_text(json.dumps(cache, indent=1))
    return hrefs


class Raster:
    _open: Dict[str, object] = {}
    _bad: set = set()

    @classmethod
    def close_all(cls):
        for ds in cls._open.values():
            try:
                ds.close()
            except Exception:  # noqa: BLE001
                pass
        cls._open.clear()

    @classmethod
    def get(cls, url: str, log=print):
        if url in cls._bad:
            return None
        if url not in cls._open:
            import rasterio
            try:
                path = url if url.startswith("/") or "://" not in url else "/vsicurl/" + url
                cls._open[url] = rasterio.open(path)
            except Exception as e:  # noqa: BLE001
                cls._bad.add(url)
                log(f"    DEM open failed: {url.split('/')[-1]} ({str(e).splitlines()[0][:120]})")
                return None
        return cls._open[url]


def _to_crs(ds, lats, lons):
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
    return tr.transform(np.asarray(lons), np.asarray(lats))


def _sample_points(ds, lats, lons, max_px=512):
    """Nearest-pixel values at points (NaN outside / nodata). Points along a line are read in chunks:
    one small window per chunk instead of one request per point."""
    from rasterio.windows import Window
    xs, ys = _to_crs(ds, lats, lons)
    out = np.full(len(xs), np.nan)
    b = ds.bounds
    inside = (xs >= b.left) & (xs <= b.right) & (ys >= b.bottom) & (ys <= b.top)
    idx = np.nonzero(inside)[0]
    if not len(idx):
        return out
    inv = ~ds.transform
    cols, rows = inv * (xs[idx], ys[idx])
    cols, rows = np.floor(cols).astype(int), np.floor(rows).astype(int)
    # greedy chunks of consecutive points whose window stays within max_px on a side (~1 COG block area)
    bounds, k, n = [], 0, len(idx)
    while k < n:
        j, cmin, cmax, rmin, rmax = k, cols[k], cols[k], rows[k], rows[k]
        while j + 1 < n:
            c, r = cols[j + 1], rows[j + 1]
            if max(cmax, c) - min(cmin, c) > max_px or max(rmax, r) - min(rmin, r) > max_px:
                break
            j += 1
            cmin, cmax, rmin, rmax = min(cmin, c), max(cmax, c), min(rmin, r), max(rmax, r)
        bounds.append((k, j + 1))
        k = j + 1
    for k0, k1 in bounds:
        sl = slice(k0, k1)
        c0, c1 = max(cols[sl].min(), 0), min(cols[sl].max(), ds.width - 1)
        r0, r1 = max(rows[sl].min(), 0), min(rows[sl].max(), ds.height - 1)
        if c1 < c0 or r1 < r0:
            continue
        arr = ds.read(1, window=Window(c0, r0, c1 - c0 + 1, r1 - r0 + 1), masked=True).astype(float).filled(np.nan)
        cc, rr = cols[sl] - c0, rows[sl] - r0
        ok = (cc >= 0) & (rr >= 0) & (cc < arr.shape[1]) & (rr < arr.shape[0])
        vals = np.full(len(cc), np.nan)
        vals[ok] = arr[rr[ok], cc[ok]]
        if ds.nodata is not None:
            vals[vals == ds.nodata] = np.nan
        out[idx[sl]] = vals
    return out


def _sample_grid(ds, lats, lons, step_m):
    """Bilinear values at many points from one decimated windowed read (for the corridor)."""
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds
    xs, ys = _to_crs(ds, lats, lons)
    b = ds.bounds
    inside = (xs >= b.left) & (xs <= b.right) & (ys >= b.bottom) & (ys <= b.top)
    out = np.full(len(xs), np.nan)
    if not inside.any():
        return out
    x0, x1 = max(xs[inside].min(), b.left), min(xs[inside].max(), b.right)
    y0, y1 = max(ys[inside].min(), b.bottom), min(ys[inside].max(), b.top)
    win = from_bounds(x0, y0, x1, y1, ds.transform).round_offsets().round_lengths()
    if win.width < 2 or win.height < 2:
        return _sample_points(ds, lats, lons)
    res_x = abs(ds.transform.a) if ds.crs.is_projected else abs(ds.transform.a) * 111320 * math.cos(math.radians(float(np.mean(lats))))
    f = max(1.0, (step_m / 2) / max(res_x, 1e-6))
    shape = (max(2, int(win.height / f)), max(2, int(win.width / f)))
    arr = ds.read(1, window=win, out_shape=shape, masked=True, resampling=Resampling.bilinear).astype(float).filled(np.nan)
    wt = ds.window_transform(win)
    sx, sy = wt.a * win.width / shape[1], wt.e * win.height / shape[0]
    col = (xs - wt.c) / sx - 0.5
    row = (ys - wt.f) / sy - 0.5
    c0, r0 = np.floor(col).astype(int), np.floor(row).astype(int)
    ok = inside & (c0 >= 0) & (r0 >= 0) & (c0 < shape[1] - 1) & (r0 < shape[0] - 1)
    i = np.nonzero(ok)[0]
    fc, fr = col[i] - c0[i], row[i] - r0[i]
    v00, v01 = arr[r0[i], c0[i]], arr[r0[i], c0[i] + 1]
    v10, v11 = arr[r0[i] + 1, c0[i]], arr[r0[i] + 1, c0[i] + 1]
    out[i] = (v00 * (1 - fc) + v01 * fc) * (1 - fr) + (v10 * (1 - fc) + v11 * fc) * fr
    return out


def sample_surface(surface: str, lats, lons, grid_step_m: Optional[float] = None, offline=False, log=print) -> Tuple[np.ndarray, List[str]]:
    """Values for 'dtm' or 'dsm' from the first source covering each point; returns (values, sources used)."""
    lats, lons = np.asarray(lats, float), np.asarray(lons, float)
    vals = np.full(len(lats), np.nan)
    used: List[str] = []
    bbox = (float(lons.min()) - 0.01, float(lats.min()) - 0.01, float(lons.max()) + 0.01, float(lats.max()) + 0.01)
    for src in SOURCES[surface]:
        if offline and not src.get("local", False):
            continue
        todo = np.isnan(vals)
        if not todo.any():
            break
        if src.get("region") and not _in_region(src["region"], lats[todo], lons[todo]):
            continue
        if src["tiling"] == "stac":
            groups = [[u] for u in _stac_items(src, bbox, log=log)]
        else:
            groups = sorted({tuple(_tile_urls(src, la, lo)) for la, lo in zip(lats[todo], lons[todo])})
        got_any = False
        for group in groups:
            ds = None
            for url in group:
                ds = Raster.get(url, log=log)
                if ds is not None:
                    break
            if ds is None:
                continue
            todo = np.isnan(vals)
            idx = np.nonzero(todo)[0]
            if not len(idx):
                break
            log(f"      {surface}: {src['name']} ← {group[0].split('/')[-1]} ({len(idx)} pts)")
            v = _sample_grid(ds, lats[idx], lons[idx], grid_step_m) if grid_step_m else _sample_points(ds, lats[idx], lons[idx])
            hit = ~np.isnan(v)
            if hit.any():
                vals[idx[hit]] = v[hit]
                got_any = True
        if got_any:
            used.append(src["name"])
    return vals, used


def _in_region(region, lats, lons):
    if region == "canada":
        return bool(np.any(lats >= 48.99))
    if region == "usa":
        return bool(np.any(lats < 49.01))
    return True


import atexit  # noqa: E402

atexit.register(_close_rasters)


# --------------------------------------------------------------------------- geometry
def line_points(lat0, lon0, azi, s_list):
    line = G.Line(lat0, lon0, azi)
    lats, lons = np.empty(len(s_list)), np.empty(len(s_list))
    for i, s in enumerate(s_list):
        p = line.Position(float(s))
        lats[i], lons[i] = p["lat2"], p["lon2"]
    return lats, lons


def offset_points(lats, lons, azis, t):
    """Points offset perpendicular (to the right for t>0) from each centreline point."""
    out_la, out_lo = np.empty(len(lats)), np.empty(len(lats))
    for i, (la, lo, az) in enumerate(zip(lats, lons, azis)):
        p = G.Direct(la, lo, az + 90.0, t)
        out_la[i], out_lo[i] = p["lat2"], p["lon2"]
    return out_la, out_lo


def frames(s, h, D, R):
    """Peak-anchored planar and globe coordinates."""
    u = np.asarray(s) - D
    h = np.asarray(h)
    return {"planar_x": u, "planar_z": h,
            "globe_x": (R + h) * np.sin(u / R), "globe_z": (R + h) * np.cos(u / R) - R}


def angles(s, h, obs_alt, R):
    """Elevation angles (deg) from the observer: planar rise-over-run, and globe geometric."""
    s, h = np.asarray(s, float), np.asarray(h, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        planar = np.degrees(np.arctan2(h - obs_alt, s))
        phi = s / R
        dx = (R + h) * np.sin(phi)
        dz = (R + h) * np.cos(phi) - (R + obs_alt)
        globe = np.degrees(np.arctan2(dz, dx))
    return planar, globe


def to_peak_frame_xy(points_obs_frame, phi_t_rad, R):
    """Convert [x, y] in the observer-anchored frame (y = r cos φ − R) to the peak-anchored frame."""
    out = []
    for x, y in points_obs_frame:
        r = math.hypot(x, y + R)
        psi = math.atan2(x, y + R) - phi_t_rad
        out.append([r * math.sin(psi), r * math.cos(psi) - R])
    return out


# --------------------------------------------------------------------------- per-event build
def _r(a, nd=1):
    return [None if (x is None or not np.isfinite(x)) else round(float(x), nd) for x in np.asarray(a, float)]


def _blocker(s, u, lats, lons, h, ang, summit_angle, near_m=300.0):
    if not np.isfinite(ang).any():
        return None
    i = int(np.nanargmax(ang))
    return {"s_m": round(float(s[i]), 1), "u_m": round(float(u[i]), 1), "lat": round(float(lats[i]), 6), "lon": round(float(lons[i]), 6),
            "h_m": round(float(h[i]), 1), "angle_deg": float(ang[i]), "near_peak": bool(abs(u[i]) <= near_m),
            "margin_vs_summit_arcsec": round(float((ang[i] - summit_angle) * 3600), 2) if summit_angle is not None else None}


def build_profile(obs, azi, D, R, ts: TerrainSettings, offline=False, log=print):
    lat0, lon0, h_obs = obs
    s = np.arange(ts.profile_step_m, D + ts.beyond_peak_m + 1e-6, ts.profile_step_m)
    lats, lons = line_points(lat0, lon0, azi, s)
    dtm, src_dtm = sample_surface("dtm", lats, lons, offline=offline, log=log)
    dsm, src_dsm = sample_surface("dsm", lats, lons, offline=offline, log=log)
    return {"s": s, "lats": lats, "lons": lons, "dtm": dtm, "dsm": dsm, "sources": {"dtm": src_dtm, "dsm": src_dsm}}


def obstruction(prof, obs, D, R, summit, atm, ts: TerrainSettings, solve, log=print):
    """Highest point along the line by planar angle, globe geometric angle and globe refracted (apparent) angle."""
    s, lats, lons = prof["s"], prof["lats"], prof["lons"]
    u = s - D
    lat0, lon0, h_obs = obs
    s_lat, s_lon, s_h = summit
    planar_summit = math.degrees(math.atan2(s_h - h_obs, D))
    _, gs = angles([D], [s_h], h_obs, R)
    globe_summit = float(gs[0])
    out = {}
    near = min(max(ts.near_field_frac * D, ts.near_field_min_m), ts.near_field_max_m)
    fg = s < near
    for surf in ("dtm", "dsm"):
        h = prof[surf]
        if not np.isfinite(h).any():
            out[surf] = None
            continue
        pa, ga = angles(s, h, h_obs, R)
        fg_info = {"near_field_m": near,
                   "max_planar_angle_deg": float(np.nanmax(pa[fg])) if np.isfinite(pa[fg]).any() else None,
                   "max_h_above_observer_m": float(np.nanmax(h[fg]) - h_obs) if np.isfinite(h[fg]).any() else None}
        pa, ga = np.where(fg, np.nan, pa), np.where(fg, np.nan, ga)
        res = {"planar": _blocker(s, u, lats, lons, h, pa, planar_summit),
               "globe_geometric": _blocker(s, u, lats, lons, h, ga, globe_summit)}
        if atm is not None and solve is not None:
            gmax = np.nanmax(ga)
            cand = np.nonzero(ga >= gmax - ts.candidate_window_arcsec / 3600.0)[0]
            cand = cand[np.argsort(-ga[cand])][: ts.max_candidates]
            app = np.full(len(s), np.nan)
            for i in cand:
                app[i] = solve(lats[i], lons[i], float(h[i]))
            summit_app = solve(s_lat, s_lon, s_h)
            res["globe_refracted"] = _blocker(s, u, lats, lons, h, app, summit_app)
            res["globe_refracted"] and res["globe_refracted"].update({"candidates_solved": int(len(cand))})
            res["summit_apparent_deg"] = summit_app
        res["summit_planar_deg"], res["summit_globe_geometric_deg"] = planar_summit, globe_summit
        res["foreground"] = fg_info
        out[surf] = res
    return out, {"planar": None}


def compute_terrain(ev: dict, atm, ts: Optional[TerrainSettings] = None, offline=False, log=print):
    """Returns (terrain_doc for out/terrain/{id}.json, summary for results.json)."""
    from .refraction import apparent_geometric_target, geometric_from_apparent_infinite, wgs84_normal_section_radius
    ts = ts or TerrainSettings()
    o, p = ev["observer"], ev["peak"]
    obs = (o["lat"], o["lon"], o["alt_m"])
    summit = (p["lat"], p["lon"], p["alt_m"])
    inv = G.Inverse(obs[0], obs[1], summit[0], summit[1])
    D, azi_peak = float(inv["s12"]), float(inv["azi1"])
    R = wgs84_normal_section_radius(obs[0], azi_peak)

    def solve(la, lo, h):
        try:
            return apparent_geometric_target(atm, obs[0], obs[1], obs[2], la, lo, h).apparent_deg
        except Exception:  # noqa: BLE001
            return float("nan")

    doc = {"id": ev["id"], "frame": {
        "origin": "peak foot; tangent plane at the peak", "u_m": "s − D (negative toward the observer)",
        "planar": "x = u, z = h", "globe": "x = (R+h)·sin(u/R), z = (R+h)·cos(u/R) − R", "R_m": R, "D_m": D},
        "observer": {"planar_xz": [-D, obs[2]],
                     "globe_xz": [(R + obs[2]) * math.sin(-D / R), (R + obs[2]) * math.cos(-D / R) - R],
                     "globe_drop_m": obs[2] - ((R + obs[2]) * math.cos(-D / R) - R)},
        "peak_azimuth_deg": azi_peak, "profiles": {}, "settings": ts.__dict__}
    summary = {"file": f"terrain/{ev['id']}.json", "peak_azimuth_deg": azi_peak}

    star_az = (ev.get("astropy") or {}).get("az_obs_deg")
    lines = [("to_peak", azi_peak)]
    if star_az is not None:
        lines.append(("star_azimuth", star_az))
        d_az = ((star_az - azi_peak + 540) % 360) - 180
        summary["star_minus_peak_azimuth_arcsec"] = round(d_az * 3600, 1)
        summary["lateral_offset_at_peak_m"] = round(math.radians(d_az) * D, 1)

    for name, azi in lines:
        log(f"    terrain profile {name} ({(D + ts.beyond_peak_m) / 1000:.1f} km)")
        prof = build_profile(obs, azi, D, R, ts, offline=offline, log=log)
        obst, _ = obstruction(prof, obs, D, R, summit, atm, ts, solve, log=log)
        s = prof["s"]
        entry = {"azimuth_deg": azi, "sources": prof["sources"], "s_m": _r(s, 1), "u_m": _r(s - D, 1),
                 "lat": _r(prof["lats"], 6), "lon": _r(prof["lons"], 6), "obstruction": obst}
        for surf in ("dtm", "dsm"):
            h = prof[surf]
            f = frames(s, h, D, R)
            pa, ga = angles(s, h, obs[2], R)
            entry[surf] = {"h_m": _r(h, 1), "globe_x_m": _r(f["globe_x"], 1), "globe_z_m": _r(f["globe_z"], 2),
                           "planar_angle_deg": _r(pa, 6), "globe_angle_deg": _r(ga, 6)}
        doc["profiles"][name] = entry
        summary[name] = {"sources": prof["sources"], "obstruction": obst}

    # corridor (peak-anchored grid) — optional
    if not ts.corridor:
        doc["corridor"] = None
    sc = None if not ts.corridor else np.arange(0.0, D + min(ts.beyond_peak_m, 5000.0) + 1e-6, ts.corridor_step_m)
    if ts.corridor:
        tc = np.arange(-ts.corridor_half_width_m, ts.corridor_half_width_m + 1e-6, ts.corridor_step_m)
        line = G.Line(obs[0], obs[1], azi_peak)
        cl, co, ca = [], [], []
        for sv in sc:
            q = line.Position(float(sv))
            cl.append(q["lat2"]), co.append(q["lon2"]), ca.append(q["azi2"])
        LAT, LON = np.empty((len(sc), len(tc))), np.empty((len(sc), len(tc)))
        for j, tv in enumerate(tc):
            if tv == 0:
                LAT[:, j], LON[:, j] = cl, co
            else:
                LAT[:, j], LON[:, j] = offset_points(cl, co, ca, tv)
        corr = {"s_m": _r(sc, 1), "u_m": _r(sc - D, 1), "t_m": _r(tc, 1)}
        for surf in ("dtm", "dsm"):
            v, used = sample_surface(surf, LAT.ravel(), LON.ravel(), grid_step_m=ts.corridor_step_m, offline=offline, log=log)
            H = v.reshape(LAT.shape)
            U = (sc - D)[:, None]
            corr[surf] = {"sources": used, "h_m": [_r(row, 1) for row in H],
                          "globe_z_m": [_r(row, 2) for row in (R + H) * np.cos(U / R) - R],
                          "globe_x_m": [_r(row, 1) for row in (R + H) * np.sin(U / R)]}
        doc["corridor"] = corr

    # DEM vs CSV heights
    chk = {}
    for surf in ("dtm", "dsm"):
        v, _ = sample_surface(surf, [summit[0], obs[0]], [summit[1], obs[1]], offline=offline, log=log)
        chk[f"{surf}_at_peak_m"] = None if not np.isfinite(v[0]) else round(float(v[0]), 1)
        chk[f"{surf}_at_observer_m"] = None if not np.isfinite(v[1]) else round(float(v[1]), 1)
    chk["csv_peak_alt_m"], chk["csv_observer_alt_m"] = summit[2], obs[2]
    doc["checks"] = summary["checks"] = chk

    # terrain-based comparison along the star's azimuth (the skyline that actually hid the star)
    a = ev.get("astropy") or {}
    comp = {}
    sa = summary.get("star_azimuth", {}).get("obstruction") if "star_azimuth" in summary else None
    if sa and a.get("true_alt_obs_deg") is not None:
        rate = a.get("alt_rate_deg_per_s")
        for surf in ("dtm", "dsm"):
            r = sa.get(surf)
            if not r:
                continue
            if r.get("globe_refracted") and math.isfinite(r["globe_refracted"]["angle_deg"]):
                star_true = geometric_from_apparent_infinite(atm, r["globe_refracted"]["angle_deg"], obs[2], R)
                d = a["true_alt_obs_deg"] - star_true
                sfx = "" if surf == "dtm" else "_dsm"
                comp[f"globe_ciddor{sfx}_arcsec"] = round(d * 3600, 2)
                comp[f"globe_ciddor{sfx}_seconds"] = round(d / rate, 2) if rate else None
                comp[f"globe_skyline{sfx}_star_true_deg"] = star_true
            if r.get("planar"):
                sfx = "" if surf == "dtm" else "_dsm"
                d = a["true_alt_obs_deg"] - r["planar"]["angle_deg"]
                comp[f"planar_geometric{sfx}_arcsec"] = round(d * 3600, 2)
                comp[f"planar_geometric{sfx}_seconds"] = round(d / rate, 2) if rate else None
                if a.get("erfa_apparent_alt_obs_deg") is not None:
                    d2 = a["erfa_apparent_alt_obs_deg"] - r["planar"]["angle_deg"]
                    comp[f"planar_refracted{sfx}_arcsec"] = round(d2 * 3600, 2)
                    comp[f"planar_refracted{sfx}_seconds"] = round(d2 / rate, 2) if rate else None
    summary["comparison"] = comp

    # canopy check: DTM (bare earth) is primary; flag lines where the DSM skyline differs enough to matter
    summary["canopy"] = {}
    for name in doc["profiles"]:
        ob = summary[name]["obstruction"]
        gt = ((ob.get("dtm") or {}).get("globe_refracted") or (ob.get("dtm") or {}).get("globe_geometric"))
        gs = ((ob.get("dsm") or {}).get("globe_refracted") or (ob.get("dsm") or {}).get("globe_geometric"))
        if gt and gs:
            d = (gs["angle_deg"] - gt["angle_deg"]) * 3600
            summary["canopy"][name] = {"dsm_minus_dtm_skyline_arcsec": round(d, 2), "dsm_minus_dtm_h_at_blocker_m": round(gs["h_m"] - gt["h_m"], 1),
                                       "flag": bool(abs(d) > ts.canopy_flag_arcsec)}
    doc["canopy"] = summary["canopy"]
    return doc, summary
