"""Public API. This is the hook for report generators.

    from celtheo.api import run, load_results
    results = run("data/occultations.csv", out_dir="out")      # compute + write files
    results = load_results("out/results.json")                 # read a previous run

``results`` is a plain dict matching schema/results.schema.json.
"""
from __future__ import annotations

import csv
import math
import json
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Optional

from . import __version__
from .compute import compute_event, configure_astropy
from .data import load_events
from .refraction import RefractionSettings

SCHEMA_VERSION = "1.1"


def run(csv_path="data/occultations.csv", out_dir: Optional[str] = "out", offline=False, only: Optional[Iterable[str]] = None,
        include_tests=True, profiles=True, refresh_stars=False, terrain=True, settings: Optional[RefractionSettings] = None, log: Callable = print) -> dict:
    settings = settings or RefractionSettings()
    configure_astropy(offline)
    events = load_events(csv_path)
    if only:
        keys = [k.lower() for k in only]
        events = [e for e in events if any(k in e.id for k in keys)]
    if not include_tests:
        events = [e for e in events if not e.test_kind]
    excl_path = Path(__file__).resolve().parent.parent / "config" / "exclusions.json"
    excluded = {k: v for k, v in (json.loads(excl_path.read_text()).items() if excl_path.exists() else []) if not k.startswith("_")}
    results = []
    t0 = time.time()
    for i, e in enumerate(events, 1):
        if e.id in excluded:
            log(f"[{i}/{len(events)}] {e.id}  (excluded: {excluded[e.id]})")
            results.append({"id": e.id, "row": e.row,
                            "flags": {"excluded": True, "excluded_reason": excluded[e.id], "warnings": [], "target_type": e.target_type,
                                      "is_test": e.test_kind is not None, "test_kind": e.test_kind},
                            "peak": {"name": e.peak}, "observer": {"name": e.observer}, "star": {"name": e.star},
                            "time": {"local_date": e.local_date, "local_time": e.local_time}})
            continue
        log(f"[{i}/{len(events)}] {e.id}")
        try:
            results.append(compute_event(e, settings, offline=offline, log=log, profiles=profiles, refresh_stars=refresh_stars, terrain=terrain))
        except Exception as ex:  # noqa: BLE001 — one bad row never stops the run
            log(f"    FAILED: {ex}")
            results.append({"id": e.id, "row": e.row, "flags": {"warnings": [f"computation failed: {ex}"], "target_type": e.target_type,
                            "is_test": e.test_kind is not None, "test_kind": e.test_kind},
                            "peak": {"name": e.peak}, "observer": {"name": e.observer}, "star": {"name": e.star}})
    if only and out_dir and (Path(out_dir) / "results.json").exists():
        prev = load_results(Path(out_dir) / "results.json")
        fresh = {r["id"]: r for r in results}
        log(f"merging {len(fresh)} recomputed event(s) into {out_dir}/results.json")
        results = [fresh.pop(r["id"], r) for r in prev["events"]] + list(fresh.values())
    finalize(results)
    doc = {
        "schema_version": SCHEMA_VERSION,
        "generator": f"celtheo {__version__}",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_csv": str(csv_path),
        "settings": asdict(settings),
        "conventions": {
            "angles": "degrees unless the key ends in _arcmin/_arcsec",
            "true_alt": "geometric (vacuum) altitude above the observer's ellipsoid horizon (astropy, pressure=0)",
            "comparison_sign": "observed-minus-model: positive means the star was truly higher than the model's occlusion angle",
            "seconds": "residual divided by the star's altitude rate; the clock offset that would close the residual",
            "frame": "peak-anchored: origin at the peak foot, tangent plane at the peak. Planar x=u, z=h; globe x=(R+h)sin(u/R), z=(R+h)cos(u/R)−R, u = s − D",
            "terrain": "dtm = bare earth (primary; unsuffixed terrain comparisons), dsm = top of canopy (suffix _dsm, plus a canopy flag); per-event arrays in terrain/{id}.json",
        },
        "runtime_s": round(time.time() - t0, 1),
        "events": results,
    }
    if out_dir:
        write(doc, out_dir)
    return doc


SUMMARY_BLOCKERS = ("CSV true altitude", "star not resolved", "CSV apparent angle at peak is below", "computation failed",
                    "astropy failed", "could not parse", "Moon:")


def finalize(results):
    """Post-processing that needs the whole set: occluder naming, terrain warnings, notes, in_summary.
    Safe to re-run on an existing results.json (celtheo finalize)."""
    for ev in results:
        ev.setdefault("flags", {}).setdefault("warnings", [])
        ev["flags"]["warnings"] = [w for w in ev["flags"]["warnings"]
                                   if not (w.startswith("terrain ") and "rises above the summit" in w) and not w.startswith("along the star's azimuth")
                                   and not w.startswith("Moon:")]
    name_occluders(results)
    for ev in results:
        f = ev["flags"]
        t = ev.get("terrain") or {}
        for line in ("star_azimuth", "to_peak"):
            b = ((((t.get(line) or {}).get("obstruction") or {}).get("dtm") or {}).get("globe_refracted")) or {}
            if not b or b.get("near_peak", True):
                continue
            n = b.get("nearest_named_peak") or {}
            if n.get("name") == (ev.get("peak") or {}).get("name"):
                b["same_peak"] = True  # same mountain: the CSV summit coordinate is just imprecise
                continue
            if n:
                f["warnings"].append(f"along the star's azimuth the skyline is {b['u_m'] / 1000:+.1f} km from {ev['peak']['name']}, "
                                     f"{n['distance_m'] / 1000:.1f} km from {n['name']}: likely the actual occluder")
            else:
                f["warnings"].append(f"terrain {b['u_m'] / 1000:+.1f} km from the peak rises above the summit on the {line.replace('_', ' ')} line")
            break
        # two-station quantities for the celestial-navigation view (O1 at the observer, O2 at the peak)
        cm, ap, cs, gm = ev.get("comparison"), ev.get("astropy") or {}, ev.get("csv") or {}, ev.get("geometry") or {}
        if cm is not None:
            num = lambda x: isinstance(x, (int, float))
            if num(ap.get("true_alt_peak_deg")) and num(ap.get("true_alt_obs_deg")):
                cm["star_shift_astropy_deg"] = ap["true_alt_peak_deg"] - ap["true_alt_obs_deg"]
                if num(gm.get("central_angle_deg")):
                    cm["star_shift_astropy_minus_central_angle_arcsec"] = round((cm["star_shift_astropy_deg"] - gm["central_angle_deg"]) * 3600, 2)
            if num(cs.get("true_peak")) and num(cs.get("true_obs")):
                cm["star_shift_csv_deg"] = cs["true_peak"] - cs["true_obs"]
            if num(cs.get("true_peak")) and num(gm.get("planar_angle_deg")):
                cm["ho_csv_deg"] = cs["true_peak"] - gm["planar_angle_deg"]
        # drop as an angle: tangent-chord (inscribed angle) theorem gives gamma/2; the frame coordinates give it directly too
        if isinstance(gm.get("central_angle_deg"), (int, float)):
            gm["drop_angle_deg"] = gm["central_angle_deg"] / 2.0
            og, op = gm.get("observer_globe_xz_m"), gm.get("observer_planar_xz_m")
            if og and op and og[0]:
                gm["drop_angle_frame_deg"] = math.degrees(math.atan2(op[1] - og[1], abs(og[0])))
        rf = ev.get("refraction") or {}
        app = rf.get("star_apparent_ciddor_obs_deg")
        if cm is not None and isinstance(app, (int, float)):
            if isinstance(rf.get("peak_apparent_deg"), (int, float)):
                cm["globe_ciddor_apparent_arcsec"] = round((app - rf["peak_apparent_deg"]) * 3600, 2)
            tc = (ev.get("terrain") or {}).get("comparison")
            sky = ((((ev.get("terrain") or {}).get("star_azimuth") or {}).get("obstruction") or {}).get("dtm") or {}).get("globe_refracted") or {}
            if tc is not None and isinstance(sky.get("angle_deg"), (int, float)):
                tc["globe_ciddor_apparent_arcsec"] = round((app - sky["angle_deg"]) * 3600, 2)
        if (ev.get("star") or {}).get("kind") == "body":
            f["warnings"].append("Moon: astropy gives the centre; the occultation was timed on a limb (about ±15′), so residuals are not comparable with stars")
        f["in_summary"] = bool(not f.get("excluded") and not f.get("is_test") and f.get("target_type") == "peak"
                               and not any(w.startswith(SUMMARY_BLOCKERS) for w in f["warnings"]))


def name_occluders(results, radius_m=2000.0):
    """Where terrain other than the named peak hid the star, name the nearest peak from the CSV (if within radius)."""
    from geographiclib.geodesic import Geodesic
    peaks = {}
    for ev in results:
        p = ev.get("peak") or {}
        if p.get("lat") is not None and (ev.get("flags") or {}).get("target_type") == "peak":
            peaks.setdefault(p["name"], (p["lat"], p["lon"], p.get("alt_m")))
    for ev in results:
        t = ev.get("terrain") or {}
        for line in ("star_azimuth", "to_peak"):
            for surf in ("dtm", "dsm"):
                ob = (((t.get(line) or {}).get("obstruction") or {}).get(surf) or {})
                for model in ("globe_refracted", "globe_geometric", "planar"):
                    b = ob.get(model)
                    if not b or b.get("near_peak", True):
                        continue
                    best = min(((Geodesic.WGS84.Inverse(b["lat"], b["lon"], la, lo)["s12"], nm, h) for nm, (la, lo, h) in peaks.items()), default=None)
                    if best and best[0] <= radius_m:
                        b["nearest_named_peak"] = {"name": best[1], "distance_m": round(best[0], 0), "csv_alt_m": best[2]}


FLAT = [
    ("id", "id"), ("peak", "peak.name"), ("star", "star.name"), ("observer", "observer.name"), ("utc", "time.utc"),
    ("target_type", "flags.target_type"), ("test_kind", "flags.test_kind"),
    ("distance_m", "geometry.geodesic_distance_m"), ("planar_angle_deg", "geometry.planar_angle_deg"),
    ("peak_apparent_deg", "refraction.peak_apparent_deg"), ("peak_geometric_deg", "refraction.peak_geometric_deg"),
    ("star_geometric_ciddor_deg", "refraction.star_geometric_deg"), ("star_uplift_ciddor_arcmin", "refraction.star_uplift_arcmin"),
    ("bennett_arcmin", "refraction.bennett_arcmin"), ("astropy_true_alt_obs_deg", "astropy.true_alt_obs_deg"),
    ("astropy_true_alt_peak_deg", "astropy.true_alt_peak_deg"), ("csv_true_obs_deg", "csv.true_obs"),
    ("csv_app_obs_deg", "csv.app_obs"), ("globe_ciddor_arcsec", "comparison.globe_ciddor_arcsec"),
    ("globe_ciddor_seconds", "comparison.globe_ciddor_seconds"), ("planar_geometric_arcsec", "comparison.planar_geometric_arcsec"),
    ("planar_refracted_arcsec", "comparison.planar_refracted_arcsec"), ("sounding_station", "refraction.sounding.station"),
    ("sounding_utc", "refraction.sounding.valid_utc"), ("sounding_source", "refraction.sounding.source"),
    ("globe_terrain_arcsec", "terrain.comparison.globe_ciddor_arcsec"),
    ("globe_terrain_seconds", "terrain.comparison.globe_ciddor_seconds"),
    ("planar_terrain_arcsec", "terrain.comparison.planar_geometric_arcsec"),
    ("planar_terrain_seconds", "terrain.comparison.planar_geometric_seconds"),
    ("canopy_flag", "terrain.canopy.star_azimuth.flag"),
    ("globe_observer_drop_m", "geometry.globe_observer_drop_m"),
    ("excluded", "flags.excluded"),
    ("in_summary", "flags.in_summary"),
    ("warnings", "flags.warnings"),
]


def _get(d, path):
    for k in path.split("."):
        if not isinstance(d, dict):
            return None
        d = d.get(k)
    return "; ".join(d) if isinstance(d, list) else d


def write(doc: dict, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for ev in doc["events"]:  # per-event terrain goes to its own file
        t = ev.pop("_terrain_doc", None)
        if t is not None:
            (out / "terrain").mkdir(exist_ok=True)
            (out / "terrain" / f"{ev['id']}.json").write_text(json.dumps(t, separators=(",", ":"), ensure_ascii=False))
    (out / "results.json").write_text(json.dumps(doc, indent=1, ensure_ascii=False))
    with open(out / "results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([k for k, _ in FLAT])
        for ev in doc["events"]:
            w.writerow([_get(ev, p) for _, p in FLAT])


def refinalize(out_dir="out"):
    """Re-apply finalize() to an existing run and rewrite results.json/.csv (no recomputation)."""
    doc = load_results(Path(out_dir) / "results.json")
    finalize(doc["events"])
    write(doc, out_dir)
    return doc


def add_star_refraction(out_dir="out", log=print):
    """Add the Ciddor (AstronomyLive) refraction of the actual star to an existing run, offline,
    using the cached soundings. Then re-finalize."""
    from .compute import star_refraction
    from .refraction import CiddorAtmosphere, RefractionSettings, wgs84_normal_section_radius
    from .refraction import EMPTY_SOUNDING
    from .soundings import SOUNDING_DIR, _to_sounding, parse_wyoming_csv
    doc = load_results(Path(out_dir) / "results.json")
    st = RefractionSettings(**doc.get("settings", {}))
    for i, ev in enumerate(doc["events"], 1):
        rf, g = ev.get("refraction") or {}, ev.get("geometry") or {}
        sd = rf.get("sounding")
        if ev["flags"].get("excluded") or not sd or not g.get("earth_radius_m") or not (ev.get("astropy") or {}).get("true_alt_obs_deg"):
            continue
        snd = EMPTY_SOUNDING
        if sd.get("source") == "radiosonde":
            f = SOUNDING_DIR / f"{sd['wmo']}_{sd['valid_utc'][:13]}Z.csv"
            if not f.exists():
                log(f"[{i}] {ev['id']}: sounding file {f.name} missing, skipped")
                continue
            snd = _to_sounding(parse_wyoming_csv(f.read_text()), {})
        atm = CiddorAtmosphere(snd, st)
        az = g.get("azimuth_deg")
        R_peak = wgs84_normal_section_radius(ev["peak"]["lat"], (az + 180.0) % 360.0) if az is not None else g["earth_radius_m"]
        for k in [k for k in rf if k.startswith(("star_apparent_ciddor", "star_refraction_", "star_ciddor_minus"))]:
            rf.pop(k)  # clear values from an earlier pass
        rf.update(star_refraction(atm, ev, g["earth_radius_m"], R_peak))
        log(f"[{i}] {ev['id']}: Ciddor {rf.get('star_refraction_ciddor_obs_arcmin')}′, Bennett {rf.get('star_refraction_bennett_obs_arcmin')}′, ERFA {rf.get('star_refraction_erfa_obs_arcmin')}′")
    finalize(doc["events"])
    write(doc, out_dir)
    return doc


def load_results(path="out/results.json") -> dict:
    return json.loads(Path(path).read_text())
