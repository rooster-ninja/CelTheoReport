"""Public API. This is the hook for report generators.

    from celtheo.api import run, load_results
    results = run("data/occultations.csv", out_dir="out")      # compute + write files
    results = load_results("out/results.json")                 # read a previous run

``results`` is a plain dict matching schema/results.schema.json.
"""
from __future__ import annotations

import csv
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
    name_occluders(results)
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
        b = ((((t.get("star_azimuth") or {}).get("obstruction") or {}).get("dtm") or {}).get("globe_refracted")) or {}
        if b.get("nearest_named_peak") and b["nearest_named_peak"]["name"] != ev["peak"]["name"]:
            n = b["nearest_named_peak"]
            ev["flags"]["warnings"].append(f"along the star's azimuth the skyline is {b['u_m'] / 1000:+.1f} km from {ev['peak']['name']}, "
                                           f"{n['distance_m'] / 1000:.1f} km from {n['name']}: likely the actual occluder")


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


def load_results(path="out/results.json") -> dict:
    return json.loads(Path(path).read_text())
