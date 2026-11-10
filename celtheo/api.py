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

SCHEMA_VERSION = "1.0"


def run(csv_path="data/occultations.csv", out_dir: Optional[str] = "out", offline=False, only: Optional[Iterable[str]] = None,
        include_tests=True, profiles=True, settings: Optional[RefractionSettings] = None, log: Callable = print) -> dict:
    settings = settings or RefractionSettings()
    configure_astropy(offline)
    events = load_events(csv_path)
    if only:
        keys = [k.lower() for k in only]
        events = [e for e in events if any(k in e.id for k in keys)]
    if not include_tests:
        events = [e for e in events if not e.test_kind]
    results = []
    t0 = time.time()
    for i, e in enumerate(events, 1):
        log(f"[{i}/{len(events)}] {e.id}")
        try:
            results.append(compute_event(e, settings, offline=offline, log=log, profiles=profiles))
        except Exception as ex:  # noqa: BLE001 — one bad row never stops the run
            log(f"    FAILED: {ex}")
            results.append({"id": e.id, "row": e.row, "flags": {"warnings": [f"computation failed: {ex}"], "target_type": e.target_type,
                            "is_test": e.test_kind is not None, "test_kind": e.test_kind},
                            "peak": {"name": e.peak}, "observer": {"name": e.observer}, "star": {"name": e.star}})
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
        },
        "runtime_s": round(time.time() - t0, 1),
        "events": results,
    }
    if out_dir:
        write(doc, out_dir)
    return doc


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
    (out / "results.json").write_text(json.dumps(doc, indent=1, ensure_ascii=False))
    with open(out / "results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([k for k, _ in FLAT])
        for ev in doc["events"]:
            w.writerow([_get(ev, p) for _, p in FLAT])


def load_results(path="out/results.json") -> dict:
    return json.loads(Path(path).read_text())
