"""Regression tests. Run with: pytest -q"""
import json
import math
from pathlib import Path

import pytest

from celtheo.refraction import (CiddorAtmosphere, RefractionSettings, apparent_geometric_target, bennett_arcmin,
                                geometric_from_apparent_infinite, wgs84_gauss_radius)
from celtheo.soundings import _to_sounding, parse_wyoming_csv

ROOT = Path(__file__).resolve().parent.parent
OBS = (38.9997, -104.49723, 2073.0)
PIKES = (38.840278, -105.044444, 4300.0)


def dms_s(deg):
    return deg * 3600.0


@pytest.fixture(scope="module")
def solved():
    s = _to_sounding(parse_wyoming_csv((ROOT / "data/soundings/72476_2025-01-28T00Z.csv").read_text()), {})
    atm = CiddorAtmosphere(s, RefractionSettings())
    sol = apparent_geometric_target(atm, *OBS, *PIKES)
    return atm, sol


def test_matches_al_refraction_analysis_readme(solved):
    """Numbers published in rooster-ninja/AL_RefractionAnalysis README (Pikes Peak / 39 Aqr)."""
    atm, sol = solved
    assert abs(dms_s(sol.apparent_deg) - dms_s(2 + 19 / 60 + 4.4 / 3600)) < 0.2
    assert abs(dms_s(sol.geometric_deg) - dms_s(2 + 17 / 60 + 20.6 / 3600)) < 0.2
    assert abs(dms_s(sol.phi_deg) - dms_s(27 / 60 + 16.6 / 3600)) < 0.2
    star = geometric_from_apparent_infinite(atm, sol.apparent_deg, OBS[2], wgs84_gauss_radius(OBS[0]))
    assert abs((sol.apparent_deg - star) * 60 - 13.33) < 0.005
    P = 1013.25 * (1 - 2.25577e-5 * OBS[2]) ** 5.25588
    assert abs(bennett_arcmin(sol.apparent_deg, P, -0.2) - 13.64) < 0.005


def test_offline_pipeline_with_gp_fixture(tmp_path, monkeypatch):
    """End-to-end offline run on the Pikes Peak event. The star position is rebuilt from the CSV's
    GP lat/long (west-positive longitude), so astropy should land within a few arcsec of the CSV."""
    import astropy.units as u
    from astropy.coordinates import TETE, SkyCoord
    from astropy.time import Time
    from celtheo import api, stars
    from celtheo.compute import configure_astropy

    configure_astropy(True)
    t = Time("2025-01-28T01:43:07")
    ra = (t.sidereal_time("apparent", longitude=0).deg - 179.92305) % 360
    c = SkyCoord(ra=ra * u.deg, dec=-14.072628 * u.deg, frame=TETE(obstime=t)).transform_to("icrs")
    cache = tmp_path / "stars.json"
    cache.write_text(json.dumps({"39 Aqr": {"query": "39 Aqr", "kind": "star", "ra_deg": c.ra.deg, "dec_deg": c.dec.deg,
                                            "pmra_masyr": 0, "pmdec_masyr": 0, "source": "fixture"}}))
    monkeypatch.setattr(stars, "CACHE", cache)
    doc = api.run(ROOT / "data/occultations.csv", out_dir=tmp_path, offline=True, only=["pikes-peak-39-aquarii-2025-01-28"],
                  include_tests=False, log=lambda *a: None)
    ev = doc["events"][0]
    assert abs(ev["astropy"]["true_alt_obs_deg"] - 2.1024) * 3600 < 10
    assert ev["refraction"]["sounding"]["station"].startswith("Grand Junction")
    assert set(ev["comparison"]) >= {"globe_ciddor_arcsec", "planar_geometric_arcsec"}
    assert (tmp_path / "results.json").exists() and (tmp_path / "results.csv").exists()
