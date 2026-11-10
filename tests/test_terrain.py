"""Terrain module on a synthetic DEM (no network): frames, obstruction, star-azimuth skyline."""
import math
from pathlib import Path

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
from rasterio.transform import from_origin

from celtheo import terrain
from celtheo.refraction import CiddorAtmosphere, RefractionSettings
from celtheo.soundings import _to_sounding, parse_wyoming_csv

ROOT = Path(__file__).resolve().parent.parent
OBS = (38.9997, -104.49723, 2073.0)
PEAK = (38.840278, -105.044444, 4300.0)


def make_dem(path, ridge=None, canopy=0.0):
    """EPSG:4326 grid: 2000 m plateau + Gaussian peak at PEAK (summit 4300 m) [+ optional ridge]."""
    res = 0.0005  # ~45 m
    west, north = -105.40, 39.10
    w, h = int(1.0 / res), int(0.6 / res)
    lon = west + (np.arange(w) + 0.5) * res
    lat = north - (np.arange(h) + 0.5) * res
    LON, LAT = np.meshgrid(lon, lat)
    kx = 111320 * math.cos(math.radians(39))
    d2 = ((LON - PEAK[1]) * kx) ** 2 + ((LAT - PEAK[0]) * 111320) ** 2
    z = 2000.0 + 2300.0 * np.exp(-d2 / (2 * 3000.0 ** 2))
    if ridge:
        rl, ro, rh = ridge
        z = np.maximum(z, 2000.0 + rh * np.exp(-(((LON - ro) * kx) ** 2 + ((LAT - rl) * 111320) ** 2) / (2 * 2500.0 ** 2)))
    z = z + canopy
    with rasterio.open(path, "w", driver="GTiff", height=h, width=w, count=1, dtype="float32", crs="EPSG:4326",
                       transform=from_origin(west, north, res, res), nodata=-9999) as ds:
        ds.write(z.astype("float32"), 1)


@pytest.fixture
def atm():
    s = _to_sounding(parse_wyoming_csv((ROOT / "data/soundings/72476_2025-01-28T00Z.csv").read_text()), {})
    return CiddorAtmosphere(s, RefractionSettings())


def use_dems(monkeypatch, dtm, dsm):
    monkeypatch.setattr(terrain, "SOURCES", {
        "dtm": [{"name": "synthetic DTM", "tiling": "single", "local": True, "url_templates": [str(dtm)]}],
        "dsm": [{"name": "synthetic DSM", "tiling": "single", "local": True, "url_templates": [str(dsm)]}]})
    terrain.Raster._open.clear()
    terrain.Raster._bad.clear()


def event(az=None):
    return {"id": "synthetic", "observer": {"lat": OBS[0], "lon": OBS[1], "alt_m": OBS[2]},
            "peak": {"lat": PEAK[0], "lon": PEAK[1], "alt_m": PEAK[2]},
            "astropy": {"az_obs_deg": az, "true_alt_obs_deg": 2.10, "alt_rate_deg_per_s": -0.003,
                        "erfa_apparent_alt_obs_deg": 2.33} if az else {}}


def test_frames_and_summit_occludes(tmp_path, monkeypatch, atm):
    make_dem(tmp_path / "dtm.tif")
    make_dem(tmp_path / "dsm.tif", canopy=20.0)
    use_dems(monkeypatch, tmp_path / "dtm.tif", tmp_path / "dsm.tif")
    ts = terrain.TerrainSettings(profile_step_m=50, beyond_peak_m=10000, corridor=True, corridor_step_m=250, corridor_half_width_m=1000)
    doc, summ = terrain.compute_terrain(event(), atm, ts, offline=True, log=lambda *a: None)
    D, R = doc["frame"]["D_m"], doc["frame"]["R_m"]
    # globe observer drops by ~D²/2R; peak is shared by both frames
    assert abs(doc["observer"]["globe_drop_m"] - D * D / (2 * R)) < 1.0
    prof = doc["profiles"]["to_peak"]
    i = int(np.argmin(np.abs(np.array(prof["u_m"]))))
    assert abs(prof["dtm"]["globe_z_m"][i] - prof["dtm"]["h_m"][i]) < 0.5      # at the peak, globe z == h
    assert abs(summ["checks"]["dtm_at_peak_m"] - 4300) < 15
    for surf in ("dtm", "dsm"):
        ob = summ["to_peak"]["obstruction"][surf]
        assert ob["planar"]["near_peak"] and ob["globe_geometric"]["near_peak"] and ob["globe_refracted"]["near_peak"]
    # canopy raises the refracted skyline by ~ 20 m / D
    a_dtm = summ["to_peak"]["obstruction"]["dtm"]["globe_refracted"]["angle_deg"]
    a_dsm = summ["to_peak"]["obstruction"]["dsm"]["globe_refracted"]["angle_deg"]
    assert abs((a_dsm - a_dtm) * 3600 - math.degrees(20 / D) * 3600) < 15
    assert len(doc["corridor"]["dtm"]["h_m"]) == len(doc["corridor"]["s_m"])
    # 20 m of canopy at ~50 km is ~80 arcsec on the skyline -> flagged for review
    assert summ["canopy"]["to_peak"]["flag"] is True


def test_ridge_beyond_peak_takes_over(tmp_path, monkeypatch, atm):
    # a 5200 m summit 12 km behind the peak on the same line (tall enough to beat the summit's angle)
    from geographiclib.geodesic import Geodesic
    inv = Geodesic.WGS84.Inverse(OBS[0], OBS[1], PEAK[0], PEAK[1])
    far = Geodesic.WGS84.Direct(OBS[0], OBS[1], inv["azi1"], inv["s12"] + 12000)
    make_dem(tmp_path / "dtm.tif", ridge=(far["lat2"], far["lon2"], 3200.0))
    make_dem(tmp_path / "dsm.tif", ridge=(far["lat2"], far["lon2"], 3200.0))
    use_dems(monkeypatch, tmp_path / "dtm.tif", tmp_path / "dsm.tif")
    ts = terrain.TerrainSettings(profile_step_m=50, beyond_peak_m=20000, corridor_step_m=500, corridor_half_width_m=500)
    doc, summ = terrain.compute_terrain(event(az=inv["azi1"] + 0.01), atm, ts, offline=True, log=lambda *a: None)
    ob = summ["to_peak"]["obstruction"]["dtm"]
    assert ob["planar"]["u_m"] > 10000 and not ob["planar"]["near_peak"]
    assert ob["globe_refracted"]["margin_vs_summit_arcsec"] > 0
    assert "globe_ciddor_arcsec" in summ["comparison"] and "planar_geometric_arcsec" in summ["comparison"]
    assert "globe_ciddor_dsm_arcsec" in summ["comparison"]
