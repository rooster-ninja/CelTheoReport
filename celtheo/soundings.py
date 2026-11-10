"""Radiosonde soundings: nearest station, nearest 00Z/12Z launch, fetch + cache.

Sources, tried in order:
  1. A local file in data/soundings/ named ``{wmo}_{YYYY-MM-DDTHHZ}.csv`` (committed or cached).
  2. Iowa Environmental Mesonet RAOB archive (CSV).
  3. University of Wyoming archive (TEXT:LIST).
If all fail, the event falls back to the ambiance standard atmosphere and is flagged.
"""
from __future__ import annotations

import io
import json
import math
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Optional

import pandas as pd

from .refraction import Sounding

ROOT = Path(__file__).resolve().parent.parent
SOUNDING_DIR = ROOT / "data" / "soundings"
STATIONS = json.loads((ROOT / "config" / "stations.json").read_text())
OVERRIDES = json.loads((ROOT / "config" / "sounding_overrides.json").read_text())


def _hav_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    a = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * math.asin(math.sqrt(a))


def nearest_station(lat, lon):
    best = min(STATIONS, key=lambda s: _hav_km(lat, lon, s["lat"], s["lon"]))
    return best, _hav_km(lat, lon, best["lat"], best["lon"])


def nearest_launch(utc: datetime) -> datetime:
    """Nearest synoptic launch (00Z or 12Z)."""
    base = utc.replace(minute=0, second=0, microsecond=0, hour=0)
    cands = [base + timedelta(hours=h) for h in (-12, 0, 12, 24)]
    return min(cands, key=lambda t: abs((t - utc).total_seconds()))


def rh_from_dewpoint(t_c, td_c):
    a, b = 17.625, 243.04  # Magnus
    return 100.0 * math.exp(a * td_c / (b + td_c) - a * t_c / (b + t_c))


def _to_sounding(df: pd.DataFrame, meta: dict) -> Sounding:
    df = df.dropna(subset=["height_m", "pres_hpa", "temp_c"]).sort_values("height_m")
    df = df.drop_duplicates("height_m")
    if "rh_pct" not in df or df["rh_pct"].isna().all():
        df["rh_pct"] = [rh_from_dewpoint(t, d) if pd.notna(d) else 0.0 for t, d in zip(df["temp_c"], df.get("dwpt_c", [float("nan")] * len(df)))]
    df["rh_pct"] = df["rh_pct"].fillna(0.0).clip(0, 100)
    return Sounding(df["height_m"].tolist(), df["pres_hpa"].tolist(), df["temp_c"].tolist(), df["rh_pct"].tolist(), **meta)


def parse_wyoming_csv(text: str) -> pd.DataFrame:
    """CSV with Wyoming-style headers ('PRES hPa','HGHT m','TEMP C','DWPT C','RELH %')."""
    df = pd.read_csv(io.StringIO(text.lstrip("﻿")))
    cols = {c: c.strip().upper() for c in df.columns}
    df = df.rename(columns=cols)
    pick = lambda *names: next((c for c in df.columns if any(c.startswith(n) for n in names)), None)
    out = pd.DataFrame({
        "pres_hpa": pd.to_numeric(df[pick("PRES")], errors="coerce"),
        "height_m": pd.to_numeric(df[pick("HGHT", "GEOPOTENTIAL", "HEIGHT")], errors="coerce"),
        "temp_c": pd.to_numeric(df[pick("TEMP")], errors="coerce"),
    })
    if pick("DWPT", "DEW"):
        out["dwpt_c"] = pd.to_numeric(df[pick("DWPT", "DEW")], errors="coerce")
    if pick("RELH"):
        out["rh_pct"] = pd.to_numeric(df[pick("RELH")], errors="coerce")
    return out


def parse_wyoming_textlist(html: str) -> pd.DataFrame:
    m = re.search(r"<PRE>(.*?)</PRE>", html, re.S | re.I)
    if not m:
        raise ValueError("no sounding table in Wyoming response")
    rows = []
    for line in m.group(1).splitlines():
        parts = line.split()
        if len(parts) >= 5 and re.match(r"^-?\d", parts[0]):
            try:
                p, h, t = float(line[0:7]), float(line[7:14]), float(line[14:21])
                td = float(line[21:28]) if line[21:28].strip() else float("nan")
                rh = float(line[28:35]) if line[28:35].strip() else float("nan")
                rows.append((p, h, t, td, rh))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["pres_hpa", "height_m", "temp_c", "dwpt_c", "rh_pct"])


def parse_iem_csv(text: str) -> pd.DataFrame:
    df = pd.read_csv(io.StringIO(text))
    df.columns = [c.strip().lower() for c in df.columns]
    return pd.DataFrame({
        "pres_hpa": pd.to_numeric(df["pressure_mb"], errors="coerce"),
        "height_m": pd.to_numeric(df["height_m"], errors="coerce"),
        "temp_c": pd.to_numeric(df["tmpc"], errors="coerce"),
        "dwpt_c": pd.to_numeric(df["dwpc"], errors="coerce"),
    })


def _cache_path(wmo, t: datetime):
    return SOUNDING_DIR / f"{wmo}_{t:%Y-%m-%dT%H}Z.csv"


def fetch(station: dict, t: datetime, offline=False, log=print) -> Optional[pd.DataFrame]:
    path = _cache_path(station["wmo"], t)
    if path.exists():
        return parse_wyoming_csv(path.read_text())
    if offline:
        return None
    import requests
    tries = [
        ("IEM", f"https://mesonet.agron.iastate.edu/cgi-bin/request/raob.py?station={station['iem']}"
                f"&sts={t:%Y-%m-%dT%H:%M}Z&ets={t:%Y-%m-%dT%H:%M}Z", parse_iem_csv),
        ("Wyoming", f"https://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST"
                    f"&YEAR={t:%Y}&MONTH={t:%m}&FROM={t:%d%H}&TO={t:%d%H}&STNM={station['wmo']}", parse_wyoming_textlist),
    ]
    for name, url, parser in tries:
        try:
            r = requests.get(url, timeout=40, headers={"User-Agent": "celtheo-pipeline"})
            r.raise_for_status()
            df = parser(r.text)
            if len(df.dropna(subset=["height_m", "temp_c"])) < 10:
                raise ValueError("too few levels")
            SOUNDING_DIR.mkdir(parents=True, exist_ok=True)
            out = df.rename(columns={"pres_hpa": "PRES hPa", "height_m": "HGHT m", "temp_c": "TEMP C",
                                     "dwpt_c": "DWPT C", "rh_pct": "RELH %"})
            out.to_csv(path, index=False)
            log(f"    sounding {station['wmo']} {t:%Y-%m-%d %HZ} from {name} ({len(df)} levels)")
            return df
        except Exception as e:  # noqa: BLE001 — any failure falls through to the next source
            log(f"    {name} sounding fetch failed for {station['wmo']} {t:%Y-%m-%d %HZ}: {e}")
    return None


def sounding_for_event(event_id: str, obs_lat, obs_lon, utc: datetime, offline=False, log=print) -> Sounding:
    ov = OVERRIDES.get(event_id, {})
    if ov.get("wmo"):
        st = next(s for s in STATIONS if s["wmo"] == ov["wmo"])
        dist = _hav_km(obs_lat, obs_lon, st["lat"], st["lon"])
    else:
        st, dist = nearest_station(obs_lat, obs_lon)
    t = datetime.fromisoformat(ov["launch_utc"]).replace(tzinfo=timezone.utc) if ov.get("launch_utc") else nearest_launch(utc)
    df = fetch(st, t, offline=offline, log=log)
    if df is None:
        return Sounding([], [], [], [], station=st["name"], wmo=st["wmo"], valid_utc=f"{t:%Y-%m-%dT%H:%MZ}",
                        source="standard-atmosphere fallback", distance_km=dist)
    return _to_sounding(df, dict(station=st["name"], wmo=st["wmo"], valid_utc=f"{t:%Y-%m-%dT%H:%MZ}",
                                 source="radiosonde", distance_km=dist))
