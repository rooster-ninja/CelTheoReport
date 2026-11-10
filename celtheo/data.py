"""Load and clean the Celestial Theodolite occultation table."""
from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import List, Optional

import pandas as pd

COLS = {
    "peak": "Peak Name", "star": "Star Name", "time": "Occlusion Time", "date": "Occlusion Date",
    "utc_offset": "UTC", "observer": "Observer Name", "peak_lat": "Peak Lat", "peak_lon": "Peak Long",
    "obs_lat": "Obs Lat", "obs_lon": "Obs Long", "peak_alt": "Peak Alt (m)", "obs_alt": "Obs Alt (m)",
    "distance": "Distance (m)", "true_peak": "Star->Peak ♈︎ (dd)", "app_peak": "Star->Peak ϵ (dd)",
    "true_obs": "Star->Obs ♈︎ (dd)", "app_obs": "Star->Obs ϵ (dd)", "canva": "Canva Links",
    "gp_lat": "Star GP Lat", "gp_lon": "Star GP Long", "gp_dist": "Star GP Dist (m)",
    "drop_rate": "Star Drop Rate (arcsec)", "ge_check": "GE Occlusion Check", "temp_c": "Temp (C)",
    "azi_peak": "Azi (Peak)",
}
AIRCRAFT_MIN_ALT_M = 8000.0


@dataclass
class Event:
    id: str
    row: int
    peak: str
    star: str
    observer: str
    local_time: str
    local_date: str
    utc_offset_h: Optional[float]
    utc: Optional[datetime]
    peak_lat: float
    peak_lon: float
    peak_alt: float
    obs_lat: float
    obs_lon: float
    obs_alt: float
    csv: dict
    target_type: str = "peak"          # peak | aircraft
    test_kind: Optional[str] = None    # None | "epsilon" | "aries"
    warnings: List[str] = field(default_factory=list)


def _num(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return float("nan")
    s = str(v).strip().replace(",", "")
    if s in ("", "#N/A", "nan"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _key(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip()).lower()


def slug(s: str) -> str:
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if c.isascii())
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def _parse_date(s: str):
    s = str(s).strip()
    for fmt in ("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            pass
    return None


def _parse_time(s: str):
    s = str(s).strip().replace(" PM", "").replace(" AM", "")
    for fmt in ("%H:%M:%S", "%H:%M"):
        try:
            return datetime.strptime(s, fmt).time()
        except ValueError:
            pass
    return None


def load_events(path) -> List[Event]:
    df = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    get = lambda r, k: r.get(COLS[k].strip(), None)
    rows = []
    for i, r in df.iterrows():
        peak, star = str(get(r, "peak") or "").strip(), str(get(r, "star") or "").strip()
        if not peak or peak == "nan":
            continue
        rows.append((i + 2, r, peak, star))  # +2 = spreadsheet row number (header is row 1)

    events: List[Event] = []
    for rownum, r, peak, star in rows:
        w = []
        test_kind = None
        if peak.startswith("ϵ"):
            test_kind, peak, star = "epsilon", peak.lstrip("ϵ").replace(" Test", ""), star.lstrip("ϵ").replace(" Test", "")
        elif peak.startswith("♈"):
            test_kind = "aries"
            peak = re.sub(r"^♈︎?", "", peak).replace(" Test", "")
            star = re.sub(r"^♈︎?", "", star).replace(" Test", "")
        off = _num(get(r, "utc_offset"))
        d, t = _parse_date(get(r, "date")), _parse_time(get(r, "time"))
        csv = {k: _num(get(r, k)) for k in ("true_peak", "app_peak", "true_obs", "app_obs", "distance", "gp_lat",
                                           "gp_lon", "gp_dist", "drop_rate", "temp_c", "azi_peak")}
        csv["ge_check"] = str(get(r, "ge_check") or "").strip() or None
        csv["canva"] = str(get(r, "canva") or "").strip() or None
        e = Event(
            id="", row=rownum, peak=peak, star=star, observer=str(get(r, "observer") or "").strip(),
            local_time=str(get(r, "time") or "").strip(), local_date=str(get(r, "date") or "").strip(),
            utc_offset_h=None if pd.isna(off) else off, utc=None,
            peak_lat=_num(get(r, "peak_lat")), peak_lon=_num(get(r, "peak_lon")), peak_alt=_num(get(r, "peak_alt")),
            obs_lat=_num(get(r, "obs_lat")), obs_lon=_num(get(r, "obs_lon")), obs_alt=_num(get(r, "obs_alt")),
            csv=csv, test_kind=test_kind, warnings=w,
        )
        e._d, e._t = d, t  # type: ignore[attr-defined]
        if e.peak_alt >= AIRCRAFT_MIN_ALT_M:
            e.target_type = "aircraft"
        if d is None or t is None:
            w.append("could not parse local date/time")
        if all(pd.isna(csv[k]) for k in ("true_peak", "app_peak", "true_obs", "app_obs")):
            w.append("no measured angles in CSV")
        tp, ap = csv["true_peak"], csv["app_peak"]
        if not pd.isna(tp) and not pd.isna(ap) and ap < tp:
            w.append("CSV apparent angle at peak is below the true angle (refraction cannot lower a star)")
        events.append(e)

    # fill missing UTC offsets from the same observer on the same date, else same observer
    for e in events:
        if e.utc_offset_h is None:
            same = [o.utc_offset_h for o in events if o.utc_offset_h is not None and o.observer == e.observer and o.local_date == e.local_date]
            same = same or [o.utc_offset_h for o in events if o.utc_offset_h is not None and o.observer == e.observer]
            if same:
                e.utc_offset_h = max(set(same), key=same.count)
                e.warnings.append(f"UTC offset missing; assumed {e.utc_offset_h:+g} from the same observer")
            else:
                e.warnings.append("UTC offset missing")
        if e._d and e._t and e.utc_offset_h is not None:
            local = datetime.combine(e._d, e._t)
            e.utc = (local - timedelta(hours=e.utc_offset_h)).replace(tzinfo=timezone.utc)

    # cross-check against the GE occlusion check time where it parses
    for e in events:
        ge = _parse_time(e.csv.get("ge_check") or "")
        if ge and e._t:
            diff = abs((datetime.combine(e._d or datetime(2000, 1, 1).date(), ge) - datetime.combine(e._d or datetime(2000, 1, 1).date(), e._t)).total_seconds())
            diff = min(diff, 86400 - diff)
            if diff > 600:
                e.warnings.append(f"GE occlusion check time differs from occlusion time by {diff / 60:.0f} min")

    seen = {}
    for e in events:
        base = slug(f"{e.peak}-{e.star}-{e.utc:%Y-%m-%d}" if e.utc else f"{e.peak}-{e.star}-row{e.row}")
        if e.test_kind:
            base += f"-test-{e.test_kind}"
        n = seen.get(base, 0) + 1
        seen[base] = n
        e.id = base if n == 1 else f"{base}-{n}"
    return events
