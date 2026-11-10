"""Resolve CSV star names to catalogue positions (SIMBAD), with a committed cache."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "stars_cache.json"
NAMES = json.loads((ROOT / "config" / "star_names.json").read_text())

BODIES = {"moon": "moon"}


def query_name(csv_name: str) -> Optional[str]:
    """SIMBAD identifier for a CSV star name, or None if it cannot be resolved by name."""
    n = csv_name.strip()
    if n in NAMES:
        return NAMES[n]
    n = re.sub(r"\s*\((\d+|Moon)\)\s*$", "", n)          # "39 Aquarii (3)", "39 Aquarii (Moon)"
    m = re.match(r"^(HD|HIP|SAO|SOA)\s*(\d+)", n, re.I)  # "HD55892 (I Puppis)" -> "HD 55892"
    if m:
        cat = "SAO" if m.group(1).upper() == "SOA" else m.group(1).upper()
        return f"{cat} {m.group(2)}"
    return re.sub(r"\s*\(.*\)\s*$", "", n) or None


def _load():
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def resolve(csv_name: str, offline=False, log=print) -> dict:
    """Return {query, kind, ra_deg, dec_deg, pmra_masyr, pmdec_masyr, source} or {kind:'unresolved'}."""
    q = query_name(csv_name)
    if q is None:
        return {"query": None, "kind": "unresolved", "source": "no identifier"}
    if q.lower() in BODIES:
        return {"query": q, "kind": "body", "body": BODIES[q.lower()], "source": "astropy ephemeris"}
    cache = _load()
    if q in cache:
        return cache[q]
    if offline:
        return {"query": q, "kind": "unresolved", "source": "offline and not in cache"}
    rec = None
    try:  # SIMBAD via astroquery: gives proper motion
        from astroquery.simbad import Simbad
        s = Simbad()
        s.add_votable_fields("pmra", "pmdec")
        t = s.query_object(q)
        if t is not None and len(t):
            cols = {c.lower(): c for c in t.colnames}
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            if "ra" in cols and t[cols["ra"]].unit is not None and str(t[cols["ra"]].unit) == "deg":
                ra, dec = float(t[cols["ra"]][0]), float(t[cols["dec"]][0])
            else:  # older astroquery: sexagesimal strings
                c = SkyCoord(str(t[cols.get("ra", "RA")][0]), str(t[cols.get("dec", "DEC")][0]), unit=(u.hourangle, u.deg))
                ra, dec = c.ra.deg, c.dec.deg
            pmra = t[cols.get("pmra", "PMRA")][0] if cols.get("pmra") else 0.0
            pmdec = t[cols.get("pmdec", "PMDEC")][0] if cols.get("pmdec") else 0.0
            pmra = 0.0 if str(pmra) in ("--", "nan") else float(pmra)
            pmdec = 0.0 if str(pmdec) in ("--", "nan") else float(pmdec)
            rec = {"query": q, "kind": "star", "ra_deg": ra, "dec_deg": dec, "pmra_masyr": pmra, "pmdec_masyr": pmdec,
                   "source": "SIMBAD (astroquery)"}
    except Exception as e:  # noqa: BLE001
        log(f"    astroquery SIMBAD failed for {q!r}: {e}")
    if rec is None:
        try:  # Sesame name resolver: position only
            from astropy.coordinates import SkyCoord
            c = SkyCoord.from_name(q)
            rec = {"query": q, "kind": "star", "ra_deg": c.ra.deg, "dec_deg": c.dec.deg, "pmra_masyr": 0.0,
                   "pmdec_masyr": 0.0, "source": "Sesame (no proper motion)"}
        except Exception as e:  # noqa: BLE001
            log(f"    could not resolve {q!r}: {e}")
            return {"query": q, "kind": "unresolved", "source": f"lookup failed: {e}"}
    cache[q] = rec
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    log(f"    resolved {csv_name!r} -> {q} ({rec['source']})")
    return rec
