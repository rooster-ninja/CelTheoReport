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


SIMBAD_SERVERS = ["simbad.cds.unistra.fr", "simbad.harvard.edu"]


def _simbad_query(q: str, server: str) -> Optional[dict]:
    from astroquery.simbad import Simbad, conf
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    conf.server = server
    s = Simbad()
    if hasattr(s, "server"):
        s.server = server
    s.add_votable_fields("pmra", "pmdec")
    t = s.query_object(q)
    if t is None or not len(t):
        return None
    cols = {c.lower(): c for c in t.colnames}
    if "ra" in cols and str(t[cols["ra"]].unit) == "deg":
        ra, dec = float(t[cols["ra"]][0]), float(t[cols["dec"]][0])
    else:  # older astroquery: sexagesimal strings
        c = SkyCoord(str(t[cols["ra"]][0]), str(t[cols["dec"]][0]), unit=(u.hourangle, u.deg))
        ra, dec = c.ra.deg, c.dec.deg
    pm = []
    for k in ("pmra", "pmdec"):
        v = t[cols[k]][0] if k in cols else 0.0
        pm.append(0.0 if str(v) in ("--", "nan") else float(v))
    return {"query": q, "kind": "star", "ra_deg": ra, "dec_deg": dec, "pmra_masyr": pm[0], "pmdec_masyr": pm[1],
            "source": f"SIMBAD {server} (astroquery)"}


def resolve(csv_name: str, offline=False, log=print, refresh_no_pm=False) -> dict:
    """Return {query, kind, ra_deg, dec_deg, pmra_masyr, pmdec_masyr, source} or {kind:'unresolved'}."""
    q = query_name(csv_name)
    if q is None:
        return {"query": None, "kind": "unresolved", "source": "no identifier"}
    if q.lower() in BODIES:
        return {"query": q, "kind": "body", "body": BODIES[q.lower()], "source": "astropy ephemeris"}
    cache = _load()
    if q in cache and not (refresh_no_pm and "no proper motion" in cache[q].get("source", "")):
        return cache[q]
    if offline:
        return {"query": q, "kind": "unresolved", "source": "offline and not in cache"}
    rec = None
    for server in SIMBAD_SERVERS:  # SIMBAD via astroquery: gives proper motion; CDS first, then the Harvard mirror
        try:
            rec = _simbad_query(q, server)
            if rec:
                break
        except Exception as e:  # noqa: BLE001
            log(f"    SIMBAD ({server}) failed for {q!r}: {str(e).splitlines()[0]}")
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
