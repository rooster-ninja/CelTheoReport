"""Build the report's data file from out/. Presentation only: values are copied, never recomputed.
Terrain arrays are thinned for drawing by keeping real samples (every sample near the observer and
the peak; elsewhere the highest sample in each bin), so skylines are preserved.

    python report/build.py      ->  report/data.js   (index.html loads it)
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
doc = json.loads((ROOT / "out" / "results.json").read_text())
terrain = {}
for ev in doc["events"]:
    f = ROOT / "out" / "terrain" / f"{ev['id']}.json"
    if ev["flags"].get("excluded") or not f.exists() or not (ev.get("terrain") or {}).get("file"):
        continue
    t = json.loads(f.read_text())
    D = t["frame"]["D_m"]
    lines = {}
    for name, p in t["profiles"].items():
        s = np.array(p["s_m"], float)
        u = np.array(p["u_m"], float)
        h = np.array([np.nan if v is None else v for v in p["dtm"]["h_m"]], float)
        keep = set(np.nonzero((s <= 1500) | (np.abs(u) <= 2500))[0].tolist())
        nb = 500
        edges = np.linspace(0, len(s), nb + 1).astype(int)
        for a, b in zip(edges[:-1], edges[1:]):
            if b > a and np.isfinite(h[a:b]).any():
                keep.add(a + int(np.nanargmax(h[a:b])))
        for surf in ("dtm",):
            for model in ("planar", "globe_refracted", "globe_geometric"):
                blk = ((p.get("obstruction") or {}).get(surf) or {}).get(model)
                if blk:
                    keep.add(int(np.argmin(np.abs(s - blk["s_m"]))))
        idx = sorted(i for i in keep if np.isfinite(h[i]))
        g = lambda arr, nd: [round(float(arr[i]), nd) if arr[i] is not None else None for i in idx]
        lines[name] = {"u": g(u, 0), "h": g(h, 1), "gx": g(p["dtm"]["globe_x_m"], 0), "gz": g(p["dtm"]["globe_z_m"], 1),
                       "az": p["azimuth_deg"]}
    terrain[ev["id"]] = lines
out = "window.CELTHEO=" + json.dumps({"doc": doc, "terrain": terrain}, separators=(",", ":"), ensure_ascii=False) + ";"
(ROOT / "report" / "data.js").write_text(out)
print(f"data.js {len(out) / 1e6:.2f} MB, {len(terrain)} terrain sets")
