# CelTheoReport

This is the pipeline for the CelTheo star-occultation observations. For each recorded
occultation of a star behind a distant peak (or aircraft), it computes:

- **where the star really was**: catalogue position from SIMBAD, with proper motion applied, then
  transformed with astropy to altitude and azimuth at the observer and at the peak;
- **where a globe Earth with a real atmosphere predicts it should be**: a Ciddor (1996) ray trace
  through the nearest radiosonde sounding on the WGS-84 ellipsoid (the `refraction.py` engine from
  [AL_RefractionAnalysis](https://github.com/rooster-ninja/AL_RefractionAnalysis));
- **where a planar model predicts it should be**: rise over run, with and without refraction;

and the residuals between them, in arc-seconds and in seconds of time.

The output is a versioned JSON data contract (`schema/results.schema.json`) built for report
generators. See **[REPORT_BRIEF.md](REPORT_BRIEF.md)**.

## Setup

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e .
```

## Run

```bash
python -m celtheo list                    # parsed events + data warnings
python -m celtheo run                     # all events -> out/results.json, out/results.csv
python -m celtheo run --only pikes-peak   # subset by event-id substring
python -m celtheo run --offline           # cached stars/soundings only
pytest -q                                 # regression tests (offline)
```

The first full run needs network access, for SIMBAD and the radiosonde archives. It caches
everything it fetches:

- **stars** in `data/stars_cache.json`;
- **soundings** in `data/soundings/{wmo}_{YYYY-MM-DDTHH}Z.csv`.

Commit both so later runs are reproducible offline.

## Inputs

| Path | What |
|---|---|
| `data/occultations.csv` | Celestial Theodolite data table (v9.1.1). Rows prefixed ϵ / ♈︎ are flagged as test rows. |
| `config/star_names.json` | CSV star name to SIMBAD identifier, where the name isn't resolvable as written. `null` = no identifier. |
| `config/stations.json` | Radiosonde stations (WMO / IEM ids). The nearest one to the observer is used. |
| `config/sounding_overrides.json` | Per-event station and launch time. Pikes Peak / 39 Aqr is pinned to Grand Junction 2025-01-28 00Z, as in AL_RefractionAnalysis. |
| `data/soundings/` | Cached soundings. Ships with Grand Junction 72476, 2025-01-28 00Z. |

Soundings come from the Iowa Environmental Mesonet RAOB archive, with the University of
Wyoming archive as a fallback. The nearest 00Z/12Z launch is used. If neither archive has
data, the event uses the standard atmosphere and is flagged.

## Validation

`tests/test_engine.py` reproduces the published AL_RefractionAnalysis results for Pikes Peak / 39 Aqr.

| Quantity | Result |
|---|---|
| Peak apparent elevation | +02°19′04.4″ |
| Peak geometric elevation | +02°17′20.6″ |
| Central angle φ | 00°27′16.6″ |
| Ciddor star uplift | 13.33′ |
| Bennett | 13.64′ |

## Changes from the original `refraction.py`

- The atmosphere is a per-event object instead of module globals, so many events can run in one pass.
- The star (infinite-distance) trace is reported with both the normal-section radius used for
  the peak (`star_geometric_deg`) and the original Gauss radius (`star_geometric_gauss_deg`).
  For Pikes Peak they differ by 0.3″.
- Ray-path sampling is added for figures. The numerics are unchanged.

## Data notes surfaced by the pipeline

- **Getaway Mtn, 2025-12-31 (2 rows):** the UTC offset is missing; −7 is assumed from the same observer.
- **HD 9061 UK:** the apparent angle at the peak is below the true angle.
- **Green Mtn:** the GE check time is 60 min off the occlusion time.
- **Spreadsheet `Star GP Long`:** this column is west-positive.
