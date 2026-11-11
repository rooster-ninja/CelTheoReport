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

It also traces the **real terrain** along the line of sight. Elevation data comes from:

- USGS 3DEP at 10 m for the US, with Copernicus GLO-30 for the canopy surface;
- NRCan HRDEM at 2 m for Canada, falling back to MRDEM-30.

It draws two profiles: one to the peak, and one along the star's azimuth at the occlusion time.
Both continue past the peak. With these profiles the pipeline:

- finds what actually hid the star in each model (planar, globe geometric, globe refracted);
- repeats the comparisons against that skyline.

Everything is expressed in a **peak-anchored 2D frame**: the peak sits at the same place in both
models, and only the globe observer drops. Bare earth is the primary surface. The canopy surface
is used only to flag events where trees could matter.

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
python -m celtheo run --only burke        # recompute matching events and merge them into out/
python -m celtheo run --offline           # cached stars/soundings only (no terrain unless DEMs are local)
python -m celtheo run --no-terrain        # skip the elevation-data profiles
python -m celtheo finalize                # re-apply warnings and in_summary to out/ without recomputing
pytest -q                                 # regression tests (offline)
```

The first full run needs network access, for SIMBAD and the radiosonde archives. It caches
everything it fetches:

- **stars** in `data/stars_cache.json`;
- **soundings** in `data/soundings/{wmo}_{YYYY-MM-DDTHH}Z.csv`.

Commit both so later runs are reproducible offline. `out/` is tracked too, so the results (including terrain profiles, which need network access to regenerate) travel with the repo.

## Inputs

| Path | What |
|---|---|
| `data/occultations.csv` | Celestial Theodolite data table (v9.1.1). Rows prefixed ϵ / ♈︎ are flagged as test rows. |
| `config/star_names.json` | CSV star name to SIMBAD identifier, where the name isn't resolvable as written. `null` = no identifier. |
| `config/stations.json` | Radiosonde stations (WMO / IEM ids). The nearest one to the observer is used. |
| `config/sounding_overrides.json` | Per-event station and launch time. Pikes Peak / 39 Aqr is pinned to Grand Junction 2025-01-28 00Z, as in AL_RefractionAnalysis. |
| `config/dem_sources.json` | Elevation sources, tried in order per point. A `local` entry is usable offline. |
| `config/time_corrections.json` | Per-event fixes to UTC offset, local time or date, applied on load. Only add an entry once you've confirmed it. |
| `config/exclusions.json` | Events left out (id → reason). They stay in results.json flagged `excluded`. |
| `data/soundings/` | Cached soundings (ships with Grand Junction 72476, 2025-01-28 00Z) and `_availability.json`, the IEM list of which stations launched at each time. |

Soundings come from the Iowa Environmental Mesonet RAOB archive (JSON API). For each event the
pipeline checks which stations IEM has data for at the nearest 00Z/12Z launch. It then uses the
nearest of those within 600 km, before trying launches 12 h either side. If nothing is found,
the event uses the standard atmosphere and is flagged. For example, Denver has no 2025-03-24 00Z
launch, so the North Peak events use Grand Junction.

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

- `profiles.ray_to_peak` is in the peak-anchored frame, as of schema 1.1.
- The atmosphere is a per-event object instead of module globals, so many events can run in one pass.
- The star (infinite-distance) trace is reported with both the normal-section radius used for
  the peak (`star_geometric_deg`) and the original Gauss radius (`star_geometric_gauss_deg`).
  For Pikes Peak they differ by 0.3″.
- Ray-path sampling is added for figures. The numerics are unchanged.

## Data checks

- **Clock check:** when the CSV's true altitude disagrees with astropy by more than 60″, the
  pipeline searches ±6 h for the clock shift that reproduces it, with the star on the peak's
  azimuth. A clean whole-hour shift is reported as a UTC-offset problem; anything else points
  to the star identity, date or time.
- **Foreground:** terrain closer than 5 % of the distance (100 m to 1 km) is the observer's
  foreground. It's reported separately and excluded from the obstruction search.
- **Aircraft:** these events skip terrain, since the occluder is the aircraft.

## Data notes surfaced by the pipeline

- **Getaway Mtn, 2025-12-31 (2 rows):** the UTC offset is missing; −7 is assumed from the same observer.
- **HD 9061 UK:** the apparent angle at the peak is below the true angle.
- **Green Mtn:** the GE check time is 60 min off the occlusion time.
- **Spreadsheet `Star GP Long`:** this column is west-positive.
