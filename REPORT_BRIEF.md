# Report brief: for the agent building the CelTheo report

You are building a report on the CelTheo star-occultation observations, using the data this
repository computes. The computation is finished and tested. **Your job is presentation:
do not recompute or "correct" any value.** Every number on the page comes from `out/results.json`.

## 1 · Getting the data

```bash
python -m celtheo run            # writes out/results.json + out/results.csv (needs network on first run)
python -m celtheo run --offline  # reuse cached stars and soundings
```

Or from Python:

```python
from celtheo.api import run, load_results
doc = load_results("out/results.json")       # dict, see schema/results.schema.json
events = [e for e in doc["events"] if not e["flags"]["is_test"]]
```

The data contract is `schema/results.schema.json` (`schema_version` 1.1). Each event's terrain
profile is in `out/terrain/{id}.json` (`schema/terrain.schema.json`). Read the
`description` fields. They say what each number means and its sign convention.

## 2 · What the report is about

Observers timed the moment a star vanished behind a distant mountain peak (and, in four
cases, behind an aircraft). For each event the pipeline answers one question:

> At the recorded time, where was the star really, and where should it have been
> for the peak to just cover it?

The pipeline produces three predictions of where the star should have been. Each is compared with the star's catalogue position (`astropy.true_alt_obs_deg`):

| Model | Prediction | Residual |
|---|---|---|
| Globe + Ciddor ray trace through the real radiosonde atmosphere | `refraction.star_geometric_deg` | `comparison.globe_ciddor_arcsec` / `_seconds` |
| Planar, no refraction (rise over run) | `geometry.planar_angle_deg` | `comparison.planar_geometric_arcsec` / `_seconds` |
| Planar with astropy's refraction | `geometry.planar_angle_deg` vs `astropy.erfa_apparent_alt_obs_deg` | `comparison.planar_refracted_arcsec` / `_seconds` |

The strongest test uses **the real skyline**. `terrain.comparison` repeats the globe and planar
comparisons against the actual terrain along the star's azimuth at the moment of occlusion,
taken from bare-earth elevation data (3DEP 10 m or HRDEM 2 m). Prefer these unsuffixed
`terrain.comparison.*` values over the summit-only ones whenever they exist. `_dsm` variants
use the canopy surface, and `terrain.canopy.*.flag` marks the events where trees could change
the answer.

**Refraction of the star.** `refraction.star_apparent_ciddor_obs_deg` is the star's apparent
altitude by the AstronomyLive Ciddor ray trace, worked from its astropy true altitude.
`star_refraction_ciddor_obs_arcmin` is the lift between them; draw it as the star's apparent place
and the wedge from true to apparent. Bennett (`star_refraction_bennett_obs_arcmin`) and astropy's
own value (`star_refraction_erfa_obs_arcmin`) are side values only. The astropy one runs low near
the horizon.

Residuals are **observed − model**. The `_seconds` form converts each residual into the clock
error that would close it, which is the most intuitive unit for readers. The spreadsheet's own
values (`csv.*`) and Bennett's formula (`refraction.bennett_arcmin`) are there for cross-checks.

Present the models side by side and let the numbers speak. Don't editorialise beyond what the
residuals show.

## 3 · Suggested structure

1. **Frontispiece and title page.** For example: *“A Memoir upon the Occultation of Fixed Stars
   by Distant Summits, with Determinations of Terrestrial Refraction.”* Include the date of
   generation, from `generated_utc`.
2. **Summary plate.** All real events (exclude `is_test`), with residuals per model per event
   as a dot or whisker chart in arc-seconds, plus the seconds-of-time equivalent.
3. **Catalogue of observations.** A ruled table with event, star, observer, UTC, distance,
   height difference, the three residuals, sounding station, and a warnings marker.
4. **One plate per event**, linkable as `#{id}`. Everything is 2D.
   - **Profile figure (the centrepiece).** Use the peak-anchored frame, as in the Desmos graph:
     the peak sits at the origin and at the same height in both models, and only the globe
     observer is lower.
     - **Planar:** the terrain is `profiles.star_azimuth.u_m` against `dtm.h_m`, with the
       observer at `geometry.observer_planar_xz_m`.
     - **Globe:** the terrain is `dtm.globe_x_m` against `dtm.globe_z_m`, with the observer at
       `geometry.observer_globe_xz_m`.
     - **Drop:** label the globe observer's drop with `geometry.globe_observer_drop_m`.
     - **Rays:** draw the refracted ray `profiles.ray_to_peak` (results.json) and the straight
       planar sight line.
     - **Occluder:** mark it with `terrain.star_azimuth.obstruction.dtm.globe_refracted` and
       `...planar`. If `near_peak` is false, say so plainly: the named peak was not what hid
       the star.
     - **Scale:** use a vertical exaggeration (e.g. ×10) and state it on the figure, so the
       curvature drop is visible.
   - **Atmosphere column:** `profiles.atmosphere`, showing n−1, T, P and RH against altitude.
   - **Angles table:** true, apparent and model angles in degrees, minutes and seconds.
   - **Notes:** `flags.warnings` verbatim, in the style of a footnote or erratum, including any
     canopy note.
5. **Aircraft occultations.** Treat the four `flags.target_type == "aircraft"` events separately.
6. **Appendix on method.** Cover Bouguer's law (b = n·r·cos α), the Ciddor (1996) refractive
   index, the WGS-84 normal-section radius, the radiosonde source per event
   (`refraction.sounding`), and the test rows (`test_kind` epsilon/aries) with their purpose.
7. **Errata.** List every warning across all events.

## 4 · Style: a nineteenth-century scientific memoir (c. 1800–1880)

Think *Philosophical Transactions of the Royal Society*, the *Nautical Almanac*, and Airy's
or Bessel's observatory reports.

- **Paper and ink.** Warm ivory or laid-paper ground, iron-gall brown-black ink, and a single
  restrained accent (oxblood or Prussian blue) for model lines. No gradients, no drop shadows.
- **Type.** Use Google Fonts, which the artifact host allows.
  - *IM Fell English* or *IM Fell DW Pica* for display and running heads.
  - *EB Garamond* or *Libre Caslon Text* for body text.
  - Old-style figures, small caps for names, roman numerals for plates
    (“Plate IV.”), section marks (§), and hanging punctuation.
- **Figures as engraved plates.** Draw hairline strokes in SVG, with hatching instead of fills
  (diagonal hatch for the Earth's body, stipple for the atmosphere). Use the copperplate
  lettering convention “Fig. 1.”, with captions set in italic below. Add a plate border with
  a double rule and plate numbers in the corners.
- **Tables.** Use full-width hairline rules above and below, and none between rows.
  Align figures on the decimal point. Use the ° ′ ″ glyphs.
- **Ornament, used sparingly.** A fleuron (❦) or rule between sections, a drop cap at the
  opening, and a colophon at the end crediting the observers by name (`observer.name`).
- **Interactivity, kept invisible until used.** Hovering a residual dot shows the event.
  Clicking a catalogue row scrolls to its plate. Stay in period: no modern UI chrome.

## 5 · Hard rules

- Use only numbers that appear in `results.json`. Format them, but never re-derive them.
  If a value is `null`, print an em dash and the relevant warning.
- Keep observed − model signs as given, and label units on every axis and column.
- Leave `flags.excluded` events out entirely, except to list them with `excluded_reason` in the errata.
- Summary statistics and the summary plate use **only** events with `flags.in_summary == true`.
  The pipeline decides this, so don't re-derive it. Every other event still gets its plate,
  with its warnings shown., but show them in the appendix.
- Always show elevation-data provenance (`terrain.*.sources`) and the DEM vs CSV height check
  (`terrain.checks`).
- Always show each event's sounding provenance: station, launch time, distance, and
  “standard-atmosphere fallback” where it applies.
- Put the source CSV row (`row`) in each plate's footnote, so readers can trace it to the
  spreadsheet.
