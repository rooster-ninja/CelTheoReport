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

The data contract is `schema/results.schema.json` (`schema_version` 1.0). Read the
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
4. **One plate per event**, linkable as `#{id}`:
   - the side-view geometry: observer, peak, Earth's curve and the refracted ray
     (`profiles.ray_to_peak`), against the straight planar line;
   - the atmosphere column: `profiles.atmosphere` (n−1, T, P, RH against altitude);
   - an angles table: true, apparent and model angles in degrees, minutes and seconds;
   - notes: `flags.warnings` verbatim, in the style of a footnote or erratum.
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
- Exclude test rows from summary statistics, but show them in the appendix.
- Always show each event's sounding provenance: station, launch time, distance, and
  “standard-atmosphere fallback” where it applies.
- Put the source CSV row (`row`) in each plate's footnote, so readers can trace it to the
  spreadsheet.
