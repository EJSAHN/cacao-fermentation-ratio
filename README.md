# Cocoa bean fermentation and methylxanthine ratio bands

This analysis asks whether cocoa beans from the same tree or material stay in the same theobromine-to-caffeine (T/C) band across fermentation stages. It applies two published sets of thresholds to measurements from Nicaragua, Colombia, Ecuador and Venezuela. Each study retains its own sampling unit.

## Run the analysis

From the repository directory, run the following command with Python 3.8 or later. The calculations use the Python standard library and the CSV inputs in `data/`.

```sh
python run_analysis.py --output results
```

The command creates a new directory with 25 CSV tables. To repeat the analysis, choose a new output directory. For measurement units and source tables, see [DATA.md](DATA.md). Rights for the software and published measurements are described separately in [LICENSE](LICENSE) and [data/LICENSE.md](data/LICENSE.md).

The boundary and sampling-day checks can be run with `python -m unittest discover -s tests`.

## Endpoint comparisons

Each comparison uses two measurements of the same source unit. A change means that the two T/C values fall in different intervals. The counts below use the reported ratios in Nicaragua, Ecuador and Venezuela and the ratios of published analyte means in Colombia.

| Study | Unit | Four-band rule | Three-band rule |
| --- | --- | ---: | ---: |
| Colombia | Clone, day 0 vs day 6 | 2/5 | 0/5 |
| Ecuador | Clone or farm material, day 0 vs last observed day | 1/14 | 0/14 |
| Nicaragua | Tree, unfermented vs fermented | 21/94 | 12/94 |
| Venezuela | Plantation treatment mean, unfermented vs fermented | 0/1 | 0/1 |

Using the four-band thresholds described by Osorio-Guarín et al., this implementation assigns T/C ≤ 2, 2 < T/C < 6, 6 ≤ T/C < 8 and T/C ≥ 8 to successive bands. The three-band rule assigns T/C < 3, 3 ≤ T/C ≤ 9 and T/C > 9 to successive bands. A ratio printed as 6.0 belongs to the third band. The three-band thresholds were described for fresh beans. This analysis applies them to dried beans at each stage. [Osorio-Guarín et al. (2017)](https://doi.org/10.3389/fpls.2017.01994) · [Peláez et al. (2016)](https://doi.org/10.17268/sci.agropecu.2016.04.01)

## Reading the results

| File | Measurement |
| --- | --- |
| `endpoint_summary.csv`, `unit_margins.csv` | Counts by study and ratio definition, with the starting and ending T/C, bands and crossings for each paired unit. |
| `external_stages.csv`, `ecuador_paths.csv` | Ratios and bands at each observed stage, including paths that leave the initial band and return. |
| `colombia_sd_detail.csv`, `colombia_sd_stress.csv` | Feasible T/C intervals and numbers of changing clones when the published component means vary within rounding bounds or within one SD plus rounding. |
| `ecuador_schedules.csv`, `ecuador_schedule_materials.csv` | Changes detected and missed for each of 16 subsets of sampling days, for all materials and for each material separately. |
| `farm_cluster_estimates.csv`, `paired_rule_contrast.csv` | Tree-level change fractions and the paired difference between the rules, with percentile intervals from 10,000 whole-farm bootstrap draws. |
| `ratio_definition_agreement.csv` | Agreement between classifications based on the published T/C columns and on the ratio of published concentration means in Nicaragua. |
| `ratio_stress.csv`, `nirs_stress.csv` | Bounds on the number of changing units when the two ratios or their component concentrations can vary within the specified limits. |

An endpoint change is forced when the feasible bands at the two stages are disjoint. It is possible when some pair of admissible values falls in different bands. These bounds describe the specified perturbations. Only the outputs from farm resampling give confidence intervals.

For Ecuador, every material has observations on days 0–5. Each schedule retains days 0 and 5 and adds a subset of days 1–4. A change counts as detected when a retained day's band differs from the day-0 band. A missed change appears on at least one observed day but on none of the retained days. The four series with day-6 observations also appear in `ecuador_day6_extension.csv`.

For both Ecuadorian series that reach the boundary, the largest printed value is 6.1. Dividing 6.1 by 61/60 gives exactly 6, which still lies in the upper band. Any larger factor permits a value below the cutoff. The stored 1.02 scenario illustrates an allowance just beyond that threshold.
