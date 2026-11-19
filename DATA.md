# Source measurements

The input files contain numerical measurements, sample identifiers and the information needed to pair observations across processing stages. The source publications document sampling and laboratory procedures. Each record represents the level of aggregation available in its source. Data-specific attribution and reuse terms are listed in [data/LICENSE.md](data/LICENSE.md).

## Nicaragua

Trognitz B, Cros E, Assemat S, Davrieux F, Forestier-Chiron N, Ayestas E, Kuant A, Scheldeman X and Hermann M. Diversity of cacao trees in Waslala, Nicaragua: associations between genotype spectra, product quality and yield potential. *PLoS ONE* 8, e54079 (2013). [Article](https://doi.org/10.1371/journal.pone.0054079) · [Table S1 and CC BY 4.0 terms](https://doi.org/10.1371/journal.pone.0054079.s001) · [License](https://creativecommons.org/licenses/by/4.0/)

`waslala_s1_numeric.csv` selects and renames fields from Table S1. Rows follow the source tree IDs, including trees with missing chemistry. A blank denotes a missing value. The prefixes `n_` and `f_` identify dried beans without fermentation and dried beans after six days of fermentation. The suffixes `theo`, `caff` and `tc` identify theobromine, caffeine and the reported T/C, respectively. Concentrations are percentages of dry matter. Of the 106 source tree IDs, 94 from 30 farms have paired stage summaries. The study spanned two harvests and had fewer observations in the second. The source does not establish two observations at each stage for every tree.

`nirs_calibration.csv` gives the component-specific standard errors of cross-validation (SECV, percentage points of concentration) and calibration sample counts from Table 4. The sensitivity calculation applies 0, 0.25, 0.5 and 1 times each SECV to the corresponding concentration, independently at each stage. [Table 4](https://doi.org/10.1371/journal.pone.0054079.t004)

## Colombia

Zapata Bustamante S, Tamayo Tenorio A and Rojano BA. Efecto de la fermentación sobre la actividad antioxidante de diferentes clones de cacao colombiano. *Revista Cubana de Plantas Medicinales* 18, 391–404 (2013). [Article, Table 2](https://www.medigraphic.com/pdfs/revcubplamed/cpm-2013/cpm133g.pdf#page=8)

`component_means.csv` gives theobromine and caffeine in mg/g for five clones from Granja Tierradura, before and after six days of fermentation. Both treatments were sun-dried for five days. `colombia_table2_sd.csv` contains the corresponding SDs. The two files are paired by clone and fermentation day. Table 2 prints both day-0 SDs for ICS 60 as 0.000. Only the rounding allowances expand that sample’s concentration bounds. Means and SDs retain the precision reported in Table 2 (two and three decimal places, respectively). Each clone contributes one pair.

## Ecuador

Espín S, Samaniego I, Wakao H and Jiménez JC. La relación teobromina/cafeína asociada a la calidad del cacao ecuatoriano. *Alimentos Ciencia e Ingeniería* 16(2), 107–109 (2007). [Article](https://repositorio.iniap.gob.ec/items/b086643d-56dd-4a1f-a674-03b274484c21) · [Digitized issue and terms](https://revistas.uta.edu.ec/erevista/index.php/aci/article/view/356) · [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)

`ecuador_table1.csv` transcribes the 88 numerical ratios in Table 1 (p. 109) and adds material type and sampling day. Days are counted from the start of fermentation. Samples were sun-dried before HPLC analysis. Ten materials end on day 5. Their day-6 cells are blank. Four materials have an additional day-6 value. The highest printed ratios for CCAT 4675 and EET 103 are each 6.1. The exact factor that reduces 6.1 to the boundary is 61/60. A slightly larger factor takes both peaks below 6. The original labels “CAT 4675” and “Vol.de Dios” are retained alongside “CCAT 4675” and “Voluntad de Dios”, the names used in the text.

## Venezuela

Álvarez CO, Pérez EE and Lares MC. Beneficio del cacao Criollo venezolano: variaciones en composición proximal, metilxantinas y polifenoles. *RIVAR* 9(27), 113–128 (2022). [Article, Table 2](https://doi.org/10.35588/rivar.v9i27.5625) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

The two Venezuelan rows in `component_means.csv` give analyte means (mg/g dry basis) and the published T/C values. These are the unroasted treatments from a plantation study with 14 sampled trees. The unfermented and three-day fermented beans were sun-dried. Pulp was removed with sawdust in the unfermented treatment. The treatment summaries contribute one paired comparison. Dividing the component means provides a second calculation on the same material.
