# NMR METABOLIC PROFILING OF MOSQUITO SPECIES TO UNDERSTAND INSECTICIDE RESISTANCE

This repository contains the custom scripts used in the analysis of the data.

## Getting Started

Both scripts were written with base packages and should be highly portable. There are two scripts provided:

* **PathTabForFisher.py**: Generates a table of metabolites per pathway for a given organism. Uses KEGG as database source.
* **Pathway_Fisher.R**: Performs qualitative pathway analysis on a set of metabolites.



### Prerequisites
* Python version 2.7

* R version 3.5.2


### Installing and usage

* **PathTabForFisher.py**: No installation required.

Creates a table of metabolites per pathway (in KEGG codes) for Anopheles gambiae (aga). Pathways are in rows and metabolites are in columns.
```
python PathTabForFisher.py aga

```
* **Pathway_Fisher.R**: Source the script in R.
PathTable: takes a data frame with pathways in columns and metabolites are in rows.
Metabolites: a vector of KEGG codes
```R
source(Pathway_Fisher.R)
pathwayFisher(PathTable, Metabolites)
```




## Authors

* **Rudi Grosman**

## License

This work is licensed under a Creative Commons Attribution-NonCommercial-Share Alike 4.0 International [CC BY-NC-SA 4.0] (https://creativecommons.org/licenses/by-nc-sa/4.0/).

See the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Dr Eva Caamano Gutierrez
* Dr Arturas Grauslys
* [KEGG](https://doi.org/10.1093/nar/gkw1092)


