PGxPOP Report
================

-   [Returned counts](#returned-counts)
-   [Calculated counts](#calculated-counts)
-   [Concordance](#concordance)
-   [Figures](#figures)

<br>

Following the pipeline described in [McInnes et al (2020).
Pharmacogenetics at Scale: An Analysis of the UK Biobank. *Clinical
Pharmacology & Therapeutics*](https://doi.org/10.1002/cpt.2122)

-   pharmacogenes are assigned a **star allele**
-   **clinical function** is assigned to the star allele

Note: Imputed data (aligned to Genome Reference Consortium Human Build
37, GRCh37 (synonym hg19)) lifted over to GRCh38 (synonym hg38) to match
exome sequencing data.

## Returned counts

McInnes et al. [returned dataset
3388](https://biobank.ndph.ox.ac.uk/ukb/dset.cgi?id=3388)

Includes the directory FFR\_33722\_1

<pre>
FFR_33722_1/
├── cpt.2122.pdf
├── pgx_calls.imputed_callset_DAedit.csv
├── pgx_calls.integrated_callset_DAedit.csv
└── README.md
</pre>

**NB**. Phenotype counts reported are for “presumptive\_phenotype”

<br>

**Returned phenotype count CYP2C19**

    ## # A tibble: 9 x 3
    ##   phenotype_presumptive           n_imputed n_integrated
    ##   <chr>                               <int>        <int>
    ## 1 Indeterminate                         114           18
    ## 2 Intermediate Metabolizer           129285        13337
    ## 3 Likely Intermediate Metabolizer       217           35
    ## 4 Likely Poor Metabolizer                46           10
    ## 5 Normal Metabolizer                 191655        19505
    ## 6 Not available                         272           76
    ## 7 Poor Metabolizer                    12895         1376
    ## 8 Rapid Metabolizer                  129505        13218
    ## 9 Ultrarapid Metabolizer              22343         2210

<br>

**Returned phenotype count CYP2D6**

    ## # A tibble: 5 x 3
    ##   phenotype_presumptive    n_imputed n_integrated
    ##   <chr>                        <int>        <int>
    ## 1 Indeterminate                83318          848
    ## 2 Intermediate Metabolizer    119303        17846
    ## 3 Normal Metabolizer          181867        27236
    ## 4 Not available                77964          857
    ## 5 Poor Metabolizer             24038         2998

<br>

## Calculated counts

CYP2C19: chromosome 10 \* *GRCh38 coordinates* 10:94,762,680-94,855,546
\* UKB exome block 24, 93795159 - 95869829

CYP2D6: chromsome 22 \* *GRCh38 coordinates* 22:42,126,498-42,130,809 \*
UKB exome block 16, 41603021 - 42642420

<br>

**Phenotype count CYP2C19**

    ## # A tibble: 7 x 3
    ##   phenotype_presumptive           n_imputed n_exome
    ##   <chr>                               <int>   <int>
    ## 1 Normal Metabolizer                 132137  129989
    ## 2 Intermediate Metabolizer            47627   43608
    ## 3 Not available                         550    6263
    ## 4 Poor Metabolizer                     4784    5001
    ## 5 Likely Intermediate Metabolizer        NA     130
    ## 6 Indeterminate                          NA      86
    ## 7 Likely Poor Metabolizer                NA      19

<br>

**Phenotype count CYP2D6**

    ## # A tibble: 5 x 3
    ##   phenotype_presumptive    n_imputed n_exome
    ##   <chr>                        <int>   <int>
    ## 1 Normal Metabolizer          134538   89461
    ## 2 Not available                17483   46571
    ## 3 Intermediate Metabolizer      1894   42156
    ## 4 Indeterminate                31176    6741
    ## 5 Poor Metabolizer                 7     167
