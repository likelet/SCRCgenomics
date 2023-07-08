# SCRCgenomics
Sourcecode for SCRC project.

# R session info 

``` r
sessioninfo::session_info()
─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────
setting  value
version  R version 4.2.1 (2022-06-23)
 os       Ubuntu 20.04.6 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Etc/UTC
 date     2023-07-08
 rstudio  2022.07.2+576 Spotted Wakerobin (server)
 pandoc   2.19.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package                     * version   date (UTC) lib source
 abind                         1.4-5     2016-07-21 [1] RSPM (R 4.2.0)
 assertthat                    0.2.1     2019-03-21 [1] RSPM (R 4.2.0)
 backports                     1.4.1     2021-12-13 [1] RSPM (R 4.2.0)
 Biobase                     * 2.58.0    2022-11-01 [1] Bioconductor
 BiocGenerics                * 0.44.0    2022-11-01 [1] Bioconductor
 BiocIO                        1.8.0     2022-11-01 [1] Bioconductor
 BiocParallel                  1.32.6    2023-03-17 [1] Bioconductor
 Biostrings                  * 2.66.0    2022-11-01 [1] Bioconductor
 bitops                        1.0-7     2021-04-24 [1] RSPM (R 4.2.0)
 boot                          1.3-28    2021-05-03 [4] CRAN (R 4.2.1)
 broom                         1.0.1     2022-08-29 [1] RSPM (R 4.2.0)
 broom.helpers                 1.9.0     2022-09-23 [1] RSPM (R 4.2.0)
 BSgenome                    * 1.66.3    2023-02-16 [1] Bioconductor
 BSgenome.Hsapiens.UCSC.hg19 * 1.4.3     2023-05-09 [1] Bioconductor
 car                           3.1-1     2022-10-19 [1] RSPM (R 4.2.0)
 carData                       3.0-5     2022-01-06 [1] RSPM (R 4.2.0)
 circlize                      0.4.15    2022-05-10 [1] RSPM (R 4.2.0)
 cli                           3.4.1     2022-09-23 [1] RSPM (R 4.2.0)
 clue                          0.3-62    2022-10-18 [1] RSPM (R 4.2.0)
 cluster                     * 2.1.3     2022-03-28 [4] CRAN (R 4.2.1)
 codetools                     0.2-18    2020-11-04 [4] CRAN (R 4.2.1)
 colorspace                    2.0-3     2022-02-21 [1] RSPM (R 4.2.0)
 ComplexHeatmap              * 2.14.0    2022-11-01 [1] Bioconductor
 cowplot                     * 1.1.1     2020-12-30 [1] RSPM (R 4.2.0)
 crayon                        1.5.2     2022-09-29 [1] RSPM (R 4.2.0)
 data.table                  * 1.14.4    2022-10-17 [1] RSPM (R 4.2.0)
 DBI                           1.1.3     2022-06-18 [1] RSPM (R 4.2.0)
 DelayedArray                  0.24.0    2022-11-01 [1] Bioconductor
 digest                        0.6.30    2022-10-18 [1] RSPM (R 4.2.0)
 doParallel                  * 1.0.17    2022-02-07 [1] RSPM (R 4.2.0)
 dplyr                       * 1.0.10    2022-09-01 [1] RSPM (R 4.2.0)
 ellipsis                      0.3.2     2021-04-29 [1] RSPM (R 4.2.0)
 evaluate                      0.17      2022-10-07 [1] RSPM (R 4.2.0)
 fansi                         1.0.3     2022-03-24 [1] RSPM (R 4.2.0)
 fastmap                       1.1.0     2021-01-25 [1] RSPM (R 4.2.0)
 finalfit                    * 1.0.5     2022-08-09 [1] RSPM (R 4.2.0)
 forcats                       0.5.2     2022-08-19 [1] RSPM (R 4.2.0)
 foreach                     * 1.5.2     2022-02-02 [1] RSPM (R 4.2.0)
 formatR                       1.12      2022-03-31 [1] RSPM (R 4.2.0)
 Formula                       1.2-4     2020-10-16 [1] RSPM (R 4.2.0)
 futile.logger               * 1.4.3     2016-07-10 [1] RSPM (R 4.2.0)
 futile.options                1.0.1     2018-04-20 [1] RSPM (R 4.2.0)
 generics                      0.1.3     2022-07-05 [1] RSPM (R 4.2.0)
 GenomeInfoDb                * 1.34.9    2023-02-02 [1] Bioconductor
 GenomeInfoDbData              1.2.9     2023-05-16 [1] Bioconductor
 GenomicAlignments             1.34.1    2023-03-09 [1] Bioconductor
 GenomicRanges               * 1.50.2    2022-12-16 [1] Bioconductor
 GetoptLong                    1.0.5     2020-12-15 [1] RSPM (R 4.0.3)
 gg.gap                      * 1.3       2019-09-30 [1] RSPM (R 4.0.0)
 ggalluvial                    0.12.3    2020-12-05 [1] RSPM (R 4.0.3)
 ggplot2                     * 3.3.6     2022-05-03 [1] RSPM (R 4.2.0)
 ggpubr                      * 0.4.0     2020-06-27 [1] RSPM (R 4.2.0)
 ggrepel                     * 0.9.1     2021-01-15 [1] RSPM (R 4.2.0)
 ggsci                       * 2.9       2018-05-14 [1] RSPM (R 4.2.0)
 ggsignif                      0.6.4     2022-10-13 [1] RSPM (R 4.2.0)
 ggthemr                     * 1.1.0     2021-12-14 [2] Github (cttobin/ggthemr@4a31e0d)
 GlobalOptions                 0.1.2     2020-06-10 [1] RSPM (R 4.0.3)
 glue                          1.6.2     2022-02-24 [1] RSPM (R 4.2.0)
 gridBase                      0.4-7     2014-02-24 [1] RSPM (R 4.0.3)
 gridExtra                   * 2.3       2017-09-09 [1] RSPM (R 4.2.0)
 gt                            0.7.0     2022-08-25 [1] RSPM (R 4.2.0)
 gtable                        0.3.1     2022-09-01 [1] RSPM (R 4.2.0)
 gtsummary                   * 1.6.2     2022-09-30 [1] RSPM (R 4.2.0)
 hms                           1.1.2     2022-08-19 [1] RSPM (R 4.2.0)
 htmltools                   * 0.5.3     2022-07-18 [1] RSPM (R 4.2.0)
 httr                          1.4.4     2022-08-17 [1] RSPM (R 4.2.0)
 IRanges                     * 2.32.0    2022-11-01 [1] Bioconductor
 iterators                   * 1.0.14    2022-02-05 [1] RSPM (R 4.2.0)
 kableExtra                    1.3.4     2021-02-20 [1] RSPM (R 4.2.0)
 km.ci                         0.5-2     2009-08-30 [2] RSPM (R 4.0.3)
 KMsurv                        0.1-5     2012-12-03 [2] RSPM (R 4.0.3)
 knitr                         1.40      2022-08-24 [1] RSPM (R 4.2.0)
 lambda.r                      1.2.4     2019-09-18 [1] RSPM (R 4.2.0)
 lattice                     * 0.20-45   2021-09-22 [4] CRAN (R 4.2.1)
 lifecycle                     1.0.3     2022-10-07 [1] RSPM (R 4.2.0)
 lpSolve                       5.6.17    2022-10-10 [1] RSPM (R 4.2.0)
 maftools                    * 2.6.05    2021-02-04 [2] Bioconductor
 magrittr                    * 2.0.3     2022-03-30 [1] RSPM (R 4.2.0)
 MASS                          7.3-57    2022-04-22 [4] CRAN (R 4.2.1)
 Matrix                        1.5-1     2022-09-13 [1] RSPM (R 4.2.0)
 MatrixGenerics                1.10.0    2022-11-01 [1] Bioconductor
 matrixStats                   0.62.0    2022-04-19 [1] RSPM (R 4.2.0)
 mice                          3.14.0    2021-11-24 [1] RSPM (R 4.2.0)
 munsell                       0.5.0     2018-06-12 [1] RSPM (R 4.2.0)
 MutationalPatterns          * 3.8.1     2023-01-26 [3] Bioconductor
 neutralitytestr             * 0.0.3     2021-02-16 [1] RSPM (R 4.0.3)
 NMF                         * 0.24.0    2022-03-29 [1] RSPM (R 4.2.1)
 pheatmap                    * 1.0.12    2019-01-04 [1] RSPM (R 4.2.0)
 pillar                        1.8.1     2022-08-19 [1] RSPM (R 4.2.0)
 pkgconfig                     2.0.3     2019-09-22 [1] RSPM (R 4.2.0)
 pkgmaker                    * 0.32.2    2020-10-20 [1] RSPM (R 4.0.3)
 plyr                        * 1.8.7     2022-03-24 [1] RSPM (R 4.2.0)
 png                           0.1-7     2013-12-03 [1] RSPM (R 4.2.0)
 pracma                        2.4.2     2022-09-22 [1] RSPM (R 4.2.0)
 purrr                         0.3.5     2022-10-06 [1] RSPM (R 4.2.0)
 R6                            2.5.1     2021-08-19 [1] RSPM (R 4.2.0)
 RColorBrewer                * 1.1-3     2022-04-03 [1] RSPM (R 4.2.0)
 Rcpp                          1.0.9     2022-07-08 [1] RSPM (R 4.2.0)
 RCurl                         1.98-1.9  2022-10-03 [1] RSPM (R 4.2.0)
 readr                         2.1.3     2022-10-01 [1] RSPM (R 4.2.0)
 registry                    * 0.5-1     2019-03-05 [1] RSPM (R 4.0.0)
 reshape2                    * 1.4.4     2020-04-09 [1] RSPM (R 4.2.0)
 restfulr                      0.0.15    2022-06-16 [1] RSPM (R 4.2.1)
 rjson                         0.2.21    2022-01-09 [1] RSPM (R 4.2.0)
 rlang                         1.0.6     2022-09-24 [1] RSPM (R 4.2.0)
 rmarkdown                     2.17      2022-10-07 [1] RSPM (R 4.2.0)
 Rmisc                       * 1.5.1     2022-05-02 [1] RSPM (R 4.2.0)
 rngtools                    * 1.5.2     2021-09-20 [1] RSPM (R 4.2.0)
 Rsamtools                     2.14.0    2022-11-01 [1] Bioconductor
 rstatix                       0.7.0     2021-02-13 [1] RSPM (R 4.2.0)
 rstudioapi                    0.14      2022-08-22 [1] RSPM (R 4.2.0)
 rtracklayer                 * 1.58.0    2022-11-01 [1] Bioconductor
 rvest                         1.0.3     2022-08-19 [1] RSPM (R 4.2.0)
 S4Vectors                   * 0.36.2    2023-02-26 [1] Bioconductor
 sampling                    * 2.9       2021-01-13 [1] RSPM (R 4.0.3)
 scales                      * 1.2.1     2022-08-20 [1] RSPM (R 4.2.0)
 sessioninfo                   1.2.2     2021-12-06 [1] RSPM (R 4.2.0)
 shape                         1.4.6     2021-05-19 [1] RSPM (R 4.2.0)
 stringi                       1.7.8     2022-07-11 [1] RSPM (R 4.2.0)
 stringr                       1.4.1     2022-08-20 [1] RSPM (R 4.2.0)
 SummarizedExperiment          1.28.0    2022-11-01 [1] Bioconductor
 survival                    * 3.3-1     2022-03-03 [4] CRAN (R 4.2.1)
 survminer                   * 0.4.9     2021-03-09 [2] RSPM (R 4.0.3)
 survMisc                      0.5.5     2018-07-05 [2] RSPM (R 4.0.3)
 svglite                       2.0.0     2021-02-20 [2] RSPM (R 4.0.4)
 systemfonts                   1.0.4     2022-02-11 [1] RSPM (R 4.2.0)
 table1                      * 1.4.2     2021-06-06 [1] RSPM (R 4.2.0)
 tibble                        3.1.8     2022-07-22 [1] RSPM (R 4.2.0)
 tidyr                       * 1.2.1     2022-09-08 [1] RSPM (R 4.2.0)
 tidyselect                    1.2.0     2022-10-10 [1] RSPM (R 4.2.0)
 tzdb                          0.3.0     2022-03-28 [1] RSPM (R 4.2.0)
 utf8                          1.2.2     2021-07-24 [1] RSPM (R 4.2.0)
 vctrs                         0.5.0     2022-10-22 [1] RSPM (R 4.2.0)
 VennDiagram                 * 1.6.20    2018-03-28 [2] RSPM (R 4.0.0)
 viridisLite                   0.4.1     2022-08-22 [1] RSPM (R 4.2.0)
 webshot                       0.5.2     2019-11-22 [2] RSPM (R 4.0.3)
 wesanderson                 * 0.3.6     2018-04-20 [1] RSPM (R 4.0.0)
 withr                         2.5.0     2022-03-03 [1] RSPM (R 4.2.0)
 xfun                          0.34      2022-10-18 [1] RSPM (R 4.2.0)
 XML                           3.99-0.11 2022-10-03 [1] RSPM (R 4.2.0)
 xml2                          1.3.3     2021-11-30 [1] RSPM (R 4.2.0)
 xtable                        1.8-4     2019-04-21 [1] RSPM (R 4.2.0)
 XVector                     * 0.38.0    2022-11-01 [1] Bioconductor
 yaml                          2.3.6     2022-10-18 [1] RSPM (R 4.2.0)
 zlibbioc                      1.44.0    2022-11-01 [1] Bioconductor
 zoo                           1.8-11    2022-09-17 [1] RSPM (R 4.2.0)
```

