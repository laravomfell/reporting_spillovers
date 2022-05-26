# (No) Spillovers in reporting domestic abuse to police

This Github repo contains the code used to produce the analysis for _(No) Spillovers in reporting domestic abuse to police_ by Lara Vomfell and Jan Povala.

This research is based on data resources provided by an unnamed English police force and were provided to us under an Information Sharing Agreement. We cannot make the data publicly available or give the precise location of where the study took place.

To nonetheless provide a working version of our code, we generate synthetic data on a circle in the unit square. The `run.R` file calls all files used in our analysis, beginning with the generation of the synthetic data.

Some parts of the inference procedure are adapted from [code](https://rss.onlinelibrary.wiley.com/hub/journal/1467985x/series-a-datasets/182_3) provided by J. Zhuang and J. Mateu for [their paper](https://doi.org/10.1111/rssa.12429) _A semiparametric spatiotemporal Hawkes-type point process model with periodic background for crime data_.

Please let me know if you find any errors!


```R
sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.1        ggplot2_3.3.3       snakecase_0.11.0    doParallel_1.0.16   iterators_1.0.13    foreach_1.5.1      
 [7] mvtnorm_1.1-1       polyCub_0.8.0       data.table_1.13.6   sf_0.9-7            fields_11.6         spam_2.5-1         
[13] dotCall64_1.0-0     Matrix_1.3-2        purrr_0.3.4         spatstat_1.64-1     rpart_4.1-15        nlme_3.1-150       
[19] spatstat.data_2.0-0

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.0     gpclib_1.5-6         splines_4.0.3        lattice_0.20-41      colorspace_1.4-1     spatstat.utils_2.0-0
 [7] vctrs_0.3.6          generics_0.0.2       mgcv_1.8-33          rlang_0.4.10         e1071_1.7-4          pillar_1.4.7        
[13] withr_2.4.1          glue_1.4.2           DBI_1.1.0            sp_1.4-4             lifecycle_0.2.0      munsell_0.5.0       
[19] gtable_0.3.0         codetools_0.2-16     class_7.3-17         Rcpp_1.0.6           KernSmooth_2.23-17   tensor_1.5          
[25] classInt_0.4-3       abind_1.4-5          deldir_0.1-29        dplyr_1.0.2          polyclip_1.10-0      tools_4.0.3         
[31] magrittr_2.0.1       maps_3.3.0           goftest_1.2-2        tibble_3.0.6         crayon_1.4.1         pkgconfig_2.0.3     
[37] ellipsis_0.3.1       praise_1.0.0         R6_2.5.0             units_0.6-7          compiler_4.0.3  
``` 