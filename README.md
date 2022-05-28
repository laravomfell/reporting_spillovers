# (No) Spillovers in reporting domestic abuse to police

This Github repo contains the code used to produce the analysis for _(No) Spillovers in reporting domestic abuse to police_ by Lara Vomfell and Jan Povala.

This research is based on data resources provided by an unnamed English police force and were provided to us under an Information Sharing Agreement. We cannot make the data publicly available or give the precise location of where the study took place.

To nonetheless provide a working version of our code, we generate synthetic data on a circle in the unit square. The `run.R` file calls all files used in our analysis, beginning with the generation of the synthetic data. You should call `run.R` from the command line to set the arguments for the synthetic data (e.g., how many events). 

Some parts of the inference procedure are adapted from [code](https://rss.onlinelibrary.wiley.com/hub/journal/1467985x/series-a-datasets/182_3) provided by J. Zhuang and J. Mateu for [their paper](https://doi.org/10.1111/rssa.12429) _A semiparametric spatiotemporal Hawkes-type point process model with periodic background for crime data_.

Please let me know if you find any errors!


```R
sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
[3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] tcltk     parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] stpp_2.0-6            splancs_2.01-43       sp_1.4-5              rpanel_1.1-5         
 [5] qqplotr_0.0.5         scales_1.2.0          ggplot2_3.3.6         snakecase_0.11.0     
 [9] doParallel_1.0.17     iterators_1.0.14      foreach_1.5.2         mvtnorm_1.1-3        
[13] polyCub_0.8.0         optparse_1.7.1        data.table_1.14.2     sf_1.0-7             
[17] fields_13.3           viridis_0.6.2         viridisLite_0.4.0     spam_2.8-0           
[21] Matrix_1.4-1          purrr_0.3.4           spatstat_2.3-4        spatstat.linnet_2.3-2
[25] spatstat.core_2.4-4   rpart_4.1-15          nlme_3.1-152          spatstat.random_2.2-0
[29] spatstat.geom_2.4-0   spatstat.data_2.2-0   tictoc_1.0.1         

loaded via a namespace (and not attached):
 [1] maps_3.4.0            jsonlite_1.7.2        splines_4.1.1         dotCall64_1.0-1      
 [5] assertthat_0.2.1      robustbase_0.95-0     pillar_1.6.4          lattice_0.20-44      
 [9] glue_1.5.1            digest_0.6.29         polyclip_1.10-0       colorspace_2.0-2     
[13] htmltools_0.5.2       spatstat.sparse_2.1-1 pkgconfig_2.0.3       misc3d_0.9-1         
[17] tensor_1.5            getopt_1.20.3         spatstat.utils_2.3-1  tibble_3.1.6         
[21] proxy_0.4-26          mgcv_1.8-36           generics_0.1.0        ellipsis_0.3.2       
[25] withr_2.4.3           cli_3.1.0             magrittr_2.0.1        crayon_1.4.2         
[29] deldir_1.0-6          fansi_0.5.0           MASS_7.3-54           class_7.3-19         
[33] tools_4.1.1           lifecycle_1.0.1       munsell_0.5.0         compiler_4.1.1       
[37] e1071_1.7-9           rlang_1.0.2           plot3D_1.4            classInt_0.4-3       
[41] units_0.7-2           grid_4.1.1            htmlwidgets_1.5.4     goftest_1.2-3        
[45] gtable_0.3.0          codetools_0.2-18      abind_1.4-5           DBI_1.1.1            
[49] R6_2.5.1              gridExtra_2.3         knitr_1.36            dplyr_1.0.7          
[53] fastmap_1.1.0         utf8_1.2.2            KernSmooth_2.23-20    Rcpp_1.0.8           
[57] vctrs_0.3.8           rgl_0.108.3.2         xfun_0.26             DEoptimR_1.0-11      
[61] tidyselect_1.1.1
``` 