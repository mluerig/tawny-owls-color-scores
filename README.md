# ACCESS INFORMATION

## 1. Licenses/restrictions placed on the data or code

Creative Commons Attribution 4.0 International

## DATA & CODE FILE OVERVIEW

This data repository consist of two data files, one code script, and this README document

## Data files

1. data_individuals.csv
-> contains records for breeding individuals in long-format 

2. data_recruits.csv: background darkness of the collected substrates
-> contains records for recruits in long-format, along with parental information


## Code scripts and workflow

1. run the script analysis.R from top to bottom to reproduce all results, including figures (except figure 1A, S1, and S3) and tables.

# SOFTWARE VERSIONS
> sessionInfo()
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8    LC_MONETARY=English_Germany.utf8 LC_NUMERIC=C                    
[5] LC_TIME=English_Germany.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] sommer_4.4.1           crayon_1.5.3           MASS_7.3-61            nadiv_2.18.0           marginaleffects_0.25.1
 [6] stringr_1.5.1          plot3Drgl_1.0.4        plot3D_1.4.1           rgl_1.3.18             YesSiR_0.0.0.9000     
[11] plotrix_3.8-4          officer_0.6.8          ggwordcloud_0.6.2      flextable_0.9.7        viridis_0.6.5         
[16] viridisLite_0.4.2      ggridges_0.5.6         cowplot_1.1.3          ggplot2_3.5.2          qgraph_1.9.8          
[21] lavaanPlot_0.8.1       lavaan_0.6-19          emmeans_1.11.0         gratia_0.10.0          lmerTest_3.1-3        
[26] lme4_1.1-37            Matrix_1.7-1           scales_1.3.0           mgcv_1.9-1             nlme_3.1-166          
[31] data.table_1.17.0      pacman_0.5.1          

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0          magrittr_2.0.3          TH.data_1.1-3          
  [6] estimability_1.5.1      farver_2.1.2            nloptr_2.2.1            rmarkdown_2.29          ragg_1.4.0             
 [11] vctrs_0.6.5             minqa_1.2.8             askpass_1.2.1           base64enc_0.1-3         htmltools_0.5.8.1      
 [16] Formula_1.2-5           htmlwidgets_1.6.4       plyr_1.8.9              sandwich_3.1-1          zoo_1.8-14             
 [21] uuid_1.2-1              misc3d_0.9-1            igraph_2.1.4            lifecycle_1.0.4         pkgconfig_2.0.3        
 [26] R6_2.6.1                fastmap_1.2.0           rbibutils_2.3           digest_0.6.37           numDeriv_2016.8-1.1    
 [31] fdrtool_1.2.18          colorspace_2.1-1        patchwork_1.3.0         textshaping_1.0.0       Hmisc_5.2-3            
 [36] abind_1.4-8             compiler_4.4.2          fontquiver_0.2.1        withr_3.0.2             glasso_1.11            
 [41] htmlTable_2.4.3         backports_1.5.0         psych_2.5.3             openssl_2.3.2           corpcor_1.6.10         
 [46] gtools_3.9.5            tools_4.4.2             pbivnorm_0.6.0          foreign_0.8-87          zip_2.3.2              
 [51] nnet_7.3-19             glue_1.8.0              quadprog_1.5-8          DiagrammeR_1.0.11       gridtext_0.1.5         
 [56] grid_4.4.2              checkmate_2.3.2         cluster_2.1.6           reshape2_1.4.4          generics_0.1.3         
 [61] gtable_0.3.6            tidyr_1.3.1             xml2_1.3.8              pillar_1.10.2           splines_4.4.2          
 [66] dplyr_1.1.4             lattice_0.22-6          survival_3.7-0          tidyselect_1.2.1        fontLiberation_0.1.0   
 [71] pbapply_1.7-2           knitr_1.50              fontBitstreamVera_0.1.1 reformulas_0.4.0        gridExtra_2.3          
 [76] stats4_4.4.2            xfun_0.52               visNetwork_2.1.2        stringi_1.8.7           boot_1.3-31            
 [81] ggokabeito_0.1.0        evaluate_1.0.3          codetools_0.2-20        tcltk_4.4.2             gdtools_0.4.2          
 [86] tibble_3.2.1            cli_3.6.4               rpart_4.1.23            xtable_1.8-4            systemfonts_1.2.2      
 [91] Rdpack_2.6.4            munsell_0.5.1           Rcpp_1.0.14             coda_0.19-4.1           png_0.1-8              
 [96] parallel_4.4.2          mvnfast_0.2.8           jpeg_0.1-11             mvtnorm_1.3-3           purrr_1.0.4            
[101] rlang_1.1.6             multcomp_1.4-28         mnormt_2.1.1  
