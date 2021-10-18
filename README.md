
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Micro4all

<!-- badges: start -->
<!-- badges: end -->

The goal of is to unify all functions used in our laboratory to get a
complete analysis of amplicon data for microbial community studies. It
includes several useful tools, as an *automatized detection of MOCK
community* and calculation of *community cut-off* as well as loops for
statistical analysis.

## Installation

You can install directly from GitHub with:

``` r
devtools::install_github("nuriamw/micro4all")
#> Downloading GitHub repo nuriamw/micro4all@HEAD
#> lifecycle    (1.0.0        -> 1.0.1       ) [CRAN]
#> pillar       (1.6.2        -> 1.6.3       ) [CRAN]
#> vroom        (1.5.4        -> 1.5.5       ) [CRAN]
#> stringi      (1.7.3        -> 1.7.5       ) [CRAN]
#> cpp11        (0.3.1        -> 0.4.0       ) [CRAN]
#> readr        (2.0.1        -> 2.0.2       ) [CRAN]
#> hms          (1.1.0        -> 1.1.1       ) [CRAN]
#> tibble       (3.1.3        -> 3.1.5       ) [CRAN]
#> data.table   (1.14.0       -> 1.14.2      ) [CRAN]
#> RcppArmad... (0.10.6.0.0   -> 0.10.7.0.0  ) [CRAN]
#> matrixStats  (0.60.1       -> 0.61.0      ) [CRAN]
#> xfun         (0.25         -> 0.26        ) [CRAN]
#> tidyr        (1.1.3        -> 1.1.4       ) [CRAN]
#> knitr        (1.33         -> 1.36        ) [CRAN]
#> R.utils      (2.10.1       -> 2.11.0      ) [CRAN]
#> bslib        (0.2.5.1      -> 0.3.1       ) [CRAN]
#> mime         (0.11         -> 0.12        ) [CRAN]
#> httpuv       (1.6.2        -> 1.6.3       ) [CRAN]
#> digest       (0.6.27       -> 0.6.28      ) [CRAN]
#> e1071        (1.7-8        -> 1.7-9       ) [CRAN]
#> htmltools    (0.5.1.1      -> 0.5.2       ) [CRAN]
#> styler       (1.5.1        -> 1.6.2       ) [CRAN]
#> shiny        (1.6.0        -> 1.7.1       ) [CRAN]
#> questionr    (0.7.4        -> 0.7.5       ) [CRAN]
#> maptools     (1.1-1        -> 1.1-2       ) [CRAN]
#> pairwiseA... (ce41196b1... -> ece560d23...) [GitHub]
#> Skipping 2 packages not available: microbiome, BiocGenerics
#> Installing 25 packages: lifecycle, pillar, vroom, stringi, cpp11, readr, hms, tibble, data.table, RcppArmadillo, matrixStats, xfun, tidyr, knitr, R.utils, bslib, mime, httpuv, digest, e1071, htmltools, styler, shiny, questionr, maptools
#> Installing packages into 'C:/Users/nuria/AppData/Local/Temp/RtmpgV2i50/temp_libpath5d074bcdc4'
#> (as 'lib' is unspecified)
#> package 'lifecycle' successfully unpacked and MD5 sums checked
#> package 'pillar' successfully unpacked and MD5 sums checked
#> package 'vroom' successfully unpacked and MD5 sums checked
#> package 'stringi' successfully unpacked and MD5 sums checked
#> package 'cpp11' successfully unpacked and MD5 sums checked
#> package 'readr' successfully unpacked and MD5 sums checked
#> package 'hms' successfully unpacked and MD5 sums checked
#> package 'tibble' successfully unpacked and MD5 sums checked
#> package 'data.table' successfully unpacked and MD5 sums checked
#> package 'RcppArmadillo' successfully unpacked and MD5 sums checked
#> package 'matrixStats' successfully unpacked and MD5 sums checked
#> package 'xfun' successfully unpacked and MD5 sums checked
#> package 'tidyr' successfully unpacked and MD5 sums checked
#> package 'knitr' successfully unpacked and MD5 sums checked
#> package 'R.utils' successfully unpacked and MD5 sums checked
#> package 'bslib' successfully unpacked and MD5 sums checked
#> package 'mime' successfully unpacked and MD5 sums checked
#> package 'httpuv' successfully unpacked and MD5 sums checked
#> package 'digest' successfully unpacked and MD5 sums checked
#> package 'e1071' successfully unpacked and MD5 sums checked
#> package 'htmltools' successfully unpacked and MD5 sums checked
#> package 'styler' successfully unpacked and MD5 sums checked
#> package 'shiny' successfully unpacked and MD5 sums checked
#> package 'questionr' successfully unpacked and MD5 sums checked
#> package 'maptools' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\nuria\AppData\Local\Temp\RtmpyEejx9\downloaded_packages
#> Downloading GitHub repo pmartinezarbizu/pairwiseAdonis@HEAD
#> 
#>          checking for file 'C:\Users\nuria\AppData\Local\Temp\RtmpyEejx9\remotes26b46d14e42\pmartinezarbizu-pairwiseAdonis-ece560d\pairwiseAdonis/DESCRIPTION' ...  v  checking for file 'C:\Users\nuria\AppData\Local\Temp\RtmpyEejx9\remotes26b46d14e42\pmartinezarbizu-pairwiseAdonis-ece560d\pairwiseAdonis/DESCRIPTION'
#>       -  preparing 'pairwiseAdonis':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
#>       -  checking for LF line-endings in source and make files and shell scripts
#>   -  checking for empty or unneeded directories
#>      Omitted 'LazyData' from DESCRIPTION
#>       -  building 'pairwiseAdonis_0.4.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/nuria/AppData/Local/Temp/RtmpgV2i50/temp_libpath5d074bcdc4'
#> (as 'lib' is unspecified)
#>          checking for file 'C:\Users\nuria\AppData\Local\Temp\RtmpyEejx9\remotes26b47481516d\nuriamw-micro4all-aee8a09/DESCRIPTION' ...     checking for file 'C:\Users\nuria\AppData\Local\Temp\RtmpyEejx9\remotes26b47481516d\nuriamw-micro4all-aee8a09/DESCRIPTION' ...   v  checking for file 'C:\Users\nuria\AppData\Local\Temp\RtmpyEejx9\remotes26b47481516d\nuriamw-micro4all-aee8a09/DESCRIPTION' (344ms)
#>       -  preparing 'micro4all': (914ms)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
#>   Warning:     Warning: C:/Users/nuria/AppData/Local/Temp/RtmpklQarz/Rbuild3ccc78127f4a/micro4all/man/ASV_table_classified_raw.Rd:16: unknown macro '\item'
#>       -  checking for LF line-endings in source and make files and shell scripts
#>       -  checking for empty or unneeded directories
#>   Removed empty directory      Removed empty directory 'micro4all/vignettes'
#>       -  building 'micro4all_0.0.0.9000.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/nuria/AppData/Local/Temp/RtmpgV2i50/temp_libpath5d074bcdc4'
#> (as 'lib' is unspecified)
```

## How to use it?

Along with this package, a full and detailed workflow has been created.
There, you can find our usual pipeline for data analysis and examples of
functions. Check the tutorial [here](.docs\tutorial)
