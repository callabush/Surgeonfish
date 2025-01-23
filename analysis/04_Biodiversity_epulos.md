---
title: "BioDiversity"
author: "Calla Bush St George"
date: "2024-09-04"
output:
  html_document: 
    code_folding: show
    theme: spacelab
    highlight: pygments
    keep_md: yes
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 3
  keep_md: true  
editor_options: 
  chunk_output_type: console
---


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.align = "center",
                      fig.path = "../figures/07_Biodiversity_epulos/",
                      dev = "png", dpi = 200) 

# send any figure output to this folder 
```

# Setting the Environment 

### Set my seed

```r
# Any number can be chose
set.seed(567890)
```

## Load Libraries 

```r
pacman::p_load(tidyverse, devtools, patchwork, iNEXT, phyloseq,rstatix, ggpubr,
               install = FALSE)
```

## Load in Data 

```r
load("data/06_PreProcessing_epulos/raw_preprocessed_physeq.RData")
raw_preprocessed_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1112 taxa and 83 samples ]
## sample_data() Sample Data:       [ 83 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 1112 taxa by 9 taxonomic ranks ]
```

```r
# Intuition Check 
min(sample_sums(raw_preprocessed_physeq))
```

```
## [1] 559
```

```r
# Setting colors for gut sections 
host_species_colors <- c(
  "Acanthurus nigros" = "dodgerblue4",
  "Acanthurus achilles" = "royalblue",
  "Acanthurus nigrofuscus" = "blue3",
  "Acanthurus olivaceus" = "blue4",
  "Acanthurus triostegus" = "slateblue3",
  "Ctenochaetus striatus" = "olivedrab",
  "Naso brevirostris" = "red2",
  "Naso lituratus" = "tomato4",
  "Naso tonganus" = "salmon2",
  "Odax pullus" = "springgreen2",
  "Zebrasoma flavescens" = "darkorange",
  "Zebrasoma scopas" = "yellow2",
  "Zebrasoma velifer" = "goldenrod")

# Make metadata dataframe
metadata_df <-
  raw_preprocessed_physeq %>%
  sample_data() %>%
  data.frame()
```



# Goals

1. Calculate the Hill Diversity of the samples. 
2. Evaluate the rarefaction curves. 
3. Evaluate the Diversity values. 
4. Makes notes of specific samples and their seq depth. 

# Diversity Calculations with iNEXT 


```r
# prepare input data 
iNEXT_input_df <- 
  raw_preprocessed_physeq %>%
  otu_table() %>%
  data.frame()
# Quick check
dim(iNEXT_input_df)
```

```
## [1] 1112   83
```

```r
# Run iNEXT: Calculate the Hill Numbers 
# Note that: Species in ROWS, Samples in COLUMNS 
# Remember to set the seed! 
iNEXT_data <- iNEXT(iNEXT_input_df, 
                   q = c(0,1,2), datatype = "abundance")

# Save the file
save(iNEXT_data, file = "data/08_Biodiversity_epulos/iNEXT_data.RData")
```

# Gut section colors

```r
gut_section_colors <- c(
  III = "darkslateblue",
  IV = "darkorchid1",
  V = "cyan3",
  feces = "cornflowerblue")
```

# Region colors

```r
region_colors <- c(
  "Cook Islands" = "olivedrab4",
  "Hawaii" = "deeppink",
  GBR = "royalblue3",
  Captivity = "brown")
```

# Diet colors

```r
diet_colors <- c(
  herbivore = "forestgreen",
  detritovore = "hotpink4",
  omnivore = "orange")
```


# Evaluate the Diversity! 

```r
load("data/08_Biodiversity_epulos/iNEXT_data.RData")
str(iNEXT_data)
```

```
## List of 3
##  $ DataInfo:'data.frame':	83 obs. of  14 variables:
##   ..$ Assemblage: chr [1:83] "A16" "A17" "A18" "AA1" ...
##   ..$ n         : num [1:83] 15903 20099 1464 2920 9125 ...
##   ..$ S.obs     : num [1:83] 73 81 28 32 43 69 20 50 92 10 ...
##   ..$ SC        : num [1:83] 1 1 1 1 1 ...
##   ..$ f1        : num [1:83] 1 2 0 0 0 0 1 1 3 0 ...
##   ..$ f2        : num [1:83] 5 7 4 2 1 5 0 5 5 2 ...
##   ..$ f3        : num [1:83] 1 2 0 0 1 3 1 3 4 0 ...
##   ..$ f4        : num [1:83] 2 2 0 2 0 2 1 0 3 0 ...
##   ..$ f5        : num [1:83] 2 0 0 1 1 3 0 0 0 0 ...
##   ..$ f6        : num [1:83] 0 1 0 1 1 2 0 1 1 0 ...
##   ..$ f7        : num [1:83] 0 1 1 0 2 4 2 0 0 0 ...
##   ..$ f8        : num [1:83] 1 1 3 0 2 1 0 2 3 0 ...
##   ..$ f9        : num [1:83] 0 1 2 0 1 2 0 1 1 0 ...
##   ..$ f10       : num [1:83] 2 2 0 1 3 0 0 0 0 0 ...
##  $ iNextEst:List of 2
##   ..$ size_based    :'data.frame':	9960 obs. of  10 variables:
##   .. ..$ Assemblage: chr [1:9960] "A16" "A16" "A16" "A16" ...
##   .. ..$ m         : num [1:9960] 1 884 1767 2651 3534 ...
##   .. ..$ Method    : chr [1:9960] "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##   .. ..$ Order.q   : num [1:9960] 0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ qD        : num [1:9960] 1 58.1 63.1 65.5 67 ...
##   .. ..$ qD.LCL    : num [1:9960] 1 57.3 62.2 64.4 65.8 ...
##   .. ..$ qD.UCL    : num [1:9960] 1 59 64.1 66.6 68.3 ...
##   .. ..$ SC        : num [1:9960] 0.051 0.99 0.997 0.998 0.999 ...
##   .. ..$ SC.LCL    : num [1:9960] 0.0495 0.9899 0.9961 0.9976 0.9982 ...
##   .. ..$ SC.UCL    : num [1:9960] 0.0525 0.991 0.9969 0.9983 0.9988 ...
##   ..$ coverage_based:'data.frame':	8130 obs. of  8 variables:
##   .. ..$ Assemblage: chr [1:8130] "A16" "A16" "A16" "A16" ...
##   .. ..$ SC        : num [1:8130] 0.051 0.99 0.997 0.998 0.999 ...
##   .. ..$ m         : num [1:8130] 1 884 1767 2651 3534 ...
##   .. ..$ Method    : chr [1:8130] "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##   .. ..$ Order.q   : num [1:8130] 0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ qD        : num [1:8130] 1 58.1 63.1 65.5 67 ...
##   .. ..$ qD.LCL    : num [1:8130] 0.98 57.14 61.88 63.96 65.28 ...
##   .. ..$ qD.UCL    : num [1:8130] 1.02 59.16 64.41 67.01 68.76 ...
##  $ AsyEst  :'data.frame':	249 obs. of  7 variables:
##   ..$ Assemblage: chr [1:249] "A16" "A16" "A16" "A17" ...
##   ..$ Diversity : chr [1:249] "Species richness" "Shannon diversity" "Simpson diversity" "Species richness" ...
##   ..$ Observed  : num [1:249] 73 30.8 19.6 81 19.4 ...
##   ..$ Estimator : num [1:249] 73.1 30.9 19.6 81.3 19.5 ...
##   ..$ s.e.      : num [1:249] 3.752 0.245 0.245 5.608 0.201 ...
##   ..$ LCL       : num [1:249] 73 30.4 19.1 81 19.1 ...
##   ..$ UCL       : num [1:249] 80.5 31.4 20.1 92.3 19.9 ...
##  - attr(*, "class")= chr "iNEXT"
```

```r
typeof(iNEXT_data)
```

```
## [1] "list"
```

# Plot Diversity 

```r
# Prepare Colors 
color_df <- 
  iNEXT_input_df %>%
  colnames() %>%
  data.frame()
# Check
head(color_df)
```

```
##     .
## 1 A16
## 2 A17
## 3 A18
## 4 AA1
## 5 AA2
## 6 AA3
```

```r
# Rename the column 
colnames(color_df)[1] <- "names"
# Check
head(color_df)
```

```
##   names
## 1   A16
## 2   A17
## 3   A18
## 4   AA1
## 5   AA2
## 6   AA3
```

```r
# Make a helper dataframe for plotting with colors 
iNEXT_color_df <- 
  color_df %>%
  # Fix the names for merging
  mutate(names = gsub(names, pattern = "[.]", replace = "-"),
         names = gsub(names, pattern = "X",  replace = "")) %>%
  # Merge with metadata
  left_join(metadata_df, by = "names") %>%
   #Merge with colors for plotting with ggiNEXT
  left_join(y=data.frame(gut_section_colors = gut_section_colors,
            gut_section = names(gut_section_colors)),
            by = "gut_section")
```

# Plot Rarefaction with `ggiNEXT`


```r
# Plot rarefaction! 
# rarefaction/extrapolation curve, type = 1 

# Order q: 
  # 0 = Richness/ Number of Total taxa
  # 1 = Exponential Shannon / Number of "Common" taxa
  # 2 = Inverse Simpson / Number of "Dominant" taxa 

ggiNEXT(iNEXT_data, type = 1, facet.var = "Order.q") + 
  facet_wrap(~Order.q, scales = "fixed") + 
  scale_color_manual(values = iNEXT_color_df$gut_section_colors, guide = FALSE) + 
  scale_fill_manual(values = iNEXT_color_df$gut_section_colors, guide = FALSE) + 
  scale_shape_manual(values = base::rep(17, nsamples(raw_preprocessed_physeq)),
                     guide = FALSE) +
  theme(legend.position = "none")
```

```
## Scale for colour is already present.
## Adding another scale for colour, which will replace the existing scale.
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

<img src="../figures/07_Biodiversity_epulos/ggiNEXT-1.png" style="display: block; margin: auto;" />


# Manually plot Diversity 

## Rarefaction

```r
iNEXT_manual_df <- 
  iNEXT_data$iNextEst$size_based %>%
  dplyr::rename(names = Assemblage) %>%
  # Fix the samples names 
  mutate(names = gsub(names, pattern = "[.]", replace = "-"),
         names = gsub(names, pattern = "X", replace = "")) %>%
  # join with metadata 
  left_join(., metadata_df, by = "names") %>%
  # Add colors to data frame
  left_join(., data.frame(gut_section_colors = gut_section_colors,
                          gut_section = names(gut_section_colors)),
            by = "gut_section") 

# Inspect 
dim(iNEXT_manual_df)
```

```
## [1] 9960   20
```

```r
str(iNEXT_manual_df)
```

```
## 'data.frame':	9960 obs. of  20 variables:
##  $ names             : chr  "A16" "A16" "A16" "A16" ...
##  $ m                 : num  1 884 1767 2651 3534 ...
##  $ Method            : chr  "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##  $ Order.q           : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ qD                : num  1 58.1 63.1 65.5 67 ...
##  $ qD.LCL            : num  1 57.3 62.2 64.4 65.8 ...
##  $ qD.UCL            : num  1 59 64.1 66.6 68.3 ...
##  $ SC                : num  0.051 0.99 0.997 0.998 0.999 ...
##  $ SC.LCL            : num  0.0495 0.9899 0.9961 0.9976 0.9982 ...
##  $ SC.UCL            : num  0.0525 0.991 0.9969 0.9983 0.9988 ...
##  $ host_species      : chr  "Acanthurus nigros" "Acanthurus nigros" "Acanthurus nigros" "Acanthurus nigros" ...
##  $ gut_section       : chr  "IV" "IV" "IV" "IV" ...
##  $ region            : chr  "Cook Islands" "Cook Islands" "Cook Islands" "Cook Islands" ...
##  $ location          : chr  "Mitiaro Island" "Mitiaro Island" "Mitiaro Island" "Mitiaro Island" ...
##  $ diet              : chr  "herbivore" "herbivore" "herbivore" "herbivore" ...
##  $ year              : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ month             : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ day               : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ sample_lab        : chr  "S521" "S521" "S521" "S521" ...
##  $ gut_section_colors: chr  "darkorchid1" "darkorchid1" "darkorchid1" "darkorchid1" ...
```

```r
# Plot it - Rarefaction Curve 
iNEXT_manual_df %>%
  # Filter out rows that are calcaulted by rarefaction from iNEXT
  dplyr::filter(Method == "Extrapolation") %>%
  # Make the actual rarefaction plot with 
  # the # of sequences on the x-axis and diversity on the y-axis
  # You can choose to pick one diversity value or plot all three 
  ggplot(aes(x = m, y= qD, color = gut_section, group = names)) + 
  # line 
  geom_line() + 
  geom_point() + 
  # Challenge: Facet with gut section
  facet_grid(Order.q~gut_section, scales = "fixed") + 
  scale_color_manual(values = gut_section_colors) + 
  theme(legend.position = "bottom")
```

<img src="../figures/07_Biodiversity_epulos/iNEXT-manual-1.png" style="display: block; margin: auto;" />


# Diversity vs Gut Section 


```r
iNEXT_manual_df %>%
  dplyr::filter(Method == "Observed") %>%
  ggplot(aes(x = gut_section, y = qD)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(.~Order.q, scales = "free") + 
  geom_point(aes(color = gut_section)) + 
  stat_smooth() + 
  labs(x = "Host species", y = "# of ASVs") + 
  scale_color_manual(values = gut_section_colors) + 
  theme(legend.position = "bottom")
```

```
## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
```

<img src="../figures/07_Biodiversity_epulos/div-vs--1.png" style="display: block; margin: auto;" />


# Session Information 

```r
# Ensure reproducibility 
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.3.2 (2023-10-31)
##  os       macOS Sonoma 14.6.1
##  system   x86_64, darwin20
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2024-09-04
##  pandoc   3.1.11 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version    date (UTC) lib source
##  abind              1.4-5      2016-07-21 [2] CRAN (R 4.3.0)
##  ade4               1.7-22     2023-02-06 [2] CRAN (R 4.3.0)
##  ape                5.7-1      2023-03-13 [2] CRAN (R 4.3.0)
##  backports          1.4.1      2021-12-13 [2] CRAN (R 4.3.0)
##  Biobase            2.62.0     2023-10-24 [2] Bioconductor
##  BiocGenerics       0.48.1     2023-11-01 [2] Bioconductor
##  biomformat         1.30.0     2023-10-24 [2] Bioconductor
##  Biostrings         2.70.2     2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
##  bitops             1.0-7      2021-04-24 [2] CRAN (R 4.3.0)
##  broom              1.0.5      2023-06-09 [2] CRAN (R 4.3.0)
##  bslib              0.6.1      2023-11-28 [2] CRAN (R 4.3.0)
##  cachem             1.0.8      2023-05-01 [2] CRAN (R 4.3.0)
##  car                3.1-2      2023-03-30 [1] CRAN (R 4.3.0)
##  carData            3.0-5      2022-01-06 [1] CRAN (R 4.3.0)
##  cli                3.6.2      2023-12-11 [2] CRAN (R 4.3.0)
##  cluster            2.1.6      2023-12-01 [2] CRAN (R 4.3.0)
##  codetools          0.2-19     2023-02-01 [2] CRAN (R 4.3.2)
##  colorspace         2.1-0      2023-01-23 [2] CRAN (R 4.3.0)
##  crayon             1.5.2      2022-09-29 [2] CRAN (R 4.3.0)
##  data.table         1.15.0     2024-01-30 [2] CRAN (R 4.3.2)
##  devtools         * 2.4.5      2022-10-11 [2] CRAN (R 4.3.0)
##  digest             0.6.34     2024-01-11 [2] CRAN (R 4.3.0)
##  dplyr            * 1.1.4      2023-11-17 [2] CRAN (R 4.3.0)
##  ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.0)
##  evaluate           0.23       2023-11-01 [2] CRAN (R 4.3.0)
##  fansi              1.0.6      2023-12-08 [2] CRAN (R 4.3.0)
##  farver             2.1.1      2022-07-06 [2] CRAN (R 4.3.0)
##  fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.0)
##  forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.3.0)
##  foreach            1.5.2      2022-02-02 [2] CRAN (R 4.3.0)
##  fs                 1.6.3      2023-07-20 [2] CRAN (R 4.3.0)
##  generics           0.1.3      2022-07-05 [2] CRAN (R 4.3.0)
##  GenomeInfoDb       1.38.6     2024-02-08 [2] Bioconductor 3.18 (R 4.3.2)
##  GenomeInfoDbData   1.2.11     2023-11-13 [2] Bioconductor
##  ggplot2          * 3.4.4      2023-10-12 [1] CRAN (R 4.3.0)
##  ggpubr           * 0.6.0.999  2024-07-30 [1] Github (kassambara/ggpubr@6aeb4f7)
##  ggsignif           0.6.4      2022-10-13 [1] CRAN (R 4.3.0)
##  glue               1.7.0      2024-01-09 [1] CRAN (R 4.3.0)
##  gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.0)
##  highr              0.10       2022-12-22 [2] CRAN (R 4.3.0)
##  hms                1.1.3      2023-03-21 [2] CRAN (R 4.3.0)
##  htmltools          0.5.7      2023-11-03 [2] CRAN (R 4.3.0)
##  htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.3.0)
##  httpuv             1.6.14     2024-01-26 [2] CRAN (R 4.3.2)
##  igraph             2.0.1.1    2024-01-30 [2] CRAN (R 4.3.2)
##  iNEXT            * 3.0.0      2022-08-29 [2] CRAN (R 4.3.0)
##  IRanges            2.36.0     2023-10-24 [2] Bioconductor
##  iterators          1.0.14     2022-02-05 [2] CRAN (R 4.3.0)
##  jquerylib          0.1.4      2021-04-26 [2] CRAN (R 4.3.0)
##  jsonlite           1.8.8      2023-12-04 [2] CRAN (R 4.3.0)
##  knitr              1.45       2023-10-30 [2] CRAN (R 4.3.0)
##  labeling           0.4.3      2023-08-29 [2] CRAN (R 4.3.0)
##  later              1.3.2      2023-12-06 [2] CRAN (R 4.3.0)
##  lattice            0.22-5     2023-10-24 [2] CRAN (R 4.3.0)
##  lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.3.0)
##  lubridate        * 1.9.3      2023-09-27 [2] CRAN (R 4.3.0)
##  magrittr           2.0.3      2022-03-30 [2] CRAN (R 4.3.0)
##  MASS               7.3-60.0.1 2024-01-13 [2] CRAN (R 4.3.0)
##  Matrix             1.6-5      2024-01-11 [2] CRAN (R 4.3.0)
##  memoise            2.0.1      2021-11-26 [2] CRAN (R 4.3.0)
##  mgcv               1.9-1      2023-12-21 [2] CRAN (R 4.3.0)
##  mime               0.12       2021-09-28 [2] CRAN (R 4.3.0)
##  miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.0)
##  multtest           2.58.0     2023-10-24 [2] Bioconductor
##  munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.0)
##  nlme               3.1-164    2023-11-27 [2] CRAN (R 4.3.0)
##  pacman             0.5.1      2019-03-11 [1] CRAN (R 4.3.0)
##  patchwork        * 1.2.0.9000 2024-05-07 [1] Github (thomasp85/patchwork@d943757)
##  permute            0.9-7      2022-01-27 [2] CRAN (R 4.3.0)
##  phyloseq         * 1.46.0     2023-10-24 [2] Bioconductor
##  pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.0)
##  pkgbuild           1.4.3      2023-12-10 [2] CRAN (R 4.3.0)
##  pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.0)
##  pkgload            1.3.4      2024-01-16 [2] CRAN (R 4.3.0)
##  plyr               1.8.9      2023-10-02 [1] CRAN (R 4.3.0)
##  profvis            0.3.8      2023-05-02 [2] CRAN (R 4.3.0)
##  promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.0)
##  purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.0)
##  R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.0)
##  Rcpp               1.0.12     2024-01-09 [2] CRAN (R 4.3.0)
##  RCurl              1.98-1.14  2024-01-09 [2] CRAN (R 4.3.0)
##  readr            * 2.1.5      2024-01-10 [1] CRAN (R 4.3.0)
##  remotes            2.4.2.1    2023-07-18 [2] CRAN (R 4.3.0)
##  reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.0)
##  rhdf5              2.46.1     2023-11-29 [2] Bioconductor
##  rhdf5filters       1.14.1     2023-11-06 [2] Bioconductor
##  Rhdf5lib           1.24.2     2024-02-07 [2] Bioconductor 3.18 (R 4.3.2)
##  rlang              1.1.3      2024-01-10 [2] CRAN (R 4.3.0)
##  rmarkdown          2.25       2023-09-18 [2] CRAN (R 4.3.0)
##  rstatix          * 0.7.2      2023-02-01 [1] CRAN (R 4.3.0)
##  rstudioapi         0.15.0     2023-07-07 [2] CRAN (R 4.3.0)
##  S4Vectors          0.40.2     2023-11-23 [2] Bioconductor
##  sass               0.4.8      2023-12-06 [2] CRAN (R 4.3.0)
##  scales             1.3.0      2023-11-28 [2] CRAN (R 4.3.0)
##  sessioninfo        1.2.2      2021-12-06 [2] CRAN (R 4.3.0)
##  shiny              1.8.0      2023-11-17 [2] CRAN (R 4.3.0)
##  stringi            1.8.3      2023-12-11 [2] CRAN (R 4.3.0)
##  stringr          * 1.5.1      2023-11-14 [2] CRAN (R 4.3.0)
##  survival           3.5-7      2023-08-14 [2] CRAN (R 4.3.0)
##  tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.0)
##  tidyr            * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
##  tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.3.0)
##  tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
##  timechange         0.3.0      2024-01-18 [2] CRAN (R 4.3.0)
##  tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.3.0)
##  urlchecker         1.0.1      2021-11-30 [2] CRAN (R 4.3.0)
##  usethis          * 2.2.2      2023-07-06 [2] CRAN (R 4.3.0)
##  utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.0)
##  vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.0)
##  vegan              2.6-4      2022-10-11 [1] CRAN (R 4.3.0)
##  withr              3.0.0      2024-01-16 [2] CRAN (R 4.3.0)
##  xfun               0.42       2024-02-08 [2] CRAN (R 4.3.2)
##  xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.0)
##  XVector            0.42.0     2023-10-24 [2] Bioconductor
##  yaml               2.3.8      2023-12-11 [2] CRAN (R 4.3.0)
##  zlibbioc           1.48.0     2023-10-24 [2] Bioconductor
## 
##  [1] /Users/cab565/Library/R/x86_64/4.3/library
##  [2] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
