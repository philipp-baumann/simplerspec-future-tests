Tests for speeding up the simplerspec processing workflow using futures
================
Philipp Baumann || <philipp.baumann@usys.ethz.ch>

Load R packages
===============

``` r
pkgs <- c("tidyverse", "simplerspec", "future", "future.apply", "doFuture")
lapply(pkgs, library, character.only = TRUE)
```

Plan the future
===============

``` r
## Define how `future()s` are resolved
# doFuture: A Universal Foreach Parallel Adaptor using the Future API of the 'future' Package
library("doFuture")
registerDoFuture()
plan(multicore)
availableCores()
```

    ## system 
    ##      8

Performance test: Reading and processing spectral data
======================================================

``` r
################################################################################
## Part 1: Read and pre-process spectra
################################################################################

## Read spectra in list ========================================================

# List of OPUS binary spectra files
lf <- dir("data/spectra", full.names = TRUE)

# Read spectra from files into R list: sequential ------------------------------

t_start <- Sys.time()
system.time(
  spc_list <- read_opus_univ(fnames = lf, extract = c("spc"))
)
```

    ##    user  system elapsed 
    ##  17.917   0.469  18.607

``` r
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 18.65672 secs

``` r
## Parallel implementation using doFuture --------------------------------------

t_start <- Sys.time()
system.time(
  spc_list <- read_opus_univ(fnames = lf, extract = c("spc"), parallel = TRUE)
)
```

    ##    user  system elapsed 
    ##  32.791   2.683   7.707

``` r
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 7.825047 secs

``` r
## Parallel implementation using doParallel ------------------------------------

# Allows to tune the models using parallel processing (e.g. use all available
# cores of a CPU); caret package automatically detects the registered backend
library("doParallel")
# Make a cluster with all possible threads (more than physical cores)
cl <- makeCluster(detectCores())
# Register backend
registerDoParallel(cl)
# Return number of parallel workers
getDoParWorkers() # 8 threads on MacBook Pro (Retina, 15-inch, Mid 2015);
```

    ## [1] 8

``` r
# Quadcore processor

t_start <- Sys.time()
system.time(
  spc_list <- read_opus_univ(fnames = lf, extract = c("spc"), parallel = TRUE)
)
```

    ##    user  system elapsed 
    ##   0.582   0.053   4.770

``` r
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 4.902125 secs

Gather spectral data from list into tibble
==========================================

``` r
# Gather spectra from nested list into tibble
t_start <- Sys.time()
system.time(
  spc_gathered <- gather_spc(spc_list)
)
```

    ##    user  system elapsed 
    ##   0.036   0.000   0.036

``` r
t_end <- Sys.time()

spc_gathered
```

    ## # A tibble: 284 x 6
    ##    unique_id        file_id    sample_id   spc       wavenumbers metadata 
    ##    <chr>            <chr>      <chr>       <list>    <list>      <list>   
    ##  1 BF_lo_01_soil_c… BF_lo_01_… BF_lo_01_s… <data.ta… <dbl [1,71… <tibble …
    ##  2 BF_lo_01_soil_c… BF_lo_01_… BF_lo_01_s… <data.ta… <dbl [1,71… <tibble …
    ##  3 BF_lo_01_soil_c… BF_lo_01_… BF_lo_01_s… <data.ta… <dbl [1,71… <tibble …
    ##  4 BF_lo_02_soil_c… BF_lo_02_… BF_lo_02_s… <data.ta… <dbl [1,71… <tibble …
    ##  5 BF_lo_02_soil_c… BF_lo_02_… BF_lo_02_s… <data.ta… <dbl [1,71… <tibble …
    ##  6 BF_lo_02_soil_c… BF_lo_02_… BF_lo_02_s… <data.ta… <dbl [1,71… <tibble …
    ##  7 BF_lo_03_soil_c… BF_lo_03_… BF_lo_03_s… <data.ta… <dbl [1,71… <tibble …
    ##  8 BF_lo_03_soil_c… BF_lo_03_… BF_lo_03_s… <data.ta… <dbl [1,71… <tibble …
    ##  9 BF_lo_03_soil_c… BF_lo_03_… BF_lo_03_s… <data.ta… <dbl [1,71… <tibble …
    ## 10 BF_lo_04_soil_c… BF_lo_04_… BF_lo_04_s… <data.ta… <dbl [1,71… <tibble …
    ## # ... with 274 more rows

Scenario 0: Use sequential processing (single-threaded R process)
-----------------------------------------------------------------

``` r
t_start <- Sys.time()
system.time(
  spc_proc <- spc_gathered %>%
    # Resample spectra to new wavenumber interval
    resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
    # Average replicate scans per sample_id
    average_spc() %>%
    # Preprocess spectra using Savitzky-Golay first derivative with
    # a window size of 21 points
    preprocess_spc(select = "sg_1_w21")
)
```

    ##    user  system elapsed 
    ##  14.098   0.079  14.192

``` r
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 14.26563 secs

Scenario 1: Chunk tibble into partitions and process partitions in parallel
---------------------------------------------------------------------------

``` r
# Inspired from multidplyr::partition();
# see https://github.com/hadley/multidplyr/blob/master/R/shard.R
partition_spc <- function(spc_tbl, 
                          groups = future::availableCores(),
                          id_nopart = sample_id) {
  id_nopart <- enquo(id_nopart)
  spc_tbl_nested <- spc_tbl %>%
    dplyr::group_by(!!id_nopart) %>%
    tidyr::nest()
  n <- nrow(spc_tbl_nested)
  m <- groups
  part_id <- sample(floor(m * (seq_len(n) - 1L) / n + 1L))
  
  spc_tbl_nested %>%
    dplyr::mutate(part_id = as.integer(part_id)) %>%
    tidyr::unnest()
}

spc_gathered %>%
    partition_spc()
```

    ## # A tibble: 284 x 7
    ##    sample_id  part_id unique_id     file_id  spc     wavenumbers metadata 
    ##    <chr>        <int> <chr>         <chr>    <list>  <list>      <list>   
    ##  1 BF_lo_01_…       6 BF_lo_01_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  2 BF_lo_01_…       6 BF_lo_01_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  3 BF_lo_01_…       6 BF_lo_01_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  4 BF_lo_02_…       3 BF_lo_02_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  5 BF_lo_02_…       3 BF_lo_02_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  6 BF_lo_02_…       3 BF_lo_02_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  7 BF_lo_03_…       6 BF_lo_03_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  8 BF_lo_03_…       6 BF_lo_03_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ##  9 BF_lo_03_…       6 BF_lo_03_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ## 10 BF_lo_04_…       2 BF_lo_04_soi… BF_lo_0… <data.… <dbl [1,71… <tibble …
    ## # ... with 274 more rows

### Scenario 1 variant a: use `furrr::future_map()`

``` r
system.time(
spc_proc <- spc_gathered %>%
    partition_spc() %>%
    split(f = .$part_id) %>%
    furrr::future_map(
      ~ .x %>%
          # Resample spectra to new wavenumber interval
          resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
          # Average replicate scans per sample_id
          average_spc() %>%
          # Preprocess spectra using Savitzky-Golay first derivative with
          # a window size of 21 points
          preprocess_spc(select = "sg_1_w21")
    )
)
```

    ##    user  system elapsed 
    ##  30.114   3.290   6.362

``` r
spc_proc
```

    ## $`1`
    ## # A tibble: 36 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_09…       1 BF_lo_09… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_09…       1 BF_lo_09… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_09…       1 BF_lo_09… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_lo_12…       1 BF_lo_12… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_lo_12…       1 BF_lo_12… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_lo_12…       1 BF_lo_12… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 BF_lo_17…       1 BF_lo_17… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 BF_lo_17…       1 BF_lo_17… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 BF_lo_17…       1 BF_lo_17… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 BF_lo_18…       1 BF_lo_18… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 26 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`2`
    ## # A tibble: 38 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_06…       2 BF_lo_06… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_06…       2 BF_lo_06… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_06…       2 BF_lo_06… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_mo_08…       2 BF_mo_08… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_mo_08…       2 BF_mo_08… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_mo_09…       2 BF_mo_09… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 BF_mo_09…       2 BF_mo_09… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 BF_mo_09…       2 BF_mo_09… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 BF_mo_17…       2 BF_mo_17… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 BF_mo_17…       2 BF_mo_17… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 28 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`3`
    ## # A tibble: 36 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_03…       3 BF_lo_03… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_03…       3 BF_lo_03… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_03…       3 BF_lo_03… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_mo_13…       3 BF_mo_13… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_mo_13…       3 BF_mo_13… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_mo_13…       3 BF_mo_13… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 CI_sb_01…       3 CI_sb_01… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 CI_sb_01…       3 CI_sb_01… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 CI_sb_01…       3 CI_sb_01… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 CI_sb_04…       3 CI_sb_04… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 26 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`4`
    ## # A tibble: 33 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_01…       4 BF_lo_01… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_01…       4 BF_lo_01… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_01…       4 BF_lo_01… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_lo_02…       4 BF_lo_02… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_lo_02…       4 BF_lo_02… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_lo_02…       4 BF_lo_02… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 BF_lo_08…       4 BF_lo_08… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 BF_lo_08…       4 BF_lo_08… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 BF_lo_08…       4 BF_lo_08… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 BF_mo_01…       4 BF_mo_01… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 23 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`5`
    ## # A tibble: 36 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_04…       5 BF_lo_04… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_04…       5 BF_lo_04… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_04…       5 BF_lo_04… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_lo_10…       5 BF_lo_10… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_lo_10…       5 BF_lo_10… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_lo_10…       5 BF_lo_10… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 BF_lo_19…       5 BF_lo_19… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 BF_lo_19…       5 BF_lo_19… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 BF_lo_19…       5 BF_lo_19… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 BF_mo_11…       5 BF_mo_11… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 26 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`6`
    ## # A tibble: 36 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_05…       6 BF_lo_05… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_05…       6 BF_lo_05… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_05…       6 BF_lo_05… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_lo_11…       6 BF_lo_11… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_lo_11…       6 BF_lo_11… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_lo_11…       6 BF_lo_11… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 BF_lo_16…       6 BF_lo_16… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 BF_lo_16…       6 BF_lo_16… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 BF_lo_16…       6 BF_lo_16… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 BF_mo_19…       6 BF_mo_19… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 26 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`7`
    ## # A tibble: 36 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_14…       7 BF_lo_14… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_14…       7 BF_lo_14… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_14…       7 BF_lo_14… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_lo_20…       7 BF_lo_20… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_lo_20…       7 BF_lo_20… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_lo_20…       7 BF_lo_20… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 CI_sb_09…       7 CI_sb_09… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 CI_sb_09…       7 CI_sb_09… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 CI_sb_09…       7 CI_sb_09… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 CI_sb_15…       7 CI_sb_15… CI_sb_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 26 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>
    ## 
    ## $`8`
    ## # A tibble: 33 x 12
    ##    sample_id part_id unique_id file_id spc   wavenumbers metadata spc_rs
    ##    <chr>       <int> <chr>     <chr>   <lis> <list>      <list>   <list>
    ##  1 BF_lo_07…       8 BF_lo_07… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  2 BF_lo_07…       8 BF_lo_07… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  3 BF_lo_07…       8 BF_lo_07… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  4 BF_lo_13…       8 BF_lo_13… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  5 BF_lo_13…       8 BF_lo_13… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  6 BF_lo_13…       8 BF_lo_13… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  7 BF_lo_15…       8 BF_lo_15… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  8 BF_lo_15…       8 BF_lo_15… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ##  9 BF_lo_15…       8 BF_lo_15… BF_lo_… <dat… <dbl [1,71… <tibble… <data…
    ## 10 BF_mo_02…       8 BF_mo_02… BF_mo_… <dat… <dbl [1,71… <tibble… <data…
    ## # ... with 23 more rows, and 4 more variables: wavenumbers_rs <list>,
    ## #   spc_mean <list>, spc_pre <list>, xvalues_pre <list>

### Scenario 1 variant b: use `future.apply::future_lapply()`

``` r
t_start <- Sys.time()
system.time(
spc_proc <- spc_gathered %>%
    partition_spc() %>%
    split(f = .$part_id) %>%
    future.apply::future_lapply(
      function(x) x %>%
          # Resample spectra to new wavenumber interval
          resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
          # Average replicate scans per sample_id
          average_spc() %>%
          # Preprocess spectra using Savitzky-Golay first derivative with
          # a window size of 21 points
          preprocess_spc(select = "sg_1_w21")
    )
)
```

    ##    user  system elapsed 
    ##  25.249   3.647   6.061

``` r
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 6.334516 secs

All in once: Partition, process in parallel, and row bind/gather results into tibble
------------------------------------------------------------------------------------

``` r
plan(multiprocess)

t_start <- Sys.time()
system.time(
  spc_proc <- spc_gathered %>%
    partition_spc() %>%
    split(f = .$part_id) %>%
    furrr::future_map(
      ~ .x %>%
      # Resample spectra to new wavenumber interval
      resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
      # Average replicate scans per sample_id
      average_spc() %>%
      # Preprocess spectra using Savitzky-Golay first derivative with
      # a window size of 21 points
      preprocess_spc(select = "sg_1_w21")
      ) %>%
      dplyr::bind_rows()
)
```

    ##    user  system elapsed 
    ##  25.234   2.667   4.854

``` r
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 5.096391 secs

Scenario 2: Partition and parallelize with foreach
--------------------------------------------------

``` r
t_start <- Sys.time()
spc_proc_l <- spc_gathered %>%
  partition_spc() %>%
  split(f = .$part_id)

  spc_proc <- foreach::foreach(i = seq_along(spc_proc_l),
    # .export = "read_opus_bin_univ",
    .errorhandling = "pass", # error object generated by task evaluation will
    # be included with the rest of the results
    # export the foreach package to the individual workers
    .packages = c("foreach", "simplerspec", "magrittr")) %dopar% {
      # Resample spectra to new wavenumber interval
      spc_proc_l[[i]] %>%
        resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
        # Average replicate scans per sample_id
        average_spc() %>%
        # Preprocess spectra using Savitzky-Golay first derivative with
        # a window size of 21 points
        preprocess_spc(select = "sg_1_w21")
    }
t_end <- Sys.time()
t_end - t_start
```

    ## Time difference of 7.108063 secs

Session info
============

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] de_CH.UTF-8/de_CH.UTF-8/de_CH.UTF-8/C/de_CH.UTF-8/de_CH.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2     doParallel_1.0.11  doFuture_0.6.0    
    ##  [4] iterators_1.0.10   future.apply_1.0.0 future_1.9.0      
    ##  [7] simplerspec_0.1.0  foreach_1.4.4      forcats_0.3.0     
    ## [10] stringr_1.3.1      dplyr_0.7.6        purrr_0.2.5       
    ## [13] readr_1.1.1        tidyr_0.8.1        tibble_1.4.2      
    ## [16] ggplot2_3.0.0      tidyverse_1.2.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.4          listenv_0.7.0            
    ##  [3] haven_1.1.2               lattice_0.20-35          
    ##  [5] colorspace_1.3-2          htmltools_0.3.6          
    ##  [7] yaml_2.2.0                utf8_1.1.4               
    ##  [9] rlang_0.2.1               pillar_1.3.0             
    ## [11] glue_1.3.0                withr_2.1.2              
    ## [13] prospectr_0.1.3           modelr_0.1.2             
    ## [15] readxl_1.1.0              bindr_0.1.1              
    ## [17] plyr_1.8.4                munsell_0.5.0            
    ## [19] gtable_0.2.0              cellranger_1.1.0         
    ## [21] rvest_0.3.2               codetools_0.2-15         
    ## [23] evaluate_0.11             RcppArmadillo_0.8.600.0.0
    ## [25] knitr_1.20                fansi_0.2.3              
    ## [27] furrr_0.1.0               broom_0.5.0              
    ## [29] Rcpp_0.12.18              scales_0.5.0             
    ## [31] backports_1.1.2           jsonlite_1.5             
    ## [33] hms_0.4.2                 digest_0.6.15            
    ## [35] stringi_1.2.4             grid_3.5.1               
    ## [37] rprojroot_1.3-2           cli_1.0.0                
    ## [39] tools_3.5.1               magrittr_1.5             
    ## [41] lazyeval_0.2.1            crayon_1.3.4             
    ## [43] pkgconfig_2.0.1           data.table_1.11.4        
    ## [45] xml2_1.2.0                lubridate_1.7.4          
    ## [47] assertthat_0.2.0          rmarkdown_1.10           
    ## [49] httr_1.3.1                rstudioapi_0.7           
    ## [51] globals_0.12.1            R6_2.2.2                 
    ## [53] hexView_0.3-3             nlme_3.1-137             
    ## [55] compiler_3.5.1
