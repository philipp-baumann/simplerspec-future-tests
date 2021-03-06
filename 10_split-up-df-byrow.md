Speed up splitting of spectral data.table by row to create list-column
================
Philipp Baumann
10 November 2018

## Motivation

The simplerspec processing pipeline functions return spectral data as
list-columns. A newly computed spectral column is represented as list of
data.table’s, one per row or single spectrum. Therefore, a single
data.table containing all spectra needs to be split by row into list of
data frames for all processing steps. This facilitates exploring the
input and output data stream at the different steps of pipeline.
However, the current implementation is rather slow because these
splitting operations rely on `base::split()`. In contrast, the opposite
operation is extremely fast, as e.g. implemented in
`data.table::rbindlist()`.

There are different ways how to split a data frame by row into a list of
data frames. Luckily, there is promising ways of splitting by row, as a
comparison by Winston Chang and updates of Jenny Bryan illustrates
[here](https://github.com/jennybc/row-oriented-workflows/blob/master/iterate-over-rows.md).
The goal is to showcase the speed gain that results from a better split
implementation using `purrr::transponse()` to do the same job, but much
faster.

## Test alternative splitting within `simplerspec::average_spc()`

``` r
# Load required packages
pkgs <- c("here", "tidyverse", "simplerspec", "doFuture")
lapply(pkgs, library, character.only = TRUE)
```

    ## Warning: package 'ggplot2' was built under R version 3.4.4

    ## Warning: package 'tibble' was built under R version 3.4.3

    ## Warning: package 'tidyr' was built under R version 3.4.4

    ## Warning: package 'purrr' was built under R version 3.4.4

    ## Warning: package 'dplyr' was built under R version 3.4.4

    ## Warning: package 'stringr' was built under R version 3.4.4

    ## Warning: package 'iterators' was built under R version 3.4.3

``` r
# Register and plan multicore processsing
registerDoFuture()
plan(multicore)

# Read the spectra
lf <- dir(here("data", "spectra"), full.names = TRUE)
spc_list <- read_opus_univ(fnames = lf, extract = c("spc"), parallel = TRUE)
```

Next, we gather the spectra into a list.

``` r
spc_gathered <- spc_list %>%
  gather_spc()
spc_gathered
```

    ## # A tibble: 284 x 6
    ##    unique_id        file_id    sample_id   metadata  spc       wavenumbers
    ##    <chr>            <chr>      <chr>       <list>    <list>    <list>     
    ##  1 BF_lo_01_soil_c… BF_lo_01_… BF_lo_01_s… <tibble … <data.ta… <dbl [1,71…
    ##  2 BF_lo_01_soil_c… BF_lo_01_… BF_lo_01_s… <tibble … <data.ta… <dbl [1,71…
    ##  3 BF_lo_01_soil_c… BF_lo_01_… BF_lo_01_s… <tibble … <data.ta… <dbl [1,71…
    ##  4 BF_lo_02_soil_c… BF_lo_02_… BF_lo_02_s… <tibble … <data.ta… <dbl [1,71…
    ##  5 BF_lo_02_soil_c… BF_lo_02_… BF_lo_02_s… <tibble … <data.ta… <dbl [1,71…
    ##  6 BF_lo_02_soil_c… BF_lo_02_… BF_lo_02_s… <tibble … <data.ta… <dbl [1,71…
    ##  7 BF_lo_03_soil_c… BF_lo_03_… BF_lo_03_s… <tibble … <data.ta… <dbl [1,71…
    ##  8 BF_lo_03_soil_c… BF_lo_03_… BF_lo_03_s… <tibble … <data.ta… <dbl [1,71…
    ##  9 BF_lo_03_soil_c… BF_lo_03_… BF_lo_03_s… <tibble … <data.ta… <dbl [1,71…
    ## 10 BF_lo_04_soil_c… BF_lo_04_… BF_lo_04_s… <tibble … <data.ta… <dbl [1,71…
    ## # ... with 274 more rows

We can now record the time needed to average spectra by `sample_id`
using the current version of `average_spc()`. The following section
realizes the split:

``` r
# Create list of averaged spectra, one spectrum is one data.table
# Method split.data.table is not yet available in data.table 1.9.7
# Wait for v2.0.0
# See https://github.com/philipp-baumann/simplerspec/blob/master/R/average-spc.R
spc_mean_list <- stats::setNames(
  split(spc_mean_noid, seq(nrow(spc_mean_noid))),
    sample_id_mean
)
```

Let’s time the current function using the above implementation:

``` r
system.time({
spc_avg <- spc_gathered %>%
    # Average replicate scans per sample_id
    average_spc(column_in = "spc")
})
```

    ##    user  system elapsed 
    ##   3.164   0.033   3.198

``` r
spc_avg
```

    ## # A tibble: 284 x 7
    ##    unique_id    file_id  sample_id  metadata  spc    wavenumbers spc_mean 
    ##    <chr>        <chr>    <chr>      <list>    <list> <list>      <list>   
    ##  1 BF_lo_01_so… BF_lo_0… BF_lo_01_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  2 BF_lo_01_so… BF_lo_0… BF_lo_01_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  3 BF_lo_01_so… BF_lo_0… BF_lo_01_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  4 BF_lo_02_so… BF_lo_0… BF_lo_02_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  5 BF_lo_02_so… BF_lo_0… BF_lo_02_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  6 BF_lo_02_so… BF_lo_0… BF_lo_02_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  7 BF_lo_03_so… BF_lo_0… BF_lo_03_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  8 BF_lo_03_so… BF_lo_0… BF_lo_03_… <tibble … <data… <dbl [1,71… <data.ta…
    ##  9 BF_lo_03_so… BF_lo_0… BF_lo_03_… <tibble … <data… <dbl [1,71… <data.ta…
    ## 10 BF_lo_04_so… BF_lo_0… BF_lo_04_… <tibble … <data… <dbl [1,71… <data.ta…
    ## # ... with 274 more rows

``` r
# Let's inspect the fist 10 columns of the first mean spectrum
spc_avg$spc_mean[[1]][, 1:10]
```

    ##       3997.4    3995.4   3993.3    3991.3    3989.2   3987.2    3985.2
    ## 1: 0.1238099 0.1232519 0.122713 0.1229396 0.1234718 0.123702 0.1239118
    ##       3983.1    3981.1      3979
    ## 1: 0.1240727 0.1236635 0.1230908

The function below is speed improved because the best performing
implementation `purrr::transpose` among the options tested in the before
mentioned
benchmark:

``` r
average_spc <- function(spc_tbl, by = "sample_id", column_in = "spc_rs") {

  # Avoid R CMD check note: `no visible binding for global variable`
  spc_rs <- sample_id <- id <- NULL

  column_in <- rlang::enquo(column_in)

  # Combine rows of all resampled spectra in one data.table
  spc <- data.table::rbindlist(dplyr::pull(spc_tbl, !!column_in))

  # Add sample_id column to resampled spectra
  spc[, id := spc_tbl[, by][[by]]]

  # Average spectra, use sample_id as index for grouping
  data.table::setkey(spc, id)
  spc_mean <- spc[, lapply(.SD, mean), by = id]

  # Create vector of sample_id from column sample_id in spc_mean
  sample_id_mean <- spc_mean[, id]
  # Delete sample_id column in data.table
  spc_mean_noid <- spc_mean[, id := NULL]

  # Create list of averaged spectra, one spectrum is one data.table
  # Use best performing alternative in
  # https://github.com/jennybc/row-oriented-workflows/blob/master/iterate-over-rows.md
  spc_mean_list <- stats::setNames(
    purrr::transpose(spc_mean_noid),
      sample_id_mean
  )

  # Quote the symbol or the string supplied to by argument
  by <- ensym(by)

  # Convert averaged spectra and sample_id to tibble
  spc_mean_tbl <- tibble::tibble(
    !! by := sample_id_mean, spc_mean = spc_mean_list
  )
  # Join mean spectra tibble spc_tbl_mean to spc_tbl
  spc_tbl <- dplyr::left_join(spc_tbl, spc_mean_tbl, by = quo_name(by))
}
```

We can now test how much faster the new function is:

``` r
system.time({
spc_proc <- spc_gathered %>%
    # Average replicate scans per sample_id
    average_spc(column_in = "spc")
})
```

    ##    user  system elapsed 
    ##   0.482   0.022   0.504

Great\! By simply replacing `base::split()` by `purrr::transpose()` a
more than 6 times speedup is achieved.

The the functions returns the same result within list-column
    `spc_mean`:

``` r
spc_avg$spc_mean[[1]][, 1:10]
```

    ##       3997.4    3995.4   3993.3    3991.3    3989.2   3987.2    3985.2
    ## 1: 0.1238099 0.1232519 0.122713 0.1229396 0.1234718 0.123702 0.1239118
    ##       3983.1    3981.1      3979
    ## 1: 0.1240727 0.1236635 0.1230908

# Test alternative splitting within `simplerspec::preprocess_spc()`

The spectral preprocessing function also splits a single data table by
row. The current version performs as follows:

``` r
system.time({
spc_proc <- spc_avg %>%
    # Preprocess spectra using Savitzky-Golay first derivative with
    # a window size of 21 points
    preprocess_spc(select = "sg_1_w21")
})
```

    ##    user  system elapsed 
    ##   8.436   0.053   8.494

The current version of `preprocess_spc()` splits spectral data to
achieve list-columns using the following lines of
code:

``` r
# Convert preprocessed spectra in data.table to list of data.table spectra
spc_pre_list <- split(spc_pre, seq(nrow(spc_pre)))
```

We replace `split()` by `purrr::transpose()`:

``` r
spc_pre_list <- purrr::transpose(spc_pre)
```

This is the entire function:

``` r
preprocess_spc <- function(spc_tbl, select, column_in = "spc_mean",
  custom_function = NULL) {

  # Convert list of spectral data.tables to one data.table
  spc_raw <- data.table::rbindlist(spc_tbl[column_in][[column_in]])

  ## Perform preprocessing =====================================================

  # Use custom function when supplying option ----------------------------------
  if (!is.null(custom_function) & select == "custom") {
    # Create full character string for parsing
    custom_fct <- paste0("custom <- ", custom_function)
    # parse converts the character string into an expression
    # eval evaluates the expression; as a result, custom object is computed
    # and saved within the current workspace
    eval(parse(text = custom_fct))
    ## x <- spc_raw
    ## custom <- eval(substitute(custom_function), envir = parent.frame())
    # -> Error in is.data.frame(X) : object 'x' not found
  }
  # -> returns error:
  # custom_function = prospectr::savitzkyGolay(X = x, m = 0, p = 3, w = 9)
  # Error in is.data.frame(X) : object 'x' not found
  # -> Maybe solution: http://stackoverflow.com/questions/30563745/non-standard-evaluation-from-another-function-in-r

  # Savitzky-Golay preprocessing
  # use different derivatives and window sizes ---------------------------------

  # Zero order Savitzky-Golay (no derivative) -> only smoothing
  if (select == "sg_0_w9") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)}
  # First derivative Savitzky-Golay
  if (select == "sg_1_w5") {
    sg_1_w5 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 5)}
  if(select == "sg_1_w9") {
    sg_1_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 9)}
  if (select == "sg_1_w11") {
    sg_1_w11 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 11)}
  if (select == "sg_1_w13") {
    sg_1_w13 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 13)}
  if (select == "sg_1_w15") {
    sg_1_w15 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 15)}
  if(select == "sg_1_w17") {
    sg_1_w17 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 17)}
  if (select == "sg_1_w19") {
    sg_1_w19 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 19)}
  # Implement window size of 21, corresponds to ICRAF standard;
  # see e.g. Terhoeven-Urselmans et al. (2010)
  if (select == "sg_1_w21") {
    sg_1_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 21)}
  if (select == "sg_1_w23") {
    sg_1_w23 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 23)}
  if(select == "sg_1_w25") {
    sg_1_w25 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 25)}
  if (select == "sg_1_w27") {
    sg_1_w27 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 27)}
  if (select == "sg_1_w35") {
    sg_1_w35 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 35)}
  if (select == "sg_1_w41") {
    sg_1_w41 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 41)}
  if (select == "sg_1_w51") {
    sg_1_w51 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 51)}
  # Second derivative Savitzky-Golay
  if (select == "sg_2_w11") {
    sg_2_w11 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 11)}
  if (select == "sg_2_w21") {
    sg_2_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 21)}
  # Savitzky-Golay (order 0) smoothing and derivative with a window size of
  # 21 points
  if (select == "sg_0_1_w21") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_1_w21 <- prospectr::savitzkyGolay(X = sg_0_w9,
      m = 1, p = 3, w = 21)}
  # Savitzky-Golay second derivative
  if (select == "sg_2_w5") {
    sg_2_w5 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 5)}
  if (select == "sg_2_w11") {
    sg_2_w11 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 11)}

  # Standard normal variate (SNV) ----------------------------------------------

  # Calculate standard normal variate (SNV) after Savitzky-Golay smoothing
  if (select == "sg_0_snv") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)}
  if (select == "sg_1_snv") {
    sg_1_snv <- prospectr::standardNormalVariate(sg_1_w5)}
  # Standard normal variate (SNV) and first gap-segment derivative
  if (select == "snv_gsd_m1_w11_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w11_s1 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 11, s = 1)}
  if (select == "snv_gsd_m1_w21_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w21_s5 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 21, s = 5)}
  if (select == "snv_gsd_m1_w31_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w31_s1 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 31, s = 5)}
  if (select == "snv_gsd_m1_w31_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w31_s5 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 31, s = 5)}
  # Standard normal variate (SNV) and second gap-segment derivative
  if (select == "snv_gsd_m2_w5_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w5_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 5, s = 1)}
  if (select == "snv_gsd_m2_w21_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w21_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 21, s = 1)}
  if (select == "snv_gsd_m2_w31_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w31_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 31, s = 5)}
  if (select == "snv_gsd_m2_w31_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w31_s5 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 31, s = 1)}
  if (select == "snv_gsd_m2_w51_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w51_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 51, s = 1)}
  if (select == "snv_gsd_m2_w51_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w51_s5 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 51, s = 5)}
  # 1rst Gap-segement derivative
  if (select == "gsd_m1_w5_s4") {
    gsd_m1_w5_s4 <- prospectr::gapDer(X = spc_raw, m = 1, w = 5, s = 4)}
  if (select == "gsd_m1_w11_s5") {
    gsd_m1_w11_s5 <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 5)}
  if (select == "gsd_m1_w11_s21") {
    gsd_m1_w11_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 21)}
  if (select == "gsd_m1_w21_s1") {
    gsd_m1_w21_s1 <- prospectr::gapDer(X = spc_raw, m = 1, w = 21, s = 1)}
  if (select == "gsd_m1_w21_s21") {
    gsd_m1_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 21, s = 21)}
  if (select == "gsd_m1_w35_s21") {
    gsd_m1_w35_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 35, s = 21)}
  if (select == "gsd_m1_w5_s21") {
    gsd_m1_w5_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 5, s = 21)}
  # 2nd Gap-segment derivative
  if (select == "gsd_m2_w21_s21") {
    gsd_m2_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 2, w = 21, s = 21)}
  # 4th Gap-segment derivative
  if (select == "gsd_m4_w21_s21") {
    gsd_m4_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 4, w = 21, s = 21)}

  # Savitzky-Golay combined with multiple scatter correction (MSC --------------
  # Savitzky-Golay with 3rd order polynomial, a window size of 21
  # and first derivative + MSC
  if (select == "sg_1_w21_msc") {
    sg_1_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 21)
    # Use msc function from the pls package; use column means of X as reference
    # spectrum
    sg_1_w21_msc <- pls::msc(X = sg_1_w21, reference = NULL)
  }
  # Savitzky-Golay combined with multiple scatter correction (MSC --------------
  # Savitzky-Golay with 4th order polynomial, a window size of 21
  # and second derivative + MSC
  if (select == "sg_2_w21_msc") {
    sg_2_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 4, w = 21)
    # Use msc function from the pls package; use column means of X as reference
    # spectrum
    sg_2_w21_msc <- pls::msc(X = sg_2_w21, reference = NULL)
  }

  # Continuum-removal ----------------------------------------------------------
  if (select == "cr") {
    cr <- prospectr::continuumRemoval(X = spc_raw,
      wav = as.numeric(colnames(spc_raw)), type = "A")}
  if (select == "cr_refl") {
    cr_refl <- prospectr::continuumRemoval(X = spc_raw,
      wav = as.numeric(colnames(spc_raw)), type = "R")}

  # Select final preprocessing based on selection argument and
  # save matrix in data.table
  pre <- select
  spc_pre <- data.table::as.data.table(get(pre))

  # Convert preprocessed spectra in data.table to list of data.table spectra
  spc_pre_list <- purrr::transpose(spc_pre)
  # Convert x-values of preprocessed spectra in list of vectors
  # prospectr only hands over new xunits in matrix colnames of type character
  xvalues_pre_list <- lapply(spc_pre_list,
    function(x) as.numeric(colnames(x)))

  # Add list of preprocessed spectra and correspoding wavenumbers to tibble
  spc_tbl <- tibble::add_column(spc_tbl,
    spc_pre = spc_pre_list, xvalues_pre = xvalues_pre_list)
}
```

``` r
system.time({
spc_proc <- spc_avg %>%
    # Preprocess spectra using Savitzky-Golay first derivative with
    # a window size of 21 points
    preprocess_spc(select = "sg_1_w21")
})
```

    ##    user  system elapsed 
    ##   0.335   0.024   0.360

With this, we achieve an about 16 times speedup :-)

Finally, we can test the speed improvement for the entire pipeline:

``` r
system.time({
spc_proc <- spc_gathered %>%
    # Average replicate scans per sample_id
    simplerspec::average_spc(column_in = "spc") %>%
    # Preprocess spectra using Savitzky-Golay first derivative with
    # a window size of 21 points
    simplerspec::preprocess_spc(select = "sg_1_w21")
})
```

    ##    user  system elapsed 
    ##  12.015   0.086  12.122

vs.

``` r
system.time({
spc_proc <- spc_gathered %>%
    # Average replicate scans per sample_id
    average_spc(column_in = "spc") %>%
    # Preprocess spectra using Savitzky-Golay first derivative with
    # a window size of 21 points
    preprocess_spc(select = "sg_1_w21")
})
```

    ##    user  system elapsed 
    ##   0.564   0.026   0.591
