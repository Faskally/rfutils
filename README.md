
# rfutils

## to install

``` r
remotes::install_github("faskally/rfutils")
```

## contributing

After making changes to the code, run:

``` r
devtools::document()
```

to update help files and the NAMESPACE

then

``` r
devtools::check()
```

to check for errors, and fix warnings and notes

Finally, the readme can be rebuilt using

``` r
rmarkdown::render("README.Rmd")
```

then changes can be commited and pushed.
