## Test environments
* local R installation, Windows 10 x64, R 3.6.2

## R CMD check results

0 errors | 0 warnings | 0 note

## Mendings

* We add output information in all .Rd files.

* We use on.exit to make sure we do not change users' par in `plot.bsrr`.
