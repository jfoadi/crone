
## Resubmission
* First resubmission: changed "description" field in DESCRIPTION file; ignored debian pre-test as there is no 'test_x.dat' file.

* Second resubmission: test_h.dat and test_x.dat actually were a product of examples execution. I have eliminated the problem by making the examples NotRun.

* Third resubmission: changed the description in file DESCRIPTION; it's now more informative. Turned all F vectors into FF. Replaced some "print" and "cat" with "message". Changed printing files in examples from "homedir" to tempdir() and made the examples runnable. Re-checked vignettes.

## Test environments
* Windows 10 - Using RStudio 1.2.1335
* Linux - Using TRAVIS (installed biber and ghostscript)
* devtools::check_rhub - failed Ubuntu platform just for vignettes because no biber found

## R CMD check results
Carried out from within RStudio.
No ERRORS.
No WARNINGS.
No NOTES.

## Downstream dependencies
Not applicable. First release.