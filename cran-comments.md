
## Resubmission
* First resubmission: changed "description" field in DESCRIPTION file; ignored debian pre-test as there is no 'test_x.dat' file.

* Second resubmission: test_h.dat and test_x.dat actually were a product of examples execution. I have eliminated the problem by making the examples NotRun.

* Third resubmission: changed the description in file DESCRIPTION; it's now more informative. Turned all F vectors into FF. Replaced some "print" and "cat" with "message". Changed printing files in examples from "homedir" to tempdir() and made the examples runnable. Re-checked vignettes.

* Fourth resubmission: added a reference in DESCRIPTION and README.

* Fifth resubmission: hopefully the DOI format in DESCRIPTION is now fine.

* Sixth resubmission: rechanged DOI format as <doi:10..1088/1361-6404/aa8188>

* Seventh resubmission: DESCRIPTION is fine. Link in README has been changed to [link]() format

* eight resubmission (now version 0.1.1): drastically changed references as "biber" was causing problems to CRAN automated building. Now references are inserted manually as there's just a handful of them.

## Test environments
* Windows 10 - Using RStudio 1.2.1335
* Linux - Using TRAVIS (installed ghostscript)
* devtools::check_rhub - Only Windows complain about missing ghostscript (not something I can control, I believe)

## R CMD check results
Carried out from within RStudio.
No ERRORS.
No WARNINGS.
No NOTES.

## Downstream dependencies
Not applicable. First release.