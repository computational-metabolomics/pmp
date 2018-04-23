# pmp
Peak Matrix Processing

Package: pmp

Type: Package

Title: Peak matrix processing

Version: 0.1.0

Author: Ralf Weber and Andris Jankevics

Maintainer: Andris Jankevics <a.jankevics@bham.ac.uk>

Description: Tools and filters for peak matrix scaling, normalisation and filtering.

License: GPL

Imports: impute, pcaMethods, missForest

Encoding: UTF-8

LazyData: true

RoxygenNote: 6.0.1

## Installation instructions

Install Bioconductor dependencies and devtools package:

`source("https://bioconductor.org/biocLite.R")`
`biocLite("impute")`
`biocLite("pcaMethods")`

`install.packages ("devtools")`

Obtain Github authorisation token fro private repositories: https://github.com/settings/tokens

Use this command from R:

`devtools::install_github("computational-metabolomics/pmp/pmp", auth_token = "YOUR_GITHUB_TOKEN")`
