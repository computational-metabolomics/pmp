==============================================
PMP: Peak Matrix Processing
==============================================

|Git| |Bioconda| |Build Status| |License| |Coverage|


------------
Install
------------

Github
------------

.. code-block:: r

  library(devtools)
  library(testthat)
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("impute")
  BiocManager::install("pcaMethods")
  install_github('computational-metabolomics/pmp')
 
Conda
------------

.. code-block:: command

   conda create -n pmp pmp -c conda-forge -c bioconda -c computational-metabolomics
   source activate pmp

------------
References
------------


.. |Build Status| image:: https://github.com/computational-metabolomics/pmp/workflows/pmp/badge.svg
   :target: https://github.com/computational-metabolomics/pmp/actions

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/pmp

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: https://bioconda.github.io/recipes/bioconductor-pmp/README.html

.. |License| image:: https://img.shields.io/badge/licence-GNU_v3-teal.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html
   
.. |Coverage| image:: https://codecov.io/gh/computational-metabolomics/pmp/branch/master/graph/badge.svg
   :target: https://codecov.io/github/computational-metabolomics/pmp?branch=master
