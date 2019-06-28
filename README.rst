==============================================
PMP: Peak Matrix Processing
==============================================

|Git| |Bioconda| |Build Status (Travis)| |Build Status (AppVeyor)| |License| |Coverage|


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


.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/pmp/master.svg?label=Travis
   :target: https://travis-ci.org/computational-metabolomics/pmp

.. |Build Status (AppVeyor)| image:: https://ci.appveyor.com/api/projects/status/github/computational-metabolomics/pmp?branch=master&svg=true
   :target: https://ci.appveyor.com/project/computational-metabolomics/pmp

.. |Build Status (AppVeyor)| image:: https://ci.appveyor.com/api/projects/status/github/computational-metabolomics/pmp?branch=master&svg=true
   :target: https://ci.appveyor.com/project/computational-metabolomcis/pmp

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/pmp

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: https://bioconda.github.io/recipes/bioconductor-pmp/README.html

.. |License| image:: https://img.shields.io/badge/licence-GNU_v3-teal.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html
   
.. |Coverage| image:: https://codecov.io/gh/computational-metabolomics/pmp/branch/master/graph/badge.svg
   :target: https://codecov.io/github/computational-metabolomics/pmp?branch=master
