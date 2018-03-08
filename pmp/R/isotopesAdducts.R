isotopesAdducts <- function (formula, adducts)
{
  require(BRAIN)
  res <- useBRAIN(aC=formula)
  masses <- res$masses
  mDifs <- read.csv(adducts)
  mDifs <- mDifs$massdiff
  out <- NULL
  for (m in masses)
  {
    out <- append(out,m+mDifs)
  }
  out
}