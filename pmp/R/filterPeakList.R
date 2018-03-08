mzRange <- function(ppm, mass)
{
  delta <- ppm*mass*10^-6
  out <- c(mass-delta,mass+delta)
  out
}

filterPeakList <- function (Data, ppm, masses)
{
  mhits <- NULL
  for (mm in masses)
  {
    Range <-mzRange(ppm, mm)
    hh <- which (Data>=Range[1] & Data<=Range[2])
    if (length(hh)>0)
    {
      mhits <- append(mhits, hh)
    }
  }
  unique(mhits)
}
