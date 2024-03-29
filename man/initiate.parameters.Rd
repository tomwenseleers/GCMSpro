% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-0.0.2.R
\name{initiate.parameters}
\alias{initiate.parameters}
\title{Initiate parameters for a new GCMS data processing run}
\usage{
initiate.parameters(base.dir, Plot = TRUE, verbose = TRUE,
  save.backup = TRUE, ncores = parallel::detectCores(), sample.file,
  compound.file, ref.samp = 1, mz.ignore = NULL,
  lib.type = "alkanes", data.filename = "data", enc = "integer",
  data.dir = NULL, bas.tau = 0.2, bas.nknots = 15, bas.degree = 3,
  bas.subsamp = 10, untargeted = TRUE, win.size = 300,
  peak.width.range = c(1, 20), peak.trim = 1e-10, cutoff = 1,
  lives = 2, stop.thres = 2e-04, tau = 0.1, cobs.tau = 0.1,
  eps = 1, maxiter = 10, tol = 0.01, basfun.u = function(x)
  return(min(x)), sel.power = 2, t.power = 2, win.ext = 25,
  scale.subsamp = 10, method.upd = "nnL0pois", max.RI.sd = 50,
  elu.win = 60)
}
\arguments{
\item{elu.win}{}
}
\value{
TODO a list
}
\description{
TODO
}
