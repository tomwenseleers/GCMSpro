

####---- DESCRIPTION ----####


# This script is template/example script for using our GC-MS processing
# pipeline.

####---- SETUP ENVIRONMENT ----####


# Clear workspace
rm(list = ls())
gc()

# Setup new workspace
# base.dir indicates the folders that contains the sample metadata table and
# eventually the target information table. Results and plots will be saved there too.
base.dir <- "C:/Users/Ento/Documents/Chris/GCMSpro/190815 Testing the pipeline/"
base.dir <- "D:/Documents/Dropbox/christophe/_tests deconvolution scripts/190815 Testing the pipeline/"
base.dir <- "H:/" # From USB stick
setwd(base.dir)
library(GCMSpro)

####---- INITIALIZE A NEW ANALYSIS ----####


# TODO make a function that will create a compatible compound file from supplied spectrum + RI libraries (get MSP functions)

# We setup some general parameters used throughout the pipeline. This
# comprises the plot and output directories, and the sample and target
# files.
params <- initiate.parameters(base.dir = base.dir,
                              sample.file = "beer sample table - sub.csv",
                              compound.file = "beer target table.csv",
                              data.dir = "D:/Documents/temp/",
                              verbose = TRUE,
                              Plot = TRUE,
                              save.backup = TRUE,
                              ncores = 20, # detectCores(),
                              untargeted.samples = c("20170201/5425006700054.cdf",
                                                     "20170129/5411081004736_170129234546.cdf",
                                                     "20170202/54004009_170203145246.cdf",
                                                     "20170523/5425026610197.cdf",
                                                     "20170213/5411858000145.cdf"),
                              # Calibration specific parameters
                              ref.samp = 1, # can be numeric (index in sample table) or character (path to sample)
                              mz.ignore = NULL,
                              lib.type = "alkanes",
                              # Preprocessing parameters
                              data.filename = "data",
                              enc = "integer", # or "double" but requires more RAM
                              bas.tau = 0.2,
                              bas.nknots = 15, # Number of knots for the B-spline basis
                              bas.degree = 3, # degree of the piecewise polynomial
                              bas.subsamp = 10, # subsampling factor when fitting the baseline
                              # Deconvolution parameters
                              win.size = 300,
                              peak.width.range = c(1,20),
                              peak.trim = 1E-10,
                              cutoff = 1,
                              lives = 2,
                              stop.thres = 0.01,
                              tau = 0.1,
                              cobs.tau = 0.1,
                              eps = 1,
                              maxiter = 10,
                              tol = 0.01,
                              basfun.u = function(x) return(min(x)),
                              sel.power = 2,
                              t.power = 2,
                              win.ext = 25,
                              scale.subsamp = 10,
                              method.upd = "nnL0pois",
                              max.RI.sd = 50,
                              elu.win = 60)


####---- CALIBRATE THE DATA ----####


# The calibration step is meant to extract the RT information of some known
# standards. The extracted information is used to build a model that maps RT
# taken on different days and/or in different technical batches to a theoretical
# RI scale. Another model is build to match this theoretical RI to a reference
# RT scale, given the reference sample supplied as 'ref.samp'. Some mz traces
# can be ignored by specifying them in 'mz.ignore'. Parallelization is possible
# and stringly advised when possible.
cal.res <- RT.calibration(params = params) # 51.8 min


####---- PREPROCESS DATA ----####


# The preprocessing step will combine all samples in a single matrix (row-
# augmentation) stored in a disk-backed file. Every sample is individually
# processed by
preprocess.res <- preprocessing(params = params, calibration.output = cal.res)


####---- DECONVOLUTION ----####


# The deconvolution step is the core step of the algorithm, it will extract from
# the preprocessed samples the elution and spectrum profiles. The deconvolution
# can occur in three steps:
#   - Untargeted deconvolution - elution profiles and spectrum profiles are
#       extracted with no prior information. This can be computionally heavy, so
#       maybe perform only on max 50 samples.
#   - Targeted deconvolution - elution profiles are extracted based on given
#       spectrum profiles. This is much faster so can be used after extracting
#       the untargeted profiles in a few samples to extract the profiles in the
#       remaining samples.
deconvolution.res <- deconvolution(params = params,
                                   calibration.output = cal.res,
                                   preprocess.output = preprocess.res)


####---- COMPOUND IDENTIFICATION ----####


# TODO make a function that takes a matrix of spectrum profiles + targets.df and returns compound names


