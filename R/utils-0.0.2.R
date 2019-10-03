# This script contains all functions for the "Development of
# a state-of-the-art pipeline for high throughput analysis
# of gas chromatography - mass spectrometry data"

# Version 2

# This script is under development

# Some required packages

library(bigmemory)
library(cobs)
library(doSNOW)
library(export)
library(ggplot2)
library(httr) # for building target library
library(quantreg)
library(L0glm)
library(Matrix)
library(nnls)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
# library(RMassBank) # for building the target table
library(rvest) # for building the target table
library(snow)
library(splines) # TODO needed?
library(tidyr) # for formating the output table
library(xcms) # TODO only used for reading in cdf files, maybe extract useful code in order to avoid all the associated dependencies
library(xml2) # for building the target table


####---- GENERAL FUNCTIONS ----####

# print.progress ####
# Function that prints the progress in a loop
# INPUT
#   current:  the current iteration
#   total:    the total number of iteration in the loop
#   ...:      additional arguments passed to format() for decimal formatting
# OUTPUT
#   void
print.progress <- function(current, total, ...){
  if(length(current)!=1 && length(total) != 1) stop("Arguments \'current\' and \'total\' must have length 1.")
  perc <- format(current/total*100, digits = 0,...)
  cat(paste0("\rProgress: ", perc, " %   "))
  if(current == total) cat ("\n")
  flush.console()
}


# get.TODOs ####
# Function to print the TODO's of a file
# INPUT
#   file:  chararacter string containing the path to the file
# OUTPUT
#   void
get.TODOs <- function(file){
  file.txt <- readLines(file)
  todo.lines <- which(grepl("TODO", file.txt))
  todos <- file.txt[todo.lines]
  todos <- sapply(todos, function(x) gsub("^.+TODO", "", x))
  for(i in 1:length(todos)){
    cat(paste0("TODO ", i, " at line ", todo.lines[i],":", todos[i], "\n"))
  }
}


# get.TOC ####
# Function to print the table of content of a file
# INPUT
#   file:  chararacter string containing the path to the file
# OUTPUT
#   void
get.TOC <- function(file){
  file.txt <- readLines(file)
  headers.lines <- which(grepl("####", file.txt))
  headers <- file.txt[headers.lines]
  headers <- sapply(headers, function(x) gsub("####.*$", "", x))
  cat(paste0(trimws(headers), "\n"))
}


# plotmat ####
# Function to plot a chromatogram with sensible parameters.
# INPUT
#   M:  the chromatogram given as a matrix
#   subsamp:  subsampling factor for the rows. A higher value will decrease
#             resolution but improve plotting speed
# gamma:    the gamma correction factor
#   xlab, ylab:   the labels for the x and y axis
#   plot.axes:    a logical indicating if axes should be plotted
#   ...:  additional arguments passed to image()
# OUTPUT
#   void
# TODO clean this, allow for user supplied scantimes (from several samples?)
plotmat <- function(M, subsamp = 10, gamma = 0.12,
                    xlab = "Retention time (mins)",
                    ylab = "Mass (m/z)",
                    plot.axes = TRUE, ...) {
  colpalette <- colorRampPalette(c("black", "black", "blue","green","yellow","orange", "red"))
  M <- M[seq(1,nrow(M),subsamp),]
  M <- as.matrix(M) # convert sparse to dense matrix
  M[is.na(M)] <- 0
  neg <- M<0
  if(any(neg)) warning("Negative values in matrix 'M'.")
  M[neg] <- 0
  image(x = 1:nrow(M), y = 1:ncol(M), z = M^gamma,
        col = colpalette(1000), useRaster = T,
        xlab = xlab, ylab = ylab, axes = FALSE, ...)
  if(plot.axes){
    if (length(rownames(M)) != 0) {
      xvals <- as.numeric(rownames(M)) # TODO adapt in case no numeric (or not monotonic)
    } else {
      xvals <- 1:nrow(M)
    }
    xrange <- range(xvals)
    xlabls <- as.character(round(xvals))
    xsel <- (as.numeric(xlabls) %% 10) == 0 & c(TRUE, diff(as.numeric(xlabls)) != 0)
    ylabls <- as.character(round(as.numeric(colnames(M))))
    ysel <- (as.numeric(ylabls) %% 100) == 0 & c(TRUE,diff(as.numeric(ylabls))!=0)
    box()
    axis(1, at = (1:nrow(M))[xsel], labels = xlabls[xsel], tick = TRUE, tcl = -0.4)
    axis(2, at = (1:ncol(M))[ysel], labels = ylabls[ysel], tick = TRUE, las = 1, tcl = -0.4)
  }
}

# read.cdf ####
# Function to read in cdf file as matrix
# INPUT
#   file:   string giving the path of the cdf file to load
#   mz.ignore:  specify the index of the ion traces to ignore
# OUTPUT
#   M:  the GC-MS data contained in file formated as a "matrix" object
read.cdf <- function(file, mz.ignore = NULL) {
  f <- suppressMessages(xcmsRaw(file, profstep = 1, profmethod="intlin")) # set profstep to 0 to read in the raw data
  mzvals <- f@mzrange[1]:f@mzrange[2]
  scantimes <- f@scantime/60 # we return scantimes in mins
  M <- t(f@env$profile) # data matrix
  M[M<0] <- 0
  colnames(M) <- mzvals
  rownames(M) <- scantimes
  if(!is.null(mz.ignore)){
    M <- M[,-which(colnames(M) %in% mz.ignore)] # Remove these ion traces
  }
  M[is.na(M)] <- 0
  return(M)
}


# check.for.backups ####
# TODO doc
check.for.backups <- function(dir, pattern = "[.]rds$"){
  bckup <- grep(list.files(path = dir), pattern = pattern, value = TRUE)
  if(length(bckup) > 0){
    cat("Backup file was found.\nDo you want to restart the run ('yes') or load the backup ('no')?")
    resp <- ""
    while(!resp %in% c("yes", "no")){
      resp <- readline(prompt = "Restart run? [yes] or [no]: ")
    }
    if(resp == "yes"){
      return(NULL)
    } else {
      if(length(bckup) > 1) bckup <- bckup[order(bckup, decreasing = TRUE)][1]
      return(readRDS(paste0(dir, bckup)))
    }
  } else {
    return(NULL)
  }
}


# cosine.sim ####
# Function that compute the cosine similarity between to matrices. The cosine
# similarity = cos(theta) = X . Y / (||X||.||Y||)
# INPUT
#   X:  a matrix of size m x n. X can be a vector if n = 1
#   Y:  a matrix of size m x p. Y can be a vector if p = 1
# OUTPUT
#   cs: the cosine similarity
# REFERENCE
# https://en.wikipedia.org/wiki/Cosine_similarity#Definition
cosine.sim <- function(X, Y=NULL){
  if(is.null(Y)) Y <- X
  # if(is.matrix(X)) warning("Find a more efficient way to compute the norms")
  cp <- crossprod(X, Y)
  norm <- diag(sqrt(crossprod(X))) %*% t(diag(sqrt(crossprod(Y)))) # TODO this is not efficient
  cs <- cp/norm
  return(cs) # similarity matrix, columns correspond to Y,, rows correspond to X
}


# normS ####
# TODO continue documentation
# utility function to normalize spectra to have base peak intensity of 1 in spectra arranged in
# columns (by=2) or rows (by=1) of a matrix
# ntype = 0: no normalization
# ntype = 1: all peaks sum to one
# ntype = 2: l2 norm, the quadratic sums to one
# ntype = 3: peaks are scaled so that the highest peak equals one.
normS <- function(S, by=2, ntype = 1){
  if(ntype == 0){
    return(S)
  } else if (ntype == 1){
    return(sweep(S, by, apply(S,by,sum), FUN="/"))
  } else if (ntype == 2){
    return(sweep(S, by, apply(S,by,function(x){sqrt(sum(x^2))}), FUN="/"))
  } else if (ntype == 3){
    return(sweep(S, by, apply(S,by,max), FUN="/"))  # renormalize spectra to have base peak intensity of 1
  }
}


####---- EXPORTED FUNCTIONS ----####


#' Initiate parameters for a new GCMS data processing run
#'
#' The function generates a list with all required arguments for controlling the
#' processing of GCMS data. Input checks on the files are performed to ensure that the
#' supplied information complies to the expected input.
#'
#' @param base.dir
#' a character string that indicates the path to the folder were to run the
#' GCMS data processing.
#' @param Plot
#' should intermediate plots be generated ? Watch out the files containing the
#' plots will take up some disk space. Make sure to have sufficients disk space.
#' #' @param verbose
#' should progress information be printed to the console ?
#' @param save.backup
#' should a backup be saved after every pipeline step ? Watch out the backup
#' files will take up some disk space. Make sure to have sufficients disk space.
#' @param ncores
#' number of cores to use for computation (using the \code{\link{snow}} package)
#' of parallelized steps.
#' @param sample.file
#' the name or path to the CSV file containing the required sample information.
#' The file should contains at least the following headers: \code{name},
#' \code{file}, \code{type}, \code{date}, \code{batch}.
#' @param compound.file
#' the name or path to the CSV file containing the candidate compounds for
#' identification. The file should contains at least the following headers:
#' \code{compound}, \code{RI}, \code{spectrum}.
#' @param ref.samp
#' the reference sample used to extract the reference time scale and m/z value
#' scale. The sample can be given as the index in \code{compound.file} or as
#' the file path to the CSV file. By default, the first sample in the supplied
#' table is chosen.
#' @param mz.ignore
#' a vector contraining the m/z values to ignore during the processing.
#' @param std.lib
#' a string indicating which library of standard compounds should be used for
#' calibration. Currently, only \code{"alkanes"} for linear alkane calibration
#' is supplied.
#' @param data.filename
#' the name of the file containing the combined data. It is advised to used a
#' different name for different runs/experiments in case the same data should be
#' analyzed with different settings.
#' @param enc
#' the type of encoding of the data. \code{integer} or \code{double}. While the
#' former takes up less storage spess, the second allows for decimal values.
#' @param data.dir
#' the path where the data should be stored. If \code{NULL}, it is saved in a
#' subdirectory called \code{disk backed data/}.
#'
#'
#' \strong{Advanced parameters.}
#'
#'
#' @param tau
#' The \code{tau} parameter for quantile regression. This serves as an
#' underapproximation parameter for fitting the signal. Should be greater
#' than 0 (complete underapproximation) and smaller than 1 (complete
#' overapproximation). \code{tau == 0.5} boils down to median regression. A
#' value of 0.2 seems to work best in most cases.
#' @param bas.tau
#' Same as \code{tau} but for fitting the baseline spectrum profiles.
#' @param bas.nknots
#' the number of knots (~ flexion points) to use to build the baseline elution
#' profiles.
#' @param bas.degree
#' the degree of the splines used to construct the baseline
#' @param bas.subsamp
#' a subsampling factor >=1 indicating by how much the data should be sampled
#' for fitting the baseline. This values means taking one scanline every
#' \code{bas.subsamp} scanlines. When \code{bas.subsamp == 1}, no subsapling is
#' performed. This allows dramatic increase in speed when preprocessing the
#' data.
#' @param untargeted.samples
#' a vector of sample names or sample indices for which an untargeted
#' deconvolution should be performed.
#' @param win.size
#' the size of the window (in number of scanlines) used for optimizing the
#' elution peak width during the untargeted deconvolution initialization.
#' @param peak.width.range
#' a range wihtin which the peak width is optimized. The peak width is the
#' number of scanlines between the inflection points#'
#' @param peak.trim
#' a threshold under which values of the initial Gaussian elution profiles are
#' clipped to 0.
#' @param cutoff
#' a threshold that filters out peak for which the ratio between peak height and
#' the baseline is smaller or equal than this threshold.
#' @param lives
#' the maximal number of times a initialized peak can be extracted during the
#' untargeted deconvolution
#' @param stop.thres
#' a threshold that stops the untargeted deconvolution once the highest
#' remaining peaks is not higher than the threshold times the height of the
#' highest initial peak.
#' @param cobs.tau
#' the tau parameters for the unimodalite constraints (see \code{\link{cobs}})
#' @param eps
#' a threshold under which measured and extracted itensities are trimmed to 0
#' @param maxiter
#' the maximum number of iterations when optimizing the spectrum profiles
#' within the untargeted deconvolution
#' @param tol
#' a threshold that stops the spectrum optimization loop once the difference in
#' spectra between two consecutive iteration is lower than that threshold.
#' @param basfun.u
#' a function to estimate a baseline given an extracted elution profiles. The
#' output of the function is then subtracted from the elution profiles.
#' @param sel.power
#' a factor that influences the weights on selective ions when estimating the
#' elution profiles from the spectrum profiles. The higher the factor, the more
#' weights on selective ions.
#' @param t.power
#' a factor that influences the weights when estimating the spectrum profiles
#' from the elution profiles. The higher the factor, the more weight is put on
#' high itensity values.
#' @param win.ext
#' the number of scanlines to extend the search window on each side when
#' extracting profiles in all samples given the optimization in the sample
#' with the highest peak height sample.
#' @param scale.subsamp
#' a subsampling factor >=1 indicating by how much the data should be sampled
#' when scaling the extracted elution profiles using the median regression.
#' @param method.upd
#' the method used to optimize simultaneously the spectrum profiles at the end
#' of the untargeted deconvolution algorithm. Available methods are weighted
#' nonnegative least squares regression (\code{WNNLS}) or nonnegative L0 poisson
#' regression (\code{nnL0pois}).
#'
#' @return
#'
#' a list containing the required parameters for running the next steps of the
#' GC-MS processing pipeline
#'
#' @export
#'
initiate.parameters <- function(base.dir, Plot = TRUE, verbose = TRUE,
                                save.backup = TRUE, ncores = parallel::detectCores(),
                                # Sample information
                                sample.file,
                                # Compound information
                                compound.file,
                                # Calibration specific parameters
                                # TODO allow that a reference sample is not supplied but that the users gives all the reference parameters: scantime, mzvalue, day, batch,...
                                ref.samp = 1, # can be numeric (index in sample table) or character (path to sample)
                                mz.ignore = NULL,
                                std.lib = "alkanes",
                                # Preprocessing parameters
                                data.filename = "data",
                                enc = "integer", # or "double" but requires more RAM
                                data.dir = NULL,
                                bas.tau = 0.2,
                                bas.nknots = 15, # Number of knots for the B-spline basis
                                bas.degree = 3, # degree of the piecewise polynomial
                                bas.subsamp = 10, # subsampling factor when fitting the baseline
                                # Deconvolution parameters
                                untargeted.samples = "all",
                                win.size = 300,
                                peak.width.range = c(1,20),
                                peak.trim = 1E-10,
                                cutoff = 1,
                                lives = 2,
                                stop.thres = 2E-4,
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
                                method.upd = "nnL0pois"){


  if(verbose) cat("====== Initializing a new GC-MS processing run ======\n\n")
  # Check the base directory
  if(!dir.exists(base.dir)) stop("Could not find the supplied 'base.dir'. This should be a string giving the path to an existing directory.")
  if(tail(strsplit(base.dir, "")[[1]], 1) != "/") base.dir <- paste0(base.dir, "/") # Make sure that the last character of the path is "/"
  # Create an output directory
  output.dir <- paste0(base.dir, "Output/")
  if(!dir.exists(output.dir) & !dir.create(output.dir, showWarnings = FALSE)) stop(paste0("Could not create the output directory: ", output.dir))
  # Create a plot directory
  if(Plot){
    plot.dir <- paste0(base.dir, "Plots/")
    if(!dir.exists(plot.dir) & !dir.create(plot.dir, showWarnings = FALSE)) stop(paste0("Could not create the plot directory: ", plot.dir))
  } else {
    plot.dir <- NULL
  }
  # Load the sample file
  if(verbose) cat("Loading the sample file...")
  if(!grepl(pattern = base.dir, x = sample.file)) sample.file <- paste0(base.dir, sample.file)
  if(!grepl(pattern = "[.]csv$", x = sample.file)) stop("Sample file is expected to be a CSV file.")
  samples.df <- read.csv(sample.file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  samples.df <- check.sample.table(samples.df)
  standards.df <- samples.df[samples.df$type == "standard",]
  samples.df <- samples.df[samples.df$type != "standard",]

  if(verbose) cat(" loaded without error.\n")
  # Load the target table
  if(!is.null(compound.file)){
    if(verbose) cat("Loading the compound file...")
    if(!grepl(pattern = base.dir, x = compound.file)) compound.file <- paste0(base.dir, compound.file)
    if(!grepl(pattern = "[.]csv$", x = compound.file)) stop("The compound file is expected to be a CSV file.")
    compounds.df <- read.csv(compound.file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
    compounds.df <- check.target.table(compounds.df)
    if(verbose) cat(" loaded without error.\n")
  } else {
    compounds.df <- NULL
  }
  # Check the reference sample
  if(is.numeric(ref.samp)) ref.samp <- samples.df$file[ref.samp]
  if(!file.exists(ref.samp)) stop(paste0("Could not find reference sample file: ", ref.samp))
  # Check the sample to undergo untargeted deconvolution
  all.samples <- (1:nrow(samples.df))[samples.df$type != "standard"]
  if(identical(untargeted.samples, "all")){
    untargeted.samples <- all.samples
  } else if(is.character(untargeted.samples)){
    untargeted.samples <- which(samples.df$name %in% untargeted.samples)
    if(length(untargeted.samples) == 0) stop("The supplied sample names for untargeted deconvolution were not found in the sample file.")
  } else {
    stop("Untargeted samples should be a vector of strings containing the sample names (at least one).")
  }
  if(length(untargeted.samples) == 0) stop("At least one sample must be processed using untargeted deconvolution. Supply 'untargeted.samples' as sample names (contained in the sample table).")
  if(max(untargeted.samples) > nrow(samples.df)) stop("Sample index exceeds number of rows of the sample table.")
  # TODO this will have to change with FAME standards
  if(any(samples.df$type[untargeted.samples] == "standard")) stop("You cannot supply samples with standards for untargeted dconvolution. This will automatically be performed in the RT calibration.")
  # Get the targeted samples (the remaining samples)
  targeted.samples <- all.samples[!all.samples %in% untargeted.samples]

  if(verbose) cat(paste0("\nBase directory: ", base.dir, "\n",
                         "Number of samples: ", nrow(samples.df), "\n",
                         "Number of library compounds: ", ifelse(is.null(compounds.df), 0, nrow(compounds.df)), "\n"))
  # Return the variable
  return(list(base.dir = base.dir,  output.dir = output.dir, plot.dir = plot.dir,
              samples.df = samples.df, standards.df = standards.df,  compounds.df = compounds.df,
              verbose = verbose, save.backup = save.backup, ncores = ncores,
              # RT calibration
              ref.samp = ref.samp,  mz.ignore = mz.ignore,  std.lib = std.lib,
              # Preprocessing parameters
              data.filename = data.filename, enc = enc, data.dir = data.dir,
              bas.tau = bas.tau, bas.nknots = bas.nknots, bas.degree = bas.degree,
              bas.subsamp = bas.subsamp,
              # Deconvolution parameters
              untargeted.samples = untargeted.samples, targeted.samples = targeted.samples,
              win.size = win.size, peak.width.range = peak.width.range,
              peak.trim = peak.trim, cutoff = cutoff, lives = lives, stop.thres = stop.thres,
              tau = tau, cobs.tau = cobs.tau, eps = eps, maxiter = maxiter, tol = tol,
              basfun.u = basfun.u, sel.power = sel.power, t.power = t.power,
              win.ext = win.ext, scale.subsamp = scale.subsamp,  method.upd = method.upd))
}


# RT.calibration ####
#' Performs retention time calibration from the standard samples
#'
#' The function extract the RT/RI information from the standard samples and
#' builds the RT calibration model given a reference sample.
#'
#' @param params
#' the ouptut of the \code{\link{initiate.parameters}} function.
#'
#' @return
#'
#' a list with the calibration results
#'
#' @export
#'
RT.calibration <- function(params){
  if(params$verbose) cat("====== RT extraction of calibration compounds ======\n\n")
  # Search for backups
  backup <- check.for.backups(dir = params$output.dir, pattern = "1-RT calibration")
  if(!is.null(backup)) return(backup)
  # Start timing
  t1 <- proc.time()[3]

  # Subfolder for plotting
  if(!is.null(params$plot.dir)){
    plot.subdir <- paste0(params$plot.dir, "1-Sample alignment/")
    dir.create(plot.subdir, showWarnings = FALSE)
  }

  if(params$verbose) cat(paste0("Number of samples: ", nrow(params$standards.df), "\n",
                                "Number of cores: ", params$ncores, "\n",
                                "Calibration compounds: ", params$std.lib, "\n\n"))

  # Extract RT information from standards
  if(params$verbose) cat("Extracting the RT information from the standard samples...")
  RT.extraction <- extract.RTs(standards.df = params$standards.df, ref.samp = params$ref.samp,
                               std.lib = params$std.lib, ncores = params$ncores,
                               mz.ignore = params$mz.ignore, plot.folder = plot.subdir,
                               # Baseline fitting
                               # TODO get the parameters out of params?
                               bas.subsamp = params$bas.subsamp, bas.tau = params$bas.tau,
                               bas.nknots = params$bas.nknots, bas.degree = params$bas.degree,
                               # Peak identification
                               win.size = params$win.size, peak.width.range = params$peak.width.range,
                               peak.trim = params$peak.trim, cutoff = params$cutoff,
                               # Deconvolution
                               tau = params$tau, cobs.tau = params$cobs.tau, eps = params$eps,
                               maxiter = params$maxiter, tol = params$tol, basfun.u = params$basfun.u,
                               sel.power = params$sel.power, t.power = params$t.power, win.ext = params$win.ext,
                               scale.subsamp = params$scale.subsamp)
  if(params$verbose) cat("done\n\n")
  RI.df <- RT.extraction$RI.df
  scantimes.ref <- RT.extraction$scantimes.ref
  mzvals.ref <- RT.extraction$mzvals.ref

  # Generate the RI to RT model based on the extracted data
  # TODO constraint fit to nonnegative predictions
  # TODO adjust dof automatically
  if(params$verbose) cat("Build the RT/RI calibration models\n\n")
  RT.to.RI.fit <- get.RT.to.RI.fit(alk.df = RI.df, rm.outliers = TRUE,
                                   formul = NULL, dof = 10, Plot = FALSE)
  RI.filt.df <- attr(RT.to.RI.fit, "data")
  # Generate a model that maps the theoretical RI to a reference RT scale
  RI.to.refRT.fit <- get.RI.to.refRT.fit(RT.to.RI.fit, ref.samp = params$ref.samp, Plot = FALSE)

  # Plot calibration models
  if(!is.null(params$plot.dir)){
    # RT to RI models
    RTpredgrid <- expand.grid(day = unique(RI.filt.df$day), RT = seq(1, max(RI.filt.df$RT, na.rm = T), length.out = 200))
    RTpredgrid$batch <- RI.filt.df$batch[match(RTpredgrid$day, RI.filt.df$day)]
    RTpredgrid$RI <- predict(RT.to.RI.fit, newdata = RTpredgrid)
    p <- qplot(data = RTpredgrid, RT, RI, geom=c("line"), col=I("red3")) + # spline linear model fits
      geom_point(data = RI.filt.df, col= RI.filt.df$is.outlier+1) + # actual measured points
      facet_wrap(facets = ~day, labeller = function(x) list(day = RI.filt.df$date[match(x$day, RI.filt.df$day)]) ) +
      xlab("Retentiom time (min)")+ylab("Retention index")
    graph2ppt(x = p, file = paste0(plot.subdir, "RT to RI calibration models per sample.ppt"))
    # RI to reference RT model
    RI.ref.df <- attr(RI.to.refRT.fit, "data")
    RI.new <- 1:(max(RI.ref.df$RI)+1000)
    RT.pred <- predict(RI.to.refRT.fit, newdata = list(RI = RI.new))
    p <- ggplot(data = RI.ref.df, x = RT, y = RI, geom=c("point"), col=I("grey20")) +
      geom_point(aes(x = RT, y = RI), col=I("grey20")) +
      geom_line(aes(x = RT, y = RI), data = data.frame(RI = RI.new, RT = RT.pred), col= "red3") +
      xlab("Reference retention time (min)") + ylab("Theoretical retention index") +
      xlim(range(scantimes.ref)) +
      labs(title = "RI to reference RT mapping") +
      theme(plot.title = element_text(hjust = 0.5))
    graph2ppt(fun = function(){ suppressWarnings(print(p)) },
              file = paste0(plot.subdir, "RI to reference RT calibration model.ppt"))
    if(params$verbose) message(paste0("Plots were saved in: ", plot.subdir))
  }

  # Return the results
  notes <- paste0("The generated list contains the RT calibration models extracted from ",
                  nrow(params$standards.df), " sample containing ", params$std.lib, " as standards.\n",
                  "Computed on: ", Sys.Date(), "\n",
                  "Timing: ",round((proc.time()[3] - t1)/60, 1) ," minutes\n",
                  "The RT to RI model is stored in \'RT.to.RI.model\' (used to convert the observed RT to theoretical RI)\n",
                  "The RI to reference RT model is stored in \'RI.to.RT.model\' (used to convert the theoretical RI to the reference RT)\n",
                  "Additionnaly, reference time scale and mz scale are stored, as well as the calibration data.\n")
  if(params$verbose) cat("\nOUTPUT\n\n", notes)
  out <- list(alk.df = RI.filt.df,
              scantimes = scantimes.ref,
              mzvals = mzvals.ref,
              RT.to.RI.model = RT.to.RI.fit,
              RI.to.RT.model = RI.to.refRT.fit,
              notes = notes)
  if(params$save.backup){
    filen <- paste0(params$output.dir, "1-RT calibration with alkane ladder-", gsub("[-]|[:]|[ ]", "", Sys.time()), ".rds")
    saveRDS(object = out, file = filen)
    message(paste0("Backup saved as: ", filen))
  }
  return(out)
}

# preprocessing ####
#' Performs preprocessing of the samples and combines data in a single matrix.
#'
#' Function to perform preprocessing on the samples and store them in a single
#' disk-backed file (\code{\link{bigmemory}} object). The preprocessing consists
#' in:
#' \itemize{
#' \item{Calibrating the RT of every sample to a reference RT}
#' \item{Interpolate time and mz axis to a common scale}
#' \item{Fit and subtract background}
#' }
#'
#' @param params
#' the ouptut of the \code{\link{initiate.parameters}} function.
#' @param calibration.output
#' the ouptut of the \code{\link{RT.calibration}} function.
#'
#' @return
#'
#' a list with the preprocessing results. Note that a disk-backed file is also
#' saved in the supplied directory (do not manually change or move the file).
#'
#' @export
#'
preprocessing <- function(params,
                          calibration.output){
  if(params$verbose) cat("======= Preprocessing samples =====\n\n")
  # Search for backups
  backup <- check.for.backups(dir = params$output.dir, pattern = "2-Preprocessing")
  if(!is.null(backup)) return(backup)
  # Start timing
  t1 <- proc.time()[3]

  # Get the reference time scans and m/z values
  mzvals <- calibration.output$mzvals
  scantimes <- calibration.output$scantimes

  # Create subfolder for plotting
  if(!is.null(params$plot.dir)){
    plot.subdir <- paste0(params$plot.dir, "2-Preprocessing/")
    dir.create(plot.subdir, showWarnings = FALSE)
  }
  # Create subfolder for storing data
  if(is.null(params$data.dir)){
    data.dir <- paste0(params$output.dir, "disk backed data/") # TODO tempdir should be outputdir but this avoids overloading Tom's Dropbox account
    dir.create(data.dir, showWarnings = FALSE)
  } else if(!dir.exists(data.dir)){
    dir.create(data.dir, showWarnings = FALSE)
  }

  # Load the reference
  if(params$verbose) cat(paste0("Loading the reference sample (", params$ref.samp, ")..."))
  M.blank <- read.cdf(params$ref.samp, mz.ignore = params$mz.ignore)
  if(params$verbose) cat("done\n\n")

  # Initialize the disk backed file
  m <- length(scantimes)
  n <- length(mzvals)
  s <- nrow(params$samples.df)
  if(params$verbose) cat(paste0("Number of samples to load: ", s, "\n",
                                "Creating the disk backed file (encoded as ", params$enc, ")\n"))
  check.bigmory.matrix(data.dir, params$data.filename)
  data <- big.matrix(nrow = m*s, ncol = n, type = params$enc,
                     backingfile = paste0(params$data.filename, ".bin"),
                     descriptorfile = paste0(params$data.filename, ".desc"),
                     backingpath = data.dir)
  if(params$verbose) cat(paste0("Data file stored in: ", data.dir, "\n\n"))

  # Preallocate the baseline components
  # TODO improve RAM usage/speed with: C.bas <- sparse.bs(x = 1:m, nknots = bas.nknots, degree = bas.degree, sparse = T) # with sparsity but not compatible with the fitting functions
  C.bas <- bs(x = 1:m, df = params$bas.nknots - params$bas.degree, degree = params$bas.degree)
  S.bas <- list()

  # Preprocess samples and store in the disk-backed file
  # TODO parallelize, but how to do this with bigmemory objects since it cannot be changed during a loop else R crashes?
  #      Solution is maybe to split samples in batches were each batch contains 'ncores' samples and parallelize within batch, save data, and serialize over beatches ?
  if(params$verbose) cat(paste0("Preprocessing samples...\n"))
  for (i in 1:s){
    if(params$samples.df$file[i] == params$ref.samp){
      # The reference sample is already in memory, requires no calibration, and is already on the reference scale
      M <- M.blank
      if(!is.null(params$plot.dir)) M0 <- M
    } else {
      # 1. Load the file
      M <- read.cdf(params$samples.df$file[i], mz.ignore = params$mz.ignore)
      if(!is.null(params$plot.dir)) M0 <- M
      scantimes.samp <- as.numeric(rownames(M))
      mzvals.samp <- as.numeric(colnames(M))
      day.samp <- params$samples.df$day[i]
      batch.samp <- params$samples.df$batch[i]

      # 2. (optional) RT recalibration
      if(!(is.null(calibration.output$RT.to.RI.model) && is.null(calibration.output$RI.to.RT.model))){
        scantimes.samp0 <- scantimes.samp # The unaligned RT scale
        # mapping RT to the theoretical RI scale
        RI.samp <- predict(calibration.output$RT.to.RI.model,
                           newdata = data.frame(RT = scantimes.samp,
                                                day = day.samp,
                                                batch = batch.samp))
        # mapping the RI to the RT scale of the reference day
        scantimes.samp <- predict(calibration.output$RI.to.RT.model,
                                  newdata = data.frame(RI = RI.samp))
      }

      # 3. Interpolate samples to a common time and mz scale
      # Interpolate scantimes of the sample to match those of reference
      if(!identical(scantimes.samp, scantimes)) {
        M <- apply(M, 2, function(yvals) approxfun(x = scantimes.samp, y = yvals, method = "linear")(scantimes)) # linearly interpolate to common x values
        rownames(M) <- scantimes
        M[is.na(M)] <- 0
      }
      # Interpolate mzvalues of the sample to match those of reference
      if(!identical(mzvals.samp, mzvals)) {
        M <- t(apply(M, 1, function(yvals) approxfun(x = mzvalues.samp, y = yvals, method = "linear")(mzvals))) # linearly interpolate to common x values
        colnames(M) <- mzvals
        M[is.na(M)] <- 0
      }
    }

    # 4. Fit and remove the baseline
    win <- seq(1, m, by = params$bas.subsamp) # Subsample the sample for computational speed up
    # Estimate (sample specific) spectrum profiles of the baseline using quantile regression
    S.bas[[i]] <- fit.S.bas(M = M, C.bas = C.bas, win = win, tau = params$bas.tau)
    M <- pmax(M - Matrix::tcrossprod(C.bas, S.bas[[i]]), 0)

    # 5. Plot the preprocessing
    if(!is.null(params$plot.dir)){
      png(file = paste0(plot.subdir, "Preprocessing - ", gsub("[.]|/", "-", params$samples.df$name[i]), ".png"), res = 300, width = 5000, height = 2500)
      par(mfrow = c(3,1), oma = c(0,0,3,0))
      zlims <- c(0, max(M0))^0.12
      plotmat(M0, subsamp = 10, plot.axes = F, zlim = zlims, main = "Raw data")
      plotmat(M, subsamp = 10, plot.axes = F, zlim = zlims, main = "Preprocessed data")
      plotmat(Matrix::tcrossprod(C.bas[win,], S.bas[[i]]), subsamp = 1, plot.axes = F, zlim = zlims, main = "Fitted baseline")
      title(main = params$samples.df$name[i], outer = T)
      dev.off()
    }

    # Store matrix in the bigmemory object
    mode(M) <- params$enc
    data[(1:m)+(i-1)*m, ] <- M
    flush(data) # Store to file
    invisible(gc()) # Free up RAM
    if(params$verbose) print.progress(i, s)
  }

  # Return results
  notes <- paste0("The list contains the preprocessed data of ", nrow(params$samples.df), " samples\n",
                  "Computed on: ", Sys.Date(), "\n",
                  "Timing: ",round((proc.time()[3] - t1)/60, 1) ," minutes\n",
                  "The bigmemory object containing the data is stored in ", data.dir, "\n",
                  "Elution profiles of the baseline are stored in C.bas and the spectrum\n",
                  "profiles of the baseline are stored in S.bas\n")
  if(params$verbose) cat(paste0("\nOUTPUT\n\n", notes, "\n"))
  out <- list(file.bin = paste0(data.dir, params$data.filename, ".bin"),
              file.desc = paste0(data.dir, params$data.filename, ".desc"),
              data.dir = data.dir,
              C.bas = C.bas,
              S.bas = S.bas,
              notes = notes)
  if(params$save.backup){
    filen <- paste0(params$output.dir, "2-Preprocessing-", gsub("[-]|[:]|[ ]", "", Sys.time()), ".rds")
    saveRDS(object = out, file = filen)
    message(paste0("Backup saved as: ", filen))
  }
  return(out)
}


# deconvolution ####
#' Performs deconvolution of the the preprocessed GC-MS data
#'
#' Function to perform deconvolution of the data. The deconvolution is split in
#' 2 steps
#' \itemize{
#' \item{The untargeted deconvolution blindly extracts the elution and spectrum
#' profiles from the set of untargeted samples.}
#' \item{The targeted deconvolution uses the profiles extracted in the
#' untargeted step as prior knowledge. The spectrum profiles areused as
#' covariates wghile the elution profiles are used to define extraction windows.}
#' }
#'
#' @param params
#' the ouptut of the \code{\link{initiate.parameters}} function.
#' @param calibration.output
#' the ouptut of the \code{\link{RT.calibration}} function.
#' @param preprocess.output
#' the ouptut of the \code{\link{preprocessing}} function.
#'
#' @return
#'
#' a list with the deconvolution results. Note that, if required, two backup
#' files are generated: one containing the untageted deconvolution results (from
#' which the prior information can be extracted), the other containing the
#' final deconvolution.
#'
#' @export
#'
deconvolution <- function(params,
                          calibration.output,
                          preprocess.output){
  if(params$verbose) cat("======= Deconvolution =====\n\n")
  # Search for backups
  deconvolution.out <- check.for.backups(dir = params$output.dir, pattern = "5-Complete deconvolution")
  if(!is.null(deconvolution.out)) return(deconvolution.out)
  # Start timing
  t1 <- proc.time()[3]

  # Get the directory where the bigmemory matrix is stored
  data.dir <- preprocess.output$data.dir

  # Get the reference scan times and m/z values
  scantimes <- calibration.output$scantimes
  mzvals <- calibration.output$mzvals
  m <- length(scantimes)
  n <- length(mzvals)
  s <- nrow(params$samples.df)

  ### A. UNTARGETED

  # Check the presence of a backup
  if(params$verbose) cat(paste0("--- Untargeted Deconvolution ---\n\n"))
  untarg.output <- check.for.backups(dir = params$output.dir, pattern = "3-Untargeted deconvolution")
  if(is.null(untarg.output)){
    # Create the plot directory for untargeted
    if(!is.null(params$plot.dir)){
      plot.dir.untarg <- paste0(params$plot.dir, "3-Untargeted deconvolution/")
      dir.create(plot.dir.untarg, showWarnings = FALSE)
    } else {
      plot.dir.untarg <- NULL
    }
    # Perform untargeted deconvolution
    if(params$verbose) cat(paste0("Performing untargeted deconvolution in ", length(params$untargeted.samples), " samples.\n"))
    untarg.output <- untargeted.deconvolution(data.filename = params$data.filename, data.dir = data.dir,
                                              m = m, n = n, s = s, samps = params$untargeted.samples,
                                              scantimes = scantimes, mzvals = mzvals,
                                              samples.df = params$samples.df,
                                              plot.dir = plot.dir.untarg,
                                              verbose = params$verbose, ncores = params$ncores,
                                              # Peak identification
                                              win.size = params$win.size, peak.width.range = params$peak.width.range,
                                              peak.trim = params$peak.trim, cutoff = params$cutoff, lives = params$lives,
                                              # Deconvoluton
                                              stop.thres = params$stop.thres, tau = params$tau, cobs.tau = params$cobs.tau,
                                              eps = params$eps, maxiter = params$maxiter, tol = params$tol,
                                              basfun.u = params$basfun.u, sel.power = params$sel.power, t.power = params$t.power,
                                              win.ext = params$win.ext, scale.subsamp = params$scale.subsamp,
                                              # Spectrum profiles update
                                              method.upd = params$method.upd)
    if(params$save.backup){
      filen <- paste0(params$output.dir, "3-Untargeted deconvolution-", gsub("[-]|[:]|[ ]", "", Sys.time()), ".rds")
      saveRDS(object = untarg.output, file = filen)
      message(paste0("Backup saved as: ", filen))
    }

    # TODO (optionally) fuse (near) identical spectrum profiles
  }

  ### B. TARGETED

  # Check the presence of a backup
  if(params$verbose) cat("\n--- Targeted Deconvolution ---\n\n")
  targ.output <- check.for.backups(dir = params$output.dir, pattern = "4-Targeted deconvolution")
  if(is.null(targ.output) && length(params$targeted.samples) > 0){
    # Create the plot directory for targeted
    if(!is.null(params$plot.dir)){
      plot.dir.targ <- paste0(params$plot.dir, "4-Targeted deconvolution/")
      dir.create(plot.dir.targ, showWarnings = FALSE)
    } else {
      plot.dir.targ <- NULL
    }
    # First fill missing elution profiles with a prior
    C.prior <- fill.C(C = untarg.output$C, m = m, n = n, s = s,
                      untargeted.samples =  params$untargeted.samples,
                      targeted.samples = params$targeted.samples,
                      win.ext = params$win.ext)
    # Perform targeted deconvolution
    out <- targeted.deconvolution(data.filename = params$data.filename, data.dir = data.dir,
                                  S = untarg.output$S, C.prior = C.prior, m = m,
                                  # Deconvolute all samples to extract profiles in targeted samples, and update profile in untargeted
                                  samps = c(params$untargeted.samples, params$targeted.samples),
                                  samples.df = params$samples.df, scantimes = scantimes, mzvals = mzvals,
                                  plot.dir = plot.dir.targ, verbose = params$verbose,
                                  ncores = params$ncores, t.power = params$t.power, cobs.tau = params$cobs.tau,
                                  smooth = FALSE, dampen = FALSE, unimodal = TRUE)
    if(save.backup){
      filen <- paste0(params$output.dir, "4-Targeted deconvolution-", gsub("[-]|[:]|[ ]", "", Sys.time()), ".rds")
      saveRDS(object = out,
              file = filen)
      message(paste0("Backup saved as: ", filen))
    }

  } else {
    out <- untarg.output
  }
  return(targ.output)
}


####---- CHECKING INPUT ----####


# check.sample.table ####
# TODO doc
check.sample.table <- function(x){
  req <- c("name", "file", "type", "date", "batch")
  misg <- which(!req %in% colnames(x))
  if(length(misg) > 0) stop(paste0("The sample table misses the following field(s): ", paste(misg, collapse = ", ")))
  x$file <- gsub("\\\\", "/", x$file)
  # Make sure there is at least 1 sample with RT calibration compounds
  if(!any(x$type == "standard")) stop("No sample of type 'standard' found. A standard sample is required in order to perform RT/RI calibration.")
  x$date <- as.Date(x$date)
  year <- as.numeric(format(x$date, "%Y"))
  x$day <- as.numeric(format(x$date, "%j")) + 365 * (year-min(year)) # avoid year ambiguities
  x$batch <- as.factor(x$batch)
  # Remove NA's and
  rem <- which(is.na(x[,req]), arr.ind = TRUE)
  if(nrow(rem) > 0){
    warning(paste0("NAs were found when checking the sample table. The associated rows (n = ", length(unique(rem[,1])) ,") are removed."))
    # TODO save removed samples in a file
    x <- x[-rem[,1],]
  }
  return(x)
}

# check.target.table ####
# TODO doc
check.target.table <- function(x){
  req <- c("compound", "RI", "spectrum")
  misg <- which(!req %in% colnames(x))
  if(length(misg) > 0) stop(paste0("The target table misses the following field(s): ", paste(misg, collapse = ", ")))
  # Remove NA's and
  rem <- which(is.na(x[,req]), arr.ind = TRUE)
  if(nrow(rem) > 0){
    warning(paste0("NAs were found when checking the target table. The associated rows (n = ", length(unique(rem[,1])) ,") are removed."))
    # TODO save removed compounds in a file
    x <- x[-rem[,1],]
  }

  # Check if values make sense
  if(any(x$RI < 0)) stop("'RI' must be nonnegative numbers.")
  # Order by expected RI
  x <- x[order(x$RI),] # sort targets with respect to RI
  return(x)
}


####---- BUILDING TARGET LIBRARY ----####



# parse.MSP ####
# TODO finish implementing MSP parsing functions
parse.MSP <- function(){


}

# generate.compound.file ###
# TODO needed ???
# TODO doc
# If compounds is a file should be read by read.table (sep = white space). The table
# should not contain headers as is converted to a character vector
#
# Function to allow user to easily load target compound.
generate.compound.file <- function(base.dir, library, compounds = "all", column = NULL){
  if(!is.character(compounds)) stop("'compounds' must be a either 'all' (fetch all compounds from the library), a vector of character strings giving compound names, or a character string giving the path to a text file containing compound names.")
  if(length(compounds) == 1 && file.exists(compounds)) {
    compounds <- as.character(as.vector(read.table(compounds, header = FALSE)))
  }
  # Load the libraries
  if(!all(file.exists(library))) stop("One the library files you supplied does not exist.")
  spectrum.df <- do.call(rbind, lapply(library, function(lib) return(readRDS(lib)) ))

}


# query.databases ####
# TODO doc
query.databases <- function(library, compounds = NULL, verbose = TRUE){
  # TODO build a function tha automatically queries all compound information
  # See D:\Documents\Dropbox\christophe\_tests deconvolution scripts\190429 Building the pipeline\Intermediate scripts\Generate the target table.R



  stop("efhroiergh")
}




# functgroupstofront ####
# TODO documentation
# function to put any functional groups in the name to the front of the name
functgroupstofront <- function(name) {
  if(grepl("<", name)) {
    inside <- gsub(", ", "", sub("^.*<(.+)>.*$", "\\1", name))
    outside <- sub("^(.*) <.*>(.*)$" , "\\1\\2", name)
    name <- paste0(inside, outside)
  }
  name <- paste(rev(strsplit(name, ", ")[[1]]), collapse = "")
  return(name)
}


# getURL2 ####
# TODO documentation
getURL2 <- function(url, type = "text", max.try = 10) {
  i <- 1
  while(i <= max.try){
    res <- tryCatch(httr::content(GET(url), type, encoding = "UTF-8"), err = function(e) return(NULL))
    if(!is.null(res)) break
    i <- i + 1
  }
  if (grepl(pattern = "Status: 404|Status: 400|Page not found", res) || is.null(res)) res <- NA
  return(res)
}


# names2inchikey ####
# TODO documentation
names2inchikey <- function(queries, verbose = FALSE){
  # Query DBs by order of proference: NIST > CACTUS > PubChem
  df <- do.call(rbind, lapply(queries, function(query){
    query0 <- query
    db <- NA
    inchikey <- NA
    # 1. Try NIST DB
    if (grepl("<", query) | grepl(", ", query) > 0){
      query <- functgroupstofront(query)
    }
    url <- paste0('https://webbook.nist.gov/cgi/cbook.cgi?Name=', URLencode(query),'&Units=SI')
    doc <- getURL2(url, type = "parsed")
    # doc <- read_html(URLencode(url))
    x <- html_nodes(doc, xpath="//li/strong/..") # all elements with a bold title
    nistname <- html_text(html_nodes(doc, "h1")[2])
    if(nistname == "Search Results"){ # In case there are mulitiple matches
      link <- html_nodes(html_nodes(doc, "main"), "a")[1] # Take only the first match
      link <- gsub(x = link, pattern = "^.*href[=]\"", "")
      link <- gsub(x = link, pattern = "\">.*$", "")
      url <- paste0("https://webbook.nist.gov", link)
      doc <- getURL2(url, type = "parsed")
      # doc <- read_html(URLencode(url))
      x <- html_nodes(doc, xpath="//li/strong/..") # all elements with a bold title
    }
    x <- x[which(grepl("IUPAC Standard InChIKey:",x))]
    inchikey <- gsub("IUPAC Standard InChIKey:\n","", html_text(x)) # std inchikey
    inchikey[identical(inchikey,character(0))] <- NA
    if(!is.na(inchikey)) db <- "NISTwb"

    # 2. Try CACTUS
    query <- query0
    if(is.na(inchikey)){
      inchikey <- getURL2(paste0("https://cactus.nci.nih.gov/chemical/structure/",URLencode(query),"/stdinchikey"))
      if (!is.na(inchikey)) inchikey <- gsub("InChIKey=", "", inchikey)
      db <- "CACTUS"
    }

    # 3. Try PubChem
    if(is.na(inchikey)){
      cid <- getURL2(paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",URLencode(query),"/cids/TXT"))
      if(is.character(cid)) cid <- strsplit(cid, "\n")[[1]][1]
      if(is.na(cid)){
        inchikey <- NA
      } else {
        inchikey <- getURL2(paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",URLencode(cid),"/property/InChIKey/TXT"))
        inchikey <- strsplit(inchikey, "\n")[[1]]
        db <- "PubChem"
      }
    }
    closeAllConnections()
    if(verbose) print.progress(which(queries == query), length(queries))
    return(data.frame(query = query, inchikey = inchikey, db.inchikey = db, stringsAsFactors = FALSE))
  }))
  return(df)
}

# inchikeys2cids ####
# TODO documentation
inchikeys2cids <- function(inchikeys, verbose = FALSE){
  cids <- sapply(inchikeys, function(inchikey){
    if(is.na(inchikey)) return(NA)
    url <- paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",URLencode(inchikey),"/cids/TXT")
    cid <- getURL2(url)
    if(!is.na(cid)) cid <- as.numeric(strsplit(cid,"\n")[[1]][1])
    if(verbose) print.progress(which(inchikeys == inchikey)[1], length(inchikeys))
    return(cid)
  })
}

# cids2RIs ####
# TODO documentation
cids2RIs <- function(cids, column_class = "Semi-standard non-polar",
                     fun = median, verbose = FALSE){
  RI.df <- do.call(rbind, lapply(cids, function(cid){
    empty.df <- data.frame(cid = cid, column_class = NA, RI = NA, RI.sd = NA, RI.min = NA, RI.max = NA, RI.n = NA)
    if(is.na(cid)) return(empty.df)
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/",cid,"/XML/?response_type=display")
    doc <- getURL2(url, type = "parsed")
    if(is.na(doc[1])) return(empty.df) # TODO there is an issue with GET when access the RI webpage for cid = 335 !!!
    doc <- read_xml(doc)
    ns <- xml_ns(doc)
    information <- xml_find_all(doc, "//d1:TOCHeading[text()='Kovats Retention Index']/following-sibling::d1:Information", ns)
    indexes <- do.call(rbind, lapply(information, function(x) {
      name <- xml_text(xml_find_first(x, "./d1:Name", ns))
      # values <- as.numeric(sapply(xml_find_all(x, "./*[contains(name(), 'NumValue')]", ns), xml_text))
      values <- as.numeric(sapply(xml_find_all(x, "./d1:Value/d1:Number", ns), xml_text))
      return(data.frame(pubchemid = cid,
                        column_class = name,
                        kovats_ri = values, stringsAsFactors = F))
    }))
    sel <- indexes$column_class %in% c(column_class, "All column types") & !is.na(indexes$kovats_ri)
    if(verbose) print.progress(which(cids == cid), length(cids))
    if(sum(sel, na.rm = T) == 0){
      return(empty.df)
    } else {
      vals <- indexes$kovats_ri[sel]
      return(data.frame(cid = cid, column_class = column_class, RI = fun(vals),
                        RI.sd = sd(vals), RI.min = min(vals), RI.max = max(vals),
                        RI.n = length(vals)))
    }
  }))
  return(RI.df)
}

# inchikeys2specs ####
# TODO dicumentation
inchikeys2specs <- function(inchikeys){
  cat("Load database...")
  load("D:/Documents/Dropbox/christophe/GCMS_libraries/nist2011DB.rda") # TODO find a way to improve this, maybe using MONA instead ?
  cat("done\n")
  lib.inchikey <- sapply(nist2011DB, function(x) x$inchikey)
  spec <- sapply(inchikeys, function(inchikey){
    if(is.na(inchikey)) return(NA)
    inchikey.isom <- strsplit(inchikey, "-")[[1]][1]
    ind <- which(grepl(x = lib.inchikey, pattern = inchikey.isom)) # Match the common isomeric inchikey in case the stereo isomer is absent from db
    if(length(ind) == 0) return(NA)
    if(length(ind) > 1){ # If there are several matches, test the match with full InChIKey, if no more matches are found, keep the first index
      ind.new <- which(lib.inchikey == inchikey)
      if(length(ind.new) == 1){
        ind <- ind.new
      } else {
        ind <- ind[1]
      }
    }
    sp <- matrtotext(nist2011DB[[ind]]$pspectrum, norm = TRUE,
                     digits = 2, sort = TRUE)
    print.progress(which(inchikeys == inchikey)[1], length(inchikeys))
    return(sp)
  })
  return(spec)
}

# matrtotext ####
# TODO documentation
# auxilliary function to convert spectra to text
matrtotext <- function(mat, norm=TRUE, digits=2, sort=TRUE) {
  if(!inherits(mat, "matrix")){ mat <- matrix(unlist(mat),ncol=2,byrow=T) }
  txt <- NA
  if (nrow(mat)>0 & !is.na(mat[[1]])) {
    if (norm){
      mat[,2] <- round(100*mat[,2]/max(mat[,2]),digits)
    } else {
      mat[,2] <- round(mat[,2],digits)
    }
    if (sort) {
      ord <- order(mat[,2],decreasing=T) # sort from high to low intensity
      mat <- mat[ord,,drop=F]
    }
    txt <- paste(apply(mat, 1, paste, collapse=" ("),")",sep="",collapse=", ")
  }
  return(txt)
}


# parse.spectrum ####
# TODO documentation
# function to convert single spectrum in txt format txtspec to single column matrix with intensities & mz vals as rownames
parse.spectrum <- function(txtspec, mzout, ntype = 1) {
  txtspec <- as.character(txtspec)
  spec <- matrix(as.numeric(sapply(strsplit(txtspec," \\(|, "),
                                   function(s) gsub(pattern = ")", replacement = "", x = s, fixed=T))),
                 ncol = 2, byrow = TRUE)
  colnames(spec) <- c("mz","int")
  S <- cbind(mz = mzout, int = 0)
  mz.match <- match(S[,"mz"], spec[,"mz"], nomatch = 0)
  S[mz.match != 0, "int"] <- spec[mz.match[mz.match != 0], "int"]
  rownames(S) <- S[,"mz"]
  S <- normS(S, ntype = ntype)
  S <- S[,2,drop=F]
  S <- Matrix(S, sparse = TRUE) # return it as sparse matrix
  return(S) }
# e.g. txt2spec(targets$spectrum[[1]], mzrange=mzrange)

# parse.spectra ####
# TODO documentation
# function to convert vector of spectra in txt format txtspec to matrix S with intensities and mz vals as rownames
parse.spectra <- function(txtspec, mzout, coln = NULL, ntype = 1){
  S <- do.call(cbind, lapply(txtspec, function(txtspec) parse.spectrum(txtspec, mzout = mzout, ntype = ntype)) )
  if(is.null(coln)) coln <- 1:ncol(S)
  colnames(S) <- coln
  return(S)
}


####---- RT CALIBRATION ----####


# load.alkane.lib ####
# TODO documentation
load.alkane.lib <- function(){
  load("./libraries/alkstandDB.rda") # TODO Put the alkane library in misc folder when building package.
  return(alkstandDB)
}


# lib.to.spec ####
# TODO documentation
lib.to.spec <- function(lib, mzvals, norm = 0){
  spec <- rep(0, length(mzvals))
  names(spec) <- mzvals
  specs <- sapply(lib, function(x){
    x <- x$pspectrum
    spec[match(x[,1], mzvals, nomatch = 0)] <- x[match(mzvals, x[,1], nomatch = 0) ,2]
    return(spec)
  })
  colnames(specs) <- sapply(lib, "[[", "Name")
  rownames(specs) <- mzvals
  # Normalization
  if(norm == 1){
    specs <- apply(specs, 2, function(x) return(x/sum(x)))
  } else if(norm == 2){
    specs <- apply(specs, 2, function(x) return(x/sqrt(sum(x^2))))
  } else if(norm == "inf"){
    specs <- apply(specs, 2, function(x) return(x/max(x)))
  }
  return(specs)
}


# find.standards ####
# TODO documentation
find.standards <- function(lib, specs, RT, RI, df = NULL, sr = NULL,
                           rm.double.match = TRUE){
  # Regress the alkane spectra against the sample spectrum profiles
  df <- do.call(rbind, lapply(1:ncol(lib), function(i){
    if(!is.null(sr)){
      if(is.null(df$RT.pred)) stop("RT.pred is missing in df")
      sel <- which(RT >= df$RT.pred[i] - sr[i] & RT <= df$RT.pred[i] + sr[i])
      if(length(sel) == 0)
        return(data.frame(alkane = colnames(lib)[i], RT = NA, RI = RI[i],
                          sim = NA, stringsAsFactors = F))
    } else {
      sel <- 1:ncol(specs)
    }
    y <- lib[,i]^(1/3)
    x <- specs[,sel,drop=F]^(1/3)
    coefs <- nnls(b = y, A = x)$x
    best <- sel[which.max(coefs)]
    sim <- cosine.sim(X = x[,sel == best], Y = y)
    return(data.frame(alkane = colnames(lib)[i], RT = RT[best], RI = RI[i],
                      sim = sim, stringsAsFactors = F))
  }))
  if(rm.double.match){
    # If some spectrum match the same alkane, keep the one with highest similarity
    for(rt.dp in unique(df$RT[duplicated(df$RT)])){
      if(is.na(rt.dp)) next()
      is.dup <- which(df$RT == rt.dp)
      sub <- df[is.dup,]
      best <- which.max(sub$sim) # TODO keeping best sim or removing all ?
      sub[-best,c("RT", "sim")] <- NA
      df[is.dup,] <- sub
    }
  }
  return(df)
}


# extract.RTs ####
# TODO documentation
extract.RTs <- function(standards.df, ref.samp, std.lib = "alkanes", ncores = 1,
                        plot.folder = NULL, mz.ignore = NULL,
                        # Baseline fitting
                        bas.subsamp = 10, bas.tau = 0.2, bas.nknots = 15,
                        bas.degree = 3,
                        # Peak identification
                        win.size = 300, peak.width.range = c(1,20),
                        peak.trim = 1E-10, cutoff = 1,
                        # Deconvoltion
                        tau = 0.1, cobs.tau = 0.1, eps = 1, maxiter = 10,
                        tol = 0.01, basfun.u = function(x) return(min(x)),
                        sel.power = 2, t.power = 2, win.ext = 25,
                        scale.subsamp = 10){
  std.lib <- match.arg(std.lib, c("alkanes", "FAMEs"))

  # Get the samples to extract RT/RI information from
  files <- standards.df$file
  if(is.null(files)) stop("No standard found")
  # Extract reference sample information
  M.ref <- read.cdf(file = ref.samp, mz.ignore = mz.ignore)
  mzvals.ref <- as.numeric(colnames(M.ref))
  scantimes.ref <- as.numeric(rownames(M.ref))
  # Get the average peak width and Gaussian banded  matrix to deconvolute TICs
  TIC.ref <- rowSums(M.ref)
  peak.widths <- fit.average.peak.width(y = TIC.ref,
                                        win.size = win.size,
                                        peak.width.range = peak.width.range)
  peak.width.fun <- peak.widths$peak.width.fun

  # Get information from the target library
  if(std.lib == "alkanes") targ.lib <- load.alkane.lib() # TODO allow also for FAMEs
  specs.lib <- lib.to.spec(lib = targ.lib, mzvals = mzvals.ref, norm = 1)
  monoMW <- round(sapply(targ.lib, "[[", "monoMW"))
  RIs <- sapply(targ.lib, "[[", "RI")
  # Prepare the cluster for parallelization
  cl <- makeCluster(ncores,type="SOCK") # on linux use type="FORK" (faster), but not available on windows
  doSNOW::registerDoSNOW(cl)
  snow::clusterExport(cl, list=c("standards.df", "specs.lib", "files", "win.size", "plot.folder", "monoMW", "RIs",
                                 "bas.subsamp","bas.tau","bas.nknots","bas.degree","win.size","peak.width.range",
                                 "peak.trim","cutoff","tau", "cobs.tau", "eps", "maxiter", "tol",
                                 "basfun.u", "sel.power", "t.power", "win.ext", "scale.subsamp",
                                 "read.cdf", "cosine.sim", "plotmat"),
                      envir = environment()) # export values that each worker needs
  # Execute parallelized extraction
  RI.list <- clusterApply(cl, files, function(file) { # Parallelize over samples
    # RI.list <- lapply(files, function(file) { # Parallelize over samples
    # Load required functions
    source("D:/Documents/Dropbox/christophe/_tests deconvolution scripts/utils-0.0.2.R") # TODO adapt when building package

    # Load the data
    M <- read.cdf(file = file, mz.ignore = mz.ignore)
    mzvals <- as.numeric(colnames(M))
    scantimes <- as.numeric(rownames(M))

    # Remove the baseline
    C.bas <- bs(x = 1:nrow(M), df = bas.nknots - bas.degree, degree = bas.degree)
    S.bas <- fit.S.bas(M = M, C.bas = C.bas, tau = bas.tau, win = seq(1,nrow(M), by = bas.subsamp))
    M <- pmax(M - Matrix::tcrossprod(C.bas, S.bas), 0)

    # Identify peaks
    # TODO try reusing get.peak.locations() by adapting the function
    TIC <- rowSums(M)
    x <- bs(1:length(TIC), df = 3, degree = 2)
    TIC.rq <- x %*% rq.fit.br(x = x, y = TIC, tau = 0.9)$coefficient # Signal used as cutoff
    # plot(TIC, log = "y", type = "l")
    # lines(TIC.rq, col = 2)
    peak.widths.preds <- predict(peak.width.fun, newdata = list(x = 1:length(TIC)))
    C.dic <- build.banded.gaus.mat(w = peak.widths.preds, tol = peak.trim)
    peak.coefs <- block.deconvolution(y = TIC, X = C.dic, win.size = win.size)
    peak.coefs <- cluster.peak.coefficients(coefs = peak.coefs)
    peak.max <- which(peak.coefs != 0)
    peak.height <- peak.coefs[which(peak.coefs != 0)]
    peak.locs <- cbind(samp = 1, max = peak.max, height = peak.height, crit = 0,
                       lives = 1, active = 0)
    peak.locs <- update.detection(Y = TIC, Dic = C.dic, peaks = peak.locs,
                                  sel = NULL, cutoff = cutoff)
    peak.locs <- peak.locs[order(peak.locs[,"active"], peak.locs[,"height"], decreasing = T),]
    peak.locs <- peak.locs[peak.locs[,"active"] == 1 & peak.locs[,"height"] > TIC.rq[peak.locs[,"max"]],]

    # Perform peak deconvolution
    deconv.res <- sequential.deconvolution(resids = M,
                                           scantimes = scantimes, mzvalues = mzvals,
                                           m = nrow(M), n = ncol(M), s = 1, samps = 1,
                                           TICs = TIC, peak.locs = peak.locs, C.dic = C.dic,
                                           tau = tau, cobs.tau = cobs.tau, eps = eps,
                                           maxiter = maxiter, tol = tol, cutoff = cutoff,
                                           basfun.u = basfun.u, sel.power = sel.power, t.power = t.power,
                                           win.ext = win.ext, scale.subsamp = scale.subsamp,
                                           Plot = FALSE, Plot.to.file = FALSE,
                                           plot.folder = NULL, verbose = FALSE, debug = FALSE)
    C.targ <- deconv.res$C
    specs.targ <- deconv.res$S
    RT.meas <- apply(C.targ, 2, which.max)

    # Search for targets
    targ.df <- find.standards(lib = specs.lib, specs = specs.targ,
                              RT = RT.meas, RI = RIs, df = NULL, sr = NULL,
                              rm.double.match = TRUE)
    # Fit a monotonic increasing model to constrain the search range
    sub <- complete.cases(targ.df)
    fit <- cobs(x = targ.df$RI[sub], y = targ.df$RT[sub], constraint = "increase",
                tau = 0.01, w = targ.df$sim[sub], nknots = 6,
                print.mesg = FALSE, print.warn = FALSE)
    targ.df$RT.pred <- predict(fit, z = targ.df$RI)[,"fit"]
    sr <- diff(targ.df$RT.pred) # define a search range
    sr <- c(sr, sr[length(sr)]) # add a last range to keep vector size
    # Refine search
    targ.df <- find.standards(lib = specs.lib, specs = specs.targ,
                              RT = RT.meas, RI = RIs, df = targ.df, sr = sr,
                              rm.double.match = TRUE)
    # Plot extraction
    if(!is.null(plot.folder)){
      filename <- paste0(plot.folder, "RT extraction ", which(files == file), ".png")
      png(filename = filename, height = 2000, width = 4000, res = 300)
      subsamp = 5
      plotmat(M, subsamp = subsamp, main = paste0("RT extraction in sample \n", file))
      text(x = targ.df$RT/subsamp,
           y = match(monoMW, rownames(specs.lib)) + 25,
           labels = colnames(specs.lib), col = "white")
      points(x = targ.df$RT/subsamp,
             y = match(monoMW, rownames(specs.lib)), col = "white")
      dev.off()
    }
    # Add required meta data for further modeling
    meta.ind <- which(standards.df$file == file)
    targ.df$sample <- file
    targ.df$day <- standards.df$day[meta.ind]
    targ.df$date <- standards.df$date[meta.ind]
    targ.df$batch <- standards.df$batch[meta.ind]
    return(targ.df)
  })
  stopCluster(cl)
  RI.df <- do.call(rbind, RI.list)
  RI.df$RT.index <- RI.df$RT
  RI.df$RT <- scantimes.ref[RI.df$RT.index]
  return(list(RI.df = RI.df,
              scantimes.ref = scantimes.ref,
              mzvals.ref = mzvals.ref))
}


# get.RT.to.RI.fit ####
# TODO documentation
get.RT.to.RI.fit <- function (alk.df, rm.outliers = TRUE,
                              formul = NULL, dof = NULL, alpha = 0.05,
                              ...){
  # TODO dirty hack, find a more clever way to pass argument to user supplied formula
  list2env(list(...), envir = environment()) # In case additional parameters are needed in the formula
  alk.df <- alk.df[complete.cases(alk.df),] # Keep only observed alkane data
  # Get the number of degrees of freedom
  if(is.null(dof)) dof <- round((max(alk.df$RT)-min(alk.df$RT))/10) # TODO is there a rational way to set this, for now 1 knot every 10 mins
  # Get the default formula depending on the available covariates
  if(is.null(formul)){
    if(length(unique(alk.df$day)) > 1){ # Include day as a covariate to take day to day bias into account
      if(length(unique(alk.df$batch)) > 1){ # Include batch as a covariate to take bias into account that occur because batch of samples have been acquired at different time periods
        dof <- max(round(dof/(length(unique(alk.df$batch))-1)), 1)
        formul <- "RI ~ batch * ns(RT, df = dof) + day"
      } else { # if only 1 batch
        formul <- "RI ~ ns(RT, df = dof) + day"
      }
    } else { # If samples were taken the same day
      formul <- "RI ~ ns(RT, df = dof)"
    }
  }
  # Fit model
  RT.to.RI.fit <- lm(formula = as.formula(formul), data = alk.df)
  # (Optional) outlier removal
  alk.df$is.outlier <- FALSE
  if(rm.outliers){
    res <- RT.to.RI.fit$residuals
    # Bonferroni outlier t-test
    out.ind <- which(2*(pt(abs(res/sd(res)), length(res)-1, lower.tail = FALSE)) < alpha )
    alk.df[out.ind,"is.outlier"] <- TRUE
    # Refit the model
    RT.to.RI.fit <- lm(formula = as.formula(formul), data = alk.df[!alk.df$is.outlier,])
  }
  if(Plot){
    plot(y = alk.df$RI, x = alk.df$RT, pch = 16,
         col = ifelse(alk.df$is.outlier, "orange2", "green4"),
         xlim = c(0, max(alk.df$RT, na.rm = T)), ylim = c(0,max(alk.df$RI)),
         main = "RT/RI calibration", xlab = "Empirical RT", ylab = "Theoretical RI")
  }
  attr(RT.to.RI.fit,"data") <- alk.df
  attr(RT.to.RI.fit,"dof") <- dof
  return(RT.to.RI.fit)
}

# get.RI.to.refRT.fit ####
# TODO documentation
get.RI.to.refRT.fit <- function(RT.to.RI.fit, ref.samp, Plot = FALSE){
  RI.ref.df <- attr(RT.to.RI.fit, "data")
  RI.ref.df <- RI.ref.df[RI.ref.df$sample == ref.samp & !RI.ref.df$is.outlier,]
  df <- 5 # TODO select automatically the dof based on BIC (or LOOCV?)
  RI.to.refRT.fit <- lm(RT ~ ns(RI, df = df), data = RI.ref.df)
  attr(RI.to.refRT.fit, "data") <- RI.ref.df
  if(Plot){
    plot(x = RI.ref.df$RT, y = RI.ref.df$RI, col = "grey20", pch = 16,
         xlim = c(0, max(RI.ref.df$RT, na.rm = T)), ylim = c(0,max(RI.ref.df$RI, na.rm = T)),
         xlab = "Reference retention time (min)", ylab = "Theoretical retention index",
         main = "RI to reference RT mapping")
    RI.new <- 1:(max(RI.ref.df$RI)+1000)
    RT.pred <- predict(RI.to.refRT.fit, newdata = list(RI = RI.new))
    lines(x = RT.pred, y = RI.new, col = "red3")
  }
  return(RI.to.refRT.fit)
}



####---- PREPROCESSING THE DATA----####


# fit.S.bas ####
# TODO documentation
fit.S.bas <- function(M, C.bas, win = 1:nrow(M), tau = 0.2){
  S.bas <- t(apply(M, 2, function(y){
    coefs <- rq.fit.fnb(x = C.bas[win,], y = y[win], tau = tau)$coefficients
    coefs <- pmax(coefs, 0) # There is no non negativity constraint, so clip negatives to 0
    return(coefs)
  }))
}

# sparse.bs ####
# TODO adapt not sparse
# Function that generates a sparse matrix with B-spline basis
# INPUT
#   x:      a numeric vector of values at which to evaluate the B-spline functions or derivatives.
#   nknots: number of knots
#   degree: degree of the piecewise polynomial-default is 3 for cubic splines.
#   intercept:  if TRUE, an intercept is included in the basis; default is FALSE.
#   sparse: if TRUE, the output is encoded in a sparse matrix of class "dgCMatrix"
# OUTPUT
#   basis:  matrix containing the spline basis
# NOTE
# The code is a simplified version of splines::bs().s
sparse.bs <- function(x, nknots, degree = 3, intercept = FALSE, sparse = FALSE){
  df <- nknots + degree
  ord <- 1L + (degree <- as.integer(degree))
  nIknots <- df - ord + (1L - intercept)
  knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                     2L)[-c(1L, nIknots + 2L)]
  knots <-   quantile(x, knots)
  Aknots <- sort(c(rep(range(x), ord), knots))
  basis <- splineDesign(Aknots, x, ord, sparse = sparse)
  if(!intercept) basis <- basis[, -1L, drop = FALSE]
  return(basis)
}


####---- UNTARGETED DECONVOLUTION ----####


# untargeted.deconvolution ####
# TODO doc
# TODO the verbose in the function is to extensive and should be reduced to a few lines instead of a recurvsive prompt
untargeted.deconvolution <- function(data.filename, data.dir,
                                     m, n, s, samps,
                                     scantimes, mzvals,
                                     plot.dir = NULL, verbose = TRUE,
                                     samples.df,
                                     # Peak identification arguments
                                     win.size = 300,
                                     peak.width.range = c(1,20),
                                     peak.trim = 1E-10, # the gaussian peak threshold under which the peak is trimmed to 0.
                                     cutoff = 1,
                                     lives = 2,
                                     ncores = 1,
                                     # Deconvolution algorithm
                                     stop.thres = 2E-4,
                                     tau = 0.1,
                                     cobs.tau = 0.1,
                                     eps = 1, # If elution values are below 1, we trim to 0 because we expect it to be integer counts
                                     maxiter = 10,
                                     tol = 0.01,
                                     basfun.u = function(x) return(min(x)),
                                     sel.power = 2,
                                     t.power = 2,
                                     win.ext = 25,
                                     scale.subsamp = 10,
                                     # Spectrum profiles update arguments
                                     method.upd = "nnL0pois"){ # or WNNLS
  t1 <- proc.time()[3]

  # Attach matrix
  data <- attach.big.matrix(dget(paste0(data.dir, data.filename, ".desc"))) # Baseline removed data
  resids.filename <- paste0(data.filename, "_resids")
  if(file.exists(paste0(data.dir, resids.filename, ".bin"))) file.remove(paste0(data.dir, resids.filename, ".bin"))
  resids <- deepcopy(data, type = "double", # TODO maybe integer is sufficient, this reduce by 2 storage space
                     backingfile = paste0(resids.filename, ".bin"),
                     descriptorfile = paste0(resids.filename, ".desc"),
                     backingpath = data.dir)
  invisible(gc())

  # Peak identification
  # Compute the TIC for all samples
  TICs <- TICs.init <- as.vector(BigRowSums(pBigMat = resids@address)) # 0.61 s for 14 samples
  ref.samp <- ceiling(which.max(TICs)/m)
  # TODO use L0glm instead of WNNLS ?
  # TODO Work on multivariate data instead of TIC?
  system.time(peak.locs <- get.peak.locations(TICs = TICs, m = m,
                                              samps = samps, ref.samp = ref.samp,
                                              win.size = win.size, peak.width.range = peak.width.range,
                                              peak.trim = peak.trim, cutoff = cutoff,
                                              lives = lives,
                                              ncores = ncores,
                                              plot.file = FALSE,
                                              verbose = verbose))
  C.dic <- peak.locs$bandC
  peak.locs.init <- peak.locs0 <- peak.locs$peak.locs

  # Deconvolution step
  # Perform sequential profile extraction
  deconv.res <- sequential.deconvolution(resids = resids,
                                         m = m, n = n, s = s, samps = samps,
                                         TICs = TICs.init,
                                         scantimes = scantimes.ref, mzvalues = mzvals.ref,
                                         peak.locs = peak.locs.init, C.dic = C.dic,
                                         stop.thres = stop.thres,
                                         tau = tau,
                                         cobs.tau = cobs.tau,
                                         eps = eps,
                                         maxiter = maxiter, tol = tol,
                                         cutoff = cutoff,
                                         basfun.u = basfun.u,
                                         sel.power = sel.power, t.power = t.power,
                                         win.ext = win.ext, scale.subsamp = scale.subsamp,
                                         Plot = FALSE, Plot.to.file = FALSE, plot.folder = NULL, # Do not plot, only for testing
                                         verbose = verbose, debug = FALSE)

  if(verbose) cat(paste0("Untargeted deconvolution finished after fitting ", deconv.res$iter, " compounds.\n\n"))
  # Plot deconvolution results
  if(!is.null(plot.dir)){
    if(verbose) cat("Plotting the deconvolution results...\n")
    test.win <- 1:m
    sub.fac <- 5
    # samp.plot <- 11
    for(samp.plot in samps){
      win.plot <- rep(test.win, length(samp.plot))+rep(samp.plot-1, each = length(test.win))*m # Plot full chromatogram of a single sample
      # Plot fitted data
      samp.name <- gsub(x = samples.df$name[samp.plot], pattern = "/|[.]", replacement = "-")
      file.name <- paste0(plot.dir, "Deconvolution-",samp.name, ".png")
      plot.deconvolution.fit(M = data, M.res = resids,
                             C = deconv.res$C, S = deconv.res$S,
                             file.name = file.name,
                             win = win.plot, m = m,
                             gam = 0.1, sub.fac = sub.fac)
      print.progress(which(samps == samp.plot), length(samps))
    }
  }

  # Simultaneous update of the spectrum
  # Order profiles by increasing elution time (ie scanline where maximum occures
  # in profile)
  if(verbose) cat("Sorting profiles by increasing retention time.\n")
  ord <- scanline.order(C = deconv.res$C, m = m) # 35 s for 1876 compounds and 50 samples
  # Update spectrum profiles
  if(verbose) cat(paste0("Updating spectrum profiles\n",
                         "Number of compounds: ", ncol(deconv.res$S), "\n",
                         "Update method: ", method.upd, "\n"))
  system.time(upd <- update.spectrum(data = data,
                                     C = deconv.res$C[,ord,drop=F],
                                     S = deconv.res$S[,ord,drop=F],
                                     m = m,
                                     method = method.upd, focus = TRUE, verbose = verbose))
  # Plot the change in spectrum profiles
  if(!is.null(plot.dir)){
    png(file = paste0(plot.dir, "Spectrum update.png"), height = 3000, width = 3000, res = 300)
    par(mfrow = c(2,1))
    plotmat(t(as.matrix(deconv.res$S[,ord,drop=F])), subsamp = 1, xlab = "Compound index", main = "Spectrum profiles before update")
    plotmat(t(as.matrix(upd$S)), subsamp = 1, xlab = "Compound index", main = "Spectrum profiles after update")
    dev.off()
  }

  # Return results
  notes <- paste0("The list contains the deconvolution results for ", s, " samples.\n",
                  "The algorithm stopped after fitting ", deconv.res$iter, " components\n",
                  "(based on a stopping threshold of ", stop.thres, ").\n",
                  "'deconvolution.results' contains the main function output, 'C'\n",
                  "and 'S' contain de deconvoluted elution and spectrum profiles,\n",
                  "respectively.\n",
                  "Computed on: ", Sys.Date(), "\n",
                  "Timing: ",round((proc.time()[3] - t1)/60, 1) ," minutes\n",
                  "The data is stored in: ", data.dir, "\n")
  if(verbose) cat(paste0("\nOUTPUT\n\n", notes, "\n"))
  out <- list(C = upd$C, S = upd$S, # Most important !
              # Disk-backed residual file
              resids.desc = paste0(data.dir, file.name, ".desc"),
              resids.bin = paste0(data.dir, file.name, ".bin"),
              # Deconvolution results
              deconvolution.results = deconv.res,
              TICs.init = TICs.init, TICs.resid = TICs,
              # Peak identification
              peak.locs.init = peak.locs.init,
              peak.locs.resid = deconv.res$peak.locs,
              C.dic = C.dic,
              # Other
              notes = notes)
  return(out)
}


# check.bigmory.matrix ####
# TODO doc
check.bigmory.matrix <- function(data.dir, data.filename){
  delete.file <- ""
  if(file.exists(paste0(data.dir, data.filename, ".bin"))){
    cat(paste0("The data file '", data.dir, data.filename, ".bin' already exists. Do you want to overwrite it [yes] or not [no] ?"))
    while(!delete.file %in% c("yes", "no")){
      delete.file <- readline(prompt = "Overwrite file? [yes] or [no]: ")
    }
    if(delete.file == "yes"){
      file.remove(paste0(data.dir, data.filename, ".bin"))
    } else {
      cat("You can either specify a new file name ('data.filename') or a new data storage directory ('data.dir')")
      stop("Aborting run.")
    }
  }
  return(NULL)
}



# fit.average.peak.width ####
# Function that optimizes the gaussian peak width such that the signal (y) is best
# decomposed by a sum of gaussian peaks of that width. The signal can be fit in windows
# separately and the width across windows is fit using natural splines with median
# regression (to limit the effect of outliers).
# INPUT
#   y:  the data (vector) to fit
#   win.size:   the size of the windows. To avoid splitting y in windows, set win.size = length(y)
#   peak.width.range:   a lower and upper bound for the peak width
# OUTPUT
#   peak.widths:    the optimized peak widths for every window
#   peak.width.fun: the median regression fit
#   peak.width.fitted: the smoothed estimated peak width at every data point
fit.average.peak.width <- function(y, win.size = 200, peak.width.range = c(1,10)){
  # 1. Split the signal in windows
  m <- length(y)
  nwins <- ceiling(m/win.size)
  y.split <- split(y, rep(1:nwins, each = win.size)[1:m])
  # 2. For every window optimize the width of
  peak.widths <- sapply(y.split, function(y.win){
    y.win <- pmax(y.win - min(y.win[y.win!=0]),0) # Remove a baseline
    # Optimize the peak width in the current window
    peak.width <- optimize(y = y.win, # fitqual.cv specific arguments
                           f = fitqual.cv, interval = peak.width.range, maximum = FALSE, tol=1E-3)$minimum
    return(peak.width)
  })
  # 3. Model the peak widths for every data point
  if(nwins > 5){
    peak.width.fun <- rq(y ~ 1 + ns(x, df = 5), # use natural spline function
                         data = data.frame(y = peak.widths,
                                           x = (1:length(peak.widths))*win.size-win.size/2), # x is the index of the center of every window
                         tau = 0.5) # Median regression is robust against outliers
  } else { # If not enough windows, take the median of the peak widths
    peak.width.fun <- list(fitted.values = median(peak.widths))
  }
  return(list(peak.widths = peak.widths,
              peak.width.fun = peak.width.fun))
}

# fitqual.cv ####
# Function giving the fit quality given a gaussian peak width that is used
# to construct a banded covariate matrix of shifted gaussian peaks which
# is fit to the data (y) using NNLS. The fit quality is given by the BIC
# based on 2-fold CV, splitting the dataset in even (training set) and odd
# (test set) points.
# INPUT
#   width:  the peak width of the gaussian model
#   y:      the data to fit the shifted gaussian model (banded matrix)
# OUTPUT
#   fitqual:  the 2-fold CV BIC of the NNLS fit
# NOTE
# The 2-fold CV is a bit less accurate than without validation since
# half of the data is used for fitting, but this avoids the spurious minimum
# that appears in the region of small peak widths.
fitqual.cv <- function(width, y) {
  train <- ((1:length(y)) %% 2) == 0
  test <- ((1:length(y)) %% 2) != 0
  X <- build.banded.gaus.mat(w = rep(width, length(y)), tol = 1E-5)
  coefs <- nnls(A = X[train,], b = y[train])$x
  fitqual <- get.BIC(y = y[test],
                     yhat = Matrix::tcrossprod(X[test,], coefs),
                     npars = sum(coefs > 0),
                     transf = function(y) sqrt(y)) # Transforming on a sqrt scale leads to better results (common tranformation to account for Poisson error structure)
  return(fitqual)
}


# build.banded.gaus.mat ####
# Function to generate a square banded matrix of shifted gaussian profiles. Different peak widths
# can be supplied, but then the matrix will not be banded strictly speaking.
# INPUT
#   w:    the widths of the gaussian peaks. The length of w will give the dimensions of the (square)
#         banded matrix. Different widths are allowed.
#   tol:  value under which the gaussian profile should be trimmed to 0.
# OUTPUT
#   mat:  the banded matrix with (trimmed) gaussian peak profiles.
# NOTE
# Inspired from https://stackoverflow.com/questions/30536062/preallocate-sparse-matrix-with-max-nonzeros-in-r
build.banded.gaus.mat <- function(w, tol = 1E-5){
  # TODO implement this in Rcpp !!!
  nrow <- length(w)
  # 1.1 Locate peak boundaries
  b.mat <- t(sapply(1:nrow, function(u){
    # Analytically solve for x: y = h * exp(-(x-u)^2 / (2*w^2)) where h = 1, y = tol, w = peak width, u = center
    # Solution: xs = u +/- w * sqrt(-2log(y/h))
    xb <- round(w[u] * sqrt(-2*log(tol))) # x border position where y == tol
    # Trim x if out of bound
    x.low <- max(u - xb, 1)
    x.up <- min(u + xb, nrow)
    return(c(x.low, x.up, x.up - x.low + 1))
  }))
  # 1.2 Preallocate vectors: el.i (row index), el.j (column index), and el (element values)
  el.i <- vector(mode = "integer", length = sum(b.mat[,3]))
  el.j <- vector(mode = "integer", length = sum(b.mat[,3]))
  el <- vector(mode = "double", length = sum(b.mat[,3]))
  start <- 1
  # 1.3 Fill the vectors
  for(b in 1:nrow(b.mat)){
    el.i[start:(start+b.mat[b,3]-1)] <- b.mat[b,1]:b.mat[b,2]
    el.j[start:(start+b.mat[b,3]-1)] <- b
    el[start:(start+b.mat[b,3]-1)] <- exp(-((b.mat[b,1]:b.mat[b,2])-b)^2/(2*w[b]^2)) # y = h * exp(-(x-u)^2 / (2*w^2))
    start <- start+b.mat[b,3]
  }
  mat <- sparseMatrix(i = el.i, j = el.j, x = el)
  return(mat)
}


# get.BIC ####
# Function that computes the BIC given fitted data (y) and test data (yhat)
# and the number of parameters of the fit. BIC can be fitted on a transformed
# scale.
# INPUT
#   y:      fitted data
#   yhat:   test data
#   npars:  number of parameters used to fit y
#   transf: a function giving the scale transformation
# OUTPUT
#   BIC:    the fit quality
get.BIC <- function (y, yhat, npars, transf = function (y) sqrt(y)) {
  n <- length(y)
  RSS <- sum((transf(y)-transf(yhat))^2) # residual SS on sqrt() scale
  min2LL <- n + n*log(2*pi) + n*log(RSS/n)
  BIC <- min2LL + log(n)*npars
  return(BIC)
}

# block.deconvolution ####
# Function to solve x: y = C x subject to x > 0, splitting the problem
# in smaller windows.
# INPUT
#   y:  the response vector
#   X:  the design matrix
#   method:   the method used for deconvolution. Available methods: "NNLS"
#   win.size: size of the window in which y is split. If NULL, no splitting.
#   ...:  parameters passed to the deconvolution method
# OUTPUT
#   x:  a vector containing the estimated non-negative coefficients
block.deconvolution <- function(y, X, method = "NNLS", win.size = NULL, ...){
  y.l <- length(y)
  if(is.null(win.size)) win.size <- y.l
  win.n <- ceiling(y.l/win.size)

  # Perform deconvolution in separate windows
  # Windows are update hierarchically, that is they are updated sequentially
  # with the next window fitted on the residuals from te previous window.
  coefs <- rep(0, ncol(X))
  for(i in 1:win.n){
    # 1. Set the current deconvolution window
    win.ind <- (1:win.size)+(i-1)*win.size
    win.ind <- win.ind[win.ind <= y.l]
    # 2. Subset the data
    A <- X[win.ind,win.ind]
    b <- y[win.ind]
    # 3. Compute residuals
    if(i > 1){ # there is no previous window for the first window
      win.prev <- (1:win.size)+(i-2)*win.size # "-2" to get the previous index
      b <- pmax(y[win.ind] - Matrix::tcrossprod(X[win.ind, win.prev], coefs[win.prev]), 0)
    }
    # 5. Perform deconvolution
    if(method == "NNLS"){
      coefs[win.ind] <- nnls(A = A, b = b)$x
    } else {
      # TODO add Tom's algorithm for best subset selection here
      stop("This deconvolution method is not implemented.")
    }
  }
  return(coefs)
}


# cluster.peak.coefficients ####
# Function to cluster contiguous coefficients obtained from the deconvolution of a signal
# with a banded matrix of shifted model peaks.
# INPUT
#   coefs:  the vector of coefficients estimated after deconvolution
# OUTPUT
#   coefs.new:  the clustered coefficients
cluster.peak.coefficients <- function(coefs){
  # TODO: allow to cluster in case of sub- or supersampling
  coef.start <- which(diff(sign(coefs)) == 1)+1
  coef.stop <- which(diff(sign(coefs)) == -1)
  if(max(coef.start) > max(coef.stop)) coef.stop <- c(coef.stop,length(coefs))
  if(min(coef.start) > min(coef.stop)) coef.start <- c(1, coef.start)
  if(length(coef.start)!= length(coef.stop)) stop("there is a bug")
  coefs.new <- rep(0,length(coefs))
  for(b in 1:length(coef.start)){
    x.c <- round(coef.start[b]+(coef.stop[b]-coef.start[b])/2) # the new peak location after clustering
    coefs.new[x.c] <- sum(coefs[coef.start[b]:coef.stop[b]]) # combine the peak heights by summing the clustered coefficients
  }
  return(coefs.new)
}


# Function that fits quantile regression along one dimension of a matrix
# given a profile from the other dimension.
# INPUT
#   method: the implementation to use: "quantreg" (adapted simplex algorithm),
#           "uqr.cpp" (univariate quantile regression optimized in C++)
#   M:  the matrix to regress
#   u:  the covariate vector
#   w:  the weights vector
#   tau:  the quantile to fit
#   rm.bas:   logical indicating whether a baseline should be removed before fitting
# OUTPUT
#   v: the estimate vector
# TODO ideally this should be accelerated in Rcpp
fit.profile <- function(method = "quantreg", ...){
  if(!method %in% c("quantreg", "uqr.cpp")) stop("Wrong method supplied")
  switch(method,
         quantreg = fit.profile.rq(...),
         uqr.cpp = fit.profile.uwqr(...))
  # Compare methods
  # quantreg = fit.profile.rq(...)
  # uqr.cpp = fit.profile.uwqr(...)
  # matplot(cbind(quantreg, uqr.cpp)+1, log = "y", type = "l", lty = 1, col = 2:1)
}


# fit.profile.rq ####
# Function that quantile regression along one dimension of a matrix
# given a profile from the other dimension. The solution is computed
# using the Barrodale and Roberts simplex algorithm (from "quantreg").
# INPUT
#   M:  the matrix to regress
#   u:  the covariate vector
#   w:  the weights vector
#   tau:  the quantile to fit
#   rm.bas:   logical indicating whether a baseline should be removed before fitting
#   ...:  ignored argument
# OUTPUT
#   v: the estimate vector
fit.profile.rq <- function(M, u, w = rep(1, length(u)), tau = 0.5, rm.bas = TRUE, ...){
  nz <- u != 0 & w != 0
  x <- as.matrix(u[nz]) * w[nz] # subset and apply weights
  v <- apply(M[nz,,drop = F], 2, function(y){
    if(rm.bas) y <- y - min(y) # subtract baseline
    y <- y * w[nz] # apply weights
    coef <- rq.fit.br(y = y, x = x, tau = tau)$coefficients
    return(coef)
  })
  return(as.vector(v))
}


# fit.profile.uwqr ####
# Function that calculates the univariate weighted quantile regression
# over the columns of the matrix M given the profile u and the weights w.
# INPUT
#   M:  the matrix to regress
#   u:  the covariate vector
#   w:  the weights vector
#   tau:  the quantile to fit
#   ...:  ignored argument
# OUTPUT
#   v: the estimate vector
# TODO ideally this should be accelerated in Rcpp
fit.profile.uwqr <- function(M, u, w = rep(1, length(u)), tau = 0.5, rm.bas = TRUE, ...){
  v <- uni_WQR_by_col(Y = M, x = u, w = w, tau = tau, rm_bas = rm.bas)
  return(as.vector(v))
}



# compute.wv ####
# Function that computes the weights associated to columns of the matrix
# M given an elution profile u. The weights are coputed based on the
# similarity betwen the columns of M and u on a log10 scale.
# INPUT
#   M:  the chromatogram matrix
#   u:  the row factor.
#   v:  the column factor. If supplied, w_i is set to 0 when v_i == 0.
#   power:  the weights are power transformed.
#   eps:  small value to avoid infinite values when applying log transf.
# OUTPUT
#   w:  the computed weights
compute.wv <- function(M, u, v = NULL, power = 1, eps = .Machine$double.eps){
  if(power == 0) return(rep(1, length(u)))
  w <- apply(M, 2, function(y){
    if(sum(log10(y + eps)) == 0) return(0)
    # TODO check the influence of eps because eps is also used for trimming u and those two eps have not the same meaning
    sim <- cosine.sim(Y = log10(y + eps), X = log10(u + eps))
    return(sim)
  })^power
  if(!is.null(v)) w[v == 0] <- 0
  w <- w/sum(w)*length(w)
  return(w)
}


# compute.wu ####
# Function that computes the weights associated to rows of the matrix
# M given an elution profile u.
# INPUT
#   u:  the row factor. If supplied, w_i is set to 0 when u_i == 0.
#   power:  the weights are power transformed.
# OUTPUT
#   w:  the computed weights
compute.wu <- function(u, power = 1){
  if(power == 0) return(rep(1, length(u)))
  w <- u^power
  w <- w/sum(w)*length(w)
  return(w)
}


#### pickLogConPeak ####
# Function that approximate the given signal (y) by a log-concave smooth signal. The
# function uses the cobs algorithm which is based on quatile regression
# INPUT
#   y:    the signal (vector) to approximate
#   tau:  the quantile of the signal to fit
#   weights:    vector of weights for the quantile regression
#   local.fit:  a logical indicating if the fitting range be bound by zero values
#   mode: a numeric giving the index of the expected mode. Used only when local.fit==TRUE.
#         The fit is performed on the part of the signal that contains the mode.
#   basfun: a function that fits a baseline that is subtracted before fitting.
#   Plot: a logical indicating if the fitting should be plotted
# OUTPUT
#   y.new:  the fitted signal
# TODO make this a generic function to fit different kind of peak shape constraints for instance: gaussian, exponetially modified gaussian, or any user defined paramteric shape
pickLogConPeak <- function(y, tau = 0.1, weights = NULL, local.fit = TRUE,
                           mode = NULL, basfun = function(x)return(0),
                           Plot = FALSE, ...){
  y <- pmax(y-basfun(y),0)
  y.new <- rep(0, length(y))
  if(sum(y) == 0) return(y.new)

  # Local Fit
  if(local.fit){
    starts <- pmax(which(diff(sign(c(0,y))) == 1),1) # before: 0, after: pos value
    stops <- pmin(which(diff(sign(c(y,0))) == -1),length(y)) # before: pos value, after: 0
    if(is.null(mode)){ # if no mode is supplied, take the chunk that has highest norm
      norm.ysub <- sapply(1:length(starts), function(l){
        return(sqrt(sum(y[starts[l]:stops[l]]^2)))
      })
      start <- starts[which.max(norm.ysub)]
      stop <- stops[which.max(norm.ysub)]

    } else {
      if(y[mode] == 0) return(y.new)
      start <- tail(starts[starts <= mode],1)
      stop <- head(stops[stops >= mode],1)
      if(length(start) * length(stop) == 0) return(y.new)
    }
  } else {
    start <- 1
    stop <- length(y)
  }
  y.sub <- y[start:stop]
  x <- start:stop

  # Weights
  if(is.null(weights)){
    offs <- max(y.sub)*1E-4
    weights <- sqrt(y.sub+offs)
  } else {
    weights <- weights[start:stop]
  }
  weights <- weights/sum(weights)*length(weights)

  # Knots
  nknots <- round(length(y.sub)/3)
  if(nknots <= 3) return(y.new)

  # Log-concave fit using cobs
  require(cobs)
  offs <- max(y.sub)*1E-6
  fit <- suppressWarnings(cobs(x=x,y=log10(y.sub+offs),
                               constraint= "concave", lambda=0,
                               knots = seq(min(x),max(x),length.out = nknots),
                               nknots=nknots, knots.add = FALSE, repeat.delete.add = FALSE,
                               w=weights,
                               tau=tau, print.warn = F, print.mesg = F,
                               ...))
  y.fit <- fit$fitted
  y.fit <- pmax(10^y.fit - offs,0)
  y.fit[is.infinite(y.fit)] <- 0 # Happens when cobs does not converge and that the fit is exponentiated
  y.new[start:stop] <- y.fit

  if(Plot){
    plot(y+1, log = "y", type = "l")
    lines(y.new+1, col = 4)
    lines(x = start:stop, y = weights*max(y.sub)/max(weights)+1, col = "grey")
    abline(v = c(start, stop), lty = 2)
    abline(v = mode, col = 2)
  }
  return(y.new)
}


# plot.compound.fit ####
# Function that plots the elution profile and spectrum profile along the input data
# and the fitted data. Profiles are shown as a line plot, where as data and fit are
# shown as heatmaps (low: black -> blue -> green -> yellow -> red: high).
# INPUT
#   M:    the input data matrix
#   u:    the elution profile
#   v:    the spectrum profile
#   scantimes:  the time scale associated to the rows of M
#   mzvalues:   the mass scale associated to the columns of M
#   iter: the iteration during which the profiles were fit
#   gam:  the gamma correction applied to the heatmaps
# OUTPUT
#   NULL
plot.compound.fit <- function(M, u, v, u.raw = NULL,
                              scantimes, mzvalues,
                              iter = 0, gam = 0.1){
  par(mfrow=c(2,2), oma = c(0,0,2,0), mar = c(3,3,2.5,2.5), mgp = c(2,1,0))
  zlims <- c(0, max(M))^gam
  plotmat(M, subsamp = 1, gamma = gam, plot.axes = F, zlim = zlims,
          main = "Input data", xlab = "Time", ylab = "Mass")
  if(!is.null(u.raw)){
    plot(x = scantimes, y = u.raw+1, log = "y", pch = 16, col = "grey", ylim = c(1, max(u.raw)),
         main = "Estimated elution profile", xlab = "Time (min)", ylab = "Intensity")
    lines(x = scantimes, y = u+1, col = "red3")
  } else {
    plot(x = scantimes, y = u+1, log = "y", type = "l", col = "red3", ylim = c(1, max(u)),
         main = "Estimated elution profile", xlab = "Time (min)", ylab = "Intensity")
  }
  plotmat(Matrix::tcrossprod(u,v), subsamp = 1, gamma = gam, plot.axes = F, zlim = zlims,
          main = "Fitted data", xlab = "Time", ylab = "Mass")
  suppressWarnings(plot(x = mzvalues, y = v, log = "y", type = "h", col = "blue4",
                        main = "Estimated spectrum profile", xlab = "Mass", ylab = "Intensity"))
  # suppressWarnings(plot(x = mzvalues, y = v, log = "y", type = "h", col = ifelse(wv > sort(wv, decreasing = T)[20], "green4", "red3"),
  #                       main = "Estimated spectrum profile", xlab = "Mass", ylab = "Intensity"))


  title(paste0("Estimation of the elution and spectrum profiles after iteration ", iter), outer = T)
}


# optimize.profiles ####
# Function that iteratively updates a given u and v profile from the data
# M.all, which can be a full chromatogram or several chromatograms in a
# row-augmented matrix. u and v are fitted using quantile regression, and
# the fitting window is temporaly constrained where u != 0, and is
# adaptively update between iterations. The u profile is constrained to be
# log-concave.
# INPUT
#   M.all:  the (row-augmented) matrix
#   u:      the initial elution profile to update (should not contain zeros.)
#   v:      the initial spectrum profile to update
#   win:    the row index of M.all from which u was fit
#   peak.mode:  the row index of the expected mode of u
#   m:      the number of rows of 1 sample in M.all.
#   tau:    the tau to use for quantile regression when fitting u and v
#   cobs.tau:   the tau to use for the log concave approximation with cobs
#   qr.method:  the algorithm used to fit the quantile regression model.
#               Simplex method ("quantreg"), univariate quantile regression ("uqr.cpp")
#   power.u:    the power for u weight transformation. Set to 0 for no weighting.
#   power.v:    the power for v weight transformation. Set to 0 for no weighting.
#   eps:    small value to avoid infinite values when computing weights.
#   basfun.u:   a function that estimates a baseline to remove after
#               fitting  u
#   max.expand: the maximum number of indices that the window can be
#               extended per iteration
#   maxiter:    the maximum number of iterations
#   tol:    a cutoff that stops the iteration when the relative difference
#           in v between two consecutive iterations is lower that the cutoff
#   verbose:    logical indicating whether to print iteration progress to console
#   debug:      logical indicating whether to pause between every iteration
#   Plot:   logical indicating whether to plot update after every iteration
# OUTPUT
#   u:    the updated u
#   v:    the updated v
#   wu:   the updated u weights
#   wu:   the updated v weights
#   win:  the updated row index
#   success:  a logical indicating whether the profiles where successfully
#             update or if one of the profiles is all zero.
optimize.profiles <- function(M.all, u, v,
                              scantimes, mzvalues,
                              win, peak.mode, m,
                              tau = 0.1, cobs.tau = 0.1, qr.method = "quantreg",
                              power.u = 0, power.v = 0,
                              eps = .Machine$double.eps,
                              basfun.u = function(x)return(0),
                              max.expand = 300,
                              maxiter = 10, tol = 1E-2,
                              verbose = FALSE, debug = FALSE, Plot = FALSE){
  samp <- ceiling(min(win)/m)
  j <- 1
  delta <- vector(length = maxiter)
  success <- TRUE
  while(j <= maxiter && (j == 1 || delta[j-1] > tol)){
    t1 <- proc.time()[3]

    if(verbose){
      cat(paste0("\rIteration ", j, " out of ",maxiter, "   "))
      flush.console()
    }
    M <- M.all[win,]
    vp <- v
    up <- u

    # 1. Update time weights
    # This is meant to put more weights on local time scans
    wu <- compute.wu(u = u, power = power.u)

    # 2. Update v
    v <- fit.profile(M = M, u = u, v.init = v, w = wu, tau = tau, method = qr.method, rm.bas = TRUE)
    if(sum(v) == 0) {
      success <- FALSE
      break()
    }
    v <- v/norm(v, type = "2")

    # 3. Update v weights
    wv <- compute.wv(M = M, u = u, power = power.v, eps = eps)

    # 4. Update u
    u <- u.raw <- fit.profile(M = t(M), u = v, v.init = u, w = wv, tau = tau, method = qr.method, rm.bas = TRUE)
    # Impose unimodality using cobs (concavity on log scale)
    u.constr <- pickLogConPeak(y = u, mode = peak.mode, tau = cobs.tau, weights = wu, basfun = basfun.u)
    if(sum(u.constr) == 0) {
      success <- FALSE
      break()
    }

    if(j <= maxiter/2){ # TODO implement an automatic criterion based on the rate change between win.old and win
      # 5. Adaptively relocate window
      win.old <- win
      u.ext.l <- u.ext.r <- c()
      if(u.constr[1] != 0 & min(win) - (samp-1) * m > 1){ # if the elution profile is interrupted by the window size, extend fit
        win.expand <- (min(win)-1):max(min(win)-max.expand, (samp-1)*m+1)
        M.exp <- M.all[win.expand,]
        for(l in 1:length(win.expand)){
          y <- M.exp[l,]
          coef <- fit.profile(M = as.matrix(y), u = v, w = wv, tau = tau, method = qr.method, rm.bas = FALSE)
          u.ext.l <- c(coef, u.ext.l)
          if(coef == 0){ break() }
        }
        u.start <- 1
        win.start <- win.expand[l]
      } else {
        u.start <- which.max(u.constr != 0)
        win.start <- min(win) + u.start - 1
      }
      if(u.constr[length(u.constr)] != 0 & max(win) - (samp-1)*m < m){ # if the cobs fit is interrupted by the window size, extend fit
        win.expand <- (max(win)+1):min(max(win) + max.expand, ceiling(max(win)/m)*m)
        M.exp <- M.all[win.expand,]
        for(l in 1:length(win.expand)){
          y <- M.exp[l,]
          coef <- fit.profile(M = as.matrix(y), u = v, w = wv, tau = tau, method = qr.method, rm.bas = FALSE)
          u.ext.r <- c(u.ext.r,coef)
          if(coef == 0){ break() } # When a 0 is included, stop expanding
        }
        u.stop <- length(u)
        win.stop <- win.expand[l]
      } else {
        u.stop <- length(u) - which.max(rev(u.constr) != 0) + 1
        win.stop <- min(win) + u.stop - 1
      }
      u <- u.raw <- c(u.ext.l, u[u.start:u.stop], u.ext.r)
      u <- pickLogConPeak(y = u, mode = which.max(u), tau = cobs.tau, weights = u^2, basfun = basfun.u) # TODO weighting not consistent with rest of code
      if(sum(u) == 0) {
        success <- FALSE
        break()
      }

      # 6. Update the variables depending on the window index
      win <- win.start:win.stop
      win <- win[u != 0] # Trim the zeros in u (happens after cobs fit when extending window)
      u.raw <- u.raw[u != 0]
      u <- u[u != 0]
      peak.mode <- which.max(u)
      up.new <- rep(0, length(u))
      up.new[win %in% win.old] <- up[win.old %in% win]
      up <- up.new
    } else {
      u <- u.constr
    }

    # 7. Plot
    if(Plot){
      plot.compound.fit(M = M.all[win,], u = u, v = v, u.raw = u.raw,
                        scantimes = scantimes[win - (ceiling(win/m)-1)*m], mzvalues = mzvalues, iter = j)
    }

    # 8. Update error
    # Convergence criterion: continue until the relative change of v with
    # previous iteration is smaller than tol
    delta[j] <- norm(vp-v, "2")/norm(vp, "2")
    j <- j+1
  }
  if(success){ # Control whether a peak could be extracted
    if(verbose){ cat(paste0(" (", ifelse(j>= maxiter, "Max iteration reached", paste0("Converged with tolerance ", tol)),").\n\n"))}
  } else {
    if(verbose){ cat("\n\n") }
  }
  wu <- compute.wu(u = u, power = power.u)

  if(debug){ if(readline(prompt = "Press [enter]") != "") stop() }
  return(list(u = u, v = v,
              wu = wu, wv = wv,
              win = win,
              success = success))
}



# fit.lin.bas ####
# Fit a linear baseline through every column (dim=2) or row (dim=1) of
# the matrix M. The fitting considers only 2 edge points.
# Baseline: y = ax + b
# INPUT:
#   M:    the matrix to fit the linear baseline
#   dim:  should the baseline be fitted column (dim = 2) or row (dim = 1) wise?
#   x1:   first edge
#   x2:   second edge
# OUTUT:
#   bas:  a matrix of size M, containing the linear baseline
fit.lin.bas <- function(M, dim = 2, ext = 0){
  x1 <- 1
  x2 <- dim(M)[-dim]
  # TODO put this in cpp
  bas <- apply(M, dim, function(y){
    y1 <- mean(y[x1:(x1+ext)]) # TODO needed ?
    y2 <- mean(y[(x2-ext):x2]) # TODO needed ?
    # y1 <- y[x1]
    # y2 <- y[x2]
    b <- rep(0, length(y))
    a <- (y2-y1)/(x2-x1)
    b[x1:x2] <- (x1:x2) * a + y[x1] # b = y1

    # x <- cbind(1, 1:length(y))
    # b <- x %*% rq.fit.br(x = x, y = y, tau = 0.1)$coefficients

    return(pmax(b,0))
  })
  return(bas)
}


# scale.factors ####
# Function that scales the u factor using median regression based on initial values of u and v.
# INPUT
#   M:  matrix from which to scale the u factor
#   u:  the initial u factor
#   v:  the initial v factor
#   U.p:  a matrix containing additional u factors. This can be supplied in order to account for
#         overlapping factors
#   V.p:  a matrix containing additional v factors. This can be supplied in order to account for
#         overlapping factors
#   win:  a vector of indices corresponding to the on to subset the data to compute the scale. If supplied,
#         the length of win must match the length of u
#   bas:  a matrix containing a baseline to ignore from M. If supplied, the dimensions of bas must
#         match the dimensions of M
#   subsamp:  a numeric factor indicating how much subsampling can be performed when scaling the factors
# OUTPUT
#   u:  the scaled factor u
# NOTE
# M can be supplied as a bigmemory object.
scale.factors <- function(M, u, v, U.p = NULL, V.p = NULL, Plot = FALSE,
                          win = NULL, bas = NULL, subsamp = 1){

  # 0. Argument checks
  if(is.null(win)){
    win <- 1:length(u)
  } else {
    if(length(u) != length(win)) stop(paste0("length of u (", length(u), ") and length of win (", length(win), ") must match."))
  }
  if(is.null(bas)){
    bas <- Matrix(0, nrow = nrow(M), ncol = ncol(M), sparse = TRUE)
  } else {
    if(!all(dim(bas) == dim(M))) stop("bas and M must have the same size")
  }
  # 1. Subset the data to reduce the problem size
  u.temp <- u
  subs <- seq(1, length(u), by = subsamp)
  u <- u[subs]
  win <- win[subs]

  # 2. Build the covariate matrix
  if(!is.null(U.p) & !is.null(V.p)){
    # Select factors that overlap with the current factor and subset elements of V that are non zero in at leat 1 factor
    comp.sel <- which(Matrix::colSums(U.p[win,]) > 0) # select the overlapping compounds
    nzero.v.ind <- which(Matrix::rowSums(cbind(v, V.p[,comp.sel])) > 0) # get the element indices where there is a non zero element in at least 1 factor
    # Create the covariate matrix
    cov.mat <- matrix(as.vector(Matrix::tcrossprod(u, v[nzero.v.ind])), ncol = 1) # The first covariate is the newly fitted factors
    cov.add <- sapply(comp.sel, function(l){ return(as.vector(Matrix::tcrossprod(U.p[win,l], V.p[nzero.v.ind,l]))) })
    cov.mat <- cbind(cov.mat, unlist(cov.add))
  } else {
    nzero.v.ind <- which(v > 0)
    cov.mat <- matrix(as.vector(Matrix::tcrossprod(u, v[nzero.v.ind])), ncol = 1) # If U.p and V.p are not supplied, scale u with no other covariate.
  }
  # 3. Scale the factors using median regression
  y <- as.vector(pmax(M[win, nzero.v.ind] - bas[win, nzero.v.ind], 0)) # work on the original data, ignore points were fit is 0 for all factors
  # TODO put weights ?
  sigmas <- rq.fit(y = y, x = cov.mat, tau = 0.5)$coefficients # Update the scale only for the current profile
  # sigmas <- nnls(b = y, A = cov.mat)$x
  sigma <- max(sigmas[1], 0)

  # Note: NNLS is too sensitive when there is data from overlapping compounds that are not yet extracted (overestimation of sigma)
  u <- u.temp * sigma

  if(Plot){ # TODO need to implement plotting ?
    stop("This is useless...")
    ext <- 50
    M.plot <- M[win,]
    bas.plot <- bas[win,]

    par(mfrow = c(3,3))
    for(k in which(v != 0)){
      plot(M.plot[,k]+1, log = "y", main = k)
      lines(bas.plot[,k]+1, col = 2)
      lines(u*v[k]+1, col = "green4")
    }
  }
  return(u)
}


# update.heights ####
# Function that updates the supplied peak heights from the univariate signal Y using a dictionnary of peak shapes.
# INPUT
#   peak.locs:  matrix containing a column "max" that is the peak location index with respect to Y, and a column
#               "height" holding the current peak heights
#   Dic:    a dictionnary of peak shapes. The number of column should equal the length of Y or m if supplied
#   Y:      the univariate signal to approximate
#   win:    a vector of index on which to subset
#   m:      the length of a sample. It can me smaller than the length of Y in case that Y contains several samples
#   method:     the method used to estimate the new peak heights. Implemented methods: "NNLS"
# OUTPUT
#   height.df:  matrix similar to peak.locs with two additional columns: "height.new" contains the updates peak
#               heights and "ratio" contains the ratio between "height" and "height.new".
# TODO clean this code, this is too messy !
update.heights <- function(peak.locs, Dic, Y,
                           win = 1:length(Y), m = length(Y), method = "NNLS"){
  # 1. Get peak location that might have changed
  peak.sel <- peak.locs[,"max"] %in% win
  if(sum(peak.sel) == 0) return(NULL)
  model.in.win <- which(Matrix::colSums(Dic[win-(ceiling(win/m)-1)*m,,drop=F]) != 0) # position index of all model peak shapes that have at least 1 non zero signal
  peak.max.samp <- peak.locs[, "max"] - (peak.locs[, "samp"]-1)*m # position index of all peak contained in the current sample
  peak.to.fit <- model.in.win[model.in.win %in% peak.max.samp]

  # 2. Update peak heights using NNLS
  y <- Y[win]
  X <- Dic[win - (ceiling(win/m)-1)*m, peak.to.fit, drop = F]
  if(method == "NNLS") {
    coefs <- nnls(b = y, A = X)$x
  } else {
    stop(paste0("Method \"", method, "\" is not supported." ))
  }
  coefs <- coefs[peak.to.fit %in% (peak.locs[peak.sel, "max"] - (peak.locs[peak.sel, "samp"]-1)*m)] # subset coefs to keep only the peaks that must be updated

  # 3. Compute ratio between old and new heights
  height.df <- peak.locs[peak.sel,, drop = F]
  height.df <- height.df[order(height.df[,"max"]),,drop = F]
  height.df <- cbind(height.df, height.new = coefs)
  height.df <- cbind(height.df, ratio = height.df[,"height.new"]/height.df[,"height"])

  return(height.df)
}


# update.detection ####
# Function that updates the peak detection based on the updated peak
# coefficients.
# TODO describe wat is the final criterion that is used
# INPUT
#   Y:      the univariate signal to approximate
#   Dic:    a dictionnary of peak shapes. The number of column should equal the
#           length of Y or m if supplied
#   peaks:  matrix containing a column "max" that is the peak location index
#           with respect to Y, a column "height" holding the current
#           peak heights, and a column "samp" indicating in which sample the
#           peak is located (in case m < length(Y))
#   cutoff:   ....
# TODO document cutoff
#   sel:    index vector indicating for which peaks should the detection be
#           updated
#   m:      the length of a sample. It can me smaller than the length of Y in
#           case that Y contains several samples
# OUTPUT
#   TODO document output
update.detection <- function(Y, Dic, peaks, cutoff, sel = NULL, m = length(Y)){
  if(is.null(sel)) sel <- 1:nrow(peaks)
  # TODO worth implementing in Rcpp? but watch out the qunatreg and Matrix functions !!!
  detect.df <- do.call(rbind, lapply(sel, function(k){
    # Extract peak information
    row <- peaks[k,]
    if(row["lives"] == 0 && row["active"] == 0) return(row)
    old.active <- row["active"]
    samp <- row["samp"]
    peak.pos <- row["max"] - (samp-1)*m
    peak.shape <- Dic[, peak.pos]
    win <- which(peak.shape > 1E-3) # Focus only on the significant part of the shape
    peak.shape <- peak.shape[win]
    # Fit a baseline
    y <- log10(Y[win + (samp-1)*m]+1) # Fit baseline on a log scale
    x.bas <- cbind(1, 1:length(y)) # The baseline is made of a slope + intercept
    bas <- x.bas %*% suppressWarnings(rq.fit.br(y = y, x = x.bas, tau = 0.5)$coefficients) # warnings:  In rq.fit.br(y = y, x = x.bas, tau = 0.1) : Solution may be nonunique
    # Compute the filtering criterion
    row["crit"] <- log10(row["height"])/bas[which.max(peak.shape)] # Compute the ratio between baseline and peak signal at max
    row["active"] <- as.numeric(row["crit"] > cutoff)
    row["lives"] <- ifelse(row["active"] - old.active == 1, row["lives"]-1, row["lives"]) # remove a "life" when peak gets reactivated
    return(row)
  }))
  return(detect.df)
}


# get.peak.locations ####
# Function to identify peaks from the total ion chromatigram (TIC) of one
# or multiple samples. First, an average peak width is optimized. The peak
# width is then used to build a banded matrix of shifted Gaussians peak
# shapes which is used with NNLS to decompose the TIC (per sample). The
# NNLS coefficients are then formated as peak location information.
# INPUT
#   TICs:   a vector containing the concatenated TIC of one or multiple samples
#   m:      the length of the TIC of one sample. m is assumed identical across
#           samples and must be a multiple or identical to the length of TICs
#   samps:  sample index from which to get the peak locations
#   ref.samp:   the sample index from which to optimize the peak width
#   win.size:   size of the window for the peak width estimation and NNLS fitting
#   peak.width.range:   a range of peak widths where the average peak width is
#                       expected to lie.
#   peak.trim:  the gaussian peak threshold under which the peak is trimmed to 0.
#   cutoff:     ...
# TODO document cutoff
#   lives:    a numeric indicating up to how many times a peak can be extracted.
#             if a peak should be extracted only once, set "lives = 1".
#   ncores:   the number of cores to use for parallelization
#   plot.file:  a character string indicating the pdf file where to plot the results
#   verbose:    logical indicating whether to print progress to console
# OUTPUT
#   peak.locs:  a matrix containing the identified peak heights and locations of the
#               peak maximums
#   bandC:      a banded matrix containing shifted Gaussians with the optimized peak
#               widths
# NOTE
# Peak width is defined as the full width at half maximum
get.peak.locations <- function(TICs, m,
                               samps = NULL,  ref.samp = samps[1],
                               win.size = 300, peak.width.range = c(1,20),
                               peak.trim = 1E-5, cutoff = 0.9,
                               lives = 3,
                               ncores = 1,
                               plot.file = NULL, verbose = FALSE){
  if(is.null(samps)) samps <- 1:(length(TICs) / m)

  # TODO also allow actual peak shape in function of RT to be estimated? or
  # allow other functional forms like exp mod gaussian?
  if(verbose) cat("Extracting peak locations\n\n")

  ## Estimate an average peak width from a reference sample
  if(verbose) cat("Estimate an average peak width from the reference sample\n\n")
  peak.widths <- fit.average.peak.width(y = TICs[(1:m)+m*(ref.samp-1)],
                                        win.size = win.size,
                                        peak.width.range = peak.width.range) # 6.6 s for win.size = 300; refsamp = 9; peak.width.range = c(1,20)
  peak.width.fun <- peak.widths$peak.width.fun
  peak.widths <- peak.widths$peak.width.fun$fitted.values

  ## Fit data using the average peak width per window
  if(verbose) cat("Deconvolution of the TIC using a banded matrix with \nshifted Gaussian peak shapes (with the optimized peak width)\n\n")
  peak.widths.preds <- predict(peak.width.fun, newdata = list(x = 1:m))
  # Generate a sparse matrix containing shifted gaussian peaks with the predicted peak widths
  bandC <- build.banded.gaus.mat(w = peak.widths.preds, tol = peak.trim)
  # Deconvolute TIC with the banded matrix
  if(ncores > 1){
    # Prepare the cluster
    cl <- makeCluster(ncores,type="SOCK") # on linux use type="FORK" (faster), but not available on windows
    doSNOW::registerDoSNOW(cl)
    snow::clusterExport(cl, list=c("block.deconvolution", "cluster.peak.coefficients","m", "TICs", "bandC", "win.size"), envir = environment()) # export values that each worker needs
    # Excute computation
    detected.gaus.peak <- clusterApply(cl, samps, function(i) {
      library(Matrix)
      library(nnls)
      peak.coefs <- block.deconvolution(y = TICs[(1:m)+m*(i-1)], X = bandC, win.size = win.size)
      # Cluster contiguous coefficients
      peak.coefs <- cluster.peak.coefficients(coefs = peak.coefs)
      return(peak.coefs)
    })
    stopCluster(cl)
  } else {
    detected.gaus.peak <- lapply(samps, function(i){
      # Deconvolution of the current sample TIC using the banded matrix of gaussians
      peak.coefs <- block.deconvolution(y = TICs[(1:m)+m*(i-1)], X = bandC, win.size = win.size)
      # Cluster contiguous coefficients
      peak.coefs <- cluster.peak.coefficients(coefs = peak.coefs)
      return(peak.coefs)
    })
  }

  ## Format the extracted peak locations
  if(verbose) cat("Format the extracted peak locations \n\n")
  peak.locs <- do.call(rbind, lapply(1:length(samps), function(i){
    peak.max <- which(detected.gaus.peak[[i]] != 0)
    peak.height <- detected.gaus.peak[[i]][which(detected.gaus.peak[[i]] != 0)]
    return(cbind(samp = samps[i],
                 max = peak.max + m*(samps[i]-1),
                 height = peak.height,
                 crit = 0,
                 lives = lives,
                 active = 0))
  }))

  # Perform peak filtering
  if(verbose) cat(paste0("Filtering peaks\n",
                         "Initial number of peaks:", nrow(peak.locs), "\n"))

  peak.locs <- update.detection(Y = TICs, Dic = bandC, peaks = peak.locs,
                                sel = NULL, cutoff = cutoff, m = m)
  if(verbose) cat(paste0("Selected peaks:", sum(peak.locs[,"active"]), "\n\n"))

  ## Plot results
  if(is.character(plot.file)){
    if(verbose) cat("Plot results\n\n")
    pdf(file = plot.file, height = 8, width = 8)
    for(samp in samps){
      par(mfrow = c(1,1))
      plot(1, col = "white", axes = F, xlab = "", ylab = "")
      text(1, 1, labels = paste0("Gaussian deconvolution \nfor sample ", samp), cex = 2.5)
      TIC <- TICs[1:m+m*(samp-1)]
      win.n <- ceiling(length(TIC)/win.size)
      par(mfrow = c(2,2))
      for(i in 1:win.n){
        # Get the current window
        win.ind <- (1:win.size)+(i-1)*win.size
        win.ind <- win.ind[win.ind <= length(TIC)]
        # Find the peaks that have signal in the window
        peak.inds <- which(peak.locs[,"max"] %in% (Matrix::which(Matrix::colSums(bandC[win.ind,]) != 0) + m * (samp- 1)))
        peak.sub <- peak.locs[peak.inds,,drop=F]
        peak.pos <- peak.sub[,"max"] - m * (samp- 1)
        # Get associated coefficients
        coefs <- peak.sub[,"height"]
        # Subset the peak shape dictionary
        X <- bandC[win.ind, peak.pos]
        X <- sweep(X, 2, coefs, "*")
        # Plot results
        plot(TIC[win.ind]+1, pch = 16, col = "grey60", log = "y", cex = 0.5,
             ylim = c(min(coefs)/10, max(TIC[win.ind])), ylab = "Intensity",
             main = paste0("Window: ", i))
        points(x = peak.pos - (i-1) * win.size + 1, y = coefs, type = "h", col = c("red2", "green4")[peak.sub[,"active"]+1])
        matlines(X, col = "orange2", lty = 1)
        lines(Matrix::rowSums(X), lty = 1, col = 2)
        text(x = peak.pos - (i-1) * win.size + 1, y = coefs, labels = round(peak.sub[,"crit"], 3), cex = 0.5)
      }
    }
    dev.off()
  }

  ## Sort peaks from large to small
  if(verbose) cat("Sorting peaks\n\n")
  peak.locs <- peak.locs[order(peak.locs[,"active"], peak.locs[,"height"], decreasing = T),]
  if(verbose) cat("===================================\n\n")
  return(list(peak.locs = peak.locs,
              peak.with.optim = peak.width.fun,
              bandC = bandC))
}


# plot.peak.identification ####
# Function that plot the identified peaks from the signal y given the banded
# matrix used for the fitting.
# INPUT
#   y:    signal vector that is decomposed
#   Dic:  the banded matrix used for the decomposition containing the peak
#         shapes
#   peaks:  the information matrix of the identified peaks
#   win:  the window index in y to focus on
#   m:    length of a sample, in case y is split in different samples
#   scantimes:  the time labels associated to the rows of Dic
#   stick.cols:   names of a column to use for coloring for the sticks (up to 2
#                 colors)
#   labels:   name of a column to use as a label to plot on top of every peak
#             shape
#   ...:  further arguments passed to plot()
plot.peak.identification <- function(y, Dic, peaks, win,
                                     m = length(y), scantimes = NULL,
                                     stick.cols = NULL, labels = NA, ...){
  if(length(y) %% nrow(Dic) != 0) stop("length(y) should be a mutliple of nrow(Dic)")
  plot(y[win], log = "y", col = "grey", pch = 16,
       ylim = c(1, max(y[win])), xaxs = "i", cex = 0.5, xaxt = "n",
       ylab = "Signal", xlab = ifelse(is.null(scantimes), "Index", "Retention time (min)"), ...)
  if(is.null(scantimes)) scantimes <- 1:nrow(Dic)
  sep.val <- c(1,which(diff(ceiling(win/m)) != 0)+1, length(win))
  sep.samp <- ceiling(win[sep.val]/m)
  axis(at = sep.val, labels = rep("", length(sep.val)), side = 1)
  for(i in 2:length(sep.val)){
    start <- sep.val[i-1]
    incr <- round((sep.val[i]-start)/4)
    at <- start + (1:3) * incr
    axis(at = at, side = 1, mgp = c(0,0.5,0),
         labels = round(scantimes[win[at] - m*(ceiling(win[at]/m)-1)], 2))
  }
  if(nrow(peaks) == 0) return()
  peaks <- peaks[peaks[,"max"] %in% win,,drop = F]
  for(i in 1:nrow(peaks)){
    row <- peaks[i,]
    peak.shape <- Dic[win[ceiling(win/m) == row["samp"]] - m*(row["samp"]-1), row["max"]-m*(row["samp"]-1)]
    if(is.null(stick.cols)){
      col <- rep("green4", nrow(peaks))
    } else {
      col <- ifelse(row[stick.cols], "green4", "red3")
    }
    lines(x = which(ceiling(win/m) == row["samp"]),
          y = peak.shape * row["height"] + 1,
          col = "orange2")
    lines(x = which(win == row["max"]),
          y = row["height"] + 1,
          col = col, type = "h")
    text(x = which(win == row["max"]), y = row["height"] + 1, labels = round(row[labels], 3), cex = 0.75)
  }
  abline(v = sep.val, col = 1)
}



# sequential.deconvolution ####
# TODO document the function description
# INPUT
# TODO order and finish documenting the inputs
#   maxiter:  Number of inner iterations for optimizing the spectrum profile
#   tol:      if the relative difference of u and v with the previous iteration
#             is lower than tol, break the inner loop
#   verbose:  Print intermediate info to console
#   Plot:     Plot intermediate info
#   Plot.to.file:   Plot the results to pdf file instead of console
#   debug:    Pause in code until user response
#   tau:      The quantile to fit in quantreg
#   cobs.tau:     The quantile to fit when approximating the elution profile
#                 with log concave profile
#   cutoff
# TODO document cutoff
#   eps:      small constant to remove numerical instabilities and unecessary
#             small values
#   basfun.u:     function that estimates a baseline to remove before
#                 constraining u. Note: for no baseline subtraction, function(x) return(0)
#   sel.power:    selectivity factor. Increasing this factor puts more weights
#                 on selective ions. Set to 0 for no weighting
#   t.power:  the power tranformation for time weights. Increasing this value
#             will put more weights on time points close to the current mode of
#             the leution profile. Set to 0 for no weighting
# TODO allow win.ext to be asymmetric
#   win.ext:  the amount of time scans to add on both sides to the window
#             (optimized in 1 sample) in order to take into account
#             concentration shifts for all samples. Note: win.ext=100, will
#             extend ~20s on both sides
#   scale.subsamp:  the amount of under sampling performed during compound
#                   scaling.
# OUTPUT
# TODO document output
# s.h
sequential.deconvolution <- function(resids, data = NULL,
                                     m, n, s, samps,
                                     scantimes, mzvalues,
                                     TICs, peak.locs, C.dic,
                                     stop.thres = 0,
                                     tau = 0.1,
                                     cobs.tau = 0.1,
                                     eps = max(TICs)*1E-10,
                                     maxiter = 10, tol = 0.01,
                                     cutoff = 1,
                                     basfun.u = function(x) return(min(x)),
                                     sel.power = 2, t.power = 2,
                                     win.ext = 25, scale.subsamp = 10,

                                     Plot = FALSE, Plot.to.file = FALSE,
                                     plot.folder = "./",
                                     verbose = FALSE, debug = FALSE){
  if(is.null(data)) data <- resids

  # Initialize C and S
  C <- Matrix(0, nrow = m*s, ncol = nrow(peak.locs), sparse = T)
  S <- Matrix(0, nrow = n, ncol = nrow(peak.locs), sparse = T)

  # Create pdf
  plot.file <- NULL
  if(Plot.to.file){
    Plot <- TRUE
    plot.file <- paste0(plot.folder, "deconvolution_", Sys.Date(), "_", gsub("^.*[ ]|[:]", "", as.character(Sys.time())), ".pdf")
    pdf(file = plot.file, width = 10, height = 10)
  }

  is.active <- as.logical(peak.locs[,"active"])
  n.active <- sum(is.active)
  max.start <- max.cur <- max(peak.locs[is.active,"height"])
  t.tot <- proc.time()[3] # Time measure
  i <- 1
  while(max(peak.locs[peak.locs[,"active"] == 1, "height"]) > max.start * stop.thres &&
        sum(peak.locs[,"active"]) != 0 && i <= ncol(C)) {

    t1 <- proc.time()[3] # Time measure
    if(verbose) cat(paste0("\n====== Compound ", i, " ======\n\n"))

    # 1. Initialization

    # Initialize peak location, deconvolution window, and elution profile
    peak.df <- peak.locs[1,]
    peak <- peak.df["max"]
    samp.cur <- peak.df["samp"]
    u <- C.dic[,peak - (samp.cur-1)*m, drop = F]
    win <- Matrix::which(u != 0) + (samp.cur-1)*m
    u <- u[Matrix::which(u!=0)]
    M <- resids[win,]

    # Initialize weights
    # Mass weights are proportional to the similarity between u and a mass trace
    wv <- compute.wv(M = M, u = u, v = NULL, power = sel.power, eps = eps)
    wu <- compute.wu(u = u, power = t.power)

    # Fit the initial spectrum profile
    v <- fit.profile(M = M, u = u, w = wu, tau = tau, method = "quantreg", rm.bas = TRUE)
    nv <- norm(v,"2")
    u <- u*nv
    v <- v/nv

    # Plot intialization
    if(Plot & debug){
      par(mfrow = c(1,1))
      plot(1, axes = F, pch = "", xlab = "", ylab = "")
      text(x = 1, y = 1, labels = paste0("Fitting compound ", i), cex = 3)
      plot.compound.fit(M = M, u = u, v = v, iter = 0, scantimes = scantimes[win - (ceiling(win/m)-1)*m], mzvalues = mzvalues)
    }
    if(verbose){ cat(paste0("Sample: ", peak.df["samp"], "\n",
                            "Peak height: ", peak.df["height"], "\n",
                            "Peak location: ", peak, " (rowaugmented index) - ", peak - (ceiling(peak/m)-1)*m, " (sample index)", "\n",
                            "\n"))}

    # 2. Select sample to focus on

    # Extend window by a bit, because the concentration shifts could lead to
    # small peak eluting outside of the highly concentrated peak window
    win.start <- max(min(win - (ceiling(win/m)-1)*m) - win.ext, 1)
    win.stop <- min(max(win - (ceiling(win/m)-1)*m) + win.ext, m)
    win <- win.start:win.stop
    samp.ind <- rep(samps, each = length(win))
    win <- win + m * (samp.ind-1) # extend window for all samples (augmented row index)

    # Fit elution profile in all samples using quantreg
    M <- resids[win,] # TODO check if this can be done in memory for >100 samples
    u <- u.raw <- fit.profile(M = t(M), u = v, w = wv, tau = tau, method = "quantreg", rm.bas = TRUE)
    # Adjust the expected peak location to the sub window index
    sub.peak <- which(win == peak)
    sub.peak <- sub.peak - (ceiling(sub.peak/length(win)*length(samps))-1)*length(win)/length(samps)
    # sub.peak <- which(win == peak) %% (length(win)/length(samps)) # Adjust the peak location to the sub window
    # Constrain elution profiles
    u <- u.unimod <- unlist(lapply(unique(samp.ind), function(l){ # for every sample
      wu <- compute.wu(u = u[samp.ind == l], power = t.power)
      # Impose logconcavity on the extracted elution profile of the current sample
      u.constr <- pickLogConPeak(u[samp.ind == l], mode = sub.peak,
                                 tau = cobs.tau, weights = wu, basfun = basfun.u)
      return(u.constr)
    }))
    if(sum(u)== 0){
      if(verbose){ cat("This compound seems to be a false positive. Proceeding to next peak.\n\n") }
      peak.locs[1, "active"] <- 0 # Inactivate peak
      peak.locs <- peak.locs[c(2:nrow(peak.locs), 1),] # Put the first peak last
      next()
    }

    # Select the sample for which the compounds is most concentrated

    sub.peaks <- sapply(split(u, f = rep(1:length(samps), each = length(win)/length(samps))), which.max)
    sub.peaks <- sub.peaks + length(win)/length(samps) * ((1:length(samps))-1)
    u.heights <- u[sub.peaks]
    sub.samp <- samps[which.max(u.heights)]
    sub.win <- win[samp.ind == sub.samp]
    u <- u[samp.ind == sub.samp]
    sub.peak <- which.max(u)

    # Plot the sample selection
    if(verbose){ cat(paste0("Highest profile found in sample ", sub.samp, ".\n\n"))}
    if(Plot & debug){
      par(mfrow = c(2,1), oma = c(0,0,2,0))
      gam <- 0.1
      zlims <- c(0,max(M))^gam
      plotmat(M, subsamp = 1, plot.axes = F, zlim = zlims, gamma = gam,
              main = "Input data", xlab = "Time", ylab = "Mass")
      abline(v = length(win)/length(samps)*(1:(length(samps)-1)), col = "white")
      abline(v = which(peak == win), col = "red", lwd = 2)
      plotmat(Matrix::tcrossprod(u.unimod,v), subsamp = 1, plot.axes = F, zlim = zlims, gamma = gam,
              main = "Extracted data", xlab = "Time", ylab = "Mass")
      abline(v = length(u.unimod)/length(samps)*(1:(length(samps)-1)), col = "white")
      abline(v = sub.peaks, col = "red", lwd = 2)
      abline(v = sub.peaks[which.max(u.heights)], col = "green4", lwd = 2)
      title("Initial elution profile estimation accross samples", outer = T)
    }

    # 3. Optimize the spectrum and elution profile

    # Optimization
    if(verbose){ cat(paste0("Optimize the spectrum in sample ", sub.samp, ".\n"))}
    out <- optimize.profiles(M.all = resids,
                             scantimes = scantimes, mzvalues = mzvalues,
                             u = u, v = v,
                             win = sub.win, peak.mode = sub.peak, m = m,
                             tau = tau, cobs.tau = cobs.tau, qr.method = "quantreg", # "uqr.cpp", "quantreg"
                             basfun.u = basfun.u,
                             power.u = t.power, power.v = sel.power,
                             max.expand = 300,
                             maxiter = maxiter, tol = tol,
                             verbose = verbose, debug = FALSE, Plot = Plot & debug)
    if(!out$success){
      if(verbose){ cat("This compound seems to be a false positive. Proceeding to next peak.\n\n") }
      # peak.locs <- peak.locs[-1,,drop=F]
      peak.locs[1, "active"] <- 0 # Inactivate peak
      peak.locs <- peak.locs[c(2:nrow(peak.locs), 1),] # Put the first peak last
      next()
    }

    # Update variables
    u <- out$u
    v <- out$v
    wv <- out$wv
    wu <- out$wu
    win.opt <- out$win

    # 4. Extract elution profiles from all samples using the optimized target spectrum

    if(verbose){ cat(paste("Fit elution profiles with optimized spectrum in all samples.", "\n"))}

    # Extend window by a bit, because the concentration shifts could lead to
    # small peak eluting outside of the highly concentrated peak window
    win.start <- max(min(win.opt - (ceiling(win.opt/m)-1)*m) - win.ext, 1)
    win.stop <- min(max(win.opt - (ceiling(win.opt/m)-1)*m) + win.ext, m)
    win <- win.start:win.stop
    win <- rep(win, length(samps)) + rep(m*(samps-1), each = length(win)) # Extend the window to all samples
    M <- resids[win,]
    u.all.fit <- lapply(unique(samp.ind), function(samp){
      # Adapt window to current sample
      sub.win <- win[ceiling(win/m) == samp]
      M <- resids[sub.win,]

      # Fit elution profile in current sample
      u.fit <- fit.profile(M = t(M), u = v, w = wv, tau = tau, method = "quantreg", rm.bas = TRUE)
      wu <- compute.wu(u = u.fit, power = t.power)
      u.fit <- pickLogConPeak(u.fit, mode = which.max(u.fit), tau = cobs.tau, weights = wu, basfun = basfun.u)

      # Scale the elution profile
      bas <- Matrix(0, nrow = m*s, ncol = n, sparse = TRUE) # Generate a 0 baseline even in case no factor was extracted; dim(bas) == dim(resids)
      if(sum(u.fit) == 0 ) return(list(u = u.fit, bas = bas[sub.win,]))
      # TODO find a way that the condition below is not required
      if(length(sub.win) < 400) bas[sub.win, ] <- fit.lin.bas(M = resids[sub.win, ], dim = 2, ext = round(length(sub.win)/10))
      u.scaled <- rep(0, length(u.fit))
      u.scaled[u.fit != 0] <- scale.factors(M = resids, u = u.fit[u.fit != 0], v = v,
                                            U.p = NULL, V.p = NULL, # scale without considering other coeluting compounds
                                            win = sub.win[u.fit != 0], bas = bas,
                                            subsamp = 1)
      bas[sub.win,v!=0] <- pmin(bas[sub.win, v !=0], M[, v != 0])  # this will avoid that data is added
      return(list(u = u.scaled, bas = bas[sub.win,]))
    })
    u.all <- unlist(lapply(u.all.fit, function(x) x$u))
    bas <- do.call(rbind, lapply(u.all.fit, function(x) x$bas))
    if(sum(u.all) == 0){
      if(verbose){ cat("\nThis compound seems to be a false positive. Proceeding to next peak.\n\n") }
      peak.locs[1, "active"] <- 0 # Inactivate peak
      peak.locs <- peak.locs[c(2:nrow(peak.locs), 1),] # Put the first peak last
      next()
    }
    if(verbose){ cat(paste0("Compound extracted in samples: ", paste0(unique(ceiling(win[u.all != 0]/m)), collapse = ", "),"\n"))}

    # 5. Update elution profile, spectrum profile, residuals, TICS, peak detection

    # Store elution and spectrum profiles
    if(verbose){ cat("Updating residuals, ")}
    S[,i] <- v
    u.all[u.all < eps] <- 0 # remove numerical instabilities
    C[win,i] <- u.all

    # Update residuals
    nzero.win <- win[u.all != 0]
    nzero.u.ind <- u.all != 0
    nzero.v.ind <- v != 0
    M.fit <- as.matrix(Matrix::tcrossprod(C[nzero.win,], S[nzero.v.ind,]))
    M.obs <- resids[nzero.win, nzero.v.ind]
    M.bas <- as.matrix(bas[nzero.u.ind, nzero.v.ind])
    M.res <- pmax(M.obs - M.fit , 0) # Classical residuals
    M.res <- pmax(M.res^2/M.obs, M.bas) # adapted Pearson residuals
    M.res[is.na(M.res) | is.infinite(M.res)] <- M.bas[is.na(M.res) | is.infinite(M.res)]
    resids[nzero.win, nzero.v.ind] <- M.res

    # Update TIC
    if(verbose){ cat("TIC, ")}
    TICs.old <- TICs
    TICs[nzero.win] <- rowSums(M.res)

    # Update the peak heights
    if(verbose){ cat("peak heights, ")}
    peak.locs.old <- peak.locs # For plotting
    peak.upds <- do.call(rbind, lapply(samps, function(samp){
      peak.upd <- update.heights(peak.locs = peak.locs[peak.locs[,"samp"] == samp,, drop = F], # Update peak heights only were the fitted u is not 0 for the current sample
                                 win = nzero.win[ceiling(nzero.win/m)==samp], # adjust window to current sample
                                 Dic = C.dic, Y = TICs, m = m)
      return(peak.upd)
    }))
    if(length(peak.upds) == 0){
      if(verbose){ cat("\nThis compound seems to be a false positive. Proceeding to next peak.\n\n") }
      peak.locs[1, "active"] <- 0 # Inactivate peak
      peak.locs <- peak.locs[c(2:nrow(peak.locs), 1),] # Put the first peak last
      next()
    }
    inds.upd <- match(peak.upds[,"max"], peak.locs[,"max"], nomatch = 0)
    peak.locs[inds.upd, "height"] <- peak.upds[, "height.new"]

    # Update the set of active peaks
    if(verbose){ cat("peak filter\n")}
    peak.upds <- update.detection(Y = TICs, Dic = C.dic, peaks = peak.locs, sel = inds.upd, cutoff = cutoff, m = m)
    inds.upd <- match(peak.upds[,"max"], peak.locs[,"max"], nomatch = 0)
    peak.locs[inds.upd, c("crit", "active")] <- peak.upds[, c("crit", "active")]
    peak.locs[1, "active"] <- 0 # Inactivate the first peak to avoid infinite loops

    # Sort peaks
    peak.locs <- peak.locs[order(peak.locs[,"active"], peak.locs[,"height"], decreasing = T),]

    # Plot the peak update
    if(Plot & debug){
      par(mfrow = c(2,1))
      plot.peak.identification(y = TICs.old, Dic = C.dic, peaks = peak.locs.old, m = m, labels = "crit", stick.cols = "active",
                               win = nzero.win, scantimes = scantimes, main = "Peak identification before update")
      plot.peak.identification(y = TICs, Dic = C.dic, peaks = peak.locs, m = m,  labels = "crit", stick.cols = "active",
                               win = nzero.win, scantimes = scantimes, main = "Peak identification after update")
    }
    if(verbose){
      cat(paste0("Number of active peaks: ",  sum(peak.locs[,"active"]), "\n"))
      peak.thresh <- peak.locs[peak.locs[,"active"] == 1, "height"] > max.start * stop.thres
      cat(paste0("Number of active peaks above threshold: ",  sum(peak.thresh), "\n\n"))
    }

    # 6. Plot end of iteration results

    if(Plot){
      win.plot <- unique(win - (ceiling(win/m)-1) * m)
      start <- max(min(win.plot) - 100, 1)
      stop <- min(max(win.plot) + 100, m)
      win.plot <- start:stop
      win.plot <- rep(win.plot, length(samps)) + rep(m * (samps-1), each = length(win.plot))

      par(mfrow = c(3,1))
      gam <- 0.1
      zlims <- c(0,max(data[win.plot,]))^gam
      # Input
      # TODO data is not imported as a function argument
      plotmat(data[win.plot,], subsamp = 1, plot.axes = F, zlim = zlims, main = "Input data", gamma = gam)
      abline(v = length(win.plot)/length(samps)*(1:(length(samps)-1)), col = "white")
      abline(v = which(peak == win.plot), col = "red", lwd = 2)
      # Fitted
      plotmat(Matrix::tcrossprod(C[win.plot,i],v), subsamp = 1, plot.axes = F, zlim = zlims, main = "Extracted data", gamma = gam)
      abline(v = length(win.plot)/length(samps)*(1:(length(samps)-1)), col = "white")
      # Residuals
      plotmat(resids[win.plot,], subsamp = 1, plot.axes = F, zlim = zlims, main = "Residual data after extraction", gamma = gam)
      abline(v = length(win.plot)/length(samps)*(1:(length(samps)-1)), col = "white")
      abline(v = which(peak == win.plot), col = "red", lwd = 2)

      Sys.sleep(0.5) # Otherwise the plots are not shown correctly
    }

    # 7. Save back up intermediate results in case of crash
    # TODO find a generic list (or make a class?) to store intermediate and final results
    # TODO recevory should be performed every X iterations
    # TODO uncomment lines below and implement the saving of recovery files in case of crash
    # if is bigmatrix -> flush(resids) # write to file !!! too slow, so do every X iterations
    # function to save backup


    t1 <- proc.time()[3]-t1
    i <- i + 1
    if(verbose){ cat(paste0("End of iteration ", i-1, ".\tTotal time: ", round(t1,1), " s.\n\n"))}
    if(debug){ if(readline(prompt = "Press [enter]") != "") break() }
  }

  # Save back up of end results
  # TODO implement backup

  # Close the pdf
  if(Plot.to.file) dev.off()

  return(list(C = C[,1:(i-1)],
              S = S[,1:(i-1)],
              iter = i-1,
              peak.locs = peak.locs,
              TICs = TICs,
              timing = proc.time()[3]-t.tot))
}



# plot.deconvolution.fit ####
# TODO write documentation
plot.deconvolution.fit <- function(M, M.res = NULL, C, S, m = NULL, win = NULL, gam = 0.1, sub.fac = 1, file.name = NULL){
  # TODO display scantimes on the x axis
  if(is.null(m)) m <- nrow(M)
  if(is.null(win)) win <- 1:m
  if(is.null(M.res)){
    M.res <- pmax(M[win,] - Matrix::tcrossprod(C[win,],S), 0) # TODO allow for adapted pearson residuals
  } else {
    M.res <- M.res[win,]
  }
  samps <- unique(ceiling(win/m))
  samp.sep <- ceiling(length(win)/length(samps)*(1:(length(samps)-1))/sub.fac)+0.5
  # Plot the heatmaps
  if(!is.null(file.name)) png(file = file.name, height = 3000, width = 3000, res = 300)
  par(mfrow = c(3,1))
  zlims <- c(0,max(M[win,]))^gam
  # Input
  plotmat(M[win,], subsamp = sub.fac, plot.axes = F, zlim = zlims, main = "Input data", gamma = gam)
  abline(v = samp.sep, col = "white")
  # Fitted
  plotmat(Matrix::tcrossprod(C[win,],S), subsamp = sub.fac, plot.axes = F, zlim = zlims, main = "Extracted data", gamma = gam)
  abline(v = samp.sep, col = "white")
  # Residuals
  plotmat(M.res, subsamp = sub.fac, plot.axes = F, zlim = zlims, main = "Residual data after extraction", gamma = gam)
  abline(v = samp.sep, col = "white")
  if(!is.null(file.name)) dev.off()
}


# plot.elution.overlay ####
# TODO write documentation
plot.elution.overlay <- function(C, scantimes, facs, samps = NULL, const = 1, ...){
  if(nrow(C) %% length(scantimes) != 0) stop("nrow(C) must be a multiple of length(scantimes)")
  if(is.null(samps)) samps <- 1:(nrow(C) / length(scantimes))
  for(k in facs){
    win <- Matrix::which(C[,k,drop=F] != 0)
    win <- range(win - (ceiling(win/m)-1) * m)
    win <- max(1,win[1]-10):min(m, win[2]+10)
    elus <- sapply(samps, function(samp) return(C[win + rep((samp-1)*m, length(win)),k]))
    matplot(x = scantimes[win], y = elus+const,
            type = "l", xlab = "Retention time (min)", ylab = "Intensity",
            main = paste0("Compound ", k), ...)
  }
}

# update.spectrum ####
# TODO make sure "data" is the baseline removed data !
# TODO write documentation
#   data:   bigmemory object
#   method:  "WNNLS" or "nnL0pois"
#   focus:  should the spectrum profiles of a compound be focused on the sample
#           where it is highest or on all samples ?
#   verbose: should the progress be printed
update.spectrum <- function(data, C, S, m = NULL, method = "nnL0pois",
                            focus = TRUE, verbose = FALSE){
  if(is.null(m)) m <- nrow(data)
  n <- ncol(data)

  # TODO parallelize this
  for(k in 1:ncol(C)){
    x <- C[,k,drop=F]
    # Get the elution window
    win.elu <- Matrix::which(x != 0)
    if(focus){
      # Get the sample to fit the spectrum
      # Get position of highest elution profile
      samp.fit <- Matrix::which(x == max(x)) # computationally more efficient than which.max()
      samp.fit <- ceiling(samp.fit/m)
      win.elu <- win.elu[ceiling(win.elu/m) == samp.fit]
    }
    # TODO is extending the elution window necessary ? It is possible that a profile
    # that is cut by the window fits perfectly the data but not at all outside (typical
    # example would be acid peak), but only the focal compound is updated to so maybe
    # there will be only limited impact...
    win.coel <- win.elu
    if(FALSE){ # currently disable extending the window
      # Extend elution window to coeluting compounds
      x.coel <- C[, which(Matrix::colSums(C[win.elu,]) != 0), drop=F]
      win.coel <- Matrix::which(Matrix::rowSums(x.coel) != 0)
      win.coel <- win.coel[ceiling(win.coel/m) == samp.fit]
    }
    # Build the covariate matrix with all overlaping compounds within the
    # coelution window
    X <- C[win.coel,]
    X <- X[,c(k, (1:ncol(X))[-k])] # Put the current elution profile first
    k.nz <- Matrix::which(Matrix::colSums(X) != 0) # non-zero compounds in the window
    X <- as.matrix(X[,k.nz]) # Subsetting removes the sparsity in the matrix
    # L2 normalize the profiles
    x.n <- apply(X, 2, norm, type = "2")
    X <- sweep(X, 2, x.n, "/")
    # Fit the spectrum profile fit fitting the covariate matrix to all m/z traces.
    # Fit data considering the Poisson structure of the data
    Y <- data[win.coel,]
    S[,k] <- sapply(1:n, function(j){
      # for(j in 1:n){
      if(sum(Y[,j]) == 0) return(0)
      if(method == "nnL0pois"){
        # Fit model with Poisson noise and non-negativity constraint and L0 penalty
        # Get the initial values
        s0 <- S[j, c(k, (1:ncol(S))[-k])]  # reoder the compounds as for X
        s0 <- s0[k.nz] * x.n # use the same compounds as for X and scale with respect to L2 norm of X
        # Fit model
        s1 <- L0glm.fit(X = X, y = Y[,j], family = poisson(identity), lambda = 1,
                        start = s0, nonnegative = TRUE, normalize = FALSE,
                        control.l0 = control.l0.gen(maxit = 25),
                        control.iwls = control.iwls.gen(maxit = 1, thresh = 1E2),
                        control.fit = control.fit.gen())$coefficients
      } else if(method == "WNNLS"){
        # Compute weights as 1/variance to account for poisson noise
        w <- 1/sqrt(Y[,j]+1)
        w <- w/sum(w)*length(w)
        # Fit model
        s1 <- nnls(A = sweep(X, 1, sqrt(w), "*"), b = Y[,j] * sqrt(w))$x # WNNLS
      } else {
        stop(paste0("Deconvolution method '", method, "' not implemented."))
      }
      #   S[j,k] <- s1[1]
      # }
      return(s1[1]) # return the first value which corresponds to the current compound
    })
    # Normalize S_k
    if(sum(S[,k]) != 0) S[,k] <- S[,k]/sqrt(sum(S[,k,drop=F]^2))
    if(verbose) print.progress(k, ncol(C))
  }
  # Remove empty profiles
  sel <- Matrix::colSums(S) != 0
  if(sum(sel) == 0) stop("All spectrum profiles were removed.")
  C <- C[,sel,drop=F]
  S <- Matrix(S[,sel,drop=F], sparse = TRUE)
  return(list(S = S, C = C))
}

# reorder.columns ####
# TODO doc
scanline.order <- function(C, m = NULL) {
  if(is.null(m)) m <- nrow(C)
  # Reorder matrices from early to late scanline
  ord <- vector(mode = "numeric", length = ncol(C))
  for(k in 1:ncol(C)){
    ord[k] <- Matrix::which(C[,k,drop=F] == max(C[,k,drop=F]))
  }
  ord <- ord - (ceiling(ord/m)-1)*m
  ord <- order(ord)
  return(ord)
}


####---- TARGETED DECONVOLUTION ----####


# targeted.deconvolution ####
#
# Function to initialize elution profiles given elution search window and
# expected spectrum profiles of targets. Elution profiles are initialized
# using L0 penalized nonnegative Poisson regression.
#
# INPUT
#   data.desc:  a path (as character string) giving the location of the descriptor
#               file of the bigmemory matrix containing the data to deconvolute
#   S:  a matrix with the expected spectrum profiles. The dimension should be
#       n by k, where n is the number of column of the data and k the number of
#       target compounds.
#   m:  the number of time scans in a single chromatogram
#   C.prior:  a matrix containing starting values for C. The dimensions should
#             be ms by k, where ms is the number of rows of the data and k the
#             number of target compounds.
#   sample.names:   a vector of names of the samples contained in the bigmemory
#                   matrix
#   samps:    an indexing of the samples to analyze. This allow to deconvolute
#             only a part of the samples.
#   ncores:   the number of cores to use for parallelization.
#   smooth:   should the profiles be smoothed
#   dampen:   should the profiles be dampen with a kernel function depending on
#             RT. The function is:  (1 - abs(dist_to_exp_RT)^3)^5
#   unimodal: should the profles be contrained to unimodality.
#   w.power:
#   plot.dir:   path indicating where to plot the initialization results. Note
#               it will create a new sub folder. If NULL, nothing is plotted.
# OUTPUT
#   C:  a sparse matrix containing the row-augmented elution profiles for all
#       samples
targeted.deconvolution <- function(data.filename, data.dir,
                                   S, C.prior,
                                   m, samps, samples.df,
                                   scantimes, mzvals,
                                   plot.dir = NULL, verbose = TRUE,
                                   ncores = 1,
                                   t.power = 2, cobs.tau = 0.2,
                                   smooth = FALSE, dampen = FALSE, unimodal = FALSE){
  t1 <- proc.time()[3]
  # Get the data description file
  data.desc <- paste0(data.dir, data.filename, ".desc")

  # Extract the elution profiles given the target spectrum profiles
  if(verbose) cat("Extracting the target elution profiles.\n")
  # Prepare cluster
  cl <- makeCluster(ncores, type = "SOCK") # TODO autodetect: on linux use type="FORK" (faster), but not available on windows
  registerDoSNOW(cl)
  snow::clusterExport(cl, list = c("S", "m", "samples.df", "C.prior",
                                   "smooth", "dampen", "unimodal", "data.desc",
                                   "cobs.tau", "pickLogConPeak", "plot.dir",
                                   "plotmat"), envir = environment()) # export values that each worker needs
  # Parallelize over samples
  C <- clusterApply(cl, samps, function(samp) {
    require(Matrix, quietly = TRUE)
    require(bigmemory, quietly = TRUE)
    require(L0glm, quietly = TRUE)
    # Load disk backed file and store current samples in memory (more efficient, but memory intensive !)
    M <- attach.big.matrix(dget(data.desc))
    win.samp <- (1:m) + (samp-1)*m
    M.samp <- M[win.samp,]
    C.samp <- Matrix(0, nrow = m, ncol = ncol(S)) # C.samp is not a large object (+/- 30MB) so more efficient to fill it as a dense matrix

    for(i in 1:nrow(M.samp)){
      # Select all compounds for which the search window overlaps the current time scan
      prior <- C.prior[win.samp[i],]
      comp.sel <- prior != 0
      # comp.sel <- which(RT.f$ind.up >= i & RT.df$ind.low <= i)
      if(length(comp.sel) == 0) next()
      # Fit an L0 penalized nonnegative linear model with identity link Poisson error structure
      coefs <- L0glm.fit(X = S[,comp.sel, drop = FALSE], y = M.samp[i,],
                         family = poisson(link = "identity"),
                         normalize = FALSE, # Already normalized when parsed
                         nonnegative = TRUE, lambda = 1, # it seems that lambda = 1 is the best lambda for any L0 penalized Poisson regression
                         start = prior[comp.sel],
                         control.l0 = control.l0.gen(maxit = 1),
                         control.iwls = control.iwls.gen(maxit = 1, thresh = 1E2), # Note, thresh = 1E2 means that we trim all coefficients < 1E2 to 0
                         control.fit = control.fit.gen(maxit = 1, block.size = NULL) # No block deconvolution, data should be small enough
      )$coefficients
      # Update the elution profiles
      C.samp[i,comp.sel] <- coefs
    }

    # Process the elution profiles
    if(any(smooth, dampen, unimodal)){
      for(k in 1:ncol(C.samp)){
        nz <- Matrix::which(C.samp[,k] != 0)
        if(length(nz) == 0) next()
        win <- min(nz):max(nz)
        y0 <- y <- C.samp[win,k]

        # Smooth profiles
        if(smooth){
          x <- 1:length(y)
          smooth.win <- 15 # TODO add argument above ?
          y <- pmax(loess(y ~ x, degree = 1, span = smooth.win/length(y))$fitted, 0)
          # y <- pmax(sgolayfilt(runmed(y, k = 15), p = 1, n = 15),0)
        }

        # Dampen profiles with kernel function
        if(dampen){
          win.c <- win[round(length(win)/2)] # center of the window
          reldev <- abs(win - win.c)/(length(win)/2)
          kern <- (1-reldev^3)^5 # kernel function = (1 - reldev^a)^b where a = 3, b = 5
          y <- y * kern
        }

        # Impose unimodality constraints
        if(unimodal){
          # TODO think about how to best impose unimodality, eventually selecting best peak based on the weighted cosine similarity, also remove empty compounds at the end of the update
          w.y <- y^2/sum(y^2)*length(y)
          peak.pos <- which.max(y)
          y <- pickLogConPeak(y = y, tau = cobs.tau, weights = w.y,
                              local.fit = TRUE, mode = peak.pos, Plot = FALSE)
        }
        C.samp[win,k] <- y
      }
    }

    # Plot elution profiles
    if(!is.null(plot.dir)){
      png(filename = paste0(plot.dir, "Sample - ", gsub(pattern = "[.]|/", replacement = "", x = samples.df$name[samp]), ".png"),
          res = 300, height = 4000, width = 6000)
      par(mfrow = c(2,1))
      gam <- 0.1
      zlims <- c(0, max(M.samp))^gam
      subf <- 10
      plotmat(M.samp[seq(1, m, by = subf),], gamma = gam, zlim = zlims,
              subsamp = 1, main = paste0("Input data - Sample ", samples.df$name[samp]))
      plotmat(Matrix::tcrossprod(C.samp[seq(1, m, by = subf),], S), gamma = gam,
              zlim = zlims, subsamp = 1, main = "Extracted elution profiles")
      dev.off()
    }
    return(Matrix(C.samp, sparse = TRUE))
  })
  stopCluster(cl)
  C <- do.call(rbind, C)

  # Return results
  notes <- paste0("Targeted deconvolution.\n\n",
                  "The results contain the extracted elution profiles for ", length(samps),
                  " samples and ", ncol(S), " compounds.\n",
                  "The elution profiles are stored in 'C' and the spectrum profiles are ",
                  "stored in 'S'.\n",
                  "Computed on: ", Sys.Date(), "\n",
                  "Timing: ", round((proc.time()[3] - t1)/60, 2), " minutes.\n")
  if(verbose) cat(paste0("\nOUTPUT\n\n", notes))
  out <- list(C = C, S = S, notes = notes)
  return(out)
}


# get.C.prior ####
# TODO doc
get.C.prior <- function(RT.df, m, s){
  elu <- Matrix(0, nrow = m, ncol = 1, sparse = TRUE)
  # create the prior for a sample
  C.samp <- do.call(cbind, lapply(1:nrow(RT.df), function(k){
    elu[RT.df$ind.low[k]:RT.df$ind.up[k], drop = F] <- 1E3 # 1E3 otherwise the all will be remove during the L0 regression (because thresh = 1E2)
    return(elu)
  }))
  # Duplicate the prior to all samples
  C <- do.call(rbind, lapply(1:s, function(i) return(C.samp)))
  return(C)
}


# fill.C ####
# TODO doc
fill.C <- function(C, m, n, s,
                   untargeted.samples, targeted.samples,
                   win.ext = 25){
  # Initialize the elution profiles for the targeted samples, this is common to all the targeted samples
  C.targ <- Matrix(0, nrow = m, ncol = ncol(C),  sparse = TRUE)
  for(k in 1:ncol(C)){
    elus <- C[rep(1:m, length(untargeted.samples)) + (rep(untargeted.samples, each = m)-1)*m, k, drop=F]
    nz <- Matrix::which(elus > 0)
    if(length(nz) == 0) next()
    win <- range(nz - (ceiling(nz/m)-1)*m)
    win <- max(win[1]-win.ext, 1):min(win[2]+win.ext, m)
    # TODO possible options: average profile, constant prior where nonzero
    method <- "window"
    if(method == "average"){
      elus <- sapply(untargeted.samples, function(samp) return(C[win + rep((samp-1)*m, length(win)),k]))
      elus <- rowMeans(elus)
      C.targ[win, k] <- elus
    } else if(method == "window"){
      C.targ[win, k] <- 1E3
    }
  }
  C <- do.call(rbind, lapply(1:s, function(samp){
    if(samp %in% untargeted.samples){
      return(C[(1:m)+(samp-1)*m,, drop = FALSE])
    } else {
      return(C.targ)
    }
  }))
  return(C)
}

