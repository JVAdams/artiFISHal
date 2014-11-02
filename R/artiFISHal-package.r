#' Pelagic Fish Community Simulator
#'
#' \pkg{artiFISHal} is a pelagic fish community simulator, which can be used to create artificial lakes and populate them with known numbers of fish 
#' (identified in species-size groups) to mimic pelagic fish communities. 
#' \pkg{artiFISHal} can then be used to sample these fish with virtual acoustic and midwater trawl surveys.
#' 
#' An example of how to use the functions in \pkg{artiFISHal} is given in this 
#' \href{https://github.com/JVAdams/artiFISHal/blob/master/Vignette.md}{vignette}.
#' Use \code{\link{SimFish}} to create an artificial lake and populate it with known numbers of fish (identified in species-size groups) 
#' to mimic pelagic fish communities.
#' Use \code{\link{SampFish}} to sample the fish population with virtual acoustic and midwater trawl surveys.
#' 
#' \href{http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ}{Yule et al. (2013)} used \pkg{artiFISHal} 
#' to evaluate several different approaches for estimating the biomass of pelagic fish species in the Great Lakes.
#' 
#' \emph{U.S. Geological Survey} (USGS) Computer Program \pkg{artiFISHal} version 2014-07.
#' Written by Jean V. Adams, \href{http://www.glsc.usgs.gov/}{USGS - Great Lakes Science Center}, Ann Arbor, Michigan, USA.
#' Written in programming language R (R Core Team, 2014, www.R-project.org), version 3.1.1 (2014-07-10).
#' Run on a PC with Intel(R) Core(TM) I7-4600m CPU, 2.90 GHz processor, 16.0 GB RAM, and 
#' Microsoft Windows 7 Enterprise operating system 2009 Service Pack 1.
#' Source code is available from Jean V. Adams on \href{https://github.com/JVAdams/artiFISHal}{GitHub}, _jvadams (at) usgs (dot) gov_.
#'
#' \emph{Disclaimer:}  Although this program has been used by the USGS, no warranty, expressed or implied, is
#' made by the USGS or the United States Government as to the accuracy and functioning of
#' the program and related program material nor shall the fact of distribution constitute
#' any such warranty, and no responsibility is assumed by the USGS in connection therewith.
#' 
#' @references Yule, DL, JV Adams, DM Warner, TR Hrabik, PM Kocovsky, BC Weidel, LG Rudstam, and PJ Sullivan.  2013.  
#' \href{http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ}{Evaluating 
#' analytical approaches for estimating pelagic fish biomass using simulated fish communities}. 
#' Canadian Journal of Fisheries and Aquatic Sciences 70:1845-1857.
#' @name artiFISHal
#' @docType package
NULL
