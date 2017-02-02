artiFISHal
==========

**artiFISHal** is a package of functions for the [R programming language](http://www.r-project.org/).  **artiFISHal** is a pelagic fish community simulator, which can be used to create artificial lakes and populate them with known numbers of fish (identified in species-size groups) to mimic pelagic fish communities. **artiFISHal** can then be used to sample these fish with virtual acoustic and midwater trawl surveys.

[Yule et al. (2013)](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ) used **artiFISHal** to evaluate several different approaches for estimating the biomass of pelagic fish species in the Great Lakes. An example of how to use the functions in **artiFISHal** is given in this [vignette](https://github.com/JVAdams/artiFISHal/blob/master/Vignette.md).

- - -

You can access the functions by installing the package from within R.

    source("https://raw.githubusercontent.com/MangoTheCat/remotes/master/install-github.R")$value("mangothecat/remotes")
    remotes::install_github("JVAdams/artiFISHal")
    library(artiFISHal)
	
- - -

*NOTE:  If you are looking for files that were used in the Great Lakes Acoustic Users Group 2014 Workshop on Trawl Performance, they are now located at the [GLAUG2014](https://github.com/JVAdams/GLAUG2014) repository instead.

- - -

_Thanks to Hilary Parker whose blog post [Writing an R package from scratch](http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/) encouraged me to create my first R package._

- - -

_U.S. Geological Survey_ (USGS) Computer Program **artiFISHal** version 0.0.0.9003. 
Written by Jean V. Adams, [USGS - Great Lakes Science Center](http://www.glsc.usgs.gov/), Ann Arbor, Michigan, USA. 
Written in programming language R (R Core Team, 2015, www.R-project.org), version 3.2.3 (2015-12-10). 
Run on a PC with Intel(R) Core(TM) I7-4600m CPU, 2.90 GHz processor, 16.0 GB RAM, and Microsoft Windows 7 Enterprise operating system 2009 Service Pack 1. 
Source code is available from Jean V. Adams on [GitHub](https://github.com/JVAdams/artiFISHal), _jvadams (at) usgs (dot) gov_.

_Disclaimer:_ This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.
