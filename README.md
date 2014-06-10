artiFISHal
==========

**artiFISHal** is a package of functions for the [R programming language](http://www.r-project.org/).  
**artiFISHal** is a pelagic fish community simulator, which can be used to create artificial lakes and 
populate them with known numbers of fish (identified in species-size groups) 
to mimic pelagic fish communities. 
**artiFISHal** can then be used to sample these fish with virtual acoustic and midwater trawl surveys.

[Yule et al. (2013)](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ) used **artiFISHal** 
to evaluate several different approaches for estimating the biomass of pelagic fish species in the Great Lakes.

You should be able to access the functions by installing them directly from within R.

	library("devtools")
	devtools::install_github("JVAdams/artiFISHal")
	library(artiFISHal)

If you don't already have `Rtools` and `devtools`, you will need to download and install Rtools 3.1 from [CRAN](http://cran.r-project.org/bin/windows/Rtools/), 
	the Comprehensive R Archive Network, then run the following lines of code before submitting the code above:

	find_rtools()
	install.packages("devtools")

_Thanks to Hilary Parker whose blog post [Writing an R package from scratch](http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/)
encouraged me to create my first R package._

- - -
*NOTE:  If you are looking for the GITHUB repository that used to be at this location,
look at the [GLAUG2014](https://github.com/JVAdams/GLAUG2014) repository instead.
I moved all of the files there that were used in the Great Lakes Acoustic Users Group 2014 Workshop on Trawl Performance.*
- - -
