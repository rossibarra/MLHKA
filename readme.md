##MLHKA 

A maximum likelihood ratio test of natural selection, using polymorphism and divergence data. 

Please cite 

	Wright, S.I. and Charlesworth, B. The HKA test revisited: A maximum likelihood ratio test of the standard neutral model. Genetics. 168: 1071-1076.

http://labs.eeb.utoronto.ca/wright/Stephen_I._Wright/Programs.html

The program can be run two ways. It can be run as described in the manual (see below) with no command-line parameters, or it can be run with command line parameters:

	MLHKA -i <inpute file name> -o <output file name> -s <random seed> -c <length of MCMC chain to run>

The program is available as the source code (MLHKA_version2.1.cpp), help file (README.pdf) and sample input file (infile.txt). Original versions are available at Stephen Wright's website:

	http://labs.eeb.utoronto.ca/wright/Stephen_I._Wright/Programs.html

The source code requires a c++ compliler, e.g. on unix, type: 'g++ MLHKA_version2.cpp –o MLHKA'.

For questions, please e-mail stephen.wright@utoronto.ca .  For problems with this version of the code, please email rossibarra at ucdavis dot edu.
