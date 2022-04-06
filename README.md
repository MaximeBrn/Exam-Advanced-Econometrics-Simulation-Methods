# Exam Advanced Econometrics Simulation Methods


### Replication files
- See https://github.com/waynejtaylor/Single-Agent-Dynamic-Choice/blob/master/AM2011Table1cols2356.R
- See https://www.econometricsociety.org/publications/econometrica/2011/11/01/conditional-choice-probability-estimation-dynamic-discrete
- ;

### First Look at the Code
<ins>**i) R**</ins>
- R code runs pretty well (might encounter some issues if you run it with Jupyter notebook due to a version of R and package Rcpp but runs well on **RStudio**)
- Only first example (_#Bus_), please run __Example1_BUS_Monte_Carlo.R__ in _Exam-Advanced-Econometrics-Simulation-Methods\Code Original Arcidiacono Miller 2011\CODE R\Main_ after having changed the directory.
- Some functions are in C++ to gain speed and converted to R language via "Rcpp" package.
- In addition you can have a look at **Example1_BUS_DGP**: the code for the data generating process.

<ins>**ii) Matlab**</ins>
- Both examples but split into numerous small files.
- For example 2, you can **sumresults.m** after having changed the directory to have a look at the papers' results for the second example.
