# MutExMatSorting

![alt text](https://github.com/AleVin1995/MutExMatSorting/blob/master/web/MExMas_logo.jpg)

**Mutual Exclusivity Sorting of Binary Matrices**

This package implements an heuristic algorithm that takes in input a sparse binary matrix and sorts its rows and columns in a way that the patterns of non-null entries have a minimal overlap across rows. This highlights possible mutual exclusive trends among these patterns.

Install
--

Install the R package [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) and then install `MutExMatSorting`

```
install.packages("devtools")
library(devtools)

install_github("AleVin1995/MutExMatSorting")
```
