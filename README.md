# The Origins of Inequality

## Getting Started

### System requirements

This software has been tested on a MacBook Pro.

### Installation

To run this code, you will need to [install R](https://www.r-project.org/) and 
the following R packages:

```r
install.packages(
  c("ape", "crew", "ggtree", "patchwork", "phangorn", 
    "rstan", "tarchetypes", "targets", "tidyverse")
)
```

You will also need to download the software BayesTraits V5.0.3 from 
[here](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html)
and copy the executable file `BayesTraitsV5` into the `/bayestraits` folder
of this repository.

### Execute code

To run the pipeline:

1. Clone this repository to your local machine
2. Open the R Project file `origins-inequality.Rproj`
3. In the console, run the full pipeline with `targets::tar_make()`

## Help

Any issues, please email scott.claessens@gmail.com.

## Authors

Scott Claessens, scott.claessens@gmail.com
