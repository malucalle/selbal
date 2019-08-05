# selbal


 `selbal` is an R package for selection of balances in microbiome compositional data. As described in Rivera-Pinto et al. 2018 _Balances:  a new perspective for microbiome analysis_ https://msystems.asm.org/content/3/4/e00053-18, `selbal` implements a forward-selection method for the identification of two groups of taxa whose relative abundance, or balance, is associated with the response variable of interest.

## Getting Started


### Installation

To get a full access to the functions implemented in `selbal` we only need to run 
the following instructions:

```
# Installing the files in the repository
  # Option 1 (non - Windows users)
    devtools::install_github(repo = "UVic-omics/selbal")
  # Option 2 (for Windows users)
    devtools::install_url(url="https://github.com/UVic-omics/selbal/archive/master.zip", 
                INSTALL_opt= "--no-multiarch")
    
# Loading the library
  library("selbal")
```

### Running `selbal`

To start using `selbal` we recomend to:

- Read the manuscript [https://msystems.asm.org/content/3/4/e00053-18]
- Use the `help()` functions for getting a detailed instructions of their
  use.
- Read the associated vignette (see  https://htmlpreview.github.io/?https://github.com/UVic-omics/selbal/blob/master/vignettes/vignette.html).






