# MCR

Code to reproduce results from Section 9.1 & 10 of [Fisher et al., 2019](https://arxiv.org/abs/1801.01489), for calculating model reliance, model class reliance, and other functions for Rashomon sets.

Before running this code, users should install the R packages below.

```{r}
library(devtools)

#For solving (possibly non-convex) quadratic programs with 1 quadratic constraint:
install_github('aaronjfisher/qp1qc', ref='9a26035e5d1fae563ab2dc01a5bf6283a708d6a3')

#For computing Model Class Reliance:
install_github('aaronjfisher/mcr', ref='5b07e12eb7545eec50f9c23aed0ca8aad8cd4900')

# These references are to the specific commits associated with our final manuscript.
```

After installing the above packages, users can run code from the directories contained here.

* *simulations* - contains code to generate the polynomial classifier example in our paper. Run time ~1-3 minutes.
* *propublica-analysis* - contains code to generate our applied results. The data file used in this analysis was downloaded from the "compas-scores-two-years-violent.csv" file, at [https://github.com/propublica/compas-analysis](https://github.com/propublica/compas-analysis/tree/bafff5da3f2e45eca6c2d5055faad269defd135a).  Run time ~2-10 hours, depending on number of bootstrap samples.


Alternatively, package dependencies for this project are also managed via [packrat](https://github.com/rstudio/packrat). A bundled version of this project is available in the packrat/bundles folder (see [packrat::unbundle()](https://www.rdocumentation.org/packages/packrat/versions/0.5.0/topics/unbundle),  [packrat::on()](https://www.rdocumentation.org/packages/packrat/versions/0.5.0/topics/packrat-mode), and general documentation on [packrat](https://github.com/rstudio/packrat)).
