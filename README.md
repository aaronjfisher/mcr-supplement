# Code Supplement for Model Class Reliance

Code to reproduce results from Section 9.1 & 10 of [Fisher et al., 2019](https://arxiv.org/abs/1801.01489), for calculating model reliance and model class reliance.

Before describing installation instructions, we include a publication disclaimer.

### Disclaimer

Following the [JMLR release form](http://www.jmlr.org/forms/release.pdf), we make the following disclaimer:

> THIS SOURCE CODE IS SUPPLIED “AS IS” WITHOUT WARRANTY OF ANY KIND, AND ITS AUTHOR AND THE JOURNAL OF MACHINE LEARNING RESEARCH (JMLR) AND JMLR’S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND NEITHER THE AUTHOR NOR JMLR, NOR JMLR’S PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality of the foregoing, neither the author, nor JMLR, nor JMLR’s publishers and distributors, warrant that the Source Code will be error-free, will operate without interruption, or will meet the needs of the user.

We would ask users who encounter difficulties to please contact afishe27@alumni.jh.edu.

### Installation


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


Alternatively, package dependencies for this project are also managed via [packrat](https://github.com/rstudio/packrat). A bundled version of this project is available in the packrat/bundles folder (see [packrat::unbundle()](https://www.rdocumentation.org/packages/packrat/versions/0.5.0/topics/unbundle),  [packrat::on()](https://www.rdocumentation.org/packages/packrat/versions/0.5.0/topics/packrat-mode), and [related documentation](https://github.com/rstudio/packrat)).
