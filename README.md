# MCR

Code to reproduce results from Section 9.1 & 10 of [Fisher et al., 2019](https://arxiv.org/abs/1801.01489), for calculating model reliance, model class reliance, and other functions for Rashomon sets.

Before running this code, users should install the R packages below.

```{r}
library(devtools)

#For solving (possibly non-convex) quadratic programs with 1 quadratic constraint:
install_github('aaronjfisher/qp1qc', ref='bd41afc2d7a29a7b602ac1d40c4522335e8455c7')

#For computing Model Class Reliance:
install_github('aaronjfisher/mcr', ref='ddff960ac0c17ba96a032003a06ed6f56d676ee9')

# These references are to the specific commits associated with our final manuscript.
```

After installing the above packages, users can run code from the directories contained here.

* *simulations* - contains code to generate the polynomial classifier example in our paper. Run time ~1-3 minutes.
* *propublica-analysis* - contains code to generate our applied results. The data file used in this analysis was downloaded from the "compas-scores-two-years.csv" file, at [https://github.com/propublica/compas-analysis](https://github.com/propublica/compas-analysis/tree/bafff5da3f2e45eca6c2d5055faad269defd135a).  Run time ~1-3 hours.



