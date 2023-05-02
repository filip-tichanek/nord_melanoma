# nord_melanoma

This Repository shows R code and source data used in the publication '[SURVIVAL IN MELANOMA IN THE NORDIC COUNTRIES INTO THE ERA OF TARGETED AND IMMUNOLOGICAL THERAPIES](https://www.ejcancer.com/article/S0959-8049(23)00157-0/fulltext)', 
published in the [European Journal of Cancer](https://www.ejcancer.com/)

All data used were downloaded from the [NORDCAN website](https://nordcan.iarc.fr/en). The NORDCAN data include the estimates of 1-year and 5-years survival and their confidence intervals, based on the Pohar Perme estimator. This code utilizes the pre-calculated NORDCAN data (including the pre-calculated confidence intervals) to get continuous survival estimates and estimation of the magnitude of the change (slope) across a continuous timescale. The code also estimates 95% credible intervals for both the survival rate and the magnitude of change. This is done via Generalized Additive Models in the Bayesian framework  with [brms package](https://cran.r-project.org/web/packages/brms/index.html)

Please, cite the [original article](https://www.ejcancer.com/article/S0959-8049(23)00157-0/fulltext) when using this R code
