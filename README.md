Multivariate Bayesian spectral inference
================
2023-05-07

# 1 gibbs_vnpc()

This is the source code of the multivariate corrected likelihood (VNPC)
approach proposed for multivariate Bayesian spectral inference. The
execution of this approach is on R, while the core functions used in the
algorithm are written on C++ to enhance the performance.

The reference of VNPC is coming soon to arXiv…

# 2 A Metropolis-within-Gibbs algorithm

# 3 The Southern Oscillation example

Here are the estimates for the famous bivariate Southern Oscillation
study. There are two series in this example, the Southern Oscillation
Index (SOI) and the associated fish recruitment series, which are
simultaneously recorded to explore the El Nino cycle (see Shumway &
Stoffer (2011)).

Both series consist of a period of 453 monthly values over the years
1950-1987, which are shown in Figure <a href="#fig:Fig1">3.1</a>. The
data is available in the R package `astsa` as datasets `soi` and `rec`
(see Stoffer (2022)).

<div class="figure" style="text-align: center">

<img src="SOI_REC.jpg" alt="The SOI and Recruitment series" width="80%" style="text-align: center" />

<p class="caption">

Figure 3.1: The SOI and Recruitment series

</p>

</div>

# 4 References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-shumway2011" class="csl-entry">

Shumway, R. H., & Stoffer, D. S. (2011). *Time series analysis and its
applications with r examples* (3rd ed.). Springer.

</div>

<div id="ref-stoffer2022" class="csl-entry">

Stoffer, D. (2022). *Astsa: Applied statistical time series analysis*. R
package version 1.16. <https://cran.r-project.org/web/packages/astsa/>

</div>

</div>
