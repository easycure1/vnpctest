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

This Section introduces the MCMC algorithm used by the VNPC approach to
sample the joint posterior distributions of the Hpd-Gamma process
parameters $\mathbf{\Phi}$, the number of Bernstein polynomials $k$ and
the coefficient matrices of the parametric VAR working model
$\underline{\boldsymbol{B}}$. In particular, $\mathbf{\Phi}$ is
hierachically parametrised by $r_{l}$’s, $x_{l}$’s and
$\boldsymbol{U}_{l}$’s; furthermore, we employ the approach proposed by
Mittelbach et al. (2012) (see its Section IV) to reparametrise
$\boldsymbol{U}_{l}$’s by some hyperspherical coordinates $\varphi_{l}$’
s. For more details, please see Section 5.2 in ??? (upcoming thesis). It
should note that there are no recognised distributions for the full
conditional posterior, so the classic Gibbs sampler is not suitable here
for updating the parameter values. As a result, we use a
Metropolis-Hastings step to achieve the updates, and a
Metropolis-within-Gibbs algorithm is proposed to be used.

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

<img src="SOI_REC.jpg" alt="The SOI and Recruitment series" width="80%" />

<p class="caption">

Figure 3.1: The SOI and Recruitment series

</p>

</div>

Both series are estimated simultaneously. Figures
<a href="#fig:Fig2">3.2</a> and <a href="#fig:Fig3">3.3</a> show the
estimates of the VNPC approach using different parametric working model
orders p. As a comparison, we also consider the parametric vector
autoregressive model (VAR) with the same orders. Figure
<a href="#fig:Fig4">3.4</a> gives the psterior credible regions for both
series. We consider the pointwise and uniform regions for the VNPC
procedure with parametric working model order $5$. For the definitions
of these two regions, see ??? (upcoming thesis). Furthermore, the
cosspectrum and the quadrature spectrum are also given.

<div class="figure" style="text-align: center">

<img src="SOI_Fit.jpg" alt="Spectral estimates for the SOI series by the VNPC(p) procedure (red line) and the VAR(p) procedure (black dashed line). The periodogram is given in grey." width="80%" />

<p class="caption">

Figure 3.2: Spectral estimates for the SOI series by the VNPC(p)
procedure (red line) and the VAR(p) procedure (black dashed line). The
periodogram is given in grey.

</p>

</div>

<div class="figure" style="text-align: center">

<img src="Rec_Fit.jpg" alt="Spectral estimates for the recruitment series by the VNPC(p) procedure (red line) and the VAR(p) procedure (black dashed line). The periodogram is given in grey." width="80%" />

<p class="caption">

Figure 3.3: Spectral estimates for the recruitment series by the VNPC(p)
procedure (red line) and the VAR(p) procedure (black dashed line). The
periodogram is given in grey.

</p>

</div>

<div class="figure" style="text-align: center">

<img src="SOI_Rec_VNPC_05.jpg" alt="Posterior credible regions for the SOI and recruitment series for the VNPC procedure with a parametric working model with order 5. Pointwise 90% region is visualised in shaded pink and uniform 90% region is in shaded blue. The posterior median is given by the black solid line and the periodogram is shown in grey." width="80%" />

<p class="caption">

Figure 3.4: Posterior credible regions for the SOI and recruitment
series for the VNPC procedure with a parametric working model with order
5. Pointwise 90% region is visualised in shaded pink and uniform 90%
region is in shaded blue. The posterior median is given by the black
solid line and the periodogram is shown in grey.

</p>

</div>

To investigate the correlation between the two series, we consider the
squared coherence. The coherence $\kappa$ is defined as
$$\kappa(\omega|\boldsymbol{f})=\frac{f_{12}(\omega)}{(f_{11}(\omega)f_{22}(\omega))^{1/2}},\quad 0\leq\omega\leq\pi.$$
Figure <a href="#fig:Fig5">3.5</a> demonstrates the squared coherence of
the VNPC procedure with the parametric working model order 5.

<div class="figure" style="text-align: center">

<img src="SOI_Rec_Squared_Coherency_VNPC_5__chol.jpg" alt="Estimated squared coherence for the SOI and recruitment series by the VNPC procedure with a parametric working model with order 5. The posterior median is given by the black line and the pointwise 90% credible region is in shaded pink." width="80%" />

<p class="caption">

Figure 3.5: Estimated squared coherence for the SOI and recruitment
series by the VNPC procedure with a parametric working model with order
5. The posterior median is given by the black line and the pointwise 90%
credible region is in shaded pink.

</p>

</div>

# 4 References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-mittelbach2012" class="csl-entry">

Mittelbach, M., Matthiesen, B., & Jorswieck, E. A. (2012). Sampling
uniformly from the set of positive definite matrices with trace
constraint. *IEEE Trans. Signal Process.*, *60(5)*, 2167–2179.
[10.1109/TSP.2012.2186447](https://10.1109/TSP.2012.2186447)

</div>

<div id="ref-shumway2011" class="csl-entry">

Shumway, R. H., & Stoffer, D. S. (2011). *Time series analysis and its
applications with r examples* (3rd ed.). Springer.

</div>

<div id="ref-stoffer2022" class="csl-entry">

Stoffer, D. (2022). *Astsa: Applied statistical time series analysis*. R
package version 1.16. <https://cran.r-project.org/web/packages/astsa/>

</div>

</div>
