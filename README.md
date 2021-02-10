# Simulated examples for learning individualized treatment rules using electronic health records data
The programs and data in this repository are used to identify latent subject subgroups and learn individualized treatment rules (ITRs) using electronic health records (EHRs) data. 
Since we were not allowed to release our EHRs, we implemented proposed algorithm to simulated datasets for the purpose of illustration. 
Please refer to the following two sections about how these files are used to illustrate our algorithm. 
The locations and summaries of the programs and data are provideds in this document. 
Also, you can find annotations and comments in R scripts.

## 1   Identifying latent subgroups
### 1.1  Flow chart
![Web Figure 1](https://github.com/jitonglou/EHR_ITR/blob/main/README_figures/TRule_simdata_flowchart_1.png?raw=true "Web Figure 1")

### 1.2  Simulation and parameter estimation settings
In [this explanatory work](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_GenerateData_seed72.R), we simulated a dataset of two health markers for 10,000 subjects. 
For the *i*th subject, we assumed *Y<sub>i1</sub>(t)* was Gaussian distributed and *Y<sub>i2</sub>(t)* was Bernoulli distributed. 
Thus, the inverse functions of canonical link functions *g<sub>1</sub><sup>-1</sup>(z)=z* and *g<sub>2</sub><sup>-1</sup>(z)=e<sup>z</sup>/(1+e<sup>z</sup>)*. 
Since the distribution of *Y<sub>i1</sub>(t)* has a dispersion parameter, we set *&phi;<sub>i1</sub>(t)=0.1*.
We generated two covariates *X<sub>i1</sub> ~ Normal(0,1/3)* and *X<sub>i2</sub> ~ Bernoulli(0.5)-0.5*. 
Thus, *__X__<sub>i</sub>=(1, X<sub>i1</sub>, X<sub>i2</sub>)<sup>T</sup>* was a 3-dimensional vector of baseline covariates.</br> 
The maximum observation time *T<sub>i</sub>* for each subject was set to 12 (days). 
The measured time points for simulated markers were generated from two Poisson processes *dN<sub>ik</sub>(t)*, *k=1, 2*.
Their intensity functions are set to
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{E}\[dN_{i1}(t)|\boldsymbol{X}_i\]=\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i11}(t)-0.1L_{i12}(t)\}dt" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{E}\[dN_{i1}(t)|\boldsymbol{X}_i\]=\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i11}(t)-0.1L_{i12}(t)\}dt" title="\mathbb{E}\[dN_{i1}(t)|\boldsymbol{X}_i\]=\exp\{0.5X_{i1}+0.25X_{i2}+0.3L_{i11}(t)-0.1L_{i12}(t)\}dt" /></a>
</p>
&nbsp;&nbsp;&nbsp;&nbsp;and
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{E}\[dN_{i2}(t)|\boldsymbol{X}_i\]=1.2\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i21}(t)-0.1L_{i22}(t)\}dt." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{E}\[dN_{i2}(t)|\boldsymbol{X}_i\]=1.2\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i21}(t)-0.1L_{i22}(t)\}dt." title="\mathbb{E}\[dN_{i2}(t)|\boldsymbol{X}_i\]=1.2\exp\{0.5X_{i1}+0.25X_{i2}+0.3L_{i21}(t)-0.1L_{i22}(t)\}dt." /></a>
</p>

Thus, *__&gamma;__<sub>1</sub>=__&gamma;__<sub>2</sub>=(0.5, 0.25)<sup>T</sup>* and *__&eta;__<sub>1</sub>=__&eta;__<sub>2</sub>=(0.3, -0.1)<sup>T</sup>*.
Additionally, in the intensity functions, we let *L<sub>ik1</sub>(t)=1* if there exists measurements of *k*th marker in \[*t*-3,*t*); otherwise, *L<sub>ik1</sub>(t)=0*.
If *L<sub>ik1</sub>(t)=1*, then *L<sub>ik2</sub>(t)* is average value of all *Y<sub>ik</sub>(t)* in \[*t*-3,*t*); otherwise, *L<sub>ik2</sub>(t)=0*.</br>

The true values of *__&beta;__<sub>k</sub>(t)* were assumed to be
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{equation}\color{Black}{&space;\begin{pmatrix}&space;\boldsymbol{\beta}_1^T(t)&space;\\&space;\boldsymbol{\beta}_2^T(t)&space;\end{pmatrix}&space;=&space;\begin{pmatrix}&space;-1.36&plus;\frac{t}{10}&space;&&space;\sin(0.76&plus;t)&space;&&space;\cos(-0.3&plus;t)&space;\\&space;\cos(-0.25&plus;t)&space;&&space;0.37&plus;\frac{t}{10}&space;&&space;\sin(-0.68&plus;t)&space;\end{pmatrix}.}&space;\end{equation*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{equation}\color{Black}{&space;\begin{pmatrix}&space;\boldsymbol{\beta}_1^T(t)&space;\\&space;\boldsymbol{\beta}_2^T(t)&space;\end{pmatrix}&space;=&space;\begin{pmatrix}&space;-1.36&plus;\frac{t}{10}&space;&&space;\sin(0.76&plus;t)&space;&&space;\cos(-0.3&plus;t)&space;\\&space;\cos(-0.25&plus;t)&space;&&space;0.37&plus;\frac{t}{10}&space;&&space;\sin(-0.68&plus;t)&space;\end{pmatrix}.}&space;\end{equation*}" title="\begin{equation}\color{Black}{ \begin{pmatrix} \boldsymbol{\beta}_1^T(t) \\ \boldsymbol{\beta}_2^T(t) \end{pmatrix} = \begin{pmatrix} -1.36+\frac{t}{10} & \sin(0.76+t) & \cos(-0.3+t) \\ \cos(-0.25+t) & 0.37+\frac{t}{10} & \sin(-0.68+t) \end{pmatrix}.} \end{equation*}" /></a>
</p>

Furthermore, we assumed the correlation structure of multivariate latent processes to be 
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{equation*}\color{Black}{&space;\boldsymbol{\Omega}(t)=&space;\begin{pmatrix}&space;0.5&space;&&space;-0.25\\&space;-0.25&space;&&space;0.5&space;\end{pmatrix}&space;&plus;&space;\frac{1}{20}&space;\begin{pmatrix}&space;\sin(t&plus;2)&space;&&space;\cos(t-0.5)\\&space;\cos(t-0.5)&space;&&space;\sin(t&plus;3)&space;\end{pmatrix}}&space;\end{equation*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{equation*}\color{Black}{&space;\boldsymbol{\Omega}(t)=&space;\begin{pmatrix}&space;0.5&space;&&space;-0.25\\&space;-0.25&space;&&space;0.5&space;\end{pmatrix}&space;&plus;&space;\frac{1}{20}&space;\begin{pmatrix}&space;\sin(t&plus;2)&space;&&space;\cos(t-0.5)\\&space;\cos(t-0.5)&space;&&space;\sin(t&plus;3)&space;\end{pmatrix}}&space;\end{equation*}" title="\begin{equation*}\color{Black}{ \boldsymbol{\Omega}(t)= \begin{pmatrix} 0.5 & -0.25\\ -0.25 & 0.5 \end{pmatrix} + \frac{1}{20} \begin{pmatrix} \sin(t+2) & \cos(t-0.5)\\ \cos(t-0.5) & \sin(t+3) \end{pmatrix}} \end{equation*}" /></a>
</p>
 &nbsp;&nbsp;&nbsp;&nbsp;and
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\operatorname{Cov}(\boldsymbol{\epsilon}(t),&space;\boldsymbol{\epsilon}(s))&space;=&space;\exp\left\{-\left\(\frac{t-s}{b}\right)^2\right\}\times&space;\frac{\boldsymbol{\Omega}(t)&plus;\boldsymbol{\Omega}(s)}{2}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\operatorname{Cov}(\boldsymbol{\epsilon}(t),&space;\boldsymbol{\epsilon}(s))&space;=&space;\exp\left\{-\left\(\frac{t-s}{b}\right)^2\right\}\times&space;\frac{\boldsymbol{\Omega}(t)&plus;\boldsymbol{\Omega}(s)}{2}," title="\operatorname{Cov}(\boldsymbol{\epsilon}(t), \boldsymbol{\epsilon}(s)) = \exp\left\{-\left\(\frac{t-s}{b}\right)^2\right\}\times \frac{\boldsymbol{\Omega}(t)+\boldsymbol{\Omega}(s)}{2}," /></a>
</p>

where *b=0.5*.</br>
In the [uploaed simulation dataset](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_GenerateData_seed72_workspace.RData), some of the 10000 subjects has 0 measurement for both health markers. Thus, we excluded these subjects and the sample size for parameter estimation is 8338 subjects. We use the scaled Epanechnikov kernel
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=K_{h_{1n}}(u)=\frac{3}{4h_{1n}}\left[1-\left\(\frac{u}{h_{1n}}\right\)^2\right\]_{&plus;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K_{h_{1n}}(u)=\frac{3}{4h_{1n}}\left[1-\left\(\frac{u}{h_{1n}}\right\)^2\right\]_{&plus;}" title="K_{h_{1n}}(u)=\frac{3}{4h_{1n}}\left[1-\left\(\frac{u}{h_{1n}}\right\)^2\right\]_{+}" /></a>
</p>

as the kernel function to estimate *__&beta;__<sub>k</sub>(t)*. Similarly, the kernel function for estimating *__&Omega;__(t)* was set to the product of two scaled univariate Epanechnikov kernels
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\widetilde{K}_{h_{2n}}(u_1,u_2)=\frac{9}{16h_{2n}^2}\left\{1-\left\(\frac{u_1}{h_{2n}}\right\)^2\right\}_{&plus;}\left\{1-\left\(\frac{u_2}{h_{2n}}\right\)^2\right\}_{&plus;}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widetilde{K}_{h_{2n}}(u_1,u_2)=\frac{9}{16h_{2n}^2}\left\{1-\left\(\frac{u_1}{h_{2n}}\right\)^2\right\}_{&plus;}\left\{1-\left\(\frac{u_2}{h_{2n}}\right\)^2\right\}_{&plus;}." title="\widetilde{K}_{h_{2n}}(u_1,u_2)=\frac{9}{16h_{2n}^2}\left\{1-\left\(\frac{u_1}{h_{2n}}\right\)^2\right\}_{+}\left\{1-\left\(\frac{u_2}{h_{2n}}\right\)^2\right\}_{+}." /></a>
</p>

We extended the data-adaptive method in [Cao et al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4643299/) and selected the optimal bandwidths among {0.1,0.2,...,0.5}. We found *h=0.3* for *__&beta;__<sub>k</sub>(t)* and *h=0.2* for *__&Omega;__(t)* were close to the optimal values of bandwidths. This set of bandwidth was used in the uploaded simulation data. Since the proposed estimation method was expected to have more stable performance at time points that not on two ends, we used the estimated *__&beta;__<sub>k</sub>(t)* and *__&Omega;__(t)* at time points *t=3,4,...,11* to identify latent subgroups.

### 1.3  Results
   - Following [TRule_simdata_GenerateData_seed72.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_GenerateData_seed72.R), one could obtain estimated *__&gamma;__<sub>k</sub>* and *__&eta;__<sub>k</sub>*, *k=1, 2*, as
   
   > | Marker | Parameter | True value | Estimator |
   > | ------ | --------- | ---------- | --------- |
   > | *Y<sub>1</sub>* | *&gamma;<sub>11</sub>* | 0.5  | 0.456  |
   > | Continuous      | *&gamma;<sub>12</sub>* | 0.25 | 0.229  |
   > |                 | *&eta;<sub>11</sub>*   | 0.3  | 0.297  |
   > |                 | *&eta;<sub>12</sub>*   | -0.1 | -0.102 |
   > | *Y<sub>2</sub>* | *&gamma;<sub>21</sub>* | 0.5  | 0.456  |
   > | Binary          | *&gamma;<sub>22</sub>* | 0.25 | 0.228  |
   > |                 | *&eta;<sub>21</sub>*   | 0.3  | 0.304  |
   > |                 | *&eta;<sub>22</sub>*   | -0.1 | -0.100 |

   - Following [TRule_simdata_est_case0_seed72_tp10.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_est_case0_seed72_tp10.R), one could obtain estimated *__&beta;__<sub>k</sub>(t)*, *k=1, 2*, and *__&Omega;__(t)* at *t=10* as
   > | Marker | Parameter | True value | Estimator |
   > | ------ | --------- | ---------- | --------- |
   > | *Y<sub>1</sub>* | *&beta;<sub>10</sub>*   | -0.360 | -0.363 |
   > | Continuous      | *&beta;<sub>11</sub>*   | -0.972 | -1.044 |
   > |                 | *&beta;<sub>12</sub>*   | -0.962 | -0.945 |
   > | *Y<sub>2</sub>* | *&beta;<sub>20</sub>*   | -0.948 | -0.923 |
   > | Binary          | *&beta;<sub>21</sub>*   | 1.370  | 1.418 |
   > |                 | *&beta;<sub>22</sub>*   | 0.105  | 0.093 |

 &nbsp;&nbsp;&nbsp;&nbsp;and
   
   > | Parameter | True value | Estimator |
   > | --------- | ---------- | --------- |
   > | *&sigma;<sub>11</sub>*  | 0.473 | 0.415 |
   > | *&sigma;<sub>22</sub>*  | 0.521 | 0.486 |
   > | *&sigma;<sub>12</sub>*  | -0.300 | -0.261 |

   - Following [TRule_simdata_DistMat_case0_seed72.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_DistMat_case0_seed72.R), one could obtain the hierarchical clustering results as Web Figure 3. We classified the 8338 subjects to 4 subgroups. However, subgroup 4 only has one subject, ID = 1814. We checked this subject and found it has an outlier measurement of health marker 1 at *t=6.744*, so we excluded this subject in the next workflow.
   ![Web Figure 3](https://github.com/jitonglou/EHR_ITR/blob/main/README_figures/TRule_simdata_DistMat_rectdend.png?raw=true "Web Figure 3")

### 1.4  Programs
   - [codes/functions_simdata.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/functions_simdata.R)\
   R functions for generating simulated data.
   - [codes/functions_estparams.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/functions_estparams.R)\
   R functions for estimating parameters (intensity parameters *&gamma;<sub>k</sub>* and *&eta;<sub>k</sub>*, regression coefficients *&beta;<sub>k</sub>(t)*, covariance matrix of latent processes *&Omega;(t)*).
   - [codes/functions_DistMat.cpp](https://github.com/jitonglou/EHR_ITR/blob/main/codes/functions_DistMat.cpp)\
   C++ functions for calculating the similarity between each pair of subjects *S<sub>ij</sub>*.
   - [codes/TRule_simdata_GenerateData_seed72.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_GenerateData_seed72.R)\
   Explanatory R scripts for generating a simulated dataset and estimating intensity parameters. R random seed=72.
   - [codes/TRule_simdata_est_case0_seed72_tp10.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_est_case0_seed72_tp10.R)\
   Explanatory R scripts for estimating regression coefficients and the covariance matrix. Time point *t*=10, bandwidth *h<sub>1n</sub>=0.3*, bandwidth *h<sub>2n</sub>=0.2*.  
   - [codes/TRule_simdata_DistMat_case0_seed72.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_DistMat_case0_seed72.R)\
   Explanatory R scripts for calculating the similarity between subjects and identifying latent subgroups.
   
### 1.5  Data
   - [rdata/TRule_simdata_GenerateData_seed72_workspace.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_GenerateData_seed72_workspace.RData)\
   R workspace for estimating intensity parameters, regression coefficients, and covariance matrix. 
   - [rdata/TRule_simdata_GenerateData_seed72_intensity_est.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_GenerateData_seed72_intensity_est.RData)\
   R workspace of estimated intensity parameters *&gamma;<sub>k</sub>* and *&eta;<sub>k</sub>*.
   - [rdata/TRule_simdata_est_case0_seed72.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_est_case0_seed72.RData)\
   R workspace of estimated regression coefficients *&beta;<sub>k</sub>(t)* and covariance matrix *&Omega;(t)* across time points.
   - [rdata/TRule_simdata_DistMat_case0_seed72_dendrogram.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_DistMat_case0_seed72_dendrogram.RData)\
   R workspace of hierarchical clustering results based the similarity metric between each pair of subjects *S<sub>ij</sub>*.
   - [TRule_simdata_DistMat_rectdend.png](https://github.com/jitonglou/EHR_ITR/blob/main/README_figures/TRule_simdata_DistMat_rectdend.png)\
   Visualization of identified latent subgroups for 8338 subjects and 4 subgroups.
   
## 2  Learning ITRs
### 2.1  Flow chart
![Web Figure 2](https://github.com/jitonglou/EHR_ITR/blob/main/README_figures/TRule_simdata_flowchart_2.png?raw=true "Web Figure 2")

### 2.2  Simulation and parameter estimation settings
Following [TRule_simdata_ITRpreprocess_case0_seed72.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_ITRpreprocess_case0_seed72.R), one could obtain a simulated dataset for learning ITRs. The dataset contains the information of pre-treatment covariates and variables reflecting recent patterns of health markers *Z<sub>i</sub>*, treatment classes *A<sub>i</sub>*, clinical outcomes *R<sub>i</sub>*, and identified subgroup memberships.</br>
In specific, we created two variables reflecting recent patterns of health markers by calculating the average value of *Y<sub>i1</sub>(t)* and *Y<sub>i2</sub>(t)* between *t=3* and *t=11*. We denoted these two variables to *V<sub>i1</sub>* and *V<sub>i2</sub>*. Next, we simulated the outcome reward
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=R_i=1&plus;X_{i1}-X_{i2}&plus;2V_{i1}-2V_{i2}&plus;e_i," target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_i=1&plus;X_{i1}-X_{i2}&plus;2V_{i1}-2V_{i2}&plus;e_i," title="R_i=1+X_{i1}-X_{i2}+2V_{i1}-2V_{i2}+e_i," /></a>
</p>

where *e<sub>i</sub>~Normal(0,1/3)*. Then, the *i*th subject was randomly assigned a treatment *A<sub>i</sub>* in {A, B, C} with equal probabilities. The resulting dataset has 8330 subjects of 3 subgroups. </br>
We let *__Z__<sub>i</sub>=(X<sub>i1</sub>, X<sub>i2</sub>, V<sub>i1</sub>, V<sub>i2</sub>)<sup>T</sup>* and denoted *u* to a certain treatment in {A, B, C}. Following [TRule_simdata_ITRlearn_case0_seed72_group3.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_ITRlearn_case0_seed72_group3.R), in subgroup 3, we estimated the propensity scores *&pi;(__Z__<sub>i</sub>)=P(A<sub>i</sub>=u|__Z__<sub>i</sub>)* by a 10-fold cross-validation random forest with 3 repeats. Next, we estimated the prognostic scores *&psi;(__Z__<sub>i</sub>)=E(R<sub>i</sub>|__Z__<sub>i</sub>)* by a gradient boosting model with 5000 trees of which maximum depth was 4. Lastly, we used *__H__<sub>i</sub>=(__Z__<sub>i</sub>, &pi;(__Z__<sub>i</sub>), &psi;(__Z__<sub>i</sub>))*, *A<sub>i</sub>*, and *R<sub>i</sub>* to estimate ITRs for this subgroup. </br>
After the estimation, a certain treatment rule *D(__H__)* can be evaluated by its empirical value function, which is defined as
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\sum_{u}\sum_{i}{I(A_i={D}(\boldsymbol{H}_i)=u)R_i}/{{P}(A_i=u|\boldsymbol{Z}_i)}}&space;{\sum_{u}\sum_{i}{I(A_i={D}(\boldsymbol{H}_i)=u)}/{{P}(A_i=u|\boldsymbol{Z}_i)}}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\sum_{u}\sum_{i}{I(A_i={D}(\boldsymbol{H}_i)=u)R_i}/{{P}(A_i=u|\boldsymbol{Z}_i)}}&space;{\sum_{u}\sum_{i}{I(A_i={D}(\boldsymbol{H}_i)=u)}/{{P}(A_i=u|\boldsymbol{Z}_i)}}." title="\frac{\sum_{u}\sum_{i}{I(A_i={D}(\boldsymbol{H}_i)=u)R_i}/{{P}(A_i=u|\boldsymbol{Z}_i)}}&space;{\sum_{u}\sum_{i}{I(A_i={D}(\boldsymbol{H}_i)=u)}/{{P}(A_i=u|\boldsymbol{Z}_i)}}." /></a>
</p>

### 2.3  Results
Following [TRule_simdata_ITRlearn_case0_seed72_group3.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_ITRlearn_case0_seed72_group3.R), one could obtain the optimal ITRs, which are estimated by the proposed approach, of 2027 subjects in subgroup 3. Within this file, running the commands after `save(ITR_OVO, df_evf, file=file.path(wsdir,paste0(casename,"_est.RData")))`, one could investigate the optimal ITRs in terms of empirical value functions and distribution of treatments. </br>

| Obs\ITR | A | B | C | Total |
|------:|-----|------|------|------|
| __A__ | 326 | 196 | 129 | 651 |
| __B__ | 117 | 328 | 235 | 680 |
| __C__ | 74 | 268 | 354 | 696 |
| __Total__ | 517 | 792 | 718 | 2027 |

![Web Figure 4](https://github.com/jitonglou/EHR_ITR/blob/main/README_figures/TRule_simdata_ITRlearn_distribution.png?raw=true "Web Figure 4")

For example, as shown in the above contingency table, treatment A is the optimal ITR for 50.08% (326/651) of the subjects who were assigned treatment A, while 30.11% (196/651) and 19.82% (129/651) of the subjects should switch to treatment B and treatment C, respectively. Similarly, for subjects that were already prescribed treatment B and treatment C, 48.24% (328/680) and 50.86% (354/696) of them are still assigned the same treatments. However, alternative treatments might have better effects on other subjects. </br>
In general, the optimal ITRs decrease the number of subjects receiving treatment A (651 -> 517). At the same time, the optimal ITRs increase the subjects receiving treatment B (680 -> 792), while that for treatment C (696 -> 718) does not change drastically. This pattern is also presented in Web Figure 4. </br>
From Web Figure 4, the optimal ITRs achieve an empirical value function of -0.459. However, the value functions are -1.362, -1.339, -1.360 for three universal rules (i.e., “one-size-fits-all” rules) which assign treatment A only, B only, and C only to all subjects. Since the proposed method aimed to maximize the value function, we could conclude the proposed method outperforms the three universal rules in this case. </br>

### 2.4  Programs
   - [codes/functions_ITRlearn.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/functions_ITRlearn.R)\
   R functions for learning ITRs using the matched learning model for multicategory treatments.
   - [codes/TRule_simdata_ITRpreprocess_case0_seed72.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_ITRpreprocess_case0_seed72.R)\
   Explanatory R scripts for creating a simulated dataset for learning ITRs. The dataset contains the information of pre-treatment covariates and variables reflecting recent patterns of health markers *Z<sub>i</sub>*, treatment classes *A<sub>i</sub>*, clinical outcomes *R<sub>i</sub>*, and identified latent subgroups. The dataset has 8330 subjects and 3 treatments.
   - [codes/TRule_simdata_ITRlearn_case0_seed72_group3.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_ITRlearn_case0_seed72_group3.R)\
   Explanatory R scripts for calculating propensity, estimating prognostic scores, and learning ITRs within a subject subgroup. At the end of this file, we showed how to calculate the empirical value function of optimal ITRs and universal rules. Also, we visualized the distribution of treatments in observed assignments and that in optimal ITRs.
   - [codes/TRule_simdata_ITRpredict_case0_seed72_group3.R](https://github.com/jitonglou/EHR_ITR/blob/main/codes/TRule_simdata_ITRpredict_case0_seed72_group3.R)\
   Explanatory R scripts for implementing estimated ITRs to predict the optimal treatments for a new simulated dataset.

### 2.5 Data
   - [rdata/TRule_simdata_ITRpreprocess_case0_seed72_df.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_ITRpreprocess_case0_seed72_df.RData)\
   R workspace of the simulated dataset used to estimate propensity and prognostic scores. Columns of this dataset include *__Z__<sub>i</sub>*, *A<sub>i</sub>*, and *R<sub>i</sub>*.
   - [rdata/TRule_simdata_ITRlearn_case0_seed72_group3_workspace.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_ITRlearn_case0_seed72_group3_workspace.RData)\
   R workspace for learning ITR. This file includes the trained random forest models for propensity scores and trained gradient boosting machines for prognostic scores. 
   - [rdata/TRule_simdata_ITRlearn_case0_seed72_group3_est.RData](https://github.com/jitonglou/EHR_ITR/blob/main/rdata/TRule_simdata_ITRlearn_case0_seed72_group3_est.RData)\
   R workspace of optimal ITRs learned by the proposed approach. 
   - [TRule_simdata_ITRlearn_distribution.png](https://github.com/jitonglou/EHR_ITR/blob/main/README_figures/TRule_simdata_ITRlearn_distribution.png)\
   A barplot comparing the distribution of treatments in observed assignments with that in optimal ITRs.
