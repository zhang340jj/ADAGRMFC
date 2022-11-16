---------------------------------------------------------------------------

*** How to run ADA-GRMFC? ***

The starting point for running ADA-GRMFC is:

RUN_Cross_Validation.m
    This is for estimating the prediction performance of ADA-GRMFC using cross 
    validation. More precisely, 5 repetitions of 10-fold cross validation 
    are performed, and then the AUPRs from the five repetitions are 
    averaged to give the final AUPR.

To run, simply access the above file in MATLAB, and press F5.


In the above file, there are many options that may be set, including:

> classifier
> use_WKNKN
> cv_setting
> cross validation parameters: m, n

---------------------------------------------------------------------------

*** Algorithm entry points ***

The main entry scripts for the different algorithms are:

> alg_blm_nii.m
> alg_rls_wnn.m
> alg_cmf.m
> alg_grmf.m
> alg_wgrmf.m
> alg_wknkn.m
> alg_SRCMF
> alg_ADAGRMFC

Notes:  alg_ADAGRMFC is ADA-GRMFC in codes.

Some parameters may be manually set to a fixed value, while others are 
estimated via nested CV. These parameters may be found in the above files, 
or if their values are determined by nested CV, they may be found in 
scripts that have the suffix "parest_":

> alg_grmf_parest.m
> alg_cmf_parest.m
> alg_rls_wnn_parest.m

However, for WKNKN in particular, the parameters may be found in 
RUN_Cross_Validation.m as WKNKN may be used as a preprocessing method by 
any of the other methods.

---------------------------------------------------------------------------


---------------------------------------------------------------------------

*** Relevant Publication ***

Graph regularized non-negative matrix factorization with prior knowledge consistency
constraint for drug-target interaction prediction
Junjun Zhang and Minzhu Xie

Drug-target interaction prediction with graph-regularized matrix factorization
Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh

Xu, Y. Y., Yin, W. T., Wen, Z. W., & Zhang, Y. (2012). 
An alternating direction algorithm for matrix completion with nonnegative factors.
 Frontiers of Mathematics in China, 7(2), 365-384. doi:10.1007/s11464-012-0194-5

---------------------------------------------------------------------------