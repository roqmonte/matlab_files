%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index.

% 00 Functions (32 functions).
% 1.  acf               : Estimate ACF/PAC functions.
% 2.  adf               : ADF test unit roor.
% 3.  calendar          : Create dates var/labels to plot time series data.
% 4.  chis_prb          : Computes chi-squared probability function.
% 5.  data_process      : Data processing functions.
% 6.  Data_xprl1        : Function explores all data with scatter plot.
% 7.  Data_xprl2        : Function explores y matrix vs x matrix.
% 8.  EstimateBlockVAR  : VARX estimation with block exogeneity.
% 9.  EstimateVAR       : VARX estimation.
% 10. EvalFore          : Forecast Evaluation of a set of models.
% 11. Fluctiation       : Fluctuation test of GR.
% 12. GS_conv           : Function analize convergence of Gibbs Sampler.
% 13. GS_plot_states    : Plot of state variables fom Gibbs Sampler.
% 14. hpfilter          : HP filter.
% 15. IRFar_build       : Function builds IRF for an AR(p) model.
% 16. iwpQ              : Draw Sigma from Inverse wihart.
% 17. KalmanCarterKohn  : Kalman filter and Carter Kohn algorithm.
% 18. KalmanCarerKohn_sv: Kalman filter and Carter Kohn algorithm for SV model.
% 19. KF_KimShephardChib: Procedure for mixture of log chi-2 dist.
% 20. kfilter_LJ        : Kalman filter from Hamilton's text.
% 21. LagN              : Generates matrix with lags.
% 22. loc               : Finds the id of a string in a vector of strings.
% 23. Mat_sel           : Selects data for given threshold.
% 24. minMaxMat         : Max/min of a matrix.
% 25. NWhac             : Computes the Newey & West HAC matrix.
% 26. Plot_tvp_sg2      : Plot Sg2(t).
% 27. shadedplot        : Produces shaded areas between two lines.
% 28. standardise       : Standardise a set of variables.
% 29. trimr             : Return a matrix (or vector) stripped of the specified rows.
% 30. VAR_FEVD_spillover: Predictive directional measurement of volatility spillovers, as in Diebold and Yilmaz (2012).
% 31. VAR_TestLagLength : Test lag order VAR models.
% 32. VAR_Plot_fit      : Plot fit and residuals VAR models.
% 33. VAR_TestBlockExo  : Test BlockExo hypothesis.
% 34. vec               : Create a matrix stacking the columns of a matrix.
% 35. VEC_johansen_CV   : Johansen test of cointegration.
% 36. VEC_johansen_test : Critical values Johansen test of cointegration.

% 01 Linear Models (12 functions).
% 1.  armax_mlike       : Log-Likelihood ARMA/ARMAX(p,q) model using Kalman filter.
% 2.  armax_model       : Estimation of ARMA/ARMAX(p,q) model using Kalman Filter.
% 3.  armax_opt         : Finds order of ARMA/ARMAX(p,q) model.
% 4.  ARxols            : Estimation ARx(p) model by OLS.
% 5.  ARxopt            : Finds lag order of ARx(p) model.
% 6.  FMOLS             : Fully modified OLS, Phillips and Hansen (1990).
% 7.  FOREgaols         : Forecasting ARx(p) model using GA.
% 8.  FOREols           : Forecasting ARx(p) model by OLS.
% 9.  OLSest            : Estimation of linear autoregressive model by OLS.
% 10. OLSgen            : Estimation of ARx(p) model by GA.
% 11. OLSrec            : Recursive/rolling estimation of ARx(p) model by OLS.
% 12. TVPmodel          : Time varying parameter model.
% 13. TVPmodel          : Time varying parameter model with stochastic volatility.

% 02 Nonlinear models (11 functions).
% 1.  FOREsetarx2       : Forecasting SETARx(2) model.
% 2.  FOREsetarx3       : Forecasting SETARx(3) model.
% 3.  MS_bayesian_mod   : Bayesian Markov Switching model.
% 4.  MSarp             : Estimation of MSAR(p) model by maximun like-lihood.
% 5.  MSarp_evalf       : Likelihood of MSAR(p) model.
% 6.  MSarp_evalf_flex  : Likelihood of MSAR(p) model.
% 7.  MSarp_flex        : Estimation of MSAR(p) model by maximun like-lihood.
% 8.  MSfore            : Forecasting MSAR(p) model.
% 9.  SETARx2           : Estimation SETARx(2) model.
% 10. SETARx3           : Estimation SETARx(3) model.
% 11. TARx2             : Estimation TARx(2) model.
% 12. TARx2             : Estimation TARx(3) model.

% 03 Neural Networks (05 functions).
% 1.  Eval_MLP          : SSR of Multilayer Perceptron Network.
% 2.  Eval_RAD          : SSR of Radial Basis Network.
% 3.  Eval_RID          : SSR of Ridgelet Network.
% 4.  NeuralNetwork     : Estimation Neural Networks by GA and numerical opt.
% 5.  NNforec           : Forecast Neural Networks.

% 04 Stochastic Search Method
% Agenetic_nnetwork     : Genetic algorithm code for Neural Networks.
/ Simanneal             : Simulating Anealing (check code).

% 05 VARS (09 libraries)
% - BSVAR_rec(lib)       : Structural Bayesian VAR model, recursive identification scheme.
% - BSVAR_rec_Bexo(lib)  : Structural VAR model, with block exogeneity, recursive identification scheme.
% - BSVAR_Sgn(lib)       : Structural VAR model, using signs and impacts restrictions.
% - BSVAR_Sgn_Bexo(lib)  : Structural VAR model, with block exogeneity and signs and impacts restrictions.
% - SVAR_AB(lib)         : Structural VAR model, AB-Model identification.
% - SVAR_BQ(lib)         : Structural VAR model, with Long-Run restrictions a la Blanchard-Quah (1989)
% - SVAR_rec(lib)        : Structural VAR model, recursive identification scheme.
% - SVEV_rec(lib)        : Structural VEC model, recursive identification scheme.
% - TVPVAR(lib)          : Time varying parameter Structural VAR model, recursive identification scheme.
% - TVPVAR_SV(lib)       : TVP SVAR model with stochastic volatility, recursive identification scheme.

% 06 Factor Models (02 libraries and functions)
% - DFM_v1               : DFM, with global and block/regional factors (1 layer).
% - DFM_v2               : DFM, with global, block/regional and groups factors (2 layers).
% 1. MIDAS_PCA_reg       : Mixed frequency model with PC factors and Almon polinomials.
% 2. PCA                 : Principal component analysis.
% 3. PCA_reg             : Regression with principal component.
% 4. PSL_reg             : Partial least squared regression.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%