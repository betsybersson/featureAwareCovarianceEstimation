# Feature Aware Covariance Estimation

Replication code for the paper:

> Bersson, E., Hoffman, K., Stapleton, H. M., and Dunson, D. B. (2025). *Feature aware covariance estimation, with application to mixtures of chemical exposures.* arXiv:2504.08220. https://arxiv.org/abs/2504.08220

Code to implement the feature aware covariance estimation prior for factor models, proposed in `Feature aware covariance estimation, with application to mixtures of chemical exposures'.

Files:
 - VIGNETTE.html: Worked example showing how to fit the FACE model to your own data. Start here.
 - SIMULATION.R: Simulation replication codes for tables in Secion 4 of the article.
 - SYNTHETIC_APPLICATION.R: Implementation of the FACE model for a synthetic mixture of exposure dataset.
 - SYNTHETIC_LOD_EXPERIMENT.R: Replication codes for the LOD experiment in Section 5.2, using the synethetic dataset.
 - functions.R: FACE gibbs sampler and accompnaying functions.
 - helpers.R: Miscallaneous functions used throughout.
