# birdsOnIGS

This repository contains the code used to conduct the analysis included in the #Birds On IGS# paper.

**Code Structure**

The code includes all necessary steps for data processing, model fitting, and prediction.

In the model fitting, three GLMMs were developed, each corresponding to one of the three ecological variables interrogated in the paper (foraging guild, foraging domain, conservation status).

This repository includes only one iteration of the final prediction step, which relates to a single ecological variable (foraging guild). The other two iterations follow the same structure but use their respective GLMMs and modify the ecological input variables accordingly.

**Reproducing Results**

To extend the analysis to all three environmental variables, update the prediction step by:

Selecting the appropriate GLMM from the included models.
Adjusting the ecological input variables used for prediction to match the ecological variable of interest.
