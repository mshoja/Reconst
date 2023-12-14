# Reconst
This repository contains MATLAB codes to make inference about the parameters of Langevin models as well as lévy-driven stochastic models using univariate time series of observations. The code LangevinReconst.m
fits a Langevin model to data while the code LevyReconst.m fits a more complex stochastic model driven by a lévy noise. You can find more details on the models, the code inputs and outputs in the upper part of 
the codes added as comments. In particular, this repository contains the codes and data for the following paper

Stochastic regimes can hide the attractors in data, reconstruction algorithms can reveal them (under review)
Babak M. S. Arani, Stephen R. Carpenter, Egbert H. van Nes, Ingrid A. van de Leemput, Chi Xu, Pedro G. Lind, and Marten Scheffer

In order to generate Figures 2-5 in the above paper you should run the main code called Allfigures.m. The repository contains an ice-core calcium dataset Cadata.csv from the following link

Fischer, Hubertus; Rasmussen, Sune Olander; Fuhrer, Katrin (2022): High-resolution impurity data from the GRIP ice core. PANGAEA, https://doi.org/10.1594/PANGAEA.942777

Additional information about ice-core data can be found in the following publications

Ditlevsen, Peter D. "Observation of α‐stable noise induced millennial climate changes from an ice‐core record." Geophysical Research Letters 26.10 (1999): 1441-1444., https://doi.org/10.1029/1999GL900252
Rasmussen, S. O., Dahl-Jensen, D., Fischer, H., Fuhrer, K., Hansen, S. B., Hansson, M., Hvidberg, C. S., Jonsell, U., Kipfstuhl, S., Ruth, U., Schwander, J., Siggaard-Andersen, M.-L., Sinnl, G., Steffensen, J. P., Svensson, A. M., and Vinther, B. M.: Ice-core data used for the construction of the Greenland Ice-Core Chronology 2005 and 2021 (GICC05 and GICC21), Earth Syst. Sci. Data, 15, 3351–3364, https://doi.org/10.5194/essd-15-3351-2023, 2023.

