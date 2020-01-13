To run the analysis on the simulated data, enter the simulations folder and run the files in the following order:

1_modelling.R
2_moransI.R

The first file fits the models, extracts the parameters and stores everything.
The second file calculates Moran's I (can be slow).

Likewise, to run the analysis on the case study data, enter the case study folder and run the files in the following order:

1_uniqueEcotopes.R
2_convertEcotopes.R
3_preprocessingEcotopes1.R
4_preprocessingEcotopes2.R
5_modellingCrossValidation.R
6_modellingLambda.R
7_AUC.R

The first file generates a list of all unique ecotopes that we then classified in suitable and unsuitable
The second file applies this classification.
The third file selects the presences of the second survey, randomly samples pseudo-absences and calculates the environmental variables for those.
The fourth file calculates the autocovariate.
The fifth file fits the models based on all but one year each (7 folds).
The sixth file fits the models for different values of lambda.
The seventh file calculates the AUC from the different folds from the cross validation.

The scripts rely on:
scenarios.json and settings.json: Files where the directories are stored. These are shared between the scripts.
modelList.csv: A file which combines scenario and run into one number (ranging from 1 to 90). This is used to be able to fit the models in parallel with command line input.
external_data/casestudy/classifiedUniqueEcotopeList.csv: A file with all ecotopes and their classification in suitable (1) and unsuitable (0).

In the seeds folder, we stored all the seeds that were generated in the fitting of the models and their evaluation. To reproduce our results, these have to be used instead of the usual generation of new random seeds.
