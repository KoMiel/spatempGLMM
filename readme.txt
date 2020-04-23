To run the analysis on the simulated data, enter the simulations folder and run the files in the following order:

1_preprocessing.R
2_modelling.R
3_modelling_singleYear.R
4_model_evaluation.R
5_moransI.R

The first file takes a subsample of all data and calculates autocovariates for different range parameters.
The second file fits the models for the different range parameters and stores predictions and slope estimates.
The third file does the same, but only for single sampling events and only for the autologistic spatial error model.
The fourth file summarizes the results of the models so that they can be easily viewed and plotted afterwards.
The fifth file calculates correlograms to track remaining residual autocorrelation.

Likewise, to run the analysis on the case study data, enter the case study folder and run the files in the following order:

1_uniqueEcotopes.R
2_convertEcotopes.R
3_preprocessing.R
4_modelling.R
5_modelling_singleYear.R
6_model_evaluation.R

The first file generates a list of all unique ecotopes that we then manually classified in suitable and unsuitable
The second file applies this classification.
The third file selects the presences of the second survey, randomly samples pseudo-absences and calculates the environmental variables for those.
The fourth file fits the models based on all but one year each (7 folds).
The fifth file does the same, but only uses data of one year for fitting and the rest for evaluation.
The sixth file summarizes the results so that they can be easily viewed and plotted afterwards.

The scripts rely on:
scenarios.json and settings.json: Files where the directories are stored. These are shared between the scripts.
modelList.csv: A file which combines scenario and run into one number (ranging from 1 to 90) for the simulation. This is used to be able to fit the models in parallel.
external_data/casestudy/classifiedUniqueEcotopeList.csv: A file with all ecotopes and their classification in suitable (1) and unsuitable (0).

In the seeds folder, we stored all the seeds that were generated in the fitting of the models and their evaluation. To reproduce our results, these have to be used instead of the usual generation of new random seeds.
