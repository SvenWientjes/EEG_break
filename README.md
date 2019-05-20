# EEG_break
This repository contains the MATLAB code used for writing the master thesis [_Unmixing Cognitive Strategies Using the EEG Signal_](https://www.overleaf.com/read/tkgqkrqbdgjr). The thesis is currently in progress of being written, and is supervised and graded by [dr. Leendert van Maanen](https://www.leendertvanmaanen.com) (Amsterdam, UvA) and [dr. Max Hinne](http://www.cs.ru.nl/~mhinne/) (Nijmegen, Radboud). The thesis contains the theoretical explanation of what happens within the functions written and used by this repository. This readme file contains additional information on how to replicate the results of the simulation study to evaluate the break detection method, hopefully aiding anyone who hopes to extend this method in the future.

## Dependencies
From scratch, this repository and its dependencies can be used for generating the EEG data under specific (differing) assumptions containing a structural break, estimating if and when any structural breaks happen, and evaluating the accuracy of this estimation.

### Simulating the EEG
For simulating the EEG data accurately, we need 2 'parts'. First are the underlying ERPs which correspond to the application of any specific cognitive strategy (or sleep stage, or anything else). Second is the typical 'EEG noise' which has a specific frequency spectrum. It is not simply white noise, so a specific generator is used for this part.

#### The ERPs
Event Related Potentials are small 'bumps' of electrical activity distributed over multiple electrodes of the EEG. The onset, amplitude, and latencies of these bumps can vary from trial to trial. The toolbox used to simulate this is [SEREEGA](https://github.com/lrkrol/SEREEGA): Simulating Event Related EEG Activity.

#### The EEG noise
The EEG is contaminated by many artefacts from many different and unknown sources. However, it is observed the frequency spectrum of the EEG is often the same. Therefore, simulating realistic EEG noise tries to replicate this frequency spectrum. A simulation procedure by Yeung, Bogacz, Holroyd, Nieuwenhuis & Cohen simulates the noise by adding 50 random-phase sinusoids together. Their procedure can be downloaded [here](https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator).

#### EEGLAB
For handling the data, the standard Matlab package EEGLAB is used in this thesis. The method for structural break detection also critically relies on EEGLAB. EEGLAB can be downloaded [here](https://sccn.ucsd.edu/eeglab/download.php).

### Estimating the break points
For estimating the break points, first we decompose the signal into Independent Component Analysis (ICA) models. However, one ICA is not sufficient since there is a structural break: The EEG over the full experiment is a mixture of two ICA models. Therefore, we use Adaptive Mixture Independent Component Analysis (AMICA), which estimates different ICA models for self-similar parts of the data. AMICA returns the ICA models, but critically also the probabilities of each estimated ICA model being active over time. Structural Break detection can be applied to this timeseries of model probabilities.

#### AMICA
The mixture ICA estimator used for this thesis is AMICA, which can be downloaded [here](https://sccn.ucsd.edu/~jason/amica_web.html). It can be run from EEGLAB directly, which is useful. For replication of the large-scale simulation study in this thesis, relying on the EEGLAB-GUI of AMICA is not recommended.

#### Univariate Break Detection
To estimate the break points in the AMICA model probabilities, univariate break detection as proposed [here]() is used. This method is implemented in a toolbox named 'm_break' which can be downloaded [here](http://people.bu.edu/perron/code/m-break-matlab.zip). For more implementations of break detection, look [here](http://people.bu.edu/perron/code.html).

### Further dependencies
Some of these toolboxes come with their own dependencies, which may or may not already be installed in your version of Matlab. Consult the documentation for each of these toolboxes if something is not working, the solution should lie there.

## Replication Tutorial
Here I can discuss the scripts used for actually running the simulation study. Within the functions, not much needs to change.

### Get_Run_Full
This

### Run AMICA through EEGLAB
Settings

### AMICA_analysis and Batch_Analysis
Evaluations
