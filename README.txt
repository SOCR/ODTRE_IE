Supplementary information / reproducible research files for the paper 

Title: Optimal dynamic treatment regime estimation using information extraction from unstructured clinical text

Authors: Nina Zhou, Robert Brook, Ivo Dinov, and Lu Wang

Code was written by Nina Zhou.
In case of questions or comments please contact zhounina@umich.edu.

The main part of the code was written/evaluated in R using RStudio (Version 1.3.959) with the following software versions:

R version 4.0.2 (2020-06-22)Platform: x86_64-apple-darwin17.0 (64-bit)Running under: macOS  10.16Matrix products: defaultLAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dyliblocale:[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8attached base packages:[1] stats     graphics  grDevices utils     datasets  methods   base     other attached packages: [1] latex2exp_0.4.0     nnet_7.3-14         rpart_4.1-15        missForest_1.4      [5] itertools_0.1-3     iterators_1.0.13    foreach_1.5.1       randomForest_4.6-14 [9] ggplot2_3.3.0       dplyr_1.0.6        loaded via a namespace (and not attached): [1] pillar_1.6.1     compiler_4.0.2   tools_4.0.2      digest_0.6.27    lifecycle_1.0.0  [6] tibble_3.1.2     gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.11     rstudioapi_0.11 [11] DBI_1.1.0        parallel_4.0.2   withr_2.2.0      stringr_1.4.0    generics_0.1.0  [16] vctrs_0.3.8      grid_4.0.2       tidyselect_1.1.1 glue_1.4.2       R6_2.5.0        [21] fansi_0.5.0      farver_2.0.3     purrr_0.3.4      blob_1.2.1       magrittr_2.0.1  [26] scales_1.1.0     codetools_0.2-16 ellipsis_0.3.2   assertthat_0.2.1 colorspace_1.4-1[31] labeling_0.3     utf8_1.2.1       stringi_1.4.6    munsell_0.5.0    crayon_1.4.1    

Part of the data preprocessing code was written in Python (Jupiter Notebook) with the following software versions:

Python version: Python 3.8.6 (v3.8.6:db455296be, Sep 23 2020, 13:31:39) [Clang 6.0 (clang-600.0.57)]
Platform: macOS-10.16-x86_64-i386-64bit

Jupyter notebook information:
  IPython 7.19.0
  jupyter_client 6.1.7
  jupyter_core 4.7.0
  jupyterlab 2.2.9
  notebook 6.1.5

Imported libraries:
  numpy 1.21.1
  pandas 1.3.1
  re 2.2.1

To access the full dataset (MIMIC III) used in the simulation and application sections in this paper, please complete the following steps:
	1. Complete CITI training course
	2. Create a PhysioNet account
	3. Submit a request to access MIMIC III data
	4. Accessing the archive (via download or using a Cloud service such as Google BigQuery or Amazon Web Services) 

More information about MIMIC III access is available on this website https://mimic.physionet.org/gettingstarted/access. More information about the data is available online https://physionet.org/content/mimiciii/1.4/. 


This folder contains the following data and files that can be used to reproduce all tables and figures of the simulation section in the paper.
It contains four subfolders containing the following files:

./main_program
	
  (1) master.R

  This code reproduces the two-stage and three-stage simulation studies and optimal DTR estimation in the application section in the paper. The code also generates Table 1, Figure 1, and Figure 2. It can be executed by command: R CMD BATCH --no-save master.R master_log.txt, after changing the working directory to ./main_program. Results will be saved in the same folder. 

  The simulation code section loads intermediate results of objects aa, bb, cc, and dd for faster computation. However, readers can also run the original code for themselves, which is commented out in master.R after load().

  The application code section loads preprocessed data application_data_day1.csv and application_data_day2.csv, which can be generated using SQL_queries_hyper_emergency_application.txt (on GCP), hyper_ie.ipynb (IE), and hyper_emergency_data_preprocessing.R. 
	
./data

This subfolder contains final text related simulation/application data and the Python and SQL codes to generate these datasets.

 ./data/data_simulation

  (1) SQL queries.txt

  SQL queries used in Google Cloud Platform (GCP) to extract data from MIMIC III. To request access,   please contact PhysioNet for review.

  (2) Simulation_text_info_extraction.ipynb

  Python notebook for simulating clinical text with smoking status and perform information extraction using simulated data. The input data for this notebook is downloaded from GCP as CSV files. The output is the final sim_data.csv for simulation.
   
  (3) sim_data.csv
  
  Simulation data that contains simulated and structuralized weight and smoking information from text. Obtained from Simulation_text_info_extraction.ipynb.

  (4) two_stage_sim.Rdata

  Intermediate results for two stage simulation (1000 replicates). Loaded by master.R for saved objects aa and bb.

  (5) three_stage_sim.Rdata

  Intermediate results for three stage simulation (1000 replicates). Loaded by master.R for saved objects cc and dd.

 ./data/data_hyper_emergency

  (1) SQL_queries_hyper_emergency_application.txt
  
  SQL queries used in Google Cloud Platform (GCP) to extract data from MIMIC III. To request access,   please contact PhysioNet for review. The queries will generate 8 CSV files for data preprocessing.

  (2) hyper_ie.ipynb

  Python notebook for clinical text information extraction. The notebook extracts patients smoking status, weight, and height information from notes.csv (generated from SQL code). The notebook outputs hyper_text_ie.csv to the same folder as structured patient information for data preprocessing.

  (3) hyper_emergency_data_preprocessing.R

  This R code preprocesses the CSV files generated by the SQL code and Python notebook. It outputs preprocessed data application_data_day1.csv and application_data_day2.csv for optimal DTR estimation in the application section. The code can be run using: R CMD BATCH --no-save hyper_emergency_data_preprocessing.R, after changing working directory to ./data/data_hyper_emergency.

  (4) application_data_day1.csv (please delete after review)

  Day one preprocessed application data generated from hyper_emergency_data_preprocessing.R.

  (5) application_data_day2.csv (please delete after review)

  Day two preprocessed application data generated from hyper_emergency_data_preprocessing.R.

./function

  (1) Functions.R

  Contains functions written in R that estimates optimal DTR using Tree-based Reinforcement Learning (T-RL) and simulation functions. All functions are applied in ./main_program/master.R.

  (2) plot_tree.R

  Contains a function for plotting optimal DTR from DTRtree objects. This function is applied in ./main_program/hyper_emergency_dtr_estimation.R. 

./results

  (1) Table1.txt
  
  Records the R console printouts for Table 1 in the paper from ./main_program/master.R.

  (2) Figure1_two_stage.pdf
  
  Outputted Figure 1 - two stage study from ./main_program/master.R (object A).

  (3) Figure1_three_stage.pdf

  Outputted Figure 1 - three stage study from ./main_program/master.R (object B).

