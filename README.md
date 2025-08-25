
#### **Replication code for manuscript:**

##### Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025. Demographic trends in gender identity among adults in the United States, 2014-2021. *International Journal of Transgender Health*. Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.

---

This manuscript uses data from the US Centers for Disease Control's *Behavioral Risk Factor Surveillance System* (BRFSS) from years 2014-2021. The following instructions explain how to replicate the data retrieval, cleaning, and analyses used in this manuscript.

To recreate this analysis, you will need access to R, ideally with RStudio. Original work done in the software and operating system environments indicated below.
- RStudio Version 2022.12.0+353 (2022.12.0+353), R Version 4.1.0 (2021-05-18) in Mac OS 13.1

---

**Step 1. Clone this GitHub repository to your local environment** 

**Step 2. Create folders in your local environment** 

To run the scripts in this project, you will need to create folders to store your data. Specifically, in the main repository folder (GenderIdentity_APC), create three folders, called "data - raw", "data - clean", and "plots".

**Step 3: Obtain and prepare BRFSS datasets**

As of this writing, the BRFSS datasets can be found here: https://www.cdc.gov/brfss/annual_data/annual_data.htm. Note that you must individually download each year of data from 2014-2021. Select the "SAS Transport Format" file, and save it to the "data - raw" folder that you created in Step 2. Unzip the files after downloading, if needed.


**Step 4. Prepare datasets and run analyses**

There are two different ways you can run the scripts to obtain the analyses:


## Below here is work in progress; should be complete within a few days (8/25/2025)

***Option A: Run all files***

- *Single Step: Run the Master Script*. Run "GI master script.R" to generate the necessary data sets and all project files, conduct analyses, and generate figures.

***Option B: Run specified files***

- *Prepare Data Sets*. Run the script "SOGI_00_packages.R", then "SOGI_01_generate_BRFSS_index_variable.R", then "SOGI_02_prepare_BRFSS.R". These  files together generate the final clean data files. You only need to run this step once, regardless of how many different analyses you chose to conduct in the nest step.

- *Run Project Analyses*. You may now run the or parts of the remaining scripts based on which specific elements are of interest. A full list can all be found in "GI master script.R." 

