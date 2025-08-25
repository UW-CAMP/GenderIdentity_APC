
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


## Below here is work in progress; should be complete within a few days (8/22/2025)

***Option A: Run all files***

- *Single Step: Run the Master Script*. Run "SoSP master script.R" to generate the necessary data sets and all project files, conduct analyses, and generate figures.

***Option B: Run specified files***

- *Prepare Data Sets*. Run the scripts “SoSP_00_packages.R” and then “SoSP_01_prepare.R”. These two files together generate the final clean data files. You only need to run this step once, regardless of how many different analyses you chose to conduct in the nest step.

- *Run Project Analyses*. You may now run the all or parts of the remaining scripts based on which specific elements  are of interest. Statistical analyses are in the scripts “SoSP_02_analyses.R”, while "SoSP_03_tables_figures.R" generates the figures and the data for the tables used in the manuscript. 

---
  
*You have now completed the steps to replicate all analyses used in this manuscript. For any questions, please email Steven Goodreau at goodreau@uw.edu.*



### Sexual orientation in the United States, 2014-2021: age, period and cohort effects

### Replication instructions

**Full citation**: Pending.

**Software**: To recreate this analysis, you will need access to SAS and RStudio. Original work was conducted in the software and operating system environments:

[Who put these in here? Are they up to date?]
- SAS version 9.4 in Windows: Edition	Windows Server 2022 Datacenter Version 21H2 Installed on ‎10/‎30/‎2022 OS build 20348.1487
- RStudio Version 2022.12.0+353 (2022.12.0+353), R Version 4.1.0 (2021-05-18) in Mac OS 13.1

**Option 1: Run all files**

***Single Step: Run the Master Script***

Run "SOGI master script.R" to generate the necessary data sets and all project files. This will generate a series of plots with sexual orientation (BRFSS & YRBS), gender identity (BRFSS only), and sex-of-sex-partners (YRBS only).

**Option 2: Run specified files**

***Step A: Prepare Data Sets***

In order, run the scripts whose index number begins with "0". Before running script 01, note that at line 21, there is a space after “.XPT” before the closing quotation mark. The original file was built on Mac OS 13.1 and the space is required for RStudio to properly open the file. You may need to remove this space if you are working in a Windows or other non-Mac OS.

You have generated the clean data files to conduct the analyses of your choosing.

***Step B: Run Project Analyses***

After completing Step A once, you will save local copies of the data sets used by the remaining project files in the project. You may now run the remaining scripts based on which specific analyses are of interest.
