# README for scripts folder
This folder stores the scripts used for the data analysis

## Script folders
Here is a descriptive list of the different script folders

**Flying** stores scripts for gcms-processing and statistical tests

**NM-recognition** stores scripts for gcms-processing and statistical tests

**Plasticity** stores scripts for processing the qPCR data along with a PDF 
describing the EasyqPCR package which requires R version 3.6 

**Queenless** stores scripts for gcms-processing

**Survival** stores scripts for processing survival data


## Scripts
Here is a descriptive list of the different scripts

GCMS automised processing scripts:
- 01_importing_data_automised
  - imports GCMS integration CSV files of samples and standards and separates the files into groups
    depending on the desired variables accompanied in the sample list and exports the data as 
    a nested list as Rdata
      - requires input variables of:
          - analysis - name of experiment folder
          - gcms_batch - gcms batch written as a range of run dates and year (i.e. 1708-2408_2022)
          - sample_list - name of file housing list of samples with grouping info
          - factors (desired for comparison) - i.e. Subspecies and Flying for creation of groups

- 02_align_MS_data_automised
  - Uses GCalignR package to align chromatogram data within groups, using Rdata generated in 
    01_importing_data_automised and exports CSV of RT data frames, dataframes as Rdata
    and generated alignment heatmaps to aid misalignment correction
      - requires input variable row_merging_threshold which is the threshold for merging rows step 