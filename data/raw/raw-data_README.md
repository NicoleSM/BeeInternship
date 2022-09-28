# README for raw data folder
This folder stores raw data (i.e. Data files that have not been processed by any script before). For example the XLSX file where the data has been annotated, and the corresponding CSV that have been exported from it for its posterior analysis

## Raw data folders list and description
Here is a descriptive list of the different raw data folders

**Beekeeping** stores hive checkup XLSX files  
- Currently storing outcome of 
  mite assessment for Ib_10, Ca_05, Ca_04, Ib_30, and Ca_29.

**Derivatisation** refers to the derivatisation experiment where a random bee
sample was run, cleaned, then derivatised 
- GCMS files 

  
**Flying** refers to the flying or not flying experiment
- ODS file for the flying or not flying sample list
- GCMS integration files
  - CSV files for the integration results from Qualitative Navigator
- tmp (temporary)
  - Aligned chromatogram results in CSV format
  - Alignment correction table for samples and peak movements in XLSX format
  - Fixed chromatogram alignment tables in XLSX, CSV format
  - Compound ID tables in XLSX, CSV format
  - Data for alignment as Rdata
  - Uncorrected alignment data as Rdata
  - master tables in CSV format

**NM-recognition** stores data corresponding to the nestmate recognition 
experiments
- GCMS files
- GCMS integration files
  - CSV files for the integration results from Qualitative Navigator
- tmp (temporary)
  - Aligned chromatogram results in CSV format
  - Fixed chromatogram alignment tables in XLSX format
  - Compound ID tables in XLSX format
  - Data for alignment as Rdata
  - Uncorrected alignment data as Rdata

**Plasticity** refers to the plasticity experiment (qPCR)
- data-analysis
  - REF_Best in CSV format
    - Average CT values for the reference genes GAPDH, RPL10, RPL19 and Tub per sample
  - run1 in CSV format
    - CT values for genes of interest and two best reference genes (RPL10 and RPL19) for group 1
  - run1_2 in CSV format
    - Groups 1 and 2 run data combined
  - run2 in CSV format
    - CT values for genes of interest and two best reference genes (RPL10 and RPL19) for group 2
  - SK in CSV format
    - CT values of standards for genes of interest and two best reference genes (RPL10 and RPL19)
- PDF of example qPCR pipetting scheme and group data
- PDF of empty pipetting scheme
- qPCR run list along with group info in XLSX format
- RNA extraction table with concentration and dilution data in XLSX format

**queen-less** stores data related to the queenless hive experiment 
- GCMS files
- GCMS integration files
  - CSV files for the integration results from Qualitative Navigator
- tmp (temporary)
  - Data for alignment as Rdata
- Queenless hive sample list as XLSX, CSV, and ODS files

**Queens** refers to the data of the three extracted queens
- GCMS files 

**Survival** refers to the survival experiment
- GCMS files for Control and samples 1-200 
- Control survival experiment data performed in July 2022 as XLSX and CSV files
- Survival run list keeping track of the seven survival boxes 
- Survival data as RDS file

  