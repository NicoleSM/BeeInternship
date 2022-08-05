# Install the necessary packages ----
NPacks <- c("tidyverse", "here", "GCalignR")

# Load required packages ----
pacman::p_load(char = NPacks)

# FUNCTIONS ----#

# Load the data ----
load(here("data", "raw", "NM-recognition", "tmp", "data2align.Rdata"))
  
# Data alignment ----
# ## Samples #### #
# blanks_list <- names(samples_data_list) %>% 
#   str_subset("Hexane")
 
# # Threshold for partial peak alignment step.
# # Skimming the chromatograms no separated alkadienes/alkenes peaks seemed
# # to have a separation below 0.049 min among samples
partial_alignment_threshold <- 0.045
# 
# # Keep the linear shift small. The half of the partial alignment criteria
# # is just an arbitrary value. The idea is to avoid that this shift compromises
# # the partial alignment or rows merging steps
linear_shift_criteria <- partial_alignment_threshold  / 3
# 
# # Threshold for merging rows step.
# # Skimming the chromatograms no compound seemed to have a separation
# # over 0.083 min between samples
row_merging_threshold <- 0.08


## Br_Ca ####
set.seed(12345)
aligned_Br_Ca_list <- 
  align_chromatograms(Br_Ca_data_list
                      , rt_col_name = "RT"
                      , max_linear_shift = linear_shift_criteria
                      , max_diff_peak2mean = partial_alignment_threshold
                      , min_diff_peak2peak = row_merging_threshold
                      # , blanks = blanks_list
                      )
print(aligned_Br_Ca_list)

## Br_Ib ####
set.seed(12345)
aligned_Br_Ib_list <- 
  align_chromatograms(Br_Ib_data_list
                      , rt_col_name = "RT"
                      , max_linear_shift = linear_shift_criteria
                      , max_diff_peak2mean = partial_alignment_threshold
                      , min_diff_peak2peak = row_merging_threshold
                      # , blanks = blanks_list
  )
print(aligned_Br_Ib_list)

## Wu_Ca ####
set.seed(12345)
aligned_Wu_Ca_list <- 
  align_chromatograms(Wu_Ca_data_list
                      , rt_col_name = "RT"
                      , max_linear_shift = linear_shift_criteria
                      , max_diff_peak2mean = partial_alignment_threshold
                      , min_diff_peak2peak = row_merging_threshold
                      # , blanks = blanks_list
  )
print(aligned_Wu_Ca_list)

## Wu_Ib_0322 ####
set.seed(12345)
aligned_Wu_Ib_0322_list <- 
  align_chromatograms(Wu_Ib_0322_data_list
                      , rt_col_name = "RT"
                      , max_linear_shift = linear_shift_criteria
                      , max_diff_peak2mean = partial_alignment_threshold
                      , min_diff_peak2peak = row_merging_threshold
                      # , blanks = blanks_list
  )
print(aligned_Wu_Ib_0322_list)

## Wu_Ib_0422 ####
set.seed(12345)
aligned_Wu_Ib_0422_list <- 
  align_chromatograms(Wu_Ib_0422_data_list
                      , rt_col_name = "RT"
                      , max_linear_shift = linear_shift_criteria
                      , max_diff_peak2mean = partial_alignment_threshold
                      , min_diff_peak2peak = row_merging_threshold
                      # , blanks = blanks_list
  )
print(aligned_Wu_Ib_0422_list)

## STD-0322 ####
set.seed(12345)
aligned_STD_0322_data_list <- align_chromatograms(standards_0322_data_list
                                                 , rt_col_name = "RT"
                                                 , max_linear_shift = 0.02
                                                 , max_diff_peak2mean = 0.07
                                                 , min_diff_peak2peak = 0.04)
print(aligned_STD_0322_data_list)

## STD-0422 ####
set.seed(12345)
aligned_STD_0422_data_list <- align_chromatograms(standards_0422_data_list
                                                  , rt_col_name = "RT"
                                                  , max_linear_shift = 0.02
                                                  , max_diff_peak2mean = 0.07
                                                  , min_diff_peak2peak = 0.04)
print(aligned_STD_0422_data_list)

# Check precision of the alignment ----
## Get aligned data frames ####
# ### Samples #### #
### Br-Ca ####
Br_Ca_area  <- aligned_Br_Ca_list$aligned$Area
Br_Ca_RT  <- aligned_Br_Ca_list$aligned$RT

#### Give a label to each peak.
#### This is useful for distinguishing peaks, or filtering operations.
rownames(Br_Ca_area) <- paste0("P", 1:nrow(Br_Ca_area))
rownames(Br_Ca_RT) <- paste0("P", 1:nrow(Br_Ca_RT))

### Br-Ib ####
Br_Ib_area  <- aligned_Br_Ib_list$aligned$Area
Br_Ib_RT  <- aligned_Br_Ib_list$aligned$RT

#### Give a label to each peak.
#### This is useful for distinguishing peaks, or filtering operations.
rownames(Br_Ib_area) <- paste0("P", 1:nrow(Br_Ib_area))
rownames(Br_Ib_RT) <- paste0("P", 1:nrow(Br_Ib_RT))

### Wu-Ca ####
Wu_Ca_area  <- aligned_Wu_Ca_list$aligned$Area
Wu_Ca_RT  <- aligned_Wu_Ca_list$aligned$RT

#### Give a label to each peak.
#### This is useful for distinguishing peaks, or filtering operations.
rownames(Wu_Ca_area) <- paste0("P", 1:nrow(Wu_Ca_area))
rownames(Wu_Ca_RT) <- paste0("P", 1:nrow(Wu_Ca_RT))

### Wu-Ib-0322 ####
Wu_Ib_0322_area  <- aligned_Wu_Ib_0322_list$aligned$Area
Wu_Ib_0322_RT  <- aligned_Wu_Ib_0322_list$aligned$RT

#### Give a label to each peak.
#### This is useful for distinguishing peaks, or filtering operations.
rownames(Wu_Ib_0322_area) <- paste0("P", 1:nrow(Wu_Ib_0322_area))
rownames(Wu_Ib_0322_RT) <- paste0("P", 1:nrow(Wu_Ib_0322_RT))

### Wu-Ib-0422 ####
Wu_Ib_0422_area  <- aligned_Wu_Ib_0422_list$aligned$Area
Wu_Ib_0422_RT  <- aligned_Wu_Ib_0422_list$aligned$RT

#### Give a label to each peak.
#### This is useful for distinguishing peaks, or filtering operations.
rownames(Wu_Ib_0422_area) <- paste0("P", 1:nrow(Wu_Ib_0422_area))
rownames(Wu_Ib_0422_RT) <- paste0("P", 1:nrow(Wu_Ib_0422_RT))

### STD-0322 ####
### The area of standards is not important, only the retention time will be used
### later for calculating the retention index of the different CHC
STD_0322_RT  <- aligned_STD_0322_data_list$aligned$RT
rownames(STD_0322_RT) <- paste0("P", 1:nrow(STD_0322_RT))

### STD-0422 ####
### The area of standards is not important, only the retention time will be used
### later for calculating the retention index of the different CHC
STD_0422_RT  <- aligned_STD_0422_data_list$aligned$RT
rownames(STD_0422_RT) <- paste0("P", 1:nrow(STD_0422_RT))

# Export aligned data frames  ####
## CSV of RT data frames ####
## This is useful to check up the alignments, while looking at the chromatograms
## in the GC-MS analysis software.
write.csv(Br_Ib_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_br-ib_"
                          , "049-271"
                          ,".csv")))

write.csv(Br_Ca_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_br-ca_"
                          , "049-271"
                          ,".csv")))

write.csv(Wu_Ib_0322_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_wu-ib-0322_"
                          , "049-271"
                          ,".csv")))

write.csv(Wu_Ib_0422_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_wu-ib-0422_"
                          , "049-271"
                          ,".csv")))
write.csv(Wu_Ca_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_wu-ca_"
                          , "049-271"
                          ,".csv")))

write.csv(STD_0322_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_STD_0322_"
                          , "049-271"
                          ,".csv")))

write.csv(STD_0422_RT
          , here("data"
                 , "raw"
                 , "NM-recognition"
                 , "tmp"
                 , paste0("aligned_RT_STD_0422_"
                          , "049-271"
                          ,".csv")))

## dataframes as Rdata file ####
save(list = c("Br_Ca_area", "Br_Ca_RT", "Br_Ib_area", "Br_Ib_RT"
              , "Wu_Ca_area", "Wu_Ca_RT"
              , "Wu_Ib_0322_area", "Wu_Ib_0322_RT"
              , "Wu_Ib_0422_area", "Wu_Ib_0422_RT"
              ,"STD_0322_RT", "STD_0422_RT")
     , file = here("data",  "raw"
                   , "NM-recognition"
                   , "tmp"
                   , paste0("uncorrected-alignment_"
                            , "049-271"
                            , ".Rdata")))
print("The aligned data frames were exported")

# Diagnostic plots ####
## Normalize aligned abundance (area) ####
### Br_Ca ####
Br_Ca_area_norm <- norm_peaks(data = aligned_Br_Ca_list
                             , rt_col_name = "RT"
                             , conc_col_name = "Area")
### Give a label to each peak
colnames(Br_Ca_area_norm) <- paste0("P", 1:length(colnames(Br_Ca_area_norm)))
str(Br_Ca_area_norm)


### Br_Ib ####
Br_Ib_area_norm <- norm_peaks(data = aligned_Br_Ib_list
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
### Give a label to each peak
colnames(Br_Ib_area_norm) <- paste0("P", 1:length(colnames(Br_Ib_area_norm)))
str(Br_Ib_area_norm)

### Wu_Ca ####
Wu_Ca_area_norm <- norm_peaks(data = aligned_Wu_Ca_list
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
### Give a label to each peak
colnames(Wu_Ca_area_norm) <- paste0("P", 1:length(colnames(Wu_Ca_area_norm)))
str(Wu_Ca_area_norm)


### Wu_Ib_0322 ####
Wu_Ib_0322_area_norm <- norm_peaks(data = aligned_Wu_Ib_0322_list
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
### Give a label to each peak
colnames(Wu_Ib_0322_area_norm) <- 
  paste0("P", 1:length(colnames(Wu_Ib_0322_area_norm)))
str(Wu_Ib_0322_area_norm)

### Wu_Ib_0422 ####
Wu_Ib_0422_area_norm <- norm_peaks(data = aligned_Wu_Ib_0422_list
                                   , rt_col_name = "RT"
                                   , conc_col_name = "Area")
### Give a label to each peak
colnames(Wu_Ib_0422_area_norm) <- 
  paste0("P", 1:length(colnames(Wu_Ib_0422_area_norm)))
str(Wu_Ib_0422_area_norm)



### STD_0322 ####
STD_0322_area_norm <- norm_peaks(data = aligned_STD_0322_data_list
                           , rt_col_name = "RT"
                           , conc_col_name = "Area")
### Give a label to each peak
colnames(STD_0322_area_norm) <- paste0("P", 1:length(colnames(STD_0322_area_norm)))
str(STD_0322_area_norm)

### STD_0422 ####
STD_0422_area_norm <- norm_peaks(data = aligned_STD_0422_data_list
                                 , rt_col_name = "RT"
                                 , conc_col_name = "Area")
### Give a label to each peak
colnames(STD_0422_area_norm) <- paste0("P", 1:length(colnames(STD_0422_area_norm)))
str(STD_0422_area_norm)

# Export diagnostic plots ####
#### Set color palette
heatmap_colors <- viridis::turbo(200)

#### Generate PDF file to contain plots
pdf(here("output"
         , "NM-Recognition"
         , paste0("uncorrected-alignment-plots_"
                  , "049-271"
                  , '.pdf'))
    , width = 30, height = 15)

## Br_Ca ####
gc_heatmap(aligned_Br_Ca_list
           , main_title = "Alignment of Br Ca honeybee workers CHC")
gc_heatmap(aligned_Br_Ca_list
           , type = "discrete"
           , main_title = "Alignment of Br Ca honeybee workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_Br_Ca_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(Br_Ca_area_norm %>% log1p()), 
                   main = "Clustering of Br Ca honeybee workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors
                   )}
## Br_Ib ####
gc_heatmap(aligned_Br_Ib_list
           , main_title = "Alignment of Br Ib honeybee workers CHC")
gc_heatmap(aligned_Br_Ib_list
           , type = "discrete"
           , main_title = "Alignment of Br Ib honeybee workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_Br_Ib_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(Br_Ib_area_norm %>% log1p()), 
                   main = "Clustering of Br Ib honeybee workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors
)}

## Wu_Ca ####
gc_heatmap(aligned_Wu_Ca_list
           , main_title = "Alignment of Wu Ca honeybee workers CHC")
gc_heatmap(aligned_Wu_Ca_list
           , type = "discrete"
           , main_title = "Alignment of Wu Ca honeybee workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_Wu_Ca_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(Wu_Ca_area_norm %>% log1p()), 
                   main = "Clustering of Wu Ca honeybee workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors
)}

## Wu_Ib_0322 ####
gc_heatmap(aligned_Wu_Ib_0322_list
           , main_title = "Alignment of Wu Ib (03/22) honeybee workers CHC")
gc_heatmap(aligned_Wu_Ib_0322_list
           , type = "discrete"
           , main_title = "Alignment of Wu Ib (03/22) honeybee workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_Wu_Ib_0322_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(Wu_Ib_0322_area_norm %>% log1p()), 
                   main = "Clustering of Wu Ib (03/22) honeybee workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors
)}

## Wu_Ib_0422 ####
gc_heatmap(aligned_Wu_Ib_0422_list
           , main_title = "Alignment of Wu Ib (04/22) honeybee workers CHC")
gc_heatmap(aligned_Wu_Ib_0422_list
           , type = "discrete"
           , main_title = "Alignment of Wu Ib (04/22) honeybee workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_Wu_Ib_0422_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(Wu_Ib_0422_area_norm %>% log1p()), 
                   main = "Clustering of Wu Ib (04/22) honeybee workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors
)}

## STD_0322 ####
gc_heatmap(aligned_STD_0322_data_list
           , main_title = "Alignment of standard alkanes from 03/22")
gc_heatmap(aligned_STD_0322_data_list
           , type = "discrete"
           , main_title = "Alignment of standard alkanes from 03/22")
# Throws an error that stops the script if sourcing
# plot(aligned_STD_0322_data_list
#      , which_plot = "all")
{gplots::heatmap.2(as.matrix(STD_0322_area_norm %>% log1p()), 
                   main = "Clustering of standard alkanes from 03/22",
                   srtCol = 90,
                   dendrogram = "row",
                   # RSTDv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%)",
                   #denscol = "grey",
                   #density.info = "density",
                   # RSTDSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors)}

## STD_0422 ####
gc_heatmap(aligned_STD_0422_data_list
           , main_title = "Alignment of standard alkanes from 04/22")
gc_heatmap(aligned_STD_0422_data_list
           , type = "discrete"
           , main_title = "Alignment of standard alkanes from 04/22")
# Throws an error that stops the script if sourcing
# plot(aligned_STD_0422_data_list
#      , which_plot = "all")
{gplots::heatmap.2(as.matrix(STD_0422_area_norm %>% log1p()), 
                   main = "Clustering of standard alkanes from 04/22",
                   srtCol = 90,
                   dendrogram = "row",
                   # RSTDv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "log1p of relative abundance (%)",
                   #denscol = "grey",
                   #density.info = "density",
                   # RSTDSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors)}

#### Close graphic device to export plots into the PDF file
dev.off()
print("Diagnostic plots for the alignments were exported")

capture.output(sessionInfo()
               , file = here("output", "NM-recognition", "SInf_Script02.txt"))
print("The sessionInfo report was exported. The script 02 finished running")
