# Install the necessary packages ----
NPacks <- c("tidyverse", "here", "GCalignR")

# Load required packages ----
pacman::p_load(char = NPacks)

# FUNCTIONS ----#

# Load the data ----
load(here("data", "raw", analysis, "tmp", "data2align.Rdata"))

# Align dataframes ####
align_df <- function(df){
  set.seed(12345)
  df <- align_chromatograms(df
                            , rt_col_name = "RT"
                            , max_linear_shift = linear_shift_criteria
                            , max_diff_peak2mean = partial_alignment_threshold
                            , min_diff_peak2peak = row_merging_threshold
                            # , blanks = blanks_list
  )
  df
}

# Check precision ####
## A function to obtain aligned dataframes within the master group list and label each peak
RT_df <- function(df){
  df_RT <- df$aligned$RT
  rownames(df_RT) <- paste0("P", 1:nrow(df_RT))
  df_RT
  
}

area_df <- function(df){
  df_area <- df$aligned$Area
  rownames(df_area) <- paste0("P", 1:nrow(df_area))
  df_area
  
}
# Normalise aligned abundance (area) ####
## A function to normalise area and give a label to each peak
area_norm_df <- function(df){
  df_area_norm <- norm_peaks(data = df
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
 
  colnames(df_area_norm) <- paste0("P", 1:length(colnames(df_area_norm)))
  str(df_area_norm) 
  df_area_norm
} 


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


## mg_list ####
aligned_mg_list <-  lapply(mg_list, align_df)
print(aligned_mg_list)

## standards (STD) ####
set.seed(12345)
aligned_STD_list <- align_chromatograms(standards_list
                                                  , rt_col_name = "RT"
                                                  , max_linear_shift = 0.02
                                                  , max_diff_peak2mean = 0.07
                                                  , min_diff_peak2peak = 0.04)
print(aligned_STD_list)

# Check precision of the alignment ----
## Get aligned data frames ####
# ### Samples #### #

mg_list_RT <-  lapply(aligned_mg_list, RT_df)
print(mg_list_RT)

mg_list_area <-  lapply(aligned_mg_list, area_df)
print(mg_list_area)


### STD ####
### The area of standards is not important, only the retention time will be used
### later for calculating the retention index of the different CHC
STD_RT  <- aligned_STD_list$aligned$RT
rownames(STD_RT) <- paste0("P", 1:nrow(STD_RT))


# Export aligned data frames  ####
## CSV of RT data frames ####
## This is useful to check up the alignments, while looking at the chromatograms
## in the GC-MS analysis software.
DF_names <- names(mg_list_RT) 
for(df in DF_names){
  name_change <- str_replace(df, "_", "-")
  write.csv(mg_list_RT[[df]]
            , here("data"
                   , "raw"
                   , analysis
                   , "tmp"
                   , paste0("aligned_RT_"
                            , name_change
                            ,"_"
                            , gcms_batch
                            ,".csv")))
}

### STD ####
write.csv(STD_RT
          , here("data"
                 , "raw"
                 , "queen-less"
                 , "tmp"
                 , paste0("aligned_RT_STD_"
                          , gcms_batch
                          ,".csv")))


## dataframes as Rdata file ####
save(list = c("mg_list_area", "mg_list_RT")
     , file = here("data",  "raw"
                   , analysis
                   , "tmp"
                   , paste0("uncorrected-alignment_"
                            , gcms_batch
                            , ".Rdata")))
print("The aligned data frames were exported")

# Diagnostic plots ####
## Normalize aligned abundance (area) ####
mg_list_area_norm <-  lapply(aligned_mg_list, area_norm_df)


### STD ####
STD_area_norm <- norm_peaks(data = aligned_STD_list
                                 , rt_col_name = "RT"
                                 , conc_col_name = "Area")
### Give a label to each peak
colnames(STD_area_norm) <- paste0("P", 1:length(colnames(STD_area_norm)))
str(STD_area_norm)


# Export diagnostic plots ####
#### Set color palette
heatmap_colors <- viridis::turbo(200)

#### Generate PDF file to contain plots
pdf(here("output"
         , analysis
         , paste0("uncorrected-alignment-plots_"
                  , gcms_batch
                  , '.pdf'))
    , width = 30, height = 15)

## mg_list ####
for(df in DF_names){
  name_plot <- df
  gc_heatmap(aligned_mg_list[[df]]
             , main_title = paste0("Alignment of ", name_plot, " CHC"))
  gc_heatmap(aligned_mg_list[[df]]
             , type = "discrete"
             , main_title = paste0("Alignment of ", name_plot, " CHC"))
  # Throws an error that stops the script if sourcing
  #plot(aligned_W1_Nu_list), which_plot = "all")
  {gplots::heatmap.2(as.matrix(mg_list_area_norm[[df]] %>% log1p()), 
                     main = paste0("Clustering of ", name_plot, " CHC"),
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
  )}}

## STD ####
gc_heatmap(aligned_STD_list
           , main_title = paste0("Alignment of ", gcms_batch, " standard alkanes"))
gc_heatmap(aligned_STD_list
           , type = "discrete"
           , main_title = paste0("Alignment of ", gcms_batch, " standard alkanes"))
# Throws an error that stops the script if sourcing
# plot(aligned_STD_0322_data_list
#      , which_plot = "all")
{gplots::heatmap.2(as.matrix(STD_area_norm %>% log1p()),
                   main =  paste0("Clustering of ", gcms_batch, " standard alkanes"),
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
               , file = here("output", analysis, "SInf_Script02.txt"))
print("The sessionInfo report was exported. The script 02 finished running")
