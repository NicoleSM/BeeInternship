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
  write.csv(mg_list_RT[df]
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


write.csv(W3_EIW_RT
          , here("data"
                 , "raw"
                 , analysis
                 , "tmp"
                 , paste0("aligned_RT_W3-EIW_"
                          , "001-051"
                          ,".csv")))

write.csv(W1_Nu_RT
          , here("data"
                 , "raw"
                 , "queen-less"
                 , "tmp"
                 , paste0("aligned_RT_W1-Nu_"
                          , "001-051"
                          ,".csv")))

write.csv(W3_NIW_RT
          , here("data"
                 , "raw"
                 , "queen-less"
                 , "tmp"
                 , paste0("aligned_RT_W3-NIW_"
                          , "001-051"
                          ,".csv")))
# 
# write.csv(STD_0322_RT
#           , here("data"
#                  , "raw"
#                  , "queen-less"
#                  , "tmp"
#                  , paste0("aligned_RT_STD_0322_"
#                           , "001-051"
#                           ,".csv")))


## dataframes as Rdata file ####
save(list = c("W1_Nu_area", "W1_Nu_RT", "W3_EIW_area", "W3_EIW_RT"
              , "W3_NIW_area", "W3_NIW_RT")
     , file = here("data",  "raw"
                   , "queen-less"
                   , "tmp"
                   , paste0("uncorrected-alignment_"
                            , "001-051"
                            , ".Rdata")))
print("The aligned data frames were exported")

# Diagnostic plots ####
## Normalize aligned abundance (area) ####
### W1_Nu ####
W1_Nu_area_norm <- norm_peaks(data = aligned_W1_Nu_list
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
### Give a label to each peak
colnames(W1_Nu_area_norm) <- paste0("P", 1:length(colnames(W1_Nu_area_norm)))
str(W1_Nu_area_norm)


### W3_EIW ####
W3_EIW_area_norm <- norm_peaks(data = aligned_W3_EIW_list
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
### Give a label to each peak
colnames(W3_EIW_area_norm) <- paste0("P", 1:length(colnames(W3_EIW_area_norm)))
str(W3_EIW_area_norm)

### W3_NIW ####
W3_NIW_area_norm <- norm_peaks(data = aligned_W3_NIW_list
                              , rt_col_name = "RT"
                              , conc_col_name = "Area")
### Give a label to each peak
colnames(W3_NIW_area_norm) <- paste0("P", 1:length(colnames(W3_NIW_area_norm)))
str(W3_NIW_area_norm)



# ### STD_0322 ####
# STD_0322_area_norm <- norm_peaks(data = aligned_STD_0322_data_list
#                                  , rt_col_name = "RT"
#                                  , conc_col_name = "Area")
# ### Give a label to each peak
# colnames(STD_0322_area_norm) <- paste0("P", 1:length(colnames(STD_0322_area_norm)))
# str(STD_0322_area_norm)


# Export diagnostic plots ####
#### Set color palette
heatmap_colors <- viridis::turbo(200)

#### Generate PDF file to contain plots
pdf(here("output"
         , "queen-less"
         , paste0("uncorrected-alignment-plots_"
                  , "001-051"
                  , '.pdf'))
    , width = 30, height = 15)

## W1_Nu ####
gc_heatmap(aligned_W1_Nu_list
           , main_title = "Alignment of week 1 nurses CHC")
gc_heatmap(aligned_W1_Nu_list
           , type = "discrete"
           , main_title = "Alignment of week 1 nurses CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_W1_Nu_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(W1_Nu_area_norm %>% log1p()), 
                   main = "Clustering of week 1 nurses CHC",
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
## W3_EIW ####
gc_heatmap(aligned_W3_EIW_list
           , main_title = "Alignment of week 3 elongated abdomen in-hive workers CHC")
gc_heatmap(aligned_W3_EIW_list
           , type = "discrete"
           , main_title = "Alignment of week 3 elongated abdomen in-hive workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_W3_EIW_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(W3_EIW_area_norm %>% log1p()), 
                   main = "Clustering of week 3 elongated abdomen in-hive workers CHC",
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

## W3_NIW ####
gc_heatmap(aligned_W3_NIW_list
           , main_title = "Alignment of week 3 normal in-hive workers CHC")
gc_heatmap(aligned_W3_NIW_list
           , type = "discrete"
           , main_title = "Alignment of week 3 normal in-hive workers CHC")
# Throws an error that stops the script if sourcing
#plot(aligned_W3_NIW_list), which_plot = "all")
{gplots::heatmap.2(as.matrix(W3_NIW_area_norm %>% log1p()), 
                   main = "Clustering of week 3 normal in-hive workers CHC",
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


# ## STD_0322 ####
# gc_heatmap(aligned_STD_0322_data_list
#            , main_title = "Alignment of standard alkanes from 03/22")
# gc_heatmap(aligned_STD_0322_data_list
#            , type = "discrete"
#            , main_title = "Alignment of standard alkanes from 03/22")
# # Throws an error that stops the script if sourcing
# # plot(aligned_STD_0322_data_list
# #      , which_plot = "all")
# {gplots::heatmap.2(as.matrix(STD_0322_area_norm %>% log1p()), 
#                    main = "Clustering of standard alkanes from 03/22",
#                    srtCol = 90,
#                    dendrogram = "row",
#                    # RSTDv = dend,
#                    Colv = "NA", # this to make sure the columns are not ordered
#                    trace = "none",          
#                    # margins =c(5,0.1),      
#                    key.xlab = "log1p of relative abundance (%)",
#                    #denscol = "grey",
#                    #density.info = "density",
#                    # RSTDSideColors = microlab, # to add nice colored strips        
#                    col = heatmap_colors)}


#### Close graphic device to export plots into the PDF file
dev.off()
print("Diagnostic plots for the alignments were exported")

capture.output(sessionInfo()
               , file = here("output", "queen-less", "SInf_Script02.txt"))
print("The sessionInfo report was exported. The script 02 finished running")