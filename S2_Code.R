# Data Filtering

# Plotting source & acquisition culture dates for clonal & plasmid data
ggplot(df.clonal.aggr) +
  aes(x = as.numeric(Source_Date_of_culture), 
      y = as.numeric(Acquisition_Date_of_culture), 
      colour = Acquisition_Species) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  labs(title = "Source and Acquisition Culture Dates",
       subtitle = "Data: Clonal", 
       x = "Source Date of Culture",
       x = "Acquisition Date of Culture") +
  theme_minimal()

# Before filtering
ggplot(df.plasmid.aggr) +
  aes(x = as.numeric(Source_Date_of_culture), 
      y = as.numeric(Acquisition_Date_of_culture), 
      colour = Acquisition_Species) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  labs(title = "Source and Acquisition Culture Dates",
       subtitle = "Data: Clonal", 
       x = "Source Date of Culture",
       x = "Acquisition Date of Culture") +
  theme_minimal()


# Count number of PTs from original dataset -> 698 PTs
pt_count <- rbind(df.clonal, df.plasmid) %>% 
  select(Acquisition_Patient_ID, Source_Patient_ID)
length(unique(append(pt_count$Acquisition_Patient_ID, pt_count$Source_Patient_ID)))

# Keeping environment tidy
rm(pt_count)


# Total dataset enrollment by hospital
ggplot(df.admission) +
  aes(x = Admission_Hospital) +
  geom_bar(fill = "#112446") +
  labs(title = "Dataset Enrollment by Hospital") + 
  xlab("Hospital") + 
  ylab("Number of Patients") + 
  theme_minimal()

# Dataset enrollment over time by hospital
ggplot(df.admission) +
  aes(x = as.numeric(Admission_Date), fill = as.factor(Admission_Hospital)) +
  geom_histogram(bins = 30L) +
  labs(title = "Hospital Enrollment Over Time") + 
  scale_fill_discrete(name = "Hospital") + 
  xlab("Admission Date") + 
  ylab("Number of Patients") + 
  theme_minimal()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sequence Type (ST) Data:

{ 
  # Keep observations from only the top 12 STs of each dataset:
  # Assign ST data as a factor
  df.clonal.aggr$Acquisition_ST <- as.factor(df.clonal.aggr$Acquisition_ST)
  df.plasmid.aggr$Acquisition_ST <- as.factor(df.plasmid.aggr$Acquisition_ST)
  
  # Count and record list of top 12 STs in each dataset
  st.clonal <- df.clonal.aggr %>% group_by(Acquisition_ST) %>% count(Acquisition_ST) %>% filter(Acquisition_ST!="-") %>% 
    arrange(desc(n)) %>% mutate(inf_per = round(n/nrow(df.clonal.aggr)*100)) %>% head(n=12)
  st.plasmid <- df.plasmid.aggr %>% group_by(Acquisition_ST) %>% count(Acquisition_ST) %>% filter(Acquisition_ST!="-") %>% 
    arrange(desc(n)) %>% mutate(inf_per = round(n/nrow(df.plasmid.aggr)*100)) %>% head(n=12)
  
  st.clonal.list <- as.character(st.clonal$Acquisition_ST)
  st.plasmid.list <- as.character(st.plasmid$Acquisition_ST)
  
  # Remove observations not in top 12 STs
  df.clonal.t12 <- df.clonal.aggr[df.clonal.aggr$Acquisition_ST %in% st.clonal.list,]
  df.plasmid.t12 <- df.plasmid.aggr[df.plasmid.aggr$Acquisition_ST %in% st.plasmid.list,]
  }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






# Find the average, min, and max in time-to-event data for each exposure metric
{#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ### Clonal
  ### Direct Ward Exposure: PT-PT transmission, same time & same location
  Clonal_DirWard <- df.clonal.aggr %>% filter(Direct_Ward_Contact_Hrs > 1)
  # - Total #: 66/1451 events
  nrow(Clonal_DirWard)
  # - Hot Site:
  Clonal_DirWard %>% count(Direct_Ward_Contact_Site, sort = TRUE) %>% head(n=3)
  # - Max Hrs: 576 hrs
  max(Clonal_DirWard$Direct_Ward_Contact_Hrs)/24
  # - Avg Hrs: 120.7 hrs
  mean(Clonal_DirWard$Direct_Ward_Contact_Hrs)/24
  # - Min Hrs: 23.93 hrs
  min(Clonal_DirWard$Direct_Ward_Contact_Hrs)/24
  
  ### Direct Hospital Exposure: PT-PT transmission, same time & same location
  Clonal_DirHosp <- df.clonal.aggr %>% filter(Direct_Hosp_Contact_Hrs > 1)
  # - Total #: 407/1451 events
  nrow(Clonal_DirHosp)
  # - Hot Site: H2 by a large margin
  Clonal_DirHosp %>% count(Direct_Hosp_Contact_Site, sort = TRUE) %>% head(n=3)
  # - Max Hrs: 2570 hrs
  max(Clonal_DirHosp$Direct_Hosp_Contact_Hrs)/24
  # - Avg Hrs: 227 hrs
  mean(Clonal_DirHosp$Direct_Hosp_Contact_Hrs)/24
  # - Min Hrs: 23.93 hrs
  min(Clonal_DirHosp$Direct_Hosp_Contact_Hrs)/24
  
  ### Indirect Ward Exposure: Indirect transmission w/ basis of location, SrcPT preceeding AcqPT 
  Clonal_IndWard <- df.clonal.aggr %>% filter(Indirect_Ward_Contact_Hrs > 1)
  # - Total #: 276/1451 events
  nrow(Clonal_IndWard)
  # - Max Hrs: 2185 hrs
  max(Clonal_IndWard$Indirect_Ward_Contact_Hrs)/24
  # - Avg Hrs: 288.6 hrs
  mean(Clonal_IndWard$Indirect_Ward_Contact_Hrs)/24
  # - Min Hrs: 24.01 hrs
  min(Clonal_IndWard$Indirect_Ward_Contact_Hrs)/24
  
  ### Indirect Hospital Exposure: Indirect transmission w/ basis of location, SrcPT preceeding AcqPT 
  Clonal_IndHosp <- df.clonal.aggr %>% filter(Indirect_Hosp_Contact_Hrs > 1)
  # - Total #: 968/1451 events
  nrow(Clonal_IndHosp)
  # - Max Hrs: 2858 hrs -> 119 days
  max(Clonal_IndHosp$Indirect_Hosp_Contact_Hrs)/24
  # - Avg Hrs: 591.6 hrs -> 24.65 days
  mean(Clonal_IndHosp$Indirect_Hosp_Contact_Hrs)/24
  # - Min Hrs: 23.95 hrs
  min(Clonal_IndHosp$Indirect_Hosp_Contact_Hrs)/24
  
  ### Indirect Unrestricted Exposure: Total AcqPT time in hospitals from ScrCultureTime to AcqCultureTime
  Clonal_IndUnres <- df.clonal.aggr %>% filter(Indirect_Unrestricted_Contact_Hrs>1)
  # - Total #: 1299/1451 events
  nrow(Clonal_IndUnres)
  # - Max Hrs: 6974 hrs -> 290.6 days
  max(Clonal_IndUnres$Indirect_Unrestricted_Contact_Hrs)/24
  # - Avg Hrs: 696.9 hrs -> 29 days
  mean(Clonal_IndUnres$Indirect_Unrestricted_Contact_Hrs)/24
  # - Min Hrs: 23.95 hrs 
  min(Clonal_IndUnres$Indirect_Unrestricted_Contact_Hrs)/24
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
}


{#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ### Plasmid
  ### Direct Ward Exposure: PT-PT transmission, same time & same location
  Plasmid_DirWard <- df.plasmid.aggr %>% filter(Direct_Ward_Contact_Hrs > 1)
  # - Total #: 504/30069 events
  nrow(Plasmid_DirWard)
  # - Hot Site:
  Plasmid_DirWard %>% count(Direct_Ward_Contact_Site, sort = TRUE) %>% head(n=3)
  # - Max Hrs: 1032 hrs
  max(Plasmid_DirWard$Direct_Ward_Contact_Hrs)/24
  # - Avg Hrs: 111.0 hrs
  mean(Plasmid_DirWard$Direct_Ward_Contact_Hrs)/24
  # - Min Hrs: 23.93 hrs
  min(Plasmid_DirWard$Direct_Ward_Contact_Hrs)/24
  
  ### Direct Hospital Exposure: PT-PT transmission, same time & same location
  Plasmid_DirHosp <- df.plasmid.aggr %>% filter(Direct_Hosp_Contact_Hrs > 1)
  # - Total #: 6847/30069 events
  nrow(Plasmid_DirHosp)
  # - Hot Site: H2 by a large margin
  Plasmid_DirHosp %>% count(Direct_Hosp_Contact_Site, sort = TRUE) %>% head(n=3)
  # - Max Hrs: 11358 hrs
  max(Plasmid_DirHosp$Direct_Hosp_Contact_Hrs)/24
  # - Avg Hrs: 318.4 hrs
  mean(Plasmid_DirHosp$Direct_Hosp_Contact_Hrs)/24
  # - Min Hrs: 23.93 hrs
  min(Plasmid_DirHosp$Direct_Hosp_Contact_Hrs)/24
  
  ### Indirect Ward Exposure: Indirect transmission w/ basis of location, SrcPT preceeding AcqPT 
  Plasmid_IndWard <- df.plasmid.aggr %>% filter(Indirect_Ward_Contact_Hrs > 1)
  # - Total #: 5562/30069 events
  nrow(Plasmid_IndWard)
  # - Max Hrs: 6100 hrs
  max(Plasmid_IndWard$Indirect_Ward_Contact_Hrs)/24
  # - Avg Hrs: 317.9 hrs
  mean(Plasmid_IndWard$Indirect_Ward_Contact_Hrs)/24
  # - Min Hrs: 23.95 hrs
  min(Plasmid_IndWard$Indirect_Ward_Contact_Hrs)/24
  
  ### Indirect Hospital Exposure: Indirect transmission w/ basis of location, SrcPT preceeding AcqPT 
  Plasmid_IndHosp <- df.plasmid.aggr %>% filter(Indirect_Hosp_Contact_Hrs > 1)
  # - Total #: 20689/30069 events
  nrow(Plasmid_IndHosp)
  # - Max Hrs: 11503 hrs
  max(Plasmid_IndHosp$Indirect_Hosp_Contact_Hrs)/24
  # - Avg Hrs: 881.8 hrs
  mean(Plasmid_IndHosp$Indirect_Hosp_Contact_Hrs)/24
  # - Min Hrs: 23.95 hrs
  min(Plasmid_IndHosp$Indirect_Hosp_Contact_Hrs)/24
  
  ### Indirect Unrestricted Exposure: Total AcqPT time in hospitals from ScrCultureTime to AcqCultureTime
  Plasmid_IndUnres <- df.plasmid.aggr %>% filter(Indirect_Unrestricted_Contact_Hrs>1)
  # - Total #: 29336/30069 events -> other 2.4% should be removed (no exposure)
  nrow(Plasmid_IndUnres)
  # - Max Hrs: 11503.24 hrs -> 479.3 days
  max(Plasmid_IndUnres$Indirect_Unrestricted_Contact_Hrs)/24
  # - Avg Hrs: 1024.477 hrs -> 42.6 days
  mean(Plasmid_IndUnres$Indirect_Unrestricted_Contact_Hrs)/24
  # - Min Hrs: 24.95 hrs 
  min(Plasmid_IndUnres$Indirect_Unrestricted_Contact_Hrs)/24
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Find true start & end of dataset for either clonal or plasmid data
# True start time indicated by earliest Clonal date of acquisition
# True end time indicated by latest Plasmid date of acquisition
# Clonal
min(df.clonal.aggr$Acquisition_Date_of_culture) # True start
max(df.clonal.aggr$Acquisition_Date_of_culture)
# Plasmid
min(df.plasmid.aggr$Acquisition_Date_of_culture)
max(df.plasmid.aggr$Acquisition_Date_of_culture) # True end

# Count number of days in study
time_start <- as.numeric(min(df.clonal.aggr$Acquisition_Date_of_culture))
time_end <- as.numeric(max(df.plasmid.aggr$Acquisition_Date_of_culture))
time_days <- ceiling((time_end - time_start) * 365.2422)
date_decimal(time_start)
date_decimal(time_end)
time_days
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


