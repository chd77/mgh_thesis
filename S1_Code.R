# Data Manipulation/Prep

# Load required packages
pacman::p_load(tidyr, dplyr, ivs, lubridate, purrr, ggplot2, progress, 
               readxl, reshape2, forcats, zoo, survival, ggfortify,
               contsurvplot, pammtools, survminer, stringr,
               cvms, # For cross-validation in survival analysis
               caret, # For k-fold cross-validation functions
               intsurv, # For concordance index functions
               rms, # For calibration plots
               plotly) # For 3d surface plots

# Read in Clonal, Plasmid, and Admission Data
{ df.clonal <- readxl::read_excel("NatComms_5year_Metadata_20240618.xlsx", 
                                  sheet = "clonal_pairs", col_types = "text")
  df.plasmid <- readxl::read_excel("NatComms_5year_Metadata_20240618.xlsx", 
                                   sheet = "plasmid_pairs", col_types = "text")
  df.admission <- readxl::read_excel("NatComms_5year_Metadata_20240618.xlsx", 
                                     sheet = "list_admission", col_types = "text")
  
  # Remove hospital admissions w/o time information (2349 observations)
  # Patients might not have even arrived at the location in these cases
  # Total rows in Admission Data = 13,504, Total rows w/ >0 time = 11,155
  df.admission <- df.admission %>% 
    mutate(diff_hrs = (as.numeric(Stop_Date) - as.numeric(Start_Date))) %>% 
    filter(diff_hrs != 0)
  
  # Check if code has already been run & saved
  clonal_path <- paste0(getwd(),"/df.clonal.aggr.rds")
  if (file.exists(clonal_path)) {
    df.clonal.aggr <- readRDS(file = clonal_path)
  } else {
    print("The aggregated clonal file does not exist. It will run and save automatically.")
  }
  
  plasmid_path <- paste0(getwd(),"/df.plasmid.aggr.rds")
  if (file.exists(plasmid_path)) {
    df.plasmid.aggr <- readRDS(file = plasmid_path)
  } else {
    print("The aggregated plasmid file does not exist. It will run and save automatically.")
  }
}





################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 1: Find time overlap for each AcqPT & SrcPT transmission event
# https://search.r-project.org/CRAN/refmans/ivs/html/iv-splits.html
{
  # Convert character to date variables
  df.splits <- df.admission %>%
    mutate(Start_Date = as.POSIXct(lubridate::date_decimal(as.numeric(Start_Date))),
           Stop_Date = as.POSIXct(lubridate::date_decimal(as.numeric(Stop_Date))))
  
  # Combine Start_Date & Stop_Date to an iv() Range variable
  # Remove times for Start/Stop that are the same (no time intervals)
  df.splits <- df.splits %>% 
    filter(Start_Date != Stop_Date) %>% 
    mutate(iv = iv(Start_Date, Stop_Date), .keep = "unused")
  
  # Print all the disjoint intervals at which people arrived/departed
  iv_splits(df.splits$iv)
  
  # Determine who was at the hospital at any given time (nested splits)
  df.splits <- mutate(df.splits, splits = iv_identify_splits(iv))
  
  # Unnest() the splits to generate disjoint intervals for each guest
  df.splits <- df.splits %>%
    unnest(splits) %>%
    arrange(iv) %>%
    select(Patient_ID, splits)
  df.splits
  
  # Tabulate who was there at any given time (~2,000 time splits)
  df.splits <- df.splits %>% 
    summarise(n = n(), who = list(Patient_ID), .by = splits)
  df.splits
  
  # Create a time difference variable for each split
  # iv_format() to convert iv date range to list of two strings
  # lapply() to pull 1st & 2nd index of list for Start & End
  # str_sub() to remove extra characters
  # as.POSIXct() to convert string to date
  # force_tz() to ensure timezones are constant
  df.splits <- df.splits %>% 
    mutate(
      Start = as.POSIXct(
        stringr::str_sub(
          lapply(
            iv_format(
              df.splits$splits) %>% strsplit(", "), '[', 1), start = 2)),
      End = as.POSIXct(
        stringr::str_sub(
          lapply(
            iv_format(
              df.splits$splits) %>% strsplit(", "), '[', 2), end = -2))
      ) %>%
    mutate(Start = force_tz(Start, tz = "UTC"),
           End = force_tz(End, tz = "UTC"),
           Time_Diff_Hrs = as.numeric(End - Start))
}


# Plot number of patients in each time-step!
ggplot(df.splits) + 
  geom_point(aes(x = 1:nrow(df.splits), y = n)) + 
  labs(title = "Number of Patients in Each Time Interval") + 
  xlab("Time Interval") + 
  ylab("Admitted Patient Total") + 
  theme_minimal()

# Time interval length over complete study -> trimming necessary
plot(df.splits$Time_Diff_Hrs)

# Example for querying the data by the desired Patient_ID
grep("A110", df.splits)
filter(df.splits, purrr::map_lgl(who, ~ all(c("A501", "A592") %in% .x)))

# Example for calculating time shared between two PTs (w/o location data)
df.splits %>%
  filter(purrr::map_lgl(who, ~ all(c("A501", "A592") %in% .x))) %>%
  select(Time_Diff_Hrs) %>%
  sum()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################





################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 2: Find the most probable location of transmission
# (Ward Level -> Hosp later)

# For-loop overview: 
# Filter df.splits by SrcPT and AcqPT to find overlap times
# Sum total time SrcPT and AcqPT were in contact
# Create a list of wards where SrcPT was during this time
# Create a list of wards where AcqPT was during this time
# List any matches by time, select the one which was before the Acq. Culture Date
# Calculate sum of splits until Acq. Culture Date

# For Direct Ward/Hospital transmission, the overlapping ward(s)/hospital(s) 
#  times are summed and the most recent (likely) location preceding the 
#  acquisition culture date is assumed to be the most likely location of 
#  transmission. 
# The site of transmission is assumed but can be swayed by surveillance levels.  

# For Indirect Ward/Hospital transmission, the similar locations between 
#  SrcPT and AcqPT are compared and the total time spent by the AcqPT in
#  these locations is considered the indirect contact (exposure) time.
# The time counted is restricted to being within SrcPT Date of Culture 
#  and the AcqPT Date of Culture.

if (!file.exists(clonal_path) || !file.exists(plasmid_path)) {
  
  df_name_list_from <- c("df.clonal", "df.plasmid")
  df_name_list_to <- c("df.clonal.aggr", "df.plasmid.aggr")
  
  for (run in 1:length(df_name_list_from)) {
    df <- get(df_name_list_from[run])
    df[ , 'Direct_Ward_Contact_Hrs'] = NA
    df[ , 'Direct_Ward_Contact_Site'] = NA
    df[ , 'Direct_Hosp_Contact_Hrs'] = NA
    df[ , 'Direct_Hosp_Contact_Site'] = NA
    df[ , 'Indirect_Ward_Contact_Hrs'] = NA
    df[ , 'Indirect_Hosp_Contact_Hrs'] = NA
    df[ , 'Indirect_Unrestricted_Contact_Hrs'] = NA
    
    st <- Sys.time()
    eta_avg <- 30
    for (i in 1:nrow(df)) {
      SrcID <- df$Source_Patient_ID[i]
      AcqID <- df$Acquisition_Patient_ID[i]
      SrcCultureTime <- df$Source_Date_of_culture[i]
      AcqCultureTime <- df$Acquisition_Date_of_culture[i]
      
      # Src/Acq PT site list, remove all dates outside of Src/AcqCultureTime
      # Keep Stop_Dates which are more than (after) SrcCultureTime
      # Keep Start_Dates which are less than (before) AcqCultureTime
      SrcSites <- df.admission %>% filter(Patient_ID == SrcID,
                                          Stop_Date > SrcCultureTime,
                                          Start_Date < AcqCultureTime)
      AcqSites <- df.admission %>% filter(Patient_ID == AcqID,
                                          Stop_Date > SrcCultureTime,
                                          Start_Date < AcqCultureTime)
      
      # Time correction for AcqSites (Indirect Exposures)
      # Time correction required if patient discharged before date of culture
      # Check for AcqCultureTime within patient admission data
      time_corr1 <- AcqSites %>%
        filter(Start_Date <= AcqCultureTime,
               Stop_Date > AcqCultureTime) %>%
        nrow()
      # If culture was taken during hospital stay, replace last AcqPT date with
      #  AcqCultureTime to remove extra hospital time after positive test
      if (time_corr1 > 0) {
        m <- which(AcqSites$Start_Date < AcqCultureTime &
                     AcqSites$Stop_Date > AcqCultureTime)
        AcqSites$Stop_Date[m] <- AcqCultureTime
        AcqSites$diff_hrs[m] <- as.numeric(AcqCultureTime) -
          as.numeric(AcqSites$Start_Date[m])
      }
      
      # Calculate direct contact locations to overlapping time-splits
      split_locations <- df.splits %>%
        filter(purrr::map_lgl(who, ~ all(c(SrcID, AcqID) %in% .x)))
      split_locations <- split_locations %>% 
        filter(End > date_decimal(as.numeric(SrcCultureTime)),
               Start < date_decimal(as.numeric(AcqCultureTime)))
      split_locations[ , 'Src_Hosp'] = NA
      split_locations[ , 'Src_Ward'] = NA
      split_locations[ , 'Acq_Hosp'] = NA
      split_locations[ , 'Acq_Ward'] = NA
      
      # Confirm there is direct transmission data, if not, skip to indirect
      if (nrow(split_locations) > 0) {
        for (i_sub in 1:nrow(split_locations)) {
          overlap_start_time <- split_locations[[i_sub,"Start"]]
          
          # Source PT
          SrcID_Site <- df.admission %>%
            filter(Patient_ID == SrcID) %>%
            # Convert to string, cut extra time off seconds, convert back to date
            mutate(
              Start_Date =
                floor_date(
                  as.POSIXct(
                    lubridate::date_decimal(
                      as.numeric(Start_Date))), 'seconds'),
              Stop_Date =
                floor_date(
                  as.POSIXct(
                    lubridate::date_decimal(
                      as.numeric(Stop_Date))), 'seconds')) %>%
            # Find PT location during the overlap start time
            filter(overlap_start_time >= Start_Date &
                     overlap_start_time < Stop_Date) %>%
            select(Admission_Hospital, Ward)
          
          # Acquisition PT
          AcqID_Site <- df.admission %>%
            filter(Patient_ID == AcqID) %>%
            # Convert to string, cut extra time off seconds, convert back to date
            mutate(
              Start_Date =
                floor_date(
                  as.POSIXct(
                    lubridate::date_decimal(
                      as.numeric(Start_Date))), 'seconds'),
              Stop_Date =
                floor_date(
                  as.POSIXct(
                    lubridate::date_decimal(
                      as.numeric(Stop_Date))), 'seconds')) %>%
            # Find PT location during the overlap start time
            filter(overlap_start_time >= Start_Date &
                     overlap_start_time < Stop_Date) %>%
            select(Admission_Hospital, Ward)
          
          split_locations[i_sub,'Src_Hosp'] <- SrcID_Site['Admission_Hospital']
          split_locations[i_sub,'Src_Ward'] <- SrcID_Site['Ward']
          split_locations[i_sub,'Acq_Hosp'] <- AcqID_Site['Admission_Hospital']
          split_locations[i_sub,'Acq_Ward'] <- AcqID_Site['Ward']
        }
        
        ### Time Corrections
        # Time corrections for split_locations (Direct Exposures)
        # Check for need of time correction
        time_corr2 <- split_locations %>%
          filter(Start < date_decimal(as.numeric(AcqCultureTime)),
                 End > date_decimal(as.numeric(AcqCultureTime))) %>%
          nrow()
        # Perform time correction
        if (time_corr2 > 0) {
          n <- which(split_locations$Start < date_decimal(as.numeric(AcqCultureTime)) &
                       split_locations$End > date_decimal(as.numeric(AcqCultureTime)))
          split_locations$End[n] <- date_decimal(as.numeric(AcqCultureTime))
          split_locations$Time_Diff_Hrs[n] <- (
            as.numeric(AcqCultureTime) - 
              decimal_date(split_locations$Start[n]))*365.2422*24
        }
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ### Direct Ward
        # Count time between two PTs in the same ward at the same time
        # 'Direct_Ward_Contact_Hrs', 'Direct_Ward_Contact_Site'
        
        df$Direct_Ward_Contact_Hrs[i] <- split_locations %>%
          filter(Src_Ward == Acq_Ward) %>%
          select(Time_Diff_Hrs) %>%
          sum() %>%
          round(digits = 2)
        
        # Find all shared wards between Src & Acq PTs
        # .x=Src, .y=Acq
        WardShared <- inner_join(SrcSites, AcqSites, by = "Ward",
                                 relationship = "many-to-many")
        # Keep only Wards which overlap during time
        # Start(Acq) before Stop(Src) & Start(Src) before Stop(Acq)
        # Start(Acq) < Stop(Src) & Start(Src) < Stop(Acq)
        # any((Start_Date.y < Stop_Date.x) & (Start_Date.x < Stop_Date.y))
        # any() used to compare across different rows after inner_join()
        WardSpecific <- WardShared %>% group_by(Ward) %>%
          filter(any((Start_Date.y < Stop_Date.x) & (Start_Date.x < Stop_Date.y))) %>%
          arrange(Start_Date.y) %>%
          tail(n=1)
        LikelyWard <- WardSpecific$Ward
        
        # Assign Likely Ward where the transmission event occurred
        if (length(LikelyWard) == 1){
          df$Direct_Ward_Contact_Site[i] <- LikelyWard
        } else if (length(LikelyWard) >1) {
          print("Error - LikelyWard: length(LikelyWard) > 1")
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ### Direct Hospital
        # Add time between two PTs in the same hospital at the same time
        # 'Direct_Hosp_Contact_Hrs', 'Direct_Hosp_Contact_Site'
        
        df$Direct_Hosp_Contact_Hrs[i] <- split_locations %>%
          filter(Src_Hosp == Acq_Hosp) %>%
          select(Time_Diff_Hrs) %>%
          sum() %>%
          round(digits = 2)
        
        # Find all shared wards between Src & Acq PTs
        HospShared <- inner_join(SrcSites, AcqSites, by = "Admission_Hospital",
                                 relationship = "many-to-many")
        
        # Combine and filter overlapping times and find likely transmission hospital
        # Most recent overlapping hospital to the culture date
        HospSpecific <- HospShared %>% group_by(Admission_Hospital) %>%
          filter(any((Start_Date.y < Stop_Date.x) & (Start_Date.x < Stop_Date.y))) %>%
          arrange(Start_Date.y) %>%
          tail(n=1)
        LikelyHosp <- HospSpecific$Admission_Hospital
        
        # Assign Likely Hospital where transmission event occurred
        if (length(LikelyHosp) == 1){
          df$Direct_Hosp_Contact_Site[i] <- LikelyHosp
        } else if (length(LikelyHosp) >1) {
          print("Error - LikelyHosp: length(LikelyHosp) > 1")
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
      } else {}
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      ### Indirect Ward
      # 'Indirect_Ward_Contact_Hrs'
      # Compare Wards, add up AcqPT time spent in that ward until Culture Date
      # Includes both total indirect and direct exposure time for each PlasmidPair
      # SrcPT must have been in the site BEFORE the AcqPT
      
      # Find all shared wards between Src & Acq PTs
      # .x=Src, .y=Acq
      # Keep all time/sites in which the SrcPT was positive for the Bacteria
      # Keep all time/sites in which the AcqPT was negative/untested for the Bacteria
      WardShared <- inner_join(filter(SrcSites, Stop_Date>=SrcCultureTime),
                               filter(AcqSites, Start_Date<=AcqCultureTime),
                               by = "Ward",
                               relationship = "many-to-many")
      ward_list <- unique(WardShared$Ward)
      
      # Ensure exposure time is counted after SrcPT has visited the site(s)
      if (length(ward_list) != 0) {
        contact_time <- 0
        for (i_sub2 in 1:length(ward_list)) {
          src_indirect <- SrcSites %>% filter(Ward == ward_list[i_sub2]) %>% slice(1)
          acq_indirect <- AcqSites %>% filter(Ward == ward_list[i_sub2])
          for (i_sub3 in 1:nrow(acq_indirect)) {
            # SrcPT arrives before AcqPT
            if (src_indirect$Start_Date <= acq_indirect$Start_Date[i_sub3]) {
              contact_time <- contact_time + (acq_indirect$diff_hrs[i_sub3] *365.2422*24)
              # SrcPT arrives after AcqPT (ensure overlap w/ SrcPT Start before AcqPT Stop)
            } else if (src_indirect$Start_Date >= acq_indirect$Start_Date[i_sub3] &
                       src_indirect$Start_Date <= acq_indirect$Stop_Date[i_sub3]) {
              late_src_admission <- as.numeric(acq_indirect$Stop_Date[i_sub3]) -
                as.numeric(src_indirect$Start_Date)
              contact_time <- contact_time + (late_src_admission *365.2422*24)
            }
          }
        }
        df$Indirect_Ward_Contact_Hrs[i] <- contact_time %>% round(digits = 2)
      } else {}
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      ### Indirect Hospital
      # 'Indirect_Hosp_Contact_Hrs'
      # Compare Hospital, add up AcqPT time spent in that Hospital until Culture Date
      # Includes both total indirect and direct exposure time for each PlasmidPair
      
      # Find all shared hospitals between Src & Acq PTs
      # .x=Src, .y=Acq
      # Keep all time/sites in which the SrcPT was positive for the Bacteria
      # Keep all time/sites in which the AcqPT was negative/untested for the Bacteria
      HospShared <- inner_join(filter(SrcSites, Stop_Date>=SrcCultureTime),
                               filter(AcqSites, Start_Date<=AcqCultureTime),
                               by = "Admission_Hospital",
                               relationship = "many-to-many")
      hosp_list <- unique(HospShared$Admission_Hospital)
      
      # Ensure exposure time is counted after SrcPT has visited the site(s)
      if (length(hosp_list) != 0) {
        contact_time <- 0
        for (i_sub2 in 1:length(hosp_list)) {
          src_indirect <- SrcSites %>% filter(Admission_Hospital == hosp_list[i_sub2]) %>% slice(1)
          acq_indirect <- AcqSites %>% filter(Admission_Hospital == hosp_list[i_sub2])
          for (i_sub3 in 1:nrow(acq_indirect)) {
            # SrcPT arrives before AcqPT
            if (src_indirect$Start_Date <= acq_indirect$Start_Date[i_sub3]) {
              contact_time <- contact_time + (acq_indirect$diff_hrs[i_sub3] *365.2422*24)
              # SrcPT arrives after AcqPT (ensure overlap w/ SrcPT Start before AcqPT Stop)
            } else if (src_indirect$Start_Date >= acq_indirect$Start_Date[i_sub3] &
                       src_indirect$Start_Date <= acq_indirect$Stop_Date[i_sub3]) {
              late_src_admission <- as.numeric(acq_indirect$Stop_Date[i_sub3]) -
                as.numeric(src_indirect$Start_Date)
              contact_time <- contact_time + (late_src_admission *365.2422*24)
            }
          }
        }
        df$Indirect_Hosp_Contact_Hrs[i] <- contact_time %>% round(digits = 2)
      } else {}
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      ### Indirect Unrestricted Total
      # 'Indirect_Unrestricted_Contact_Hrs'
      # Count all the AcqPT hospital exposure time between SrcPT Date of Culture
      #  and AcqPT Date of Culture
      
      # Filter through and sum all AcqPT time that meet criteria
      df$Indirect_Unrestricted_Contact_Hrs[i] <- AcqSites %>%
        mutate(diff_hrs = (diff_hrs*365.2422*24)) %>%
        select(diff_hrs) %>%
        sum() %>%
        round(digits = 2)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      # Progress check every 50 iterations
      if (i%%50==0){
        # Time for 50 iterations * number of groups of 50 left
        eta <- as.numeric(difftime(Sys.time(),st,units="mins"))
        eta <- eta*((nrow(df)-i)/50)
        eta_avg <- (eta_avg + eta)/2
        print(paste0(i,'/',nrow(df),' | ETA: ',round(eta_avg),'mins | Finish: ',
                     format(Sys.time()+eta_avg*60, format = "%H:%M"), 
                     ' | Data: ', df_name_list_to[run] ))
        st <- Sys.time()
      }
    }
    
    # Change all 0 contact hours to NA
    df[df == 0] <- NA
    
    # Save for-loop output, one for each clonal and plasmid
    path <- paste0(getwd(), "/", df_name_list_to[run], ".rds")
    saveRDS(df, file = path)
    
    # Assign for-loop output to R environment
    assign(df_name_list_to[run], df)
    
    # Clean environment
    rm(i, i_sub, i_sub2, i_sub3, SrcID, AcqID, SrcCultureTime, AcqCultureTime,
       SrcID_Site, AcqID_Site, SrcSites, AcqSites, src_indirect, acq_indirect,
       overlap_start_time, WardShared, WardSpecific, LikelyWard, ward_list,
       HospShared, HospSpecific, LikelyHosp, hosp_list, split_locations,
       contact_time, time_corr1, time_corr2, late_src_admission, n, m, eta,
       eta_avg, st, df, path)
  }
  
  # Clean environment
  rm(df_name_list_from, df_name_list_to, run, clonal_path, plasmid_path)
  
} else {
  rm(clonal_path, plasmid_path)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################


