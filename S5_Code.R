# Predictive Modelling

# Bootstrap Predictive Models to Assess Model Performance
# Failures of the model are expected if the k-folds do not contain possible species
# Check each trained model against other data
{
  path1 <- paste0(getwd(), "/df.c.index.matrix.corr10x.rds")
  path2 <- paste0(getwd(), "/df.c.index.matrix.all3x.rds")
  to_run <- c()
  
  # Check paths for previously computed and saved data  
  if (file.exists(path1)) {
    c_index_matrix_corr10x <- readRDS(file = path1)
  } else {
    to_run <- append(to_run, "c_index_matrix_corr10x")
    print("c_index_matrix_corr10x will be run and saved.")
  }
  
  if (file.exists(path2)) {
    c_index_matrix_all3x <- readRDS(file = path2)
  } else {
    to_run <- append(to_run, "c_index_matrix_all3x")
    print("c_index_matrix_all3x will be run and saved.")
  }
  
  # Model and exposure metric dataset combinations
  models <- c("cox_mo_cdw","cox_mo_cdh","cox_mo_ciu","cox_mo_pdw","cox_mo_pdh","cox_mo_piu",
              "cox_qt_cdw","cox_qt_cdh","cox_qt_ciu","cox_qt_pdw","cox_qt_pdh","cox_qt_piu",
              "cox_yr_cdw","cox_yr_cdh","cox_yr_ciu","cox_yr_pdw","cox_yr_pdh","cox_yr_piu",
              "cox_sp_cdw","cox_sp_cdh","cox_sp_ciu","cox_sp_pdw","cox_sp_pdh","cox_sp_piu",
              "cox_st_cdw","cox_st_cdh","cox_st_ciu","cox_st_pdw","cox_st_pdh","cox_st_piu",
              "cox_sp_mo_cdw","cox_sp_mo_cdh","cox_sp_mo_ciu","cox_sp_mo_pdw","cox_sp_mo_pdh","cox_sp_mo_piu",
              "cox_sp_qt_cdw","cox_sp_qt_cdh","cox_sp_qt_ciu","cox_sp_qt_pdw","cox_sp_qt_pdh","cox_sp_qt_piu",
              "cox_sp_yr_cdw","cox_sp_yr_cdh","cox_sp_yr_ciu","cox_sp_yr_pdw","cox_sp_yr_pdh","cox_sp_yr_piu",
              "cox_st_mo_cdw","cox_st_mo_cdh","cox_st_mo_ciu","cox_st_mo_pdw","cox_st_mo_pdh","cox_st_mo_piu",
              "cox_st_qt_cdw","cox_st_qt_cdh","cox_st_qt_ciu","cox_st_qt_pdw","cox_st_qt_pdh","cox_st_qt_piu",
              "cox_st_yr_cdw","cox_st_yr_cdh","cox_st_yr_ciu","cox_st_yr_pdw","cox_st_yr_pdh","cox_st_yr_piu",
              "cox_all_mo_cdw","cox_all_mo_cdh","cox_all_mo_ciu","cox_all_mo_pdw","cox_all_mo_pdh","cox_all_mo_piu",
              "cox_all_qt_cdw","cox_all_qt_cdh","cox_all_qt_ciu","cox_all_qt_pdw","cox_all_qt_pdh","cox_all_qt_piu",
              "cox_all_yr_cdw","cox_all_yr_cdh","cox_all_yr_ciu","cox_all_yr_pdw","cox_all_yr_pdh","cox_all_yr_piu")
  datasets <- c("df.cdw", "df.cdh", "df.ciu", "df.pdw", "df.pdh", "df.piu")
  
  suppressWarnings(
    if ("c_index_matrix_corr10x" %in% to_run) {
      # Create the matrix of model and dataset combinations
      c_index_matrix_corr10x <- matrix(NA, nrow = length(models), ncol = length(rep(datasets,10)))
      rownames(c_index_matrix_corr10x) <- models
      colnames(c_index_matrix_corr10x) <- rep(datasets,10)
      
      i <- 1
      iter <- nrow(c_index_matrix_corr10x)*ncol(c_index_matrix_corr10x)/6
      eta_avg <- c(iter/10)
      fold_num <- 10
      st_full <- Sys.time()
      st <- Sys.time()
      
      for (col in 1:ncol(c_index_matrix_corr10x)) {
        for (row in 1:nrow(c_index_matrix_corr10x)) {
          df <- colnames(c_index_matrix_corr10x)[col]
          model <- rownames(c_index_matrix_corr10x)[row]
          
          # Only build and compute models on their corresponding transmission types
          if (str_sub(df, -3, -1) == str_sub(model, -3, -1)) {
            
            # Define the Dependent Variable for later use 
            if (str_sub(df, -2, -1) == "dw") {
              DV <- "Direct_Ward_Contact_Hrs"
            } else if (str_sub(df, -2, -1) == "dh") {
              DV <- "Direct_Hosp_Contact_Hrs"
            } else if (str_sub(df, -2, -1) == "iu") {
              DV <- "Indirect_Unrestricted_Contact_Hrs"
            } else {stop("Error with dependent variable selection")}
            df <- get(df)
            
            # Define the k-fold cross-validation
            cv_folds <- createFolds(as.numeric(unlist(df[DV])), k=fold_num, list=TRUE)
            
            # Continue even if certain models fail k-fold CV process
            try(results <- sapply(cv_folds, function(fold) {
              train_data <- df[-fold, ]
              test_data <- df[fold, ]
              
              # Fit the model on training data
              if (str_sub(model, 5, -5) == "mo") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month, data=train_data)
              } else if (str_sub(model, 5, -5) == "qt") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar, data=train_data)
              } else if (str_sub(model, 5, -5) == "yr") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year, data=train_data)
              } else if (str_sub(model, 5, -5) == "sp") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ Acquisition_Species, data=train_data)
              } else if (str_sub(model, 5, -5) == "st") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ Acquisition_ST, data=train_data)
              } else if (str_sub(model, 5, -5) == "sp_mo") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month + Acquisition_Species, data=train_data)
              } else if (str_sub(model, 5, -5) == "sp_qt") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar + Acquisition_Species, data=train_data)
              } else if (str_sub(model, 5, -5) == "sp_yr") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year + Acquisition_Species, data=train_data)
              } else if (str_sub(model, 5, -5) == "st_mo") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month + Acquisition_ST, data=train_data)
              } else if (str_sub(model, 5, -5) == "st_qt") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar + Acquisition_ST, data=train_data)
              } else if (str_sub(model, 5, -5) == "st_yr") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year + Acquisition_ST, data=train_data)
              } else if (str_sub(model, 5, -5) == "all_mo") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month + Acquisition_Species + Acquisition_ST, data=train_data)
              } else if (str_sub(model, 5, -5) == "all_qt") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar + Acquisition_Species + Acquisition_ST, data=train_data)
              } else if (str_sub(model, 5, -5) == "all_yr") {
                cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year + Acquisition_Species + Acquisition_ST, data=train_data)
              } else {stop("Error with model selection")}
              
              # Predict on test data
              test_pred <- predict(cox_model_cv, newdata = test_data, type = "risk")
              
              # Calculate concordance index (C-index)
              concordance_index <- cIndex(time = as.numeric(unlist(test_data[DV])),
                                          event = test_data$status, risk_score = test_pred)
              return(concordance_index[1])
            }), silent=TRUE)
            
            # Average C-index across folds
            c_index_matrix_corr10x[row,col] <- mean(results)
            
            # Estimated Times
            eta <- as.numeric(difftime(Sys.time(),st,units="mins"))*(iter-(i))
            eta_avg <- append(eta_avg, eta)
            print(paste0((i),'/',iter,' | ETA: ',round(mean(eta_avg)),'mins | Finish: ',
                         format(Sys.time()+mean(eta_avg)*60, format = "%H:%M")))
            st <- Sys.time()
            i <- i + 1
          }
        }
      }
      
      saveRDS(c_index_matrix_corr10x, file = path1)
      
      print(paste0("Completed: ", (i-1), " | Compute Time: ", 
                   round(difftime(Sys.time(), st_full, units="mins")), 
                   " mins | Finished On: ", format(Sys.time(), format = "%H:%M")))
    })
  
  
  suppressWarnings(
    if ("c_index_matrix_all3x" %in% to_run) {
      # Create the matrix of model and dataset combinations
      c_index_matrix_all3x <- matrix(NA, nrow = length(models), ncol = length(rep(datasets,3)))
      rownames(c_index_matrix_all3x) <- models
      colnames(c_index_matrix_all3x) <- rep(datasets,3)
      
      i <- 1
      iter <- nrow(c_index_matrix_all3x)*ncol(c_index_matrix_all3x)
      eta_avg <- c(iter/10)
      fold_num <- 10
      st_full <- Sys.time()
      st <- Sys.time()
      
      for (col in 1:ncol(c_index_matrix_all3x)) {
        for (row in 1:nrow(c_index_matrix_all3x)) {
          df <- colnames(c_index_matrix_all3x)[col]
          model <- rownames(c_index_matrix_all3x)[row]
          
          # Define the Dependent Variable for later use 
          if (str_sub(df, -2, -1) == "dw") {
            DV <- "Direct_Ward_Contact_Hrs"
          } else if (str_sub(df, -2, -1) == "dh") {
            DV <- "Direct_Hosp_Contact_Hrs"
          } else if (str_sub(df, -2, -1) == "iu") {
            DV <- "Indirect_Unrestricted_Contact_Hrs"
          } else {stop("Error with dependent variable selection")}
          df <- get(df)
          
          # Define the k-fold cross-validation
          cv_folds <- createFolds(as.numeric(unlist(df[DV])), k=fold_num, list=TRUE)
          
          # Continue even if certain models fail k-fold CV process
          try(results <- sapply(cv_folds, function(fold) {
            train_data <- df[-fold, ]
            test_data <- df[fold, ]
            
            # Fit the model on training data
            if (str_sub(model, 5, -5) == "mo") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month, data=train_data)
            } else if (str_sub(model, 5, -5) == "qt") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar, data=train_data)
            } else if (str_sub(model, 5, -5) == "yr") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year, data=train_data)
            } else if (str_sub(model, 5, -5) == "sp") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ Acquisition_Species, data=train_data)
            } else if (str_sub(model, 5, -5) == "st") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ Acquisition_ST, data=train_data)
            } else if (str_sub(model, 5, -5) == "sp_mo") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month + Acquisition_Species, data=train_data)
            } else if (str_sub(model, 5, -5) == "sp_qt") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar + Acquisition_Species, data=train_data)
            } else if (str_sub(model, 5, -5) == "sp_yr") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year + Acquisition_Species, data=train_data)
            } else if (str_sub(model, 5, -5) == "st_mo") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month + Acquisition_ST, data=train_data)
            } else if (str_sub(model, 5, -5) == "st_qt") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar + Acquisition_ST, data=train_data)
            } else if (str_sub(model, 5, -5) == "st_yr") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year + Acquisition_ST, data=train_data)
            } else if (str_sub(model, 5, -5) == "all_mo") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_month + Acquisition_Species + Acquisition_ST, data=train_data)
            } else if (str_sub(model, 5, -5) == "all_qt") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_quar + Acquisition_Species + Acquisition_ST, data=train_data)
            } else if (str_sub(model, 5, -5) == "all_yr") {
              cox_model_cv <- coxph(formula=Surv(get(DV), status) ~ time_year + Acquisition_Species + Acquisition_ST, data=train_data)
            } else {stop("Error with model selection")}
            
            # Predict on test data
            test_pred <- predict(cox_model_cv, newdata = test_data, type = "risk")
            
            # Calculate concordance index (C-index)
            concordance_index <- cIndex(time = as.numeric(unlist(test_data[DV])),
                                        event = test_data$status, risk_score = test_pred)
            return(concordance_index[1])
          }), silent=TRUE)
          
          # Average C-index across folds
          c_index_matrix_all3x[row,col] <- mean(results)
          
          # Estimated Times
          eta <- as.numeric(difftime(Sys.time(),st,units="mins"))*(iter-(i))
          eta_avg <- append(eta_avg, eta)
          print(paste0((i),'/',iter,' | ETA: ',round(mean(eta_avg)),'mins | Finish: ',
                       format(Sys.time()+mean(eta_avg)*60, format = "%H:%M")))
          st <- Sys.time()
          i <- i + 1
        }
      }
      
      saveRDS(c_index_matrix_all3x, file = path2)
      
      print(paste0("Completed: ", (i-1), " | Compute Time: ", 
                   round(difftime(Sys.time(), st_full, units="mins")), 
                   " mins | Finished On: ", format(Sys.time(), format = "%H:%M")))
    })
  
  # Keep environment tidy
  rm(path1,path2,to_run,models,datasets,col,DV,eta,eta_avg,fold_num,i,iter,
     results,row,st,st_full)
  
  # Convert to dataframe
  c_index_df_corr10x <- as.data.frame(c_index_matrix_corr10x)
  c_index_df_corr10x$model <- rownames(c_index_df_corr10x)
  rownames(c_index_df_corr10x) <- NULL
  c_index_df_all3x <- as.data.frame(c_index_matrix_all3x)
  c_index_df_all3x$model <- rownames(c_index_df_all3x)
  rownames(c_index_df_all3x) <- NULL
  
  # Pivot longer and Calculate average model performance
  c_index_df_corr10x <- c_index_df_corr10x %>% 
    pivot_longer(names_to = "dataset", values_to = "c_index", cols = starts_with("df")) %>% 
    drop_na(c_index) %>%
    group_by(model) %>% 
    mutate(avg_c_index = mean(c_index)) %>% 
    arrange(desc(avg_c_index))
  c_index_df_all3x <- c_index_df_all3x %>% 
    pivot_longer(names_to = "dataset", values_to = "c_index", cols = starts_with("df")) %>% 
    drop_na(c_index) %>%
    group_by(model) %>% 
    mutate(avg_c_index = mean(c_index)) %>% 
    arrange(desc(avg_c_index))
}

# Plot the C-index (performance) for each model
ggplot(c_index_df_corr10x) + 
  geom_violin(aes(x = fct_reorder(model, desc(avg_c_index)), y = c_index)) +
  geom_point(aes(x = fct_reorder(model, desc(avg_c_index)), y = avg_c_index), color = "red") +
  labs(title = "Model Performance w/ K-Fold Cross-Validation",
       subtitle = "K-fold CV with models built from corresponding Clonal & Plasmid datasets (10x replicates)",
       x = "Model", y = "Average Concordance Index") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90))

ggplot(c_index_df_all3x) + 
  geom_violin(aes(x = fct_reorder(model, desc(avg_c_index)), y = c_index)) +
  geom_point(aes(x = fct_reorder(model, desc(avg_c_index)), y = avg_c_index), color = "red") +
  labs(title = "Model Performance w/ K-Fold Cross-Validation",
       subtitle = "K-fold CV with models built from all Clonal & Plasmid datasets (3x replicates)",
       x = "Model", y = "Average Concordance Index") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90))

# Best model(s): cox_all_mo_* -> time_month + Acquisition_Species + Acquisition_ST
# Models which use time_quarter are a close second, Species and ST appear to be important






################################################################################
## Exposure Metric Performance
# Check to see AIC, C-index, and LRT for different exposure metrics
# Indirect Unrestricted has the best AIC score compared to other dependent variables
extractAIC(coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
extractAIC(coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
extractAIC(coxph(formula=Surv(Indirect_Ward_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
extractAIC(coxph(formula=Surv(Indirect_Hosp_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
extractAIC(coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))

summary(coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
summary(coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
summary(coxph(formula=Surv(Indirect_Ward_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
summary(coxph(formula=Surv(Indirect_Hosp_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
summary(coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
################################################################################






################################################################################
# Best model: cox_all_mo_piu
# Time by Month, Bacterial Species, and Strain Type model using Indirect Unrestricted data as DV!

df <- df.pdh
dataset <- "Indirect Unrestricted Hospital"
model <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.piu)

# Get predicted survival probabilities
predicted_probs <- predict(cox_all_mo_piu, newdata = df.pdh, type = "surv")

# Create survival data with predictions
calibration_data <- data.frame(predicted_probs = predicted_probs, 
                               actual_event = df$Indirect_Unrestricted_Contact_Hrs)

# Define bins for calibration plot
calibration_data_mean <- calibration_data %>%
  mutate(bin = cut(predicted_probs, 
                   breaks = quantile(predicted_probs, probs = 0:10/10, na.rm = TRUE), 
                   include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(predicted_probs), mean_actual = mean(actual_event))

# Plot calibration
ggplot() +
  geom_point(data = calibration_data, aes(x = predicted_probs, y = actual_event)) +
  geom_point(data = calibration_data_mean, aes(x = mean_pred, y = mean_actual), color = 'red') +
  geom_smooth(data = calibration_data_mean, aes(x = mean_pred, y = mean_actual), se = FALSE) +
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
  labs(x = "Predicted Probability", y = "Observed Time-To-Event", 
       title = paste0("Calibration Plot - ", dataset), subtitle = model['formula']) +
  theme_minimal()

ggforest(cox_all_mo_ciu, data=as.data.frame(df.ciu), main="Best Model Clonal Hazard Ratios", fontsize=0.4)
# ggforest(cox_all_mo_piu, data=as.data.frame(df.piu), main="Best Model Plasmid Hazard Ratios", fontsize=0.4)
round(as.data.frame(summary(cox_all_mo_ciu)$conf.int),2)
round(as.data.frame(summary(cox_all_mo_piu)$conf.int),2)
################################################################################













################################################################################
# Sensitivity Analysis - Increasing Time-to-Infection (Lower FOI)
# Calculate differences in survival from baseline with various reductions 
#  in exposure hours and at various times.

# Check for previously saved sensitivity data
if(file.exists(paste0(getwd(), "/df.km.sens.results.rds"))) {
  km_sens <- readRDS(file = paste0(getwd(), "/df.km.sens.results.rds"))
} else {
  
  # Set the range of percentages to used to increase the time-to-infection
  per <- c(seq(from=0.10, to=2.50, by=0.005))
  km_sens <- data.frame(per_incr_tti=double(),day_03=double(),
                        day_04=double(),day_05=double(),day_06=double(),
                        day_07=double(),day_08=double(),day_09=double(),
                        day_10=double(),day_11=double(),day_12=double(),
                        day_13=double(),day_14=double(),day_15=double(),
                        day_16=double(),day_17=double(),day_18=double(),
                        day_19=double(),day_20=double(),day_21=double(),
                        day_22=double(),day_23=double(),day_24=double(),
                        day_25=double(),day_26=double(),day_27=double(),
                        day_28=double(),day_29=double(),day_30=double())
  km_fit_baseline <- survfit(coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df.pdh))
  
  for (i in 1:length(per)) {
    # Record percentages to be calculated
    km_sens[i,'per_incr_tti'] <- (per[i]-1)*100
    
    # Modify exposure metric by selected percentage
    df <- df.pdh %>% mutate(Indirect_Unrestricted_Contact_Hrs = Indirect_Unrestricted_Contact_Hrs*per[i])
    
    # Build the model with modified exposure metrics
    km_fit_incr <- survfit(coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month + Acquisition_Species + Acquisition_ST, df))
    
    # Calculate deviation of modified population survival from baseline from day 3-30
    for (k in 1:(ncol(km_sens)-1)) {
      day <- as.numeric(str_sub(colnames(km_sens)[k+1], -2, -1))
      surv_incr <- approx(km_fit_incr$time, km_fit_incr$surv, xout = day*24)$y
      surv_base <- approx(km_fit_baseline$time, km_fit_baseline$surv, xout = day*24)$y
      km_sens[i,k+1] <- round((surv_incr - surv_base)*100, 3)
    }
    print(paste0('Completed: ',(per[i]*100),'% / ',(max(per)*100),'%'))
  }
  
  # Save the Sensitivity Data
  saveRDS(km_sens, file = paste0(getwd(), "/df.km.sens.results.rds"))
}

km_sens

# Convert for ggplot
km_sens_long <- km_sens %>% 
  pivot_longer(names_to = 'days', values_to = 'values', cols = contains('day'))

# Plot
ggplot(km_sens_long) + 
  geom_smooth(aes(x = per_incr_tti, y = values, color = days)) + 
  labs(title = 'Averted Infections at Daily Timepoints by Modifying Time-to-Infection',
    x = 'Percent Change in Time-to-Infection',
    y = 'Percent Infections Averted') 

# Surface Topology of Varying FOI Over Time
plot_ly(y = km_sens[,1], z = ~as.matrix(km_sens[,2:ncol(km_sens)])) %>% 
  add_surface(
    contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        highlightcolor="#ff0000",
        project=list(z=TRUE)))) %>% 
  layout(title = 'Infection Percentage Averted at Daily Timepoints by Increasing Time-to-Infection', 
         scene = list(xaxis = list(title = 'Day'), 
                      yaxis = list(title = '% Increase in Time-to-Infection'), 
                      zaxis = list(title = '% Averted Infections'))) %>%
  hide_colorbar()





# Singapore Hospitals:
# Singapore General Hospital -> 1,939 beds (avg discharge rate in our dataset = 6.4/pt)
# Changi General Hospital ->    1,043 beds
# Sengkang General Hospital ->  799 beds
# KK Women's and Children's hospital -> 848 beds
# Khoo Teck Puat Hospital ->    795 beds
# Woodlands Health Campus ->    1,800 beds
# Tan Tock Seng Hospital ->     1,700 beds

# Total beds = 8924
total_beds <- sum(1939, 1043, 799, 848, 795, 1800, 1700)
# Assume 100% Occupancy??
total_pt_days_per_year <- total_beds * 365
# 7.73 to 10.32 per 100000 patient-days (Marimuthu)
# 9 = CRE infections per 100,000 PT days
total_cre_per_year <- total_pt_days_per_year * 9/100000
# 5 years in this study
total_cre <- total_cre_per_year * 5

# Graph Effect Sizes with Singapore Populations and 7 Hospitals
km_sens_long_pt <- km_sens_long %>% mutate(values = values*total_cre/100)

ggplot(km_sens_long_pt) + 
  geom_smooth(aes(x = per_incr_tti, y = values, color = days)) + 
  labs(title = 'Averted Infections at Daily Timepoints by Modifying Time-to-Infection',
       subtitle = "Assuming 1466 PTs Over 5 Years",
       x = 'Percent Change in Time-to-Infection',
       y = 'Percent Infections Averted')

################################################################################



