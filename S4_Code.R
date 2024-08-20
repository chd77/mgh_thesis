# Cox Regression Survival Analysis 
# - Estimates a Hazard Ratio of a given endpoint associated with a risk factor
# - Exposure -> for every 1hr increase in hospital exposure, x% increase in hazard
# - Cox assumption: effects are constant over time and additive

# Cox Proportional Hazard (PH) Assumptions: 
# - The hazard ratio associated with each covariate is constant over time
# - The risk of an event for one individual does not affect the risk of another

# Split into time-based cohorts of Month, Quarter, and Year
# Remove observations before 2017
{
  df.cdw <- Clonal_DirWard %>% 
    mutate(status = 1,
           Acquisition_Date_of_culture = as.numeric(Acquisition_Date_of_culture),
           Direct_Ward_Contact_Site = as.factor(Direct_Ward_Contact_Site),
           Direct_Hosp_Contact_Site = as.factor(Direct_Hosp_Contact_Site),
           Acquisition_Species = as.factor(Acquisition_Species),
           time_year = as.factor(case_when(
             floor(Acquisition_Date_of_culture) <= 2017 ~ 2017,
             floor(Acquisition_Date_of_culture) == 2018 ~ 2018,
             floor(Acquisition_Date_of_culture) == 2019 ~ 2019,
             floor(Acquisition_Date_of_culture) >= 2020 ~ 2020)),
           time_month = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 0 ~ paste0(time_year,'_M01'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 1 ~ paste0(time_year,'_M02'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 2 ~ paste0(time_year,'_M03'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 3 ~ paste0(time_year,'_M04'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 4 ~ paste0(time_year,'_M05'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 5 ~ paste0(time_year,'_M06'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 6 ~ paste0(time_year,'_M07'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 7 ~ paste0(time_year,'_M08'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 8 ~ paste0(time_year,'_M09'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 9 ~ paste0(time_year,'_M10'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 10 ~ paste0(time_year,'_M11'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 11 ~ paste0(time_year,'_M12'))),
           time_quar = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 0 ~ paste0(time_year, '_Q1'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 1 ~ paste0(time_year, '_Q2'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 2 ~ paste0(time_year, '_Q3'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 3 ~ paste0(time_year, '_Q4'))),
           ) %>% 
    filter(as.numeric(Acquisition_Date_of_culture) >= 2017)
  
  df.cdh <- Clonal_DirHosp %>%
    mutate(status = 1,
           Acquisition_Date_of_culture = as.numeric(Acquisition_Date_of_culture),
           Direct_Ward_Contact_Site = as.factor(Direct_Ward_Contact_Site),
           Direct_Hosp_Contact_Site = as.factor(Direct_Hosp_Contact_Site),
           Acquisition_Species = as.factor(Acquisition_Species),
           time_year = as.factor(case_when(
             floor(Acquisition_Date_of_culture) <= 2017 ~ 2017,
             floor(Acquisition_Date_of_culture) == 2018 ~ 2018,
             floor(Acquisition_Date_of_culture) == 2019 ~ 2019,
             floor(Acquisition_Date_of_culture) >= 2020 ~ 2020)),
           time_month = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 0 ~ paste0(time_year,'_M01'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 1 ~ paste0(time_year,'_M02'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 2 ~ paste0(time_year,'_M03'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 3 ~ paste0(time_year,'_M04'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 4 ~ paste0(time_year,'_M05'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 5 ~ paste0(time_year,'_M06'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 6 ~ paste0(time_year,'_M07'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 7 ~ paste0(time_year,'_M08'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 8 ~ paste0(time_year,'_M09'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 9 ~ paste0(time_year,'_M10'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 10 ~ paste0(time_year,'_M11'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 11 ~ paste0(time_year,'_M12'))),
           time_quar = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 0 ~ paste0(time_year, '_Q1'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 1 ~ paste0(time_year, '_Q2'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 2 ~ paste0(time_year, '_Q3'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 3 ~ paste0(time_year, '_Q4')))
           ) %>% 
    filter(as.numeric(Acquisition_Date_of_culture) >= 2017)
  
  df.ciu <- Clonal_IndUnres %>% 
    mutate(status = 1,
           Acquisition_Date_of_culture = as.numeric(Acquisition_Date_of_culture),
           Direct_Ward_Contact_Site = as.factor(Direct_Ward_Contact_Site),
           Direct_Hosp_Contact_Site = as.factor(Direct_Hosp_Contact_Site),
           Acquisition_Species = as.factor(Acquisition_Species),
           time_year = as.factor(case_when(
             floor(Acquisition_Date_of_culture) <= 2017 ~ 2017,
             floor(Acquisition_Date_of_culture) == 2018 ~ 2018,
             floor(Acquisition_Date_of_culture) == 2019 ~ 2019,
             floor(Acquisition_Date_of_culture) >= 2020 ~ 2020)),
           time_month = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 0 ~ paste0(time_year,'_M01'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 1 ~ paste0(time_year,'_M02'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 2 ~ paste0(time_year,'_M03'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 3 ~ paste0(time_year,'_M04'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 4 ~ paste0(time_year,'_M05'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 5 ~ paste0(time_year,'_M06'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 6 ~ paste0(time_year,'_M07'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 7 ~ paste0(time_year,'_M08'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 8 ~ paste0(time_year,'_M09'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 9 ~ paste0(time_year,'_M10'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 10 ~ paste0(time_year,'_M11'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 11 ~ paste0(time_year,'_M12'))),
           time_quar = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 0 ~ paste0(time_year, '_Q1'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 1 ~ paste0(time_year, '_Q2'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 2 ~ paste0(time_year, '_Q3'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 3 ~ paste0(time_year, '_Q4')))
           ) %>% 
    filter(as.numeric(Acquisition_Date_of_culture) >= 2017)
  
  df.pdw <- Plasmid_DirWard %>% 
    mutate(status = 1,
           Acquisition_Date_of_culture = as.numeric(Acquisition_Date_of_culture),
           Direct_Ward_Contact_Site = as.factor(Direct_Ward_Contact_Site),
           Direct_Hosp_Contact_Site = as.factor(Direct_Hosp_Contact_Site),
           Acquisition_Species = as.factor(Acquisition_Species),
           time_year = as.factor(case_when(
             floor(Acquisition_Date_of_culture) <= 2017 ~ 2017,
             floor(Acquisition_Date_of_culture) == 2018 ~ 2018,
             floor(Acquisition_Date_of_culture) == 2019 ~ 2019,
             floor(Acquisition_Date_of_culture) >= 2020 ~ 2020)),
           time_month = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 0 ~ paste0(time_year,'_M01'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 1 ~ paste0(time_year,'_M02'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 2 ~ paste0(time_year,'_M03'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 3 ~ paste0(time_year,'_M04'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 4 ~ paste0(time_year,'_M05'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 5 ~ paste0(time_year,'_M06'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 6 ~ paste0(time_year,'_M07'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 7 ~ paste0(time_year,'_M08'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 8 ~ paste0(time_year,'_M09'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 9 ~ paste0(time_year,'_M10'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 10 ~ paste0(time_year,'_M11'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 11 ~ paste0(time_year,'_M12'))),
           time_quar = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 0 ~ paste0(time_year, '_Q1'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 1 ~ paste0(time_year, '_Q2'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 2 ~ paste0(time_year, '_Q3'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 3 ~ paste0(time_year, '_Q4')))
           ) %>% 
    filter(as.numeric(Acquisition_Date_of_culture) >= 2017)
  
  df.pdh <- Plasmid_DirHosp %>% 
    mutate(status = 1,
           Acquisition_Date_of_culture = as.numeric(Acquisition_Date_of_culture),
           Direct_Ward_Contact_Site = as.factor(Direct_Ward_Contact_Site),
           Direct_Hosp_Contact_Site = as.factor(Direct_Hosp_Contact_Site),
           Acquisition_Species = as.factor(Acquisition_Species),
           time_year = as.factor(case_when(
             floor(Acquisition_Date_of_culture) <= 2017 ~ 2017,
             floor(Acquisition_Date_of_culture) == 2018 ~ 2018,
             floor(Acquisition_Date_of_culture) == 2019 ~ 2019,
             floor(Acquisition_Date_of_culture) >= 2020 ~ 2020)),
           time_month = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 0 ~ paste0(time_year,'_M01'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 1 ~ paste0(time_year,'_M02'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 2 ~ paste0(time_year,'_M03'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 3 ~ paste0(time_year,'_M04'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 4 ~ paste0(time_year,'_M05'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 5 ~ paste0(time_year,'_M06'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 6 ~ paste0(time_year,'_M07'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 7 ~ paste0(time_year,'_M08'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 8 ~ paste0(time_year,'_M09'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 9 ~ paste0(time_year,'_M10'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 10 ~ paste0(time_year,'_M11'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 11 ~ paste0(time_year,'_M12'))),
           time_quar = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 0 ~ paste0(time_year, '_Q1'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 1 ~ paste0(time_year, '_Q2'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 2 ~ paste0(time_year, '_Q3'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 3 ~ paste0(time_year, '_Q4')))
           ) %>% 
    filter(as.numeric(Acquisition_Date_of_culture) >= 2017)
  
  df.piu <- Plasmid_IndUnres %>% 
    mutate(status = 1,
           Acquisition_Date_of_culture = as.numeric(Acquisition_Date_of_culture),
           Direct_Ward_Contact_Site = as.factor(Direct_Ward_Contact_Site),
           Direct_Hosp_Contact_Site = as.factor(Direct_Hosp_Contact_Site),
           Acquisition_Species = as.factor(Acquisition_Species),
           time_year = as.factor(case_when(
             floor(Acquisition_Date_of_culture) <= 2017 ~ 2017,
             floor(Acquisition_Date_of_culture) == 2018 ~ 2018,
             floor(Acquisition_Date_of_culture) == 2019 ~ 2019,
             floor(Acquisition_Date_of_culture) >= 2020 ~ 2020)),
           time_month = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 0 ~ paste0(time_year,'_M01'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 1 ~ paste0(time_year,'_M02'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 2 ~ paste0(time_year,'_M03'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 3 ~ paste0(time_year,'_M04'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 4 ~ paste0(time_year,'_M05'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 5 ~ paste0(time_year,'_M06'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 6 ~ paste0(time_year,'_M07'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 7 ~ paste0(time_year,'_M08'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 8 ~ paste0(time_year,'_M09'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 9 ~ paste0(time_year,'_M10'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 10 ~ paste0(time_year,'_M11'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/12))) == 11 ~ paste0(time_year,'_M12'))),
           time_quar = as.factor(case_when(
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 0 ~ paste0(time_year, '_Q1'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 1 ~ paste0(time_year, '_Q2'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 2 ~ paste0(time_year, '_Q3'),
             floor((as.numeric(str_sub(Acquisition_Date_of_culture, 5, -1)) / (1/4))) == 3 ~ paste0(time_year, '_Q4')))
           ) %>% 
    filter(as.numeric(Acquisition_Date_of_culture) >= 2017)
}







################################################################################
# Predictor Variable 1.1: Time Cohort (by Month)
################################################################################
# Fit the Cox PH models for all transmission types:
cox_mo_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_month, data=df.cdw) # lacking data from 2020
cox_mo_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_month, data=df.cdh)
cox_mo_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month, data=df.ciu)
cox_mo_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_month, data=df.pdw)
cox_mo_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_month, data=df.pdh)
cox_mo_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_month, data=df.piu)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
# cox.zph(cox_mo_cdw) # -> Error
# cox.zph(cox_mo_cdh) # -> Error
cox.zph(cox_mo_ciu)
cox.zph(cox_mo_pdw)
# cox.zph(cox_mo_pdh) # -> Error
# cox.zph(cox_mo_piu) # -> Error

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look good(ish) -> Pass!
# plot(cox.zph(cox_mo_cdw)) # -> Error
# plot(cox.zph(cox_mo_cdh)) # -> Error
plot(cox.zph(cox_mo_ciu))
plot(cox.zph(cox_mo_pdw))
# plot(cox.zph(cox_mo_pdh)) # -> Error
# plot(cox.zph(cox_mo_piu)) # -> Error

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_mo_cdw, type = "martingale"))
plot(residuals(cox_mo_cdh, type = "martingale"))
plot(residuals(cox_mo_ciu, type = "martingale"))
plot(residuals(cox_mo_pdw, type = "martingale"))
plot(residuals(cox_mo_pdh, type = "martingale"))
plot(residuals(cox_mo_piu, type = "martingale"))

# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset

# Testing Linearity of Continuous Covariates Assumption:
# Statistical Tests:
# Compute & plot martingale residuals against time_month
plot(df.cdw$time_month, residuals(cox_mo_cdw, type = "martingale"))
plot(df.cdh$time_month, residuals(cox_mo_cdh, type = "martingale"))
plot(df.ciu$time_month, residuals(cox_mo_ciu, type = "martingale"))
# plot(df.pdw$time_month, residuals(cox_mo_pdw, type = "martingale")) # -> Error
# plot(df.pdh$time_month, residuals(cox_mo_pdh, type = "martingale")) # -> Error
# plot(df.piu$time_month, residuals(cox_mo_piu, type = "martingale")) # -> Error

# Visualize Output
round(as.data.frame(summary(cox_mo_cdw)$conf.int),2)
round(as.data.frame(summary(cox_mo_cdh)$conf.int),2)
round(as.data.frame(summary(cox_mo_ciu)$conf.int),2)
round(as.data.frame(summary(cox_mo_pdw)$conf.int),2)
round(as.data.frame(summary(cox_mo_pdh)$conf.int),2)
round(as.data.frame(summary(cox_mo_piu)$conf.int),2)

summary(cox_mo_cdw)
summary(cox_mo_cdh)
summary(cox_mo_ciu)
summary(cox_mo_pdw)
summary(cox_mo_pdh)
summary(cox_mo_piu)

ggforest(cox_mo_cdw, data = df.cdw)
ggforest(cox_mo_cdh, data = df.cdh)
ggforest(cox_mo_ciu, data = df.ciu)
ggforest(cox_mo_pdw, data = df.pdw)
ggforest(cox_mo_pdh, data = df.pdh)
ggforest(cox_mo_piu, data = df.piu)








################################################################################
# Predictor Variable 1.2: Time Cohort (by Quarter)
################################################################################
# Fit the Cox PH models for all transmission types:
cox_qt_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_quar, data=df.cdw) # lacking data from 2020
cox_qt_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_quar, data=df.cdh)
cox_qt_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_quar, data=df.ciu)
cox_qt_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_quar, data=df.pdw)
cox_qt_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_quar, data=df.pdh)
cox_qt_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_quar, data=df.piu)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
# cox.zph(cox_qt_cdw) # -> Error
cox.zph(cox_qt_cdh)
cox.zph(cox_qt_ciu)
cox.zph(cox_qt_pdw)
cox.zph(cox_qt_pdh)
cox.zph(cox_qt_piu)

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look good(ish) -> Pass!
# plot(cox.zph(cox_qt_cdw)) # -> Error
plot(cox.zph(cox_qt_cdh))
plot(cox.zph(cox_qt_ciu))
plot(cox.zph(cox_qt_pdw))
plot(cox.zph(cox_qt_pdh))
plot(cox.zph(cox_qt_piu))

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_qt_cdw, type = "martingale"))
plot(residuals(cox_qt_cdh, type = "martingale"))
plot(residuals(cox_qt_ciu, type = "martingale"))
plot(residuals(cox_qt_pdw, type = "martingale"))
plot(residuals(cox_qt_pdh, type = "martingale"))
plot(residuals(cox_qt_piu, type = "martingale"))

# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset

# Testing Linearity of Continuous Covariates Assumption:
# Statistical Tests:
# Compute & plot martingale residuals against time_quar
plot(df.cdw$time_quar, residuals(cox_qt_cdw, type = "martingale"))
plot(df.cdh$time_quar, residuals(cox_qt_cdh, type = "martingale"))
plot(df.ciu$time_quar, residuals(cox_qt_ciu, type = "martingale"))
# plot(df.pdw$time_quar, residuals(cox_qt_pdw, type = "martingale")) # -> Error
# plot(df.pdh$time_quar, residuals(cox_qt_pdh, type = "martingale")) # -> Error
# plot(df.piu$time_quar, residuals(cox_qt_piu, type = "martingale")) # -> Error


# Visualize Output
round(as.data.frame(summary(cox_qt_cdw)$conf.int),2)
round(as.data.frame(summary(cox_qt_cdh)$conf.int),2)
round(as.data.frame(summary(cox_qt_ciu)$conf.int),2)
round(as.data.frame(summary(cox_qt_pdw)$conf.int),2)
round(as.data.frame(summary(cox_qt_pdh)$conf.int),2)
round(as.data.frame(summary(cox_qt_piu)$conf.int),2)

summary(cox_qt_cdw)
summary(cox_qt_cdh)
summary(cox_qt_ciu)
summary(cox_qt_pdw)
summary(cox_qt_pdh)
summary(cox_qt_piu)

ggforest(cox_qt_cdw, data=df.cdw, main="Clonal Direct Ward Hazard Ratios")
ggforest(cox_qt_cdh, data=df.cdh, main="Clonal Direct Hospital Hazard Ratios")
ggforest(cox_qt_ciu, data=df.ciu, main="Clonal Indirect Unrestricted Hazard Ratios")
ggforest(cox_qt_pdw, data=df.pdw, main="Plasmid Direct Ward Hazard Ratios")
ggforest(cox_qt_pdh, data=df.pdh, main="Plasmid Direct Hospital Hazard Ratios")
ggforest(cox_qt_piu, data=df.piu, main="Plasmid Indirect Unrestricted Hazard Ratios")









################################################################################
# Predictor Variable 1.3: Time Cohort (by Year)
################################################################################
# Fit the Cox PH models for all transmission types:
cox_yr_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_year, data=df.cdw) # lacking data from 2020
cox_yr_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_year, data=df.cdh)
cox_yr_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_year, data=df.ciu)
cox_yr_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ time_year, data=df.pdw)
cox_yr_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ time_year, data=df.pdh)
cox_yr_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ time_year, data=df.piu)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
cox.zph(cox_yr_cdw)
cox.zph(cox_yr_cdh)
cox.zph(cox_yr_ciu)
cox.zph(cox_yr_pdw)
cox.zph(cox_yr_pdh)
cox.zph(cox_yr_piu)

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look good(ish) -> Pass!
plot(cox.zph(cox_yr_cdw))
plot(cox.zph(cox_yr_cdh))
plot(cox.zph(cox_yr_ciu))
plot(cox.zph(cox_yr_pdw))
plot(cox.zph(cox_yr_pdh))
plot(cox.zph(cox_yr_piu))

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_yr_cdw, type = "martingale"))
plot(residuals(cox_yr_cdh, type = "martingale"))
plot(residuals(cox_yr_ciu, type = "martingale"))
plot(residuals(cox_yr_pdw, type = "martingale"))
plot(residuals(cox_yr_pdh, type = "martingale"))
plot(residuals(cox_yr_piu, type = "martingale"))


# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset


# Testing Linearity of Continuous Covariates Assumption:
# Statistical Tests:
# Compute & plot martingale residuals against time_year
plot(df.cdw$time_year, residuals(cox_yr_cdw, type = "martingale"))
plot(df.cdh$time_year, residuals(cox_yr_cdh, type = "martingale"))
plot(df.ciu$time_year, residuals(cox_yr_ciu, type = "martingale"))
plot(df.pdw$time_year, residuals(cox_yr_pdw, type = "martingale"))
plot(df.pdh$time_year, residuals(cox_yr_pdh, type = "martingale"))
plot(df.piu$time_year, residuals(cox_yr_piu, type = "martingale"))

# Visualize Output
round(as.data.frame(summary(cox_yr_cdw)$conf.int),2)
round(as.data.frame(summary(cox_yr_cdh)$conf.int),2)
round(as.data.frame(summary(cox_yr_ciu)$conf.int),2)
round(as.data.frame(summary(cox_yr_pdw)$conf.int),2)
round(as.data.frame(summary(cox_yr_pdh)$conf.int),2)
round(as.data.frame(summary(cox_yr_piu)$conf.int),2)

summary(cox_yr_cdw)
summary(cox_yr_cdh)
summary(cox_yr_ciu)
summary(cox_yr_pdw)
summary(cox_yr_pdh)
summary(cox_yr_piu)

# Forest plot visualizations
ggforest(cox_yr_cdw, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios")
ggforest(cox_yr_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios")
ggforest(cox_yr_ciu, data=as.data.frame(df.ciu), main="Clonal Indirect Unrestricted Hazard Ratios")
ggforest(cox_yr_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios")
ggforest(cox_yr_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios")
ggforest(cox_yr_piu, data=as.data.frame(df.piu), main="Plasmid Indirect Unrestricted Hazard Ratios")









################################################################################
# Predictor Variable 2: Bacterial Species
################################################################################
summary(as.factor(df.cdw$Acquisition_Species))
summary(as.factor(df.cdh$Acquisition_Species))
summary(as.factor(df.ciu$Acquisition_Species))
summary(as.factor(df.pdw$Acquisition_Species))
summary(as.factor(df.pdh$Acquisition_Species))
summary(as.factor(df.piu$Acquisition_Species))

# Fit the Cox PH models for all transmission types:
cox_sp_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ Acquisition_Species, df.cdw)
cox_sp_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ Acquisition_Species, df.cdh)
cox_sp_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ Acquisition_Species, df.ciu)
cox_sp_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ Acquisition_Species, df.pdw)
cox_sp_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ Acquisition_Species, df.pdh)
cox_sp_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ Acquisition_Species, df.piu)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
cox.zph(cox_sp_cdw)
cox.zph(cox_sp_cdh)
cox.zph(cox_sp_ciu)
cox.zph(cox_sp_pdw)
cox.zph(cox_sp_pdh)
cox.zph(cox_sp_piu)

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look ok(ish) -> Pass!
plot(cox.zph(cox_sp_cdw))
plot(cox.zph(cox_sp_cdh))
plot(cox.zph(cox_sp_ciu))
plot(cox.zph(cox_sp_pdw))
plot(cox.zph(cox_sp_pdh))
plot(cox.zph(cox_sp_piu))

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_sp_cdw, type = "martingale"))
plot(residuals(cox_sp_cdh, type = "martingale"))
plot(residuals(cox_sp_ciu, type = "martingale"))
plot(residuals(cox_sp_pdw, type = "martingale"))
plot(residuals(cox_sp_pdh, type = "martingale"))
plot(residuals(cox_sp_piu, type = "martingale"))


# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset

# Testing Linearity of Continuous Covariates Assumption:

# Statistical Tests:
# Compute & plot martingale residuals against Acquisition_Species
plot(df.cdw$Acquisition_Species, residuals(cox_sp_cdw, type = "martingale"))
plot(df.cdh$Acquisition_Species, residuals(cox_sp_cdh, type = "martingale"))
plot(df.ciu$Acquisition_Species, residuals(cox_sp_ciu, type = "martingale"))
plot(df.pdw$Acquisition_Species, residuals(cox_sp_pdw, type = "martingale"))
plot(df.pdh$Acquisition_Species, residuals(cox_sp_pdh, type = "martingale"))
plot(df.piu$Acquisition_Species, residuals(cox_sp_piu, type = "martingale"))

# Visualize Output
cox_sp_cdw
cox_sp_cdh
cox_sp_ciu
cox_sp_pdw 
cox_sp_pdh 
cox_sp_piu

round(as.data.frame(summary(cox_sp_cdw)$conf.int), 2)
round(as.data.frame(summary(cox_sp_cdh)$conf.int), 2)
round(as.data.frame(summary(cox_sp_ciu)$conf.int), 2)
round(as.data.frame(summary(cox_sp_pdw)$conf.int), 2)
round(as.data.frame(summary(cox_sp_pdh)$conf.int), 2)
round(as.data.frame(summary(cox_sp_piu)$conf.int), 2)

# Forest plot visualizations
ggforest(cox_sp_cdw, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios")
ggforest(cox_sp_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios")
ggforest(cox_sp_ciu, data=as.data.frame(df.ciu), main="Clonal Indirect Unrestricted Hazard Ratios")
ggforest(cox_sp_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios")
ggforest(cox_sp_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios")
ggforest(cox_sp_piu, data=as.data.frame(df.piu), main="Plasmid Indirect Unrestricted Hazard Ratios")









################################################################################
# Predictor Variable 3: Strain Type (for associated bacterial species)
################################################################################
# Nested variable within Bacterial Species
summary(as.factor(df.cdw$Acquisition_ST))
summary(as.factor(df.cdh$Acquisition_ST))
summary(as.factor(df.ciu$Acquisition_ST))
summary(as.factor(df.pdw$Acquisition_ST))
summary(as.factor(df.pdh$Acquisition_ST))
summary(as.factor(df.piu$Acquisition_ST))

# Fit the Cox PH models for all transmission types:
cox_st_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ Acquisition_ST, df.cdw)
cox_st_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ Acquisition_ST, df.cdh)
cox_st_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ Acquisition_ST, df.ciu)
cox_st_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ Acquisition_ST, df.pdw)
cox_st_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ Acquisition_ST, df.pdh)
cox_st_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ Acquisition_ST, df.piu)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
# cox.zph(cox_st_cdw) # -> Error
# cox.zph(cox_st_cdh) # -> Error
# cox.zph(cox_st_ciu) # -> Error
# cox.zph(cox_st_pdw) # -> Error
# cox.zph(cox_st_pdh) # -> Error
# cox.zph(cox_st_piu) # -> Error

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look ok(ish) -> Pass!
# plot(cox.zph(cox_st_cdw)) # -> Error
# plot(cox.zph(cox_st_cdh)) # -> Error
# plot(cox.zph(cox_st_ciu)) # -> Error
# plot(cox.zph(cox_st_pdw)) # -> Error
# plot(cox.zph(cox_st_pdh)) # -> Error
# plot(cox.zph(cox_st_piu)) # -> Error

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_st_cdw, type = "martingale"))
plot(residuals(cox_st_cdh, type = "martingale"))
plot(residuals(cox_st_ciu, type = "martingale"))
plot(residuals(cox_st_pdw, type = "martingale"))
plot(residuals(cox_st_pdh, type = "martingale"))
plot(residuals(cox_st_piu, type = "martingale"))

# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset

# Testing Linearity of Continuous Covariates Assumption:
# Statistical Tests:
# Compute & plot martingale residuals against Acquisition_ST
plot(df.cdw$Acquisition_ST, residuals(cox_st_cdw, type = "martingale"))
plot(df.cdh$Acquisition_ST, residuals(cox_st_cdh, type = "martingale"))
plot(df.ciu$Acquisition_ST, residuals(cox_st_ciu, type = "martingale"))
plot(df.pdw$Acquisition_ST, residuals(cox_st_pdw, type = "martingale"))
plot(df.pdh$Acquisition_ST, residuals(cox_st_pdh, type = "martingale"))
plot(df.piu$Acquisition_ST, residuals(cox_st_piu, type = "martingale"))

# Visualize Output
cox_st_cdw
cox_st_cdh
cox_st_ciu
cox_st_pdw 
cox_st_pdh 
cox_st_piu

round(as.data.frame(summary(cox_st_cdw)$conf.int), 2)
round(as.data.frame(summary(cox_st_cdh)$conf.int), 2)
round(as.data.frame(summary(cox_st_ciu)$conf.int), 2)
round(as.data.frame(summary(cox_st_pdw)$conf.int), 2)
round(as.data.frame(summary(cox_st_pdh)$conf.int), 2)
round(as.data.frame(summary(cox_st_piu)$conf.int), 2)

# Forest plot visualizations
# If ST has 1 obs it gets converted to a reference -> SOMETIMES 
ggforest(cox_st_cdw, data=as.data.frame(df.cdh), main="Clonal Direct Ward Hazard Ratios")
ggforest(cox_st_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios")
ggforest(cox_st_ciu, data=as.data.frame(df.ciu), main="Clonal Indirect Unrestricted Hazard Ratios", fontsize=0.3)
# ggforest(cox_st_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios") # -> Error
ggforest(cox_st_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios")
ggforest(cox_st_piu, data=as.data.frame(df.piu), main="Plasmid Indirect Unrestricted Hazard Ratios", fontsize=0.3)












################################################################################
# Predictor Variable 4.1: Ward
################################################################################
summary(as.factor(df.cdw$Direct_Ward_Contact_Site))
summary(as.factor(df.pdw$Direct_Ward_Contact_Site))

# Fit the Cox PH models for all transmission types:
cox_wa_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ Direct_Ward_Contact_Site, df.cdw)
cox_wa_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ Direct_Ward_Contact_Site, df.pdw)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
# cox.zph(cox_wa_cdw) # -> Error
cox.zph(cox_wa_pdw)

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look bad(ish) -> Fail!
# plot(cox.zph(cox_wa_cdw)) # -> Error
plot(cox.zph(cox_wa_pdw))

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_wa_cdw, type = "martingale"))
plot(residuals(cox_wa_pdw, type = "martingale"))

# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset

# Testing Linearity of Continuous Covariates Assumption:
# Statistical Tests:
# Compute & plot martingale residuals against Hospital Site
plot(df.cdw$Direct_Ward_Contact_Site, residuals(cox_wa_cdw, type = "martingale"))
plot(df.pdw$Direct_Ward_Contact_Site, residuals(cox_wa_pdw, type = "martingale"))

# Visualize Output
cox_wa_cdw
cox_wa_pdw 

round(as.data.frame(summary(cox_wa_cdw)$conf.int), 2)
round(as.data.frame(summary(cox_wa_pdw)$conf.int), 2)

ggforest(cox_wa_cdw, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios", fontsize=0.6)
ggforest(cox_wa_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios", fontsize=0.5)








################################################################################
# Predictor Variable 4.2: Hospital
################################################################################
summary(as.factor(df.cdh$Direct_Hosp_Contact_Site))
summary(as.factor(df.pdh$Direct_Hosp_Contact_Site))

# Fit the Cox PH models for all transmission types:
cox_ho_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ Direct_Hosp_Contact_Site, df.cdh)
cox_ho_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ Direct_Hosp_Contact_Site, df.pdh)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
cox.zph(cox_ho_cdh)
cox.zph(cox_ho_pdh)

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look ok(ish) -> Pass!
plot(cox.zph(cox_ho_cdh))
plot(cox.zph(cox_ho_pdh))

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_ho_cdh, type = "martingale"))
plot(residuals(cox_ho_pdh, type = "martingale"))


# Most of these look okay, I worry about increasing the bins of time from yearly 
#  to monthly because of a lack of viable data
# Maybe I could increase to monthly time-varying steps for the Plasmid dataset


# Testing Linearity of Continuous Covariates Assumption:
# Statistical Tests:
# Compute & plot martingale residuals against Hospital Site
plot(df.cdh$Direct_Hosp_Contact_Site, residuals(cox_ho_cdh, type = "martingale"))
plot(df.pdh$Direct_Hosp_Contact_Site, residuals(cox_ho_pdh, type = "martingale"))

# Visualize Output
cox_ho_cdh
cox_ho_pdh 

round(as.data.frame(summary(cox_ho_cdh)$conf.int), 2)
round(as.data.frame(summary(cox_ho_pdh)$conf.int), 2)

ggforest(cox_ho_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios")
ggforest(cox_ho_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios")










################################################################################
# ALL Predictor Variables: Time Cohort (Month), Bacterial Species, and Sequence Types
################################################################################
# Fit the Cox PH models for all transmission types (Month):
cox_all_mo_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                          time_month + Acquisition_Species + Acquisition_ST, df.cdw)
cox_all_mo_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                          time_month + Acquisition_Species + Acquisition_ST, df.cdh)
cox_all_mo_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ 
                          time_month + Acquisition_Species + Acquisition_ST, df.ciu)
cox_all_mo_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                          time_month + Acquisition_Species + Acquisition_ST, df.pdw)
cox_all_mo_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                          time_month + Acquisition_Species + Acquisition_ST, df.pdh)
cox_all_mo_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ 
                          time_month + Acquisition_Species + Acquisition_ST, df.piu)

# Fit the Cox PH models for all transmission types (Season):
cox_all_qt_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                          time_quar + Acquisition_Species + Acquisition_ST, df.cdw)
cox_all_qt_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                          time_quar + Acquisition_Species + Acquisition_ST, df.cdh)
cox_all_qt_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ 
                          time_quar + Acquisition_Species + Acquisition_ST, df.ciu)
cox_all_qt_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                          time_quar + Acquisition_Species + Acquisition_ST, df.pdw)
cox_all_qt_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                          time_quar + Acquisition_Species + Acquisition_ST, df.pdh)
cox_all_qt_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ 
                          time_quar + Acquisition_Species + Acquisition_ST, df.piu)

# Fit the Cox PH models for all transmission types (Year):
cox_all_yr_cdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                          time_year + Acquisition_Species + Acquisition_ST, df.cdw)
cox_all_yr_cdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                          time_year + Acquisition_Species + Acquisition_ST, df.cdh)
cox_all_yr_ciu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ 
                          time_year + Acquisition_Species + Acquisition_ST, df.ciu)
cox_all_yr_pdw <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                          time_year + Acquisition_Species + Acquisition_ST, df.pdw)
cox_all_yr_pdh <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                          time_year + Acquisition_Species + Acquisition_ST, df.pdh)
cox_all_yr_piu <- coxph(formula=Surv(Indirect_Unrestricted_Contact_Hrs, status) ~ 
                          time_year + Acquisition_Species + Acquisition_ST, df.piu)

# Testing Proportional Hazards Assumption:
# This method assesses the correlation between the scaled Schoenfeld residuals 
#   and time (or other covariates) to check if they are independent over time
# Global Test (Grambsch-Therneau Test) (Formal Statistical Tests)
# cox.zph(cox_all_mo_cdw) # -> Error
# cox.zph(cox_all_mo_cdh) # -> Error
cox.zph(cox_all_mo_ciu)
cox.zph(cox_all_mo_pdw)
# cox.zph(cox_all_mo_pdh) # -> Error
# cox.zph(cox_all_mo_piu) # -> Error

# cox.zph(cox_all_qt_cdw) # -> Error
# cox.zph(cox_all_qt_cdh) # -> Error
cox.zph(cox_all_qt_ciu)
cox.zph(cox_all_qt_pdw)
# cox.zph(cox_all_qt_pdh) # -> Error
cox.zph(cox_all_qt_piu)

# cox.zph(cox_all_yr_cdw) # -> Error
cox.zph(cox_all_yr_cdh)
cox.zph(cox_all_yr_ciu)
cox.zph(cox_all_yr_pdw)
# cox.zph(cox_all_yr_pdh) # -> Error
cox.zph(cox_all_yr_piu)

# Plot Schoenfeld residuals against time
# Lines (covariate) should be approximately horizontal
# These all look ok(ish) -> Pass!
# plot(cox.zph(cox_all_mo_cdw)) # -> Error
# plot(cox.zph(cox_all_mo_cdh)) # -> Error
plot(cox.zph(cox_all_mo_ciu))
plot(cox.zph(cox_all_mo_pdw))
# plot(cox.zph(cox_all_mo_pdh)) # -> Error
# plot(cox.zph(cox_all_mo_piu)) # -> Error

# plot(cox.zph(cox_all_qt_cdw)) # -> Error
# plot(cox.zph(cox_all_qt_cdh)) # -> Error
plot(cox.zph(cox_all_qt_ciu))
plot(cox.zph(cox_all_qt_pdw))
# plot(cox.zph(cox_all_qt_pdh)) # -> Error
plot(cox.zph(cox_all_qt_piu))

# plot(cox.zph(cox_all_yr_cdw)) # -> Error
plot(cox.zph(cox_all_yr_cdh))
plot(cox.zph(cox_all_yr_ciu))
plot(cox.zph(cox_all_yr_pdw))
# plot(cox.zph(cox_all_yr_pdh)) # -> Error
plot(cox.zph(cox_all_yr_piu))

# Compute & plot martingale residuals
# These all look good -> Pass! (no systematic patters)
plot(residuals(cox_all_mo_cdw, type = "martingale"))
plot(residuals(cox_all_mo_cdh, type = "martingale"))
plot(residuals(cox_all_mo_ciu, type = "martingale"))
plot(residuals(cox_all_mo_pdw, type = "martingale"))
plot(residuals(cox_all_mo_pdh, type = "martingale"))
plot(residuals(cox_all_mo_piu, type = "martingale"))

plot(residuals(cox_all_qt_cdw, type = "martingale"))
plot(residuals(cox_all_qt_cdh, type = "martingale"))
plot(residuals(cox_all_qt_ciu, type = "martingale"))
plot(residuals(cox_all_qt_pdw, type = "martingale"))
plot(residuals(cox_all_qt_pdh, type = "martingale"))
plot(residuals(cox_all_qt_piu, type = "martingale"))

plot(residuals(cox_all_yr_cdw, type = "martingale"))
plot(residuals(cox_all_yr_cdh, type = "martingale"))
plot(residuals(cox_all_yr_ciu, type = "martingale"))
plot(residuals(cox_all_yr_pdw, type = "martingale"))
plot(residuals(cox_all_yr_pdh, type = "martingale"))
plot(residuals(cox_all_yr_piu, type = "martingale"))


# Forest plot visualizations
# ggforest(cox_all_mo_cdw, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios") # -> Error
ggforest(cox_all_mo_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios", fontsize=0.55)
ggforest(cox_all_mo_ciu, data=as.data.frame(df.ciu), main="Clonal Indirect Unrestricted Hazard Ratios", fontsize=0.55)
# ggforest(cox_all_mo_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios") # -> Error
ggforest(cox_all_mo_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios", fontsize=0.3)
# ggforest(cox_all_mo_piu, data=as.data.frame(df.piu), main="Plasmid Indirect Unrestricted Hazard Ratios") # -> Error

# ggforest(cox_all_qt_cdw, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios") # -> Error
ggforest(cox_all_qt_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios", fontsize=0.55)
ggforest(cox_all_qt_ciu, data=as.data.frame(df.ciu), main="Clonal Indirect Unrestricted Hazard Ratios", fontsize=0.55)
# ggforest(cox_all_qt_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios") # -> Error
ggforest(cox_all_qt_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios", fontsize=0.3)
ggforest(cox_all_qt_piu, data=as.data.frame(df.piu), main="Plasmid Indirect Unrestricted Hazard Ratios", fontsize=0.3)

ggforest(cox_all_yr_cdw, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios")
ggforest(cox_all_yr_cdh, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios")
ggforest(cox_all_yr_ciu, data=as.data.frame(df.ciu), main="Clonal Indirect Unrestricted Hazard Ratios")
# ggforest(cox_all_yr_pdw, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios") # -> Error
ggforest(cox_all_yr_pdh, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios")
ggforest(cox_all_yr_piu, data=as.data.frame(df.piu), main="Plasmid Indirect Unrestricted Hazard Ratios")









################################################################################
# Assessing Overall Model Fit
# Concordance index: 
#   0.5 is random guesses from the model, 1 is perfect predictive ability
# Likelihood ratio test: 
#   Ratios >1 show association with disease
#   Ratios <1 show association with lack of disease
# AIC scores: 
#   Weights predictive ability with number of variables. Lowest scores indicate a better model
################################################################################

# Comparing time cohort, sequence type, species type, and combined models
summary(cox_mo_cdw)
summary(cox_qt_cdw)
summary(cox_yr_cdw)
summary(cox_sp_cdw)
summary(cox_st_cdw)
summary(cox_wa_cdw)
summary(cox_all_mo_cdw)
summary(cox_all_qt_cdw)
summary(cox_all_yr_cdw)
extractAIC(cox_mo_cdw)
extractAIC(cox_qt_cdw) # Best performing model
extractAIC(cox_yr_cdw)
extractAIC(cox_sp_cdw)
extractAIC(cox_st_cdw)
extractAIC(cox_wa_cdw)
extractAIC(cox_all_mo_cdw)
extractAIC(cox_all_qt_cdw)
extractAIC(cox_all_yr_cdw)

summary(cox_mo_pdw)
summary(cox_qt_pdw)
summary(cox_yr_pdw)
summary(cox_sp_pdw)
summary(cox_st_pdw)
summary(cox_wa_pdw)
summary(cox_all_mo_pdw)
summary(cox_all_qt_pdw)
summary(cox_all_yr_pdw)
extractAIC(cox_mo_pdw) # Best performing model
extractAIC(cox_qt_pdw)
extractAIC(cox_yr_pdw)
extractAIC(cox_sp_pdw)
extractAIC(cox_st_pdw)
extractAIC(cox_wa_pdw)
extractAIC(cox_all_mo_pdw)
extractAIC(cox_all_qt_pdw)
extractAIC(cox_all_yr_pdw)


summary(cox_mo_cdh)
summary(cox_qt_cdh)
summary(cox_yr_cdh)
summary(cox_sp_cdh)
summary(cox_st_cdh)
summary(cox_ho_cdh)
summary(cox_all_mo_cdh)
summary(cox_all_qt_cdh)
summary(cox_all_yr_cdh)
extractAIC(cox_mo_cdh)
extractAIC(cox_qt_cdh)
extractAIC(cox_yr_cdh)
extractAIC(cox_sp_cdh)
extractAIC(cox_st_cdh)
extractAIC(cox_ho_cdh)
extractAIC(cox_all_mo_cdh)
extractAIC(cox_all_qt_cdh) # Best performing model
extractAIC(cox_all_yr_cdh) 

summary(cox_mo_pdh)
summary(cox_qt_pdh)
summary(cox_yr_pdh)
summary(cox_sp_pdh)
summary(cox_st_pdh)
summary(cox_ho_pdh)
summary(cox_all_mo_pdh)
summary(cox_all_qt_pdh)
summary(cox_all_yr_pdh)
extractAIC(cox_mo_pdh)
extractAIC(cox_qt_pdh)
extractAIC(cox_yr_pdh)
extractAIC(cox_sp_pdh)
extractAIC(cox_st_pdh)
extractAIC(cox_ho_pdh)
extractAIC(cox_all_mo_pdh) # Best performing model
extractAIC(cox_all_qt_pdh)
extractAIC(cox_all_yr_pdh)



summary(cox_mo_ciu)
summary(cox_qt_ciu)
summary(cox_yr_ciu)
summary(cox_sp_ciu)
summary(cox_st_ciu)
summary(cox_all_mo_ciu)
summary(cox_all_qt_ciu)
summary(cox_all_yr_ciu)
extractAIC(cox_mo_ciu)
extractAIC(cox_yr_ciu)
extractAIC(cox_sp_ciu)
extractAIC(cox_st_ciu)
extractAIC(cox_all_mo_ciu) # Best performing model
extractAIC(cox_all_qt_ciu)
extractAIC(cox_all_yr_ciu)

summary(cox_mo_piu)
summary(cox_qt_piu)
summary(cox_yr_piu)
summary(cox_sp_piu)
summary(cox_st_piu)
summary(cox_all_mo_piu)
summary(cox_all_qt_piu)
summary(cox_all_yr_piu)
extractAIC(cox_mo_piu)
extractAIC(cox_qt_piu)
extractAIC(cox_yr_piu)
extractAIC(cox_sp_piu)
extractAIC(cox_st_piu)
extractAIC(cox_all_mo_piu) # Best performing model
extractAIC(cox_all_qt_piu)
extractAIC(cox_all_yr_piu)


# AIC score -> best outcome from Monthly all-inclusive models
################################################################################







# Just For Fun! -> Include Location (to Month Models)
cox_all_mo_cdw_loc <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                              time_month + Acquisition_Species + Direct_Ward_Contact_Site, df.cdw)
cox_all_mo_cdh_loc <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                              time_month + Acquisition_Species + Direct_Hosp_Contact_Site, df.cdh)

cox_all_mo_pdw_loc <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                              time_month + Acquisition_Species + Direct_Ward_Contact_Site, df.pdw)
cox_all_mo_pdh_loc <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                              time_month + Acquisition_Species + Direct_Hosp_Contact_Site, df.pdh)

# Include Location (to Year Models)
cox_all_yr_cdw_loc <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                              time_year + Acquisition_Species + Direct_Ward_Contact_Site, df.cdw)
cox_all_yr_cdh_loc <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                              time_year + Acquisition_Species + Direct_Hosp_Contact_Site, df.cdh)

cox_all_yr_pdw_loc <- coxph(formula=Surv(Direct_Ward_Contact_Hrs, status) ~ 
                              time_year + Acquisition_Species + Direct_Ward_Contact_Site, df.pdw)
cox_all_yr_pdh_loc <- coxph(formula=Surv(Direct_Hosp_Contact_Hrs, status) ~ 
                              time_year + Acquisition_Species + Direct_Hosp_Contact_Site, df.pdh)

summary(cox_all_mo_cdw_loc)
summary(cox_all_mo_cdh_loc)
summary(cox_all_mo_pdw_loc)
summary(cox_all_mo_pdh_loc)

summary(cox_all_yr_cdw_loc)
summary(cox_all_yr_cdh_loc)
summary(cox_all_yr_pdw_loc)
summary(cox_all_yr_pdh_loc)

# Clonal Direct Ward model fit w/ and w/o location data
extractAIC(cox_all_mo_cdw)
extractAIC(cox_all_yr_cdw)
extractAIC(cox_all_mo_cdw_loc)
extractAIC(cox_all_yr_cdw_loc) # Best performing model

# Clonal Direct Hospital model fit w/ and w/o location data
extractAIC(cox_all_mo_cdh) # Best performing model
extractAIC(cox_all_yr_cdh)
extractAIC(cox_all_mo_cdh_loc)
extractAIC(cox_all_yr_cdh_loc)

# Plasmid Direct Ward model fit w/ and w/o location data
extractAIC(cox_all_mo_pdw)
extractAIC(cox_all_yr_pdw)
extractAIC(cox_all_mo_pdw_loc) # Best performing model
extractAIC(cox_all_yr_pdw_loc)

# Plasmid Direct Hospital model fit w/ and w/o location data
extractAIC(cox_all_mo_pdh) # Best performing model
extractAIC(cox_all_yr_pdh)
extractAIC(cox_all_mo_pdh_loc)
extractAIC(cox_all_yr_pdh_loc)

# Forest plot visualizations
# ggforest(cox_all_mo_cdw_loc, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios", fontsize=0.55)
ggforest(cox_all_yr_cdw_loc, data=as.data.frame(df.cdw), main="Clonal Direct Ward Hazard Ratios", fontsize=0.55)
ggforest(cox_all_mo_cdh_loc, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios", fontsize=0.55)
ggforest(cox_all_yr_cdh_loc, data=as.data.frame(df.cdh), main="Clonal Direct Hospital Hazard Ratios", fontsize=0.55)

ggforest(cox_all_mo_pdw_loc, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios", fontsize=0.55)
ggforest(cox_all_yr_pdw_loc, data=as.data.frame(df.pdw), main="Plasmid Direct Ward Hazard Ratios", fontsize=0.55)
ggforest(cox_all_mo_pdh_loc, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios", fontsize=0.55)
ggforest(cox_all_yr_pdh_loc, data=as.data.frame(df.pdh), main="Plasmid Direct Hospital Hazard Ratios", fontsize=0.55)



################################################################################
################################################################################




