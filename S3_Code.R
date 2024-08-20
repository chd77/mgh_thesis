# Survival Curves 
# Built manually from survival percentages and exposure metrics


{
  # Clonal Data:
  # Sorting and calculating survival percentage from time-to-infection data by day
  # Survival Curve - Direct Ward
  SC_Clonal_DirWard <- Clonal_DirWard %>% 
    select(Direct_Ward_Contact_Hrs) %>% 
    arrange(Direct_Ward_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Direct_Ward_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Direct Hospital
  SC_Clonal_DirHosp <- Clonal_DirHosp %>% 
    select(Direct_Hosp_Contact_Hrs) %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Indirect Ward
  SC_Clonal_IndWard <- Clonal_IndWard %>% 
    select(Indirect_Ward_Contact_Hrs) %>% 
    arrange(Indirect_Ward_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Indirect_Ward_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Indirect Hospital
  SC_Clonal_IndHosp <- Clonal_IndHosp %>% 
    select(Indirect_Hosp_Contact_Hrs) %>% 
    arrange(Indirect_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Indirect_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Indirect Unrestricted
  SC_Clonal_IndUnres <- Clonal_IndUnres %>% 
    select(Indirect_Unrestricted_Contact_Hrs) %>% 
    arrange(Indirect_Unrestricted_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Indirect_Unrestricted_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  
  
  
  
  
  # Plasmid Data:
  # Sorting and calculating survival percentage from time-to-infection data by day
  # Survival Curve - Direct Ward
  SC_Plasmid_DirWard <- Plasmid_DirWard %>% 
    select(Direct_Ward_Contact_Hrs) %>% 
    arrange(Direct_Ward_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),
           survival_per = survival/nrow(.),
           day = round(Direct_Ward_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Direct Hospital
  SC_Plasmid_DirHosp <- Plasmid_DirHosp %>% 
    select(Direct_Hosp_Contact_Hrs) %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Indirect Ward
  SC_Plasmid_IndWard <- Plasmid_IndWard %>% 
    select(Indirect_Ward_Contact_Hrs) %>% 
    arrange(Indirect_Ward_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Indirect_Ward_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Indirect Hospital
  SC_Plasmid_IndHosp <- Plasmid_IndHosp %>% 
    select(Indirect_Hosp_Contact_Hrs) %>% 
    arrange(Indirect_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),          
           survival_per = survival/nrow(.),
           day = round(Indirect_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Survival Curve - Indirect Unrestricted
  SC_Plasmid_IndUnres <- Plasmid_IndUnres %>% 
    select(Indirect_Unrestricted_Contact_Hrs) %>% 
    arrange(Indirect_Unrestricted_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Indirect_Unrestricted_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n)
}





# Combined Survivor Plots:
# Clonal survival percentage over time
ggplot() +
  geom_step(SC_Clonal_DirWard, mapping=aes(x=day, y=survival_per, color="red")) + 
  geom_step(SC_Clonal_DirHosp, mapping=aes(x=day, y=survival_per, color="orange")) + 
  geom_step(SC_Clonal_IndWard, mapping=aes(x=day, y=survival_per, color="yellow")) + 
  geom_step(SC_Clonal_IndHosp, mapping=aes(x=day, y=survival_per, color="lightblue")) +
  geom_step(SC_Clonal_IndUnres, mapping=aes(x=day, y=survival_per, color="forestgreen")) +
  xlab('Days Since Admission') + ylab('Survival Percentage') + 
  ggtitle("Clonal Exposure Metric Risk") + 
  scale_colour_manual(name = 'Exposure Type', 
                      values =c('red'='red','orange'='orange',
                                'yellow'='yellow','lightblue'='lightblue',
                                'forestgreen'='forestgreen'), 
                      labels = c('Indir. Unres','Indir. Hosp','Direct Hosp',
                                 'Direct Ward','Indir. Ward'))

# Plasmid survival percentage over time
ggplot() +
  geom_step(SC_Plasmid_DirWard, mapping=aes(x=day, y=survival_per, color="red")) + 
  geom_step(SC_Plasmid_DirHosp, mapping=aes(x=day, y=survival_per, color="orange")) + 
  geom_step(SC_Plasmid_IndWard, mapping=aes(x=day, y=survival_per, color="yellow")) + 
  geom_step(SC_Plasmid_IndHosp, mapping=aes(x=day, y=survival_per, color="lightblue")) +
  geom_step(SC_Plasmid_IndUnres, mapping=aes(x=day, y=survival_per, color="forestgreen")) +
  xlab('Days Since Admission') + ylab('Survival Percentage') + 
  ggtitle("Plasmid Exposure Metric Risk") + 
  scale_colour_manual(name = 'Exposure Type', 
                      values =c('red'='red','orange'='orange',
                                'yellow'='yellow','lightblue'='lightblue',
                                'forestgreen'='forestgreen'), 
                      labels = c('Indir. Unres','Indir. Hosp','Direct Hosp',
                                 'Direct Ward','Indir. Ward'))


# Clonal vs. Plasmid survival percentage over time
ggplot() +
  geom_step(SC_Clonal_IndUnres, mapping=aes(x=day, y=survival_per, color="blue")) +
  geom_step(SC_Plasmid_IndUnres, mapping=aes(x=day, y=survival_per, color = "red")) +
  xlab('Days Since Admission') + ylab('Survival Percentage') + 
  labs(title = "Clonal vs. Plasmid Risk",
       subtitle = "Indirect Unrestricted Exposure") + 
  scale_colour_manual(name = 'CRE Type', 
                      values =c('blue'='blue','red'='red'), 
                      labels = c('Clonal','Plasmid'))




# Survival Curves by Hospital - Direct Hospital
# Clonal
{
  C_H1 <- Clonal_DirHosp %>% 
    rename(Hospital = Direct_Hosp_Contact_Site) %>% 
    filter(!is.na(Hospital),
           Hospital == "H1") %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n) %>%
    select(day, survival, survival_per, Hospital)
  
  C_H2 <- Clonal_DirHosp %>% 
    rename(Hospital = Direct_Hosp_Contact_Site) %>% 
    filter(!is.na(Hospital),
           Hospital == "H2") %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n) %>%
    select(day, survival, survival_per, Hospital)
  
  C_H3 <- Clonal_DirHosp %>% 
    rename(Hospital = Direct_Hosp_Contact_Site) %>% 
    filter(!is.na(Hospital),
           Hospital == "H3") %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n) %>%
    select(day, survival, survival_per, Hospital)
  
  C_H4 <- Clonal_DirHosp %>% 
    rename(Hospital = Direct_Hosp_Contact_Site) %>% 
    filter(!is.na(Hospital),
           Hospital == "H4") %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n) %>%
    select(day, survival, survival_per, Hospital)
  
  C_H5 <- Clonal_DirHosp %>% 
    rename(Hospital = Direct_Hosp_Contact_Site) %>% 
    filter(!is.na(Hospital),
           Hospital == "H5") %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n) %>%
    select(day, survival, survival_per, Hospital)
  
  C_H6 <- Clonal_DirHosp %>%
    rename(Hospital = Direct_Hosp_Contact_Site) %>% 
    filter(!is.na(Hospital),
           Hospital == "H6") %>% 
    arrange(Direct_Hosp_Contact_Hrs) %>% 
    mutate(survival = nrow(.) - row_number(),         
           survival_per = survival/nrow(.),
           day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
    add_count(day) %>% rename(rate = n) %>%
    select(day, survival, survival_per, Hospital)
  
  C_H_comb <- bind_rows(C_H1,C_H2)
  C_H_comb <- bind_rows(C_H_comb,C_H3)
  C_H_comb <- bind_rows(C_H_comb,C_H4)
  C_H_comb <- bind_rows(C_H_comb,C_H5)
  C_H_comb <- bind_rows(C_H_comb,C_H6)
  }

ggplot(C_H_comb) +
  geom_step(aes(x = day, y = survival_per, group = Hospital, color = Hospital)) + 
  xlab('Days Since Admission') + ylab('Survival Percentage') + 
  ggtitle("Survival Curve by Hospital - Clonal Direct Hospital")

rm(C_H1,C_H2,C_H3,C_H4,C_H5,C_H6)



# Plasmid
{
  P_H1 <- Plasmid_DirHosp %>% 
  rename(Hospital = Direct_Hosp_Contact_Site) %>% 
  filter(!is.na(Hospital),
         Hospital == "H1") %>% 
  arrange(Direct_Hosp_Contact_Hrs) %>% 
  mutate(survival = nrow(.) - row_number(),         
         survival_per = survival/nrow(.),
         day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
  add_count(day) %>% rename(rate = n) %>%
  select(day, survival, survival_per, Hospital)

P_H2 <- Plasmid_DirHosp %>% 
  rename(Hospital = Direct_Hosp_Contact_Site) %>% 
  filter(!is.na(Hospital),
         Hospital == "H2") %>% 
  arrange(Direct_Hosp_Contact_Hrs) %>% 
  mutate(survival = nrow(.) - row_number(),         
         survival_per = survival/nrow(.),
         day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
  add_count(day) %>% rename(rate = n) %>%
  select(day, survival, survival_per, Hospital)

P_H3 <- Plasmid_DirHosp %>% 
  rename(Hospital = Direct_Hosp_Contact_Site) %>% 
  filter(!is.na(Hospital),
         Hospital == "H3") %>% 
  arrange(Direct_Hosp_Contact_Hrs) %>% 
  mutate(survival = nrow(.) - row_number(),         
         survival_per = survival/nrow(.),
         day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
  add_count(day) %>% rename(rate = n) %>%
  select(day, survival, survival_per, Hospital)

P_H4 <- Plasmid_DirHosp %>% 
  rename(Hospital = Direct_Hosp_Contact_Site) %>% 
  filter(!is.na(Hospital),
         Hospital == "H4") %>% 
  arrange(Direct_Hosp_Contact_Hrs) %>% 
  mutate(survival = nrow(.) - row_number(),         
         survival_per = survival/nrow(.),
         day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
  add_count(day) %>% rename(rate = n) %>%
  select(day, survival, survival_per, Hospital)

P_H5 <- Plasmid_DirHosp %>% 
  rename(Hospital = Direct_Hosp_Contact_Site) %>% 
  filter(!is.na(Hospital),
         Hospital == "H5") %>% 
  arrange(Direct_Hosp_Contact_Hrs) %>% 
  mutate(survival = nrow(.) - row_number(),         
         survival_per = survival/nrow(.),
         day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
  add_count(day) %>% rename(rate = n) %>%
  select(day, survival, survival_per, Hospital)

P_H6 <- Plasmid_DirHosp %>%
  rename(Hospital = Direct_Hosp_Contact_Site) %>% 
  filter(!is.na(Hospital),
         Hospital == "H6") %>% 
  arrange(Direct_Hosp_Contact_Hrs) %>% 
  mutate(survival = nrow(.) - row_number(),         
         survival_per = survival/nrow(.),
         day = round(Direct_Hosp_Contact_Hrs/24)) %>% 
  add_count(day) %>% rename(rate = n) %>%
  select(day, survival, survival_per, Hospital)

P_H_comb <- bind_rows(P_H1,P_H2)
P_H_comb <- bind_rows(P_H_comb,P_H3)
P_H_comb <- bind_rows(P_H_comb,P_H4)
P_H_comb <- bind_rows(P_H_comb,P_H5)
P_H_comb <- bind_rows(P_H_comb,P_H6)
}

ggplot(P_H_comb) +
  geom_step(aes(x = day, y = survival_per, group = Hospital, color = Hospital)) + 
  xlab('Days Since Admission') + ylab('Survival Percentage') + 
  ggtitle("Survival Curve by Hospital - Plasmid Direct Hospital")

rm(P_H1,P_H2,P_H3,P_H4,P_H5,P_H6)

