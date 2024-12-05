##This is for the WHO SAGE project (Benefit-risk analysis of maternal RSV vaccine)
##NB this is provisional and confidential
##For any questions, feel free to contact Ayaka Monoi(ayaka.monoi@lshtm.ac.uk)

######VE#########from ve_fit.R
##posterior samples of VE
ve_post <-resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  group_by(group,t)

##select 1000 samples of VE vs severe RSV disease at age 0 to 11(month)
# prepare list
data_list <- list()

# from age 0 to 11(repeat from t = 15 for 12 times by 30)
for (i in seq(15, 365, by = 30)) {
  # sample VE vs severe at t = i
  ve <- ve_post %>%
    filter(group == "severe") %>%
    filter(t == i) %>%
    mutate(age = (i+15)/30 - 1)
  # return the sampled VE
  data_list[[(i+15)/30]] <- ve
}

data_list
#output here is data_list, not ve.
# data_list is 1000 draws of the vaccine efficacy against SEVERE disease
# in the calculation on Slide 7 it uses vaccine efficacy against SEVERE disease, not death

##################################################
#################RSV-associated deaths (RSV_burden.R)#########
#SA
death_sa

#Kenya
death_kenya

death_comb <- rbind(death_sa, death_kenya)
# country
countries <- unique(death_comb$country)

################ main analysis if VE follows Erlang-2 distribution###########
#Calculate deaths averted under 1 year (0to11 month)
dat_list <- list()

benef_by_country <- list()

# calculate averted death by country
# wait what? this isnt vaccine efficacy against preventing death??
# yep, using vaccine efficacy against severe disease in this calculation
# as found on slide 7

for (country in countries) {
  # filter by country
  death_comb_country <- death_comb %>% filter(country == !!country)

  # averted death by month
  dat_list <- list()
  for (i in seq(1:12)) {
    ve_age <- data_list[[i]] %>% as_tibble() %>% dplyr::select(VE_t)
    death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
    death_avert_age <- ve_age$VE_t * death_age$death 
    dat_list[[i]] <- death_avert_age
  }

#under 1 year deaths averted by sample
  sum_samples <- numeric(1000)
  for (j in 1:1000) {
    sample_age <- numeric(12)
    for (i in 1:12) {
      samp_age <- sample(unlist(dat_list[[i]]), 1)
      sample_age[i] <- samp_age # drawing one sample from every GA
    }
    sum_samples[j] <- sum(sample_age) # and repeating this 1000 times
  }

  benef_median <- quantile(sum_samples, .5)
  benef_low <- quantile(sum_samples, .025)
  benef_high <- quantile(sum_samples, .975)

  benef_by_country[[country]] <- list(
    median = benef_median,
    low = benef_low,
    high = benef_high
  )
}

#Country-specific RSV-associated deaths averted
benef_by_country



####posterior samples for benefits in ZAF in main scenario####
# Define the target country
target_country <- "ZAF"

# Filter data for the target country
death_comb_country <- death_comb %>% filter(country == target_country)

# Averted death by month
dat_list <- list()
for (i in seq(1:12)) {
  ve_age <- data_list[[i]] %>% as_tibble() %>% dplyr::select(VE_t)
  death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
  death_avert_age <- ve_age$VE_t * death_age$death
  dat_list[[i]] <- death_avert_age
}

# Under 1 year deaths averted by sample
post_benef <- numeric(1000)
for (j in 1:1000) {
  sample_age <- numeric(12)
  for (i in 1:12) {
    samp_age <- sample(unlist(dat_list[[i]]), 1)
    sample_age[i] <- samp_age
  }
  post_benef[j] <- sum(sample_age)
}
post_benef%>%median()




###########Sensitivity analysis if VE = 80% for 1 year###############
#Calculate deaths averted under 1 year (0to11 month)
da_list <- list()

benef_sens <- list()

# calculate averted death by country
for (country in countries) {
  # filter by country
  death_comb_country <- death_comb %>% filter(country == !!country)

  # averted death by month
  da_list <- list()
  for (i in seq(1:12)) {
    ve <- 0.8
    death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
    death_avert_age <- ve * death_age$death
    da_list[[i]] <- death_avert_age
  }

  #under 1 year deaths averted by sample
  sum_sample <- numeric(1000)
  for (j in 1:1000) {
    sample_age <- numeric(12)
    for (i in 1:12) {
      samp_age <- sample(unlist(da_list[[i]]), 1)
      sample_age[i] <- samp_age
    }
    sum_sample[j] <- sum(sample_age)
  }

  bene_median <- quantile(sum_sample, .5)
  bene_low <- quantile(sum_sample, .025)
  bene_high <- quantile(sum_sample, .975)

  benef_sens[[country]] <- list(
    median = bene_median,
    low = bene_low,
    high = bene_high
  )
}

#Country-specific RSV-associated deaths averted
benef_sens



