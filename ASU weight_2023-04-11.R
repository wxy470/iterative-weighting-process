library(car)
library(questionr)


rm(list = ls())
setwd("~/working directory")

# ==================== ASU dataset preparation ==================== 
ASU_data <- read.csv("ASU dataset.csv") # employed only

# gender
ASU_data$gender_cat2 <- ifelse(ASU_data$gender=='Female', 1, 0)
table(ASU_data$gender_cat2)

# age
ASU_data$age_cat4 <- car::recode(ASU_data$age, "lo:29=1; 30:44=2; 45:64=3; 65:hi=4")
table(ASU_data$age_cat4)

# HH income
ASU_data$inc_cat3 <- car::recode(ASU_data$inc, "'Less than $35,000'=1; '$35,000 to $74,999'=2; '$75,000 to $99,999'=2; else=3")
table(ASU_data$inc_cat3)

# ==================== ACS wt margin ==================== 
gender_age <- matrix(c(19221912,	27925607,	30791876,	5384125,
                       18195594,	24254513,	27562872,	4380535), nrow = 2, byrow = TRUE)
hhincome <- matrix(c(32146157, 52084490, 43314083))


# ==================== iterative weighting ==================== 
wt1 <- matrix(NA, nrow = nrow(ASU_data), ncol = 1000)
wt2 <- matrix(NA, nrow = nrow(ASU_data), ncol = 1000)

wt2[,1] <- 1 # set initial weight
(n_US <- nrow(ASU_data)) # sample size

for (k in 2:200) {
  N_US <- sum(gender_age) # population size
  ratio <- n_US/N_US # sample:population ratio
  wt_tab <- gender_age/wtd.table(x=ASU_data$gender_cat2, y=ASU_data$age_cat4, weights=wt2[,k-1]) * ratio # cell weight factor
  wt_max1 <- abs(max(wt_tab-1)) # identify the largest cell weight factor
  for (i in 1:n_US) {wt1[i,k] <- wt2[i,k-1] * wt_tab[ASU_data$gender_cat2[i]+1, ASU_data$age_cat4[i]]} # update weights with the cell weight factors
  
  N_US <- sum(hhincome)
  ratio <- n_US/N_US
  wt_tab <- hhincome/t(t(wtd.table(x=ASU_data$inc_cat3, weights=wt1[,k]))) * ratio
  wt_max2 <- abs(max(wt_tab-1))
  for (i in 1:n_US) {wt2[i,k] <- wt1[i,k] * wt_tab[ASU_data$inc_cat3[i],1]}
  
  print(k)
  k_conv <- k
  if (max(c(wt_max1, wt_max2)) < 0.000001) break # converge check
  
}

ASU_data$wt <- wt2[,k_conv] # assgin the final weight to the dataset
summary(ASU_data$wt)

# plot weight distributions
ggplot(ASU_data, aes(x=wt)) + 
  geom_histogram(aes(y=..count..), 
                 color="#00a0ff", fill="#74bbf7",
                 binwidth = 1, center = 0.5) +
  stat_bin(binwidth = 1, center = 0.5, 
           aes(y=..count.., label=..count..), geom="text", vjust=-0.5, size = 3) +
  labs(title="US sample weight", x="Sample weight", y="Count") +
  theme_bw() + 
  theme(panel.grid.minor=element_blank(), 
        plot.title = element_text(face="bold", hjust = 0.5)) + 
  theme(panel.grid.minor=element_blank())
ggsave("Figure/US_weight.png", width = 9, height = 6)


# ==================== trim weight ==================== 
(wt_threshold <- median(ASU_data$wt) + 6*IQR(ASU_data$wt)) # calculate trimming threshold
table(ASU_data$wt > wt_threshold) # check the number of cases with extreme weights (i.e., above the threshold)

ASU_data$wt_trim <- ifelse(ASU_data$wt > wt_threshold, wt_threshold, ASU_data$wt) # assign the threshold for all extreme weights
ASU_data$wt_trim <- ASU_data$wt_trim * n_US/sum(ASU_data$wt_trim) # rescale the weights to make the weight sum equal to sample size
summary(ASU_data$wt_trim)

# plot trimmed weights
ggplot(ASU_data, aes(x=wt_trim)) + 
  geom_histogram(aes(y=..count..), 
                 color="#00a0ff", fill="#74bbf7",
                 binwidth = 1, center = 0.5) +
  stat_bin(binwidth = 1, center = 0.5, 
           aes(y=..count.., label=..count..), geom="text", vjust=-0.5, size = 3) +
  labs(title="US sample weight (trimmed)", x="Sample weight", y="Count") +
  theme_bw() + 
  theme(panel.grid.minor=element_blank(), 
        plot.title = element_text(face="bold", hjust = 0.5)) + 
  theme(panel.grid.minor=element_blank())
ggsave("Figure/US_weight_trim.png", width = 9, height = 6)


# ==================== check SED distribution ==================== 
SED_cnt <- function(var) {
  SED_cnt <- cbind(table(ASU_data[[var]]),
                   wtd.table(ASU_data[[var]], weights = ASU_data$wt), 
                   wtd.table(ASU_data[[var]], weights = ASU_data$wt_trim))
  return(SED_cnt)
}

SED_wt <- rbind(SED_cnt('gender_cat2'), SED_cnt('age_cat4'), SED_cnt('inc_cat3'))
write.csv(SED_wt, "US_weighted_SED.csv", row.names=FALSE)



#
