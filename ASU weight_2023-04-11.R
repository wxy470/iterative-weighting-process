library(car)
library(questionr)


rm(list = ls())
setwd("~/OneDrive - Georgia Institute of Technology/01-Research/00-Working/US-Germany comparison/ASU")

wtd.table <- questionr::wtd.table
wtd.mean <- questionr::wtd.mean


# ==================== ASU dataset preparation ==================== 
ASU_data <- read.csv("USAtwFULLTIME1.csv")

# employment check
table(ASU_data$worker_now) # N=57 no
table(ASU_data$w3_empl_now) # N=928 are working full-time, there are also part-time worker & unemployed & furloughed without pay

# gender
table(ASU_data$gender)
ASU_data$gender_cat2 <- ifelse(ASU_data$gender=='Female', 1, 0)
table(ASU_data$gender_cat2)

# age
table(ASU_data$ageGr)
table(ASU_data$age)
ASU_data$age_cat4 <- car::recode(ASU_data$age, "lo:29=1; 30:44=2; 45:64=3; 65:hi=4")
table(ASU_data$age_cat4)

table(ASU_data$gender_cat2, ASU_data$age_cat4)

# HH income
table(ASU_data$EconStat)
table(ASU_data$inc)

ASU_data$inc_cat3 <- car::recode(ASU_data$inc, "'Less than $35,000'=1; '$35,000 to $74,999'=2; '$75,000 to $99,999'=2; else=3")
table(ASU_data$inc_cat3)

# ==================== ACS wt margin ==================== 
gender_age <- matrix(c(19221912,	27925607,	30791876,	5384125,
                       18195594,	24254513,	27562872,	4380535), nrow = 2, byrow = TRUE)
hhincome <- matrix(c(32146157, 52084490, 43314083))


# ==================== weighting ==================== 
wt1 <- matrix(NA, nrow = nrow(ASU_data), ncol = 1000)
wt2 <- matrix(NA, nrow = nrow(ASU_data), ncol = 1000)

wt2[,1] <- 1
(n_US <- nrow(ASU_data))


for (k in 2:200) {
  N_US <- sum(gender_age)
  ratio <- n_US/N_US
  wt_tab <- gender_age/wtd.table(x=ASU_data$gender_cat2, y=ASU_data$age_cat4, weights=wt2[,k-1]) * ratio
  wt_max1 <- abs(max(wt_tab-1))
  for (i in 1:n_US) {wt1[i,k] <- wt2[i,k-1] * wt_tab[ASU_data$gender_cat2[i]+1, ASU_data$age_cat4[i]]}
  
  N_US <- sum(hhincome)
  ratio <- n_US/N_US
  wt_tab <- hhincome/t(t(wtd.table(x=ASU_data$inc_cat3, weights=wt1[,k]))) * ratio
  wt_max2 <- abs(max(wt_tab-1))
  for (i in 1:n_US) {wt2[i,k] <- wt1[i,k] * wt_tab[ASU_data$inc_cat3[i],1]}
  
  print(k)
  k_conv <- k
  if (max(c(wt_max1, wt_max2)) < 0.000001) break
  
}

ASU_data$wt <- wt2[,k_conv]
summary(ASU_data$wt)


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
median(ASU_data$wt) + 6*IQR(ASU_data$wt)
table(ASU_data$wt > median(ASU_data$wt) + 6*IQR(ASU_data$wt))
wt_threshold <- median(ASU_data$wt) + 6*IQR(ASU_data$wt)

ASU_data$wt_trim <- ifelse(ASU_data$wt > wt_threshold, wt_threshold, ASU_data$wt)
ASU_data$wt_trim <- ASU_data$wt_trim * n_US/sum(ASU_data$wt_trim)
summary(ASU_data$wt_trim)

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
