######## ESTIMATING VDPV DATE OF SEEDING ######## 
# all the analysis presented in the paper

rm(list=ls())	
#install.packages(pkgs=c("dplyr","data.table","epiDisplay","wesanderson","patchwork"))
library(stringr)
library(dplyr)
library(magrittr)
library(lubridate)
library(tidyverse)
library(janitor)
library(data.table)
library(epiDisplay)
library(wesanderson)
library(MASS)
library(cowplot)
library(sjPlot)

setwd("~/Dropbox/PolioPreventingOutbreaks/14_VDPV")
datab <- read_csv("seeding_final_emergence_MA.csv") %>% clean_names()

datab$date_first <- dmy(datab$date_first)
datab$date_last <- dmy(datab$date_last)
VDPV_events <- datab

### REMOVE NON EMERGENT EVENTS ###
names(VDPV_events)
table(datab$emergence)
VDPV_events$emergence[is.na(VDPV_events$emergence)] <- "presumed emergence"
table(VDPV_events$emergence)
# prevent double counting of events that are international spread
VDPV_events <- subset(VDPV_events, emergence!="international spread" & country!="PHILIPPINES")

############################################################
# *** first section - summary stats and regression model ***
############################################################

#### Stats for Manuscript ####
length(table(VDPV_events$country))
table(VDPV_events$serotype)/dim(VDPV_events)[1]
table(VDPV_events$start_cat)

## Detection categories - less than a year to detect? ## 
VDPV_events$detection_cat <- ""
VDPV_events$detection_cat[VDPV_events$days_to_detect<365]<- "<1 year"
VDPV_events$detection_cat[VDPV_events$days_to_detect>365 & VDPV_events$days_to_detect<1461]<- "1-4 years"
VDPV_events$detection_cat[VDPV_events$days_to_detect>1461]<- ">4 years"
prop.table(table(VDPV_events$detection_cat))

## days to contain outbreak 
VDPV_events$duration_cat <- ""
VDPV_events$duration_cat[VDPV_events$outbreak_days<365]<- "<1 year"
VDPV_events$duration_cat[VDPV_events$outbreak_days<120]<- "days120"

## Summary table ## 
type1 <- subset(VDPV_events, VDPV_events$serotype==1 & VDPV_events$nucleotide_first<30 &
                  !is.na(VDPV_events$npafp_first))
type1a <- data.table(type1)
x <- type1a %>%
  group_by(who_region,country) %>%
  dplyr::summarize(Frequency = n(),
                   total_cases = sum(cases),
                   median_duration = round(median(outbreak_days, na.rm=TRUE),0), 
                   lwr_duration= round(min(outbreak_days, na.rm=TRUE),0),
                   upr_duration= round(max(outbreak_days, na.rm=TRUE),0),
                   median_nucleotidediff= round(median(nucleotide_first, na.rm=TRUE),0),
                   lwr_nucleotidediff= round(min(nucleotide_first, na.rm=TRUE),0),
                   upr_nucleotidediff= round(max(nucleotide_first, na.rm=TRUE),0),
                   mean_npafp= round(mean(npafp_first, na.rm=TRUE),1))
median(type1a$nucleotide_first)
range(type1a$nucleotide_first)
type1a %$% ci(npafp_first, ci=0.95)

type2 <- subset(VDPV_events, VDPV_events$serotype==2 & VDPV_events$year_start>=2010 & 
                  VDPV_events$nucleotide_first<38)
type2a <- data.table(type2)
x <- type2a %>%
  group_by(who_region,country) %>%
  dplyr::summarize(Frequency = n(),
                   total_cases = sum(cases),
                   median_duration = round(median(outbreak_days, na.rm=TRUE),0), 
                   lwr_duration= round(min(outbreak_days, na.rm=TRUE),0),
                   upr_duration= round(max(outbreak_days, na.rm=TRUE),0),
                   median_nucleotidediff= round(median(nucleotide_first, na.rm=TRUE),0),
                   lwr_nucleotidediff= round(min(nucleotide_first, na.rm=TRUE),0),
                   upr_nucleotidediff= round(max(nucleotide_first, na.rm=TRUE),0),
                   mean_npafp= round(mean(npafp_first, na.rm=TRUE),1))
median(type2a$nucleotide_first)
range(type2a$nucleotide_first)
type2a %$% ci(npafp_first, ci=0.95)

type3 <- subset(VDPV_events, VDPV_events$serotype==3 & VDPV_events$nucleotide_first<30)
type3a <- data.table(type3)
x <- type3a %>%
  group_by(who_region,country) %>%
  dplyr::summarize(Frequency = n(),
                   total_cases = sum(cases),
                   median_duration = round(median(outbreak_days, na.rm=TRUE),0), 
                   lwr_duration= round(min(outbreak_days, na.rm=TRUE),0),
                   upr_duration= round(max(outbreak_days, na.rm=TRUE),0),
                   median_nucleotidediff= round(median(nucleotide_first, na.rm=TRUE),0),
                   lwr_nucleotidediff= round(min(nucleotide_first, na.rm=TRUE),0),
                   upr_nucleotidediff= round(max(nucleotide_first, na.rm=TRUE),0),
                   mean_npafp= round(mean(npafp_first, na.rm=TRUE),1))
median(type3a$nucleotide_first)
range(type3a$nucleotide_first)
type3a %$% ci(npafp_first, ci=0.95)

# overall totals (excluding type 2 prior to 2010)
tmp <- rbind(type1, type2, type3)
prop.table(table(tmp$serotype))

sum(tmp$cases)
median(tmp$outbreak_days)
range(tmp$outbreak_days)
median(tmp$nucleotide_first)
range(tmp$nucleotide_first)
tmp %$% ci(npafp_first, ci=0.95)

# results section in text stats
prop.table(table(type1$detection_cat))
prop.table(table(type2$detection_cat))
prop.table(table(type3$detection_cat))

prop.table(table(type1$duration_cat))
prop.table(table(type2$duration_cat))
prop.table(table(type3$duration_cat))

# regression modelling of the number of nc changes: count data, Negative binomial model
## variance > mean, data is over-dispersed

# not all obs had npafp_first
rdat <- VDPV_events[!is.na(VDPV_events$npafp_first),]
# to hopeully improve convergence re-scale the npafp_first parameter
rdat$npafp_first2 <- rdat$npafp_first-2

rdat$preswitch <- 1
rdat$preswitch[rdat$year_start>=2016] <- 0
rdat$st <- factor(rdat$serotype,levels=c("2","1","3"))

# First, lets test if random effects account for over-dispersion
rdat$id <- as.numeric(as.factor(rdat$iso_3_code))  # used to test random effects
table(rdat$id)
# sample 1 ID for a subset, report mean and var
tmp <- rdat[rdat$serotype==2 & rdat$country != "PHILIPPINES",]
oo <- match(c(1:26),tmp$id)
tmp$id[oo]
tmp$country[oo]
mean(tmp$nucleotide_first[oo],na.rm=T)
var(tmp$nucleotide_first[oo],na.rm=T)
hist(tmp$nucleotide_first[oo],na.rm=T)

# look at obs with lots of repeats...
mean(tmp$nucleotide_first[tmp$id==5])
var(tmp$nucleotide_first[tmp$id==5])
# variance is still much greater than the mean... looks like a negative binomial is better than mixed effects


#################################
######## TYPE 2 ANALYSIS ######## 
#################################

# Remove outiers and exclude type 2 data before definition changes
dat2 <- rdat[rdat$serotype==2 & rdat$nucleotide_first<38,]
dat2 <- subset(dat2, dat2$year_start>=2010)
range(dat2$nucleotide_first) #RE-SCALE: subtract 6, so range starts at 0
dat2$nc_adjusted <- dat2$nucleotide_first - 6
mean(dat2$nc_adjusted)
var(dat2$nc_adjusted) 
range(dat2$nc_adjusted)

# best model for NC adjusted 
nbmod2a <- glm.nb(nc_adjusted ~ npafp_first2*stool_start + classification, 
                 link=log, data=dat2, control=glm.control(maxit=1000))
summary(nbmod2a) 
exp(cbind(coef(nbmod2a), confint(nbmod2a)))

# diagnostic plots type 2 model
dat2$id <- as.numeric(as.factor(dat2$iso_3_code))
frequencies <- table(dat2$nc_adjusted)
freq <- data.frame(frequencies)
names(freq)<- c("num", "Freq")
mean(dat2$nc_adjusted, na.rm=TRUE) 
var(dat2$nc_adjusted)/mean(dat2$nc_adjusted) 
mean(dat2$nc_adjusted)^2/(var(dat2$nc_adjusted)-mean(dat2$nc_adjusted)) #k parameter
exp<-dnbinom(0:16,1, mu=3.15)*59
num <- 0:16
dat <- data.frame(num, exp)
dat <- merge(dat, freq, by="num", all.x = TRUE)
dat$Freq[is.na(dat$Freq)]<- 0
names(dat)<- c("nc", "Expected", "Observed")
dat_long <- gather(dat, "outcome", "freq", 2:3)

par(mfrow=c(1,3))
barplot(frequencies,ylab="Frequency",xlab="NC Mutations",col="gray", main = "Observed",
        ylim=c(0,20), xlim=c(0,16))
barplot(dpois(0:16,3.15)*59,names=as.character(0:16), ylab="Frequency", xlab="NC Mutations", 
        col="gray", main = "Expected Poisson", ylim=c(0,20),  xlim=c(0,16))
barplot(dnbinom(0:16,1, mu=3.15)*59,names=as.character(0:16), ylab="Frequency", xlab="NC Mutations", 
        col="gray", main = "Expected Negative Binomial", ylim=c(0,20),  xlim=c(0,16))

#observed vs. expected
p3 <- ggplot(data=dat_long, aes(x=nc, y=freq, fill=outcome)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Nucleotide Mutations") +
  ylab("Frequency of outbreaks") + 
  theme_bw() + theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) +
  scale_fill_manual("", values = c("Observed" = "gray", "Expected" = "black")) +
  scale_x_continuous(breaks=c(0:32), labels=c(6:38)) +
  ggtitle("Serotype 2: Expected vs. Observed Nucleotide Mutations")

#residuals
p1 <- ggplot(nbmod2a, aes(.fitted, .resid))+geom_point() +
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values") + ylab("Residuals") +
  ggtitle("Serotype 2: Residual vs. Fitted Plot") +
  theme_bw()

#q-q plot
p2 <- ggplot(nbmod2a, aes(qqnorm(.stdresid)[[1]], .stdresid))+geom_point(na.rm = TRUE) +
  geom_abline(aes(qqline(.stdresid))) +
  xlab("Theoretical Quantiles") + ylab("Standardized Residuals") +
  ggtitle("Serotype 2: Normal Q-Q") + 
  theme_bw()



#################################
######## TYPES 1/3 ANALYSIS #####
#################################
dat13 <- rdat[rdat$serotype!=2 & rdat$nucleotide_first<30,]
mean(dat13$nucleotide_first)
var(dat13$nucleotide_first) 
range(dat13$nucleotide_first) 
dat13$nc_adjusted <- dat13$nucleotide_first - 9 #RE-SCALE so range starts at zero
mean(dat13$nc_adjusted)
var(dat13$nc_adjusted) 
range(dat13$nc_adjusted) 

# nb model
nbmod13 <- glm.nb(nc_adjusted ~ npafp_first2*stool_start + classification,
                  link=log, data=dat13, control=glm.control(maxit=1000))
summary(nbmod13) 
exp(cbind(coef(nbmod13), confint(nbmod13)))


# diagnostic plots types 1/3 model
frequencies <- table(dat13$nc_adjusted)
freq <- data.frame(frequencies)
names(freq)<- c("num", "Freq")
mean(dat13$nc_adjusted) 
exp<-dnbinom(0:18,1, mu=7.467)*15
num <- 0:18
dat <- data.frame(num, exp)
dat <- merge(dat, freq, by="num", all.x = TRUE)
dat$Freq[is.na(dat$Freq)]<- 0
names(dat)<- c("nc", "Expected", "Observed")
dat_long <- gather(dat, "outcome", "freq", 2:3)

ggplot(data=dat_long, aes(x=nc, y=freq, fill=outcome)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Nucleotide Mutations") +
  ylim(0,4) +
  ylab("Frequency of outbreaks") + 
  theme_bw() + theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) +
  scale_fill_manual("", values = c("Observed" = "gray", "Expected" = "black")) +
  scale_x_continuous(breaks=c(0:18), labels=c(9:27)) +
  ggtitle("Serotypes 1 & 3: Expected vs. Observed Nucleotide Mutations")

#residuals
p1 <- ggplot(nbmod13, aes(.fitted, .resid))+geom_point() +
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values") + ylab("Residuals") +
  ggtitle("Serotypes 1 & 3: Residual vs. Fitted Plot") +
  theme_bw()

#q-q plot
p2 <- ggplot(nbmod13, aes(qqnorm(.stdresid)[[1]], .stdresid))+geom_point(na.rm = TRUE) +
  geom_abline(aes(qqline(.stdresid))) +
  xlab("Theoretical Quantiles") + ylab("Standardized Residuals") +
  ggtitle("Serotypes 1 & 3: Normal Q-Q") + 
  theme_bw()



##############################
## DETECTION TIME ESTIMATES ##
##############################

# Assume 1st mutation 1 happens instantaneously, then calculate based on n-1 mutations and mutation rate 
# 1) create dataframe for each serotype 2 and 1/3, rows = obs, col = 1000 samples
# 2) code below will estimate rate, as before but also generate 1000 samples
# 3) at the end have a histogram and CI for the population estimate

# Calculate estimated seeding date
# Set mutation rate (1.14*10^-2/365.25) * number of nucleotides in VP1 genome (365.25 days in the year)
Nick_rate <- (1.14*10^-2*906)/365.25

# Estimate seeding date for each isolate
VDPV_events$detection_date <- VDPV_events$date_first
VDPV_events$emergence_date1 <- rep(as.Date("0000-01-01"), nrow(VDPV_events))
VDPV_events$emergence_date_max1 <- rep(as.Date("0000-01-01"), nrow(VDPV_events))
VDPV_events$emergence_date_min1 <- rep(as.Date("0000-01-01"), nrow(VDPV_events))
VDPV_events$probability_after4 <- rep(0)
VDPV_events <- as.data.frame(VDPV_events)
type13 <- dat13
type2 <- dat2

# VDPV2 first 
iter <- 1000
VDPV_eventssample <- matrix(rep(0,iter*dim(VDPV_events)[1]),nrow=dim(VDPV_events)[1],ncol=iter)
for (i in 1:nrow(VDPV_events)) {
  n <- VDPV_events$nucleotide_first[i]
  if(is.na(n)==FALSE) { 
    k <- as.Date(VDPV_events$detection_date[i] - qgamma(p=0.5,shape=(n-1),rate=Nick_rate)) # date of detection - estimated age of virus (days) (median)
    VDPV_events$emergence_date1[i] <- k
    
    k2 <- as.Date(VDPV_events$detection_date[i] - qgamma(0.025, (n-1), Nick_rate)) # lower
    VDPV_events$emergence_date_max1[i] <- k2
    
    k3 <- as.Date(VDPV_events$detection_date[i] - qgamma(0.975, (n-1), Nick_rate)) # upper
    VDPV_events$emergence_date_min1[i] <- k3

    l <- as.numeric((as.Date(VDPV_events$emergence_date1[i] + 1461) - as.Date(VDPV_events$detection_date[i]))) 
    #1461 = 4 years; >4 years from seeding
    
    # create sample of detection time from data
    VDPV_eventssample[i,] <- rgamma(iter,shape=(n-1),rate=Nick_rate)
  } 
  else {
    VDPV_events$emergence_date1[i] <- NA
    VDPV_events$emergence_date_min1[i] <- NA
    VDPV_events$emergence_date_max1[i] <- NA
    VDPV_events$probability_after4[i] <- NA }
}

hist(VDPV_eventssample)  # full set of data

# take sample 10...
# the time of emergence varies from 150-600 days, ie there's a lot of uncertainty here
# which means that it's important to report the mean and 95% UI 
hist(VDPV_eventssample[10,])

# provide an estimate of 'general emergence' based on these samples
# type 2
# if we resample the events 1000 times and report the distribution of means
tmp2 <- subset(VDPV_eventssample,VDPV_events$serotype==2)
hist(rowMeans(tmp2),xlim=c(0,3000))
# sampling
nsamples <- 10000
ss <- sample.int(nrow(tmp2),nsamples,replace = T)
tmps2 <- tmp2[ss,]
tmpss2 <- matrix(rep(0,nsamples*nsamples),nrow=nsamples)
# and then resample
for(i in 1:nsamples){
  tt <- sample.int(1000,nsamples,replace = T) # resample
  tmpss2[i,] <- tmps2[i,tt]
}
dim(tmpss2)
# create an ECDF
mydis <- ecdf(tmpss2) 
# calculate the pr(x<4years)
y <- mydis(seq(0,4,0.5)*365)
plot(seq(0,4,0.5),y)
mydis(4*365)
hist(rowMeans(tmpss2),xlim=c(0,3000))
quantile(rowMeans(tmpss2),probs=c(0.025,0.5,0.975))

# type 1
tmp1 <- subset(VDPV_eventssample,VDPV_events$serotype==1)
nsamples <- 10000
ss <- sample.int(nrow(tmp1),nsamples,replace = T)
tmps1 <- tmp1[ss,]
tmpss1 <- matrix(rep(0,nsamples*nsamples),nrow=nsamples)
for(i in 1:nsamples){
  tt <- sample.int(1000,nsamples,replace = T) # resample
  tmpss1[i,] <- tmps1[i,tt]
}
quantile(rowMeans(tmpss1),probs=c(0.025,0.5,0.975))
mydis <- ecdf(tmpss1) 
# calculate the pr(x<4years)
mydis(4*365)

# type 3
tmp3 <- subset(VDPV_eventssample,VDPV_events$serotype==3)
nsamples <- 10000
ss <- sample.int(nrow(tmp3),nsamples,replace = T)
tmps3 <- tmp3[ss,]
tmpss3 <- matrix(rep(0,nsamples*nsamples),nrow=nsamples)
for(i in 1:nsamples){
  tt <- sample.int(1000,nsamples,replace = T) # resample
  tmpss3[i,] <- tmps3[i,tt]
}
quantile(rowMeans(tmpss3),probs=c(0.025,0.5,0.975))
mydis <- ecdf(tmpss3) 
# calculate the pr(x<4years)
mydis(4*365)


##############################
#### Graphs and Figures #####
#############################

# figure 1a - differences between different levels of stool
# to generate the predictions plot predict directly from the model using new data
range(dat2$npafp_first2)
range(dat2$stool_start)
table(dat2$classification)
new_data <- expand.grid(npafp_first2=seq(-1,15,1),stool_start=c(70,80,90),
                        classification=c("afp"))  # assume all are afp ##
new_data$nc_adjusted <- NA
new_data$nc_adjusted <- predict(nbmod2a,newdata=new_data,type=c("response"))
tmp <-  predict(nbmod2a,newdata=new_data,type=c("response"),se.fit=T)
new_data$nc_adjusted_lwr <- tmp$fit-1.96*tmp$se.fit
new_data$nc_adjusted_upr <- tmp$fit+1.96*tmp$se.fit

p1 <- ggplot(new_data,aes(x=npafp_first2+2,y=nc_adjusted+6, group=as.factor(stool_start),
                          colour=as.factor(stool_start))) + geom_line(size=0.5)

p1a <- p1 + geom_smooth(alpha=0.3,aes(ymin=new_data$nc_adjusted_lwr+6, ymax=new_data$nc_adjusted_upr+6,
                                      fill=as.factor(stool_start), colour=as.factor(stool_start)), stat = "identity") + 
  coord_cartesian(ylim = c(2,60), xlim=c(2,15)) + 
  labs(x="Non-polio AFP rate (per 100,000 children <15 years)", y="Nucleotide differences",
       colour="Percentage of\nadequate stool") +
  scale_color_manual(values = c("#F1BB7B", "#FD6467", "#5B1A18","#F1BB7B", "#FD6467", "#5B1A18")) +
  scale_fill_manual(values = c("#F1BB7B", "#FD6467", "#5B1A18","#F1BB7B", "#FD6467", "#5B1A18")) +
  guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "gray"),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

stool_legend <- cowplot::get_legend(p1a)
# check that this corresponds with the data
dat2$stool_start_cat <- cut(dat2$stool_start,breaks=c(0,70,80,90,100),levels=c(70,80,90))
p1b <- p1a + geom_point(data=dat2[dat2$classification=="afp",], aes(x=npafp_first2+2, y=nc_adjusted+6,
                                                                    group=as.factor(stool_start_cat), colour=as.factor(stool_start_cat))) + 
  xlim(0,15) +
  scale_color_manual(values = c("#F1BB7B", "#FD6467", "#5B1A18","#F1BB7B", "#FD6467", "#5B1A18")) +
  scale_fill_manual(values = c("#F1BB7B", "#FD6467", "#5B1A18","#F1BB7B", "#FD6467", "#5B1A18")) +
  theme(legend.position = "none") +
  geom_hline(yintercept=6, linetype="dashed", color = "black") 

x <- ggdraw(plot_grid(plot_grid(p1b, ncol=1),
                      plot_grid(stool_legend, NULL, ncol=3),
                      rel_widths=c(1, 0.2)))

# figure 1b - difference between ES and AFP
new_datab <- expand.grid(npafp_first2=seq(-1,15,1),stool_start=c(80),
                         classification=c("afp","es"))  # assume all are afp ##

new_datab$nc_adjusted <- NA
new_datab$nc_adjusted <- predict(nbmod2a,newdata=new_datab,type=c("response"))
tmp <-  predict(nbmod2a,newdata=new_datab,type=c("response"),se.fit=T)
new_datab$nc_adjusted_lwr <- tmp$fit-1.96*tmp$se.fit
new_datab$nc_adjusted_upr <- tmp$fit+1.96*tmp$se.fit

p2 <- ggplot(new_datab,aes(x=npafp_first2+2, y=nc_adjusted+6, group=as.factor(classification),
                           colour=as.factor(classification))) + geom_line(size=0.5)

p2a <- p2 + geom_smooth(alpha=0.3,aes(ymin=new_datab$nc_adjusted_lwr+6, ymax=new_datab$nc_adjusted_upr+6,
                                      fill=as.factor(classification), colour=as.factor(classification)), stat = "identity") + 
  coord_cartesian(ylim = c(2,60), xlim=c(2,15)) + 
  labs(x="Non-polio AFP rate (per 100,000 children <15 years)", y="Nucleotide differences",
       colour="Surveillance\nclassification") +
  scale_color_manual(values = c("#FD6467", "#5B1A18", "#FD6467", "#5B1A18")) +
  scale_fill_manual(values = c("#FD6467", "#5B1A18","#FD6467", "#5B1A18")) +
  guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "gray"),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

class_legend <- cowplot::get_legend(p2a)
# check that this corresponds with the data
p2b <- p2a + geom_point(data=dat2, 
                        aes(x=npafp_first2+2, y=nc_adjusted+6,
                            group=as.factor(classification), colour=as.factor(classification))) + 
  xlim(0,15) +
  scale_color_manual(values = c("#FD6467", "#5B1A18","#FD6467", "#5B1A18")) +
  scale_fill_manual(values = c("#FD6467", "#5B1A18","#FD6467", "#5B1A18")) +
  theme(legend.position = "none") + 
  geom_hline(yintercept=6, linetype="dashed", color = "black") 

y <- ggdraw(plot_grid(plot_grid(p2b, ncol=1),
                      plot_grid(class_legend, NULL, ncol=3),
                      rel_widths=c(1, 0.2)))

title <- ggdraw() + draw_label("        Serotype 2: predicted nucleotide differences", fontface = 'bold', x = 0, hjust = 0)
xy <- ggdraw(plot_grid(plot_grid(x, y, ncol=2))) + draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), size=15) 

plot_grid(title, xy, ncol=1, rel_heights = c(0.1, 1))


### distribution of detection times plot (figure 3) ###
type13 <- dat13
type13 <- type13[order(type13$days_to_detect), ]
type13$order <- c(1:15)

iso_3_code <- paste0(type13$iso_3_code,"-",type13$year_start)
sp <- type13 %>% ggplot(aes(y=days_to_detect, x=order,  
                            ymin=days_detect_min, ymax=days_detect_max)) +
  geom_pointrange(aes(colour=who_region)) +
  scale_x_continuous(breaks=c(1:15), labels=iso_3_code)

sp2 <- sp + scale_colour_manual(values = wes_palette("Darjeeling1",5, type = "continuous")) +
  labs(title="Estimates of time to VDPV Detection: serotypes 1 & 3", 
       x="", 
       y="Time to detection (days)", colour="WHO region") + ylim(0,2700) 

f1 <- sp2 + theme(
  axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
  axis.line = element_line(colour = "gray"),
  panel.background = element_blank(), 
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) + 
  geom_hline(yintercept=365, linetype="dashed", color = "blue") +
  geom_hline(yintercept=1461, linetype="dashed", color = "red")

# % detected within a year
sum(type13$days_detect_max<365)/nrow(type13)

type2 <- dat2
type2 <- type2[order(type2$days_to_detect), ]
type2$order <- c(1:59)
iso_3_code2 <- paste0(type2$iso_3_code,"-",type2$year_start)

sp <- type2 %>% ggplot(aes(y=days_to_detect, x=order, 
                           ymin=days_detect_min, ymax=days_detect_max)) +
  geom_pointrange(aes(colour=who_region))  +
  scale_x_continuous(breaks=c(1:59), labels=iso_3_code2)

f2 <- sp + scale_colour_manual(values = c("#FF0000","#00A08A","#F98400","#5BBCD6")) +
  labs(title="Estimates of time to VDPV Detection: serotype 2", 
       x="", 
       y="Time to detection (days)", colour="WHO region") + ylim(0,2700) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_line(colour = "gray"),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  geom_hline(yintercept=365, linetype="dashed", color = "blue") +
  geom_hline(yintercept=1461, linetype="dashed", color = "red")

# % detected within a year
sum(type2$days_detect_max<365)/nrow(type2)

# save to file
library(cowplot)
pdf("DetectionEstimates_09Feb2021.pdf",height=8,width=10)
plot_grid(f2, f1, labels = c('A', 'B'),
          ncol = 1,label_size = 10)
dev.off()

#end