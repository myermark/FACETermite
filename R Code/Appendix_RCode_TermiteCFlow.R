####
# Termite_CFlow.R  
####
# Purpose: To consolidate Picarro output into a single data object and analyze termite carbon flow experiments
####
# Written by: Mark Myer
# R Version: 3.5.1 Color Spray
# Date: 8/20/20
####
library(plyr)
library(tidyverse)
library(readr)
library(lubridate)
library(ggplot2)
library(patchwork)

#Gas data-----
#Import gas data
dat_522 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190522_151440.csv")[,-c(26:30)] 
dat_523 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190523_134634.csv")[,-c(26:30)] 
dat_524 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190524_133523.csv")[,-c(26:30)] 
dat_526 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190526_211418.csv")[,-c(26:30)] 
dat_527 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190527_132826.csv")[,-c(26:30)] 
dat_528 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190528_134031.csv")[,-c(26:30)] 
dat_530 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190530_135716.csv")[,-c(26:30)] 
dat_531 <- read.csv("gasdatasheets/CFIDS2066_SSIM_20190531_130450.csv")[,-c(26:30)] 
gaskey <- read.csv("May2019_GasAnalysis_CleanedForR.csv")

#Rowbind the data into one dataframe
gas_merged <- rbind(dat_522, dat_523, dat_524, dat_526, dat_527, dat_528, dat_530, dat_531)

rm(dat_522, dat_523, dat_524, dat_526, dat_527, dat_528, dat_530, dat_531)

write.csv(gas_merged, file="carb_gasdata_merged.csv")

#Format the date to get the day
gas_merged$Date.Time <- mdy_hms(gas_merged$Date.Time, tz="EST")

#Create a merge code to put the two files together
gas_merged$merge_code <- as.numeric(paste0(day(gas_merged$Date.Time), gas_merged$Run.Num))                                
gaskey$merge_code <- as.numeric(paste0(gaskey$Day, gaskey$Inj..))

#Merge the datasets by merge code
gas_merged <- merge(gas_merged, gaskey, by="merge_code")

#Separate out the 100ppm CH4 standards to paste on later
ch4_stds <- filter(gas_merged, Sample. == "CH4 std")

#Keep the highest re-run data from each sample for CO2
gas_merged <- arrange(gas_merged, Sample., X12CO2.Mean)

gas_merged <- gas_merged[!duplicated(gas_merged[,"Sample."], fromLast=T),] %>%

#Make a run-vessel identifier
mutate(run_v = paste0("R",Run.,"V",Vessel)) %>%    
  
#Remove the 100ppm CH4 standards
filter(Sample. != "CH4 std") %>%

#Sort by Treatment, run_v, then time..hr.
arrange(Treatment, run_v, Time..hr.) 

#Stick the CH4 standards back on the bottom
ch4_stds <- mutate(ch4_stds, run_v = paste0("R",Run.,"V",Vessel))
gas_merged <- rbind(gas_merged, ch4_stds) 

#Replace all CO2 values with "-1" if 12 CO2 mean is less than 100
gas_merged <- mutate(gas_merged, X12CO2.Mean = replace(X12CO2.Mean, X12CO2.Mean <100, -1)) %>%
mutate_if(names(gas_merged) %in% (names(select(gas_merged, contains("CO2")))), funs(ifelse(X12CO2.Mean ==-1, -1, .))) %>%
  
#Replace all High Precision CH4 values with "-1" if 12 CH4 mean is less than 1.2 and with "99999" if greater than 15
mutate(HP.12CH4.Mean = replace(HP.12CH4.Mean, HP.12CH4.Mean <1.2, -1)) %>%
mutate_if(names(gas_merged) %in% (names(select(gas_merged, contains("HP")))), funs(ifelse(HP.12CH4.Mean ==-1, -1, .))) %>%
  
mutate(HP.12CH4.Mean = replace(HP.12CH4.Mean, HP.12CH4.Mean >15 , 99999)) %>%
mutate_if(names(gas_merged) %in% (names(select(gas_merged, contains("HP")))), funs(ifelse(HP.12CH4.Mean ==99999, 99999, .))) 
  
#Specify syringe*vial and machine dilution factors, add columns for total CO2 and CH4
syr.dil = 8
mach.dil = 10
gas_merged = gas_merged %>% mutate(Total_CO2 = rep(NA, times=nrow(gas_merged)), Total_CH4 = rep(NA, times=nrow(gas_merged))) %>% 
  
#Get total CH4 for CH4 standards only
mutate(Total_CH4 = replace(Total_CH4, Sample. == "CH4 std", tail(HR.12CH4.Mean * mach.dil,sum(gas_merged$Sample.=="CH4 std")))) %>%

#Get total CO2 for nonstandards and non-BD
mutate(Total_CO2 = replace(Total_CO2, Sample. == "CH4 std" | is.na(Total_CO2),-1)) 
gas_merged$Total_CO2[which(gas_merged$Sample. != "CH4 std" & gas_merged$X12CO2.Mean != -1)] = gas_merged$X12CO2.Mean[which(gas_merged$Sample. != "CH4 std" & gas_merged$X12CO2.Mean != -1)] * syr.dil * mach.dil + gas_merged$X13CO2.Mean[which(gas_merged$Sample. != "CH4 std" & gas_merged$X12CO2.Mean != -1)] * syr.dil * mach.dil

#Import corrected CH4 data and use it instead of the CH4 data in gas_merged
ch4_merged <- read.csv("FACE_CH4_Cleaned.csv") %>%

#Make a run-vessel identifier
mutate(run_v = paste0("R",Run.,"V",Vessel))

#Get total CH4 for termite treatments
ch4_merged$Total_CH4[which(ch4_merged$Treatment %in% c("F+T", "L+T"))] <- ch4_merged$HR.12CH4.Mean[which(ch4_merged$Treatment %in% c("F+T", "L+T"))] * syr.dil * mach.dil + ch4_merged$HR.13CH4.Mean[which(ch4_merged$Treatment %in% c("F+T", "L+T"))] * syr.dil * mach.dil

#Get total CH4 for nontermite treatments
ch4_merged$Total_CH4[which(ch4_merged$Treatment %in% c("F", "L"))] <- ch4_merged$HP.12CH4.Mean[which(ch4_merged$Treatment %in% c("F", "L"))] * syr.dil * mach.dil + ch4_merged$HP.13CH4.Mean[which(ch4_merged$Treatment %in% c("F", "L"))] * syr.dil * mach.dil

#Plot total CO2 for termite treatments
tiff(filename="./Final Figures/FigureD1.tiff", res=300, width=6, height=8, units="in", pointsize=12, compression="lzw", type="cairo")
#tiff(filename="totalCO2.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
FigD1A <- ggplot(data=gas_merged[which(gas_merged$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.))+24, y=Total_CO2, group=run_v, color=Treatment, pch=Colony)) +
            geom_line()+
            geom_point() +
            scale_shape_manual(values=c(16, 17, 15, 3, 18, 8)) +
            scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112)+24) +
            scale_color_discrete(labels=c("FACE","Ambient")) +
            labs(title=expression("Colony Vessel Total CO"[2]), x="Time (Hours)", y=expression("Total CO"[2]~"(PPM)"), tag = "A") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"),
                           plot.title = element_text(hjust = 0.5))
#dev.off()

#Plot total CH4 for termite treatments
#tiff(filename="totalCH4.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
FigD1B<- ggplot(data=ch4_merged[which(ch4_merged$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.))+24, y=Total_CH4, group=run_v, color=Treatment, pch=Colony)) +
          geom_line()+
          geom_point() +
          scale_shape_manual(values=c(16, 17, 15, 3, 18, 8)) +
          scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112)+24) +
          scale_color_discrete(labels=c("FACE","Ambient")) +
          labs(title=expression("Colony Vessel Total CH"[4]), x="Time (Hours)", y=expression("Total CH"[4]~"(PPM)"), tag = "B") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5))

FigD1A / FigD1B

dev.off()


#Plot delta CO2 for termite treatments
tiff(filename="./Final Figures/Figure3.tiff", res=300, width=6, height=8, units="in", pointsize=12, compression="lzw", type="cairo")
#tiff(filename="deltaCO2.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
Fig3a <-ggplot(data=gas_merged[which(gas_merged$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.)) + 24, y=Delta.13CO2.Mean, group=run_v, color=Treatment, pch=Colony)) +
          geom_line()+
          geom_point() +
          scale_shape_manual(values=c(16, 17, 15, 3, 18, 8)) +
          scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112) + 24) +
          scale_color_discrete(labels=c("FACE","Ambient")) +
          labs(title=expression("Colony Vessel CO"[2]~delta^13~"C"), x="Time (Hours)", y=expression("CO"[2]~delta^13~"C (\u2030)"), tag = "A") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5)) 
#dev.off()

#Plot delta CH4 for termite treatments
#tiff(filename="deltaCH4.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
Fig3b <- ggplot(data=ch4_merged[!is.na(ch4_merged$HR.Delta.iCH4.Mean) & which(ch4_merged$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.)) + 24, y=HR.Delta.iCH4.Mean, group=run_v, color=Treatment, pch=Colony)) +
          geom_line()+
          geom_point() +
          scale_shape_manual(values=c(17, 15, 3, 18, 8)) +
          scale_x_continuous(breaks=c(16, 24, 40, 48, 64, 72, 88, 96, 112) + 24, limits = c(24, 136)) +
          scale_color_discrete(labels=c("FACE","Ambient")) +
          labs(title=expression("Colony Vessel CH"[4]~delta^13~"C"), x="Time (Hours)", y=expression("CH"[4]~delta^13~"C (\u2030)"), tag = "B") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5))

Fig3a / Fig3b

dev.off()

#Solids data-----
#Import solids data
solidkey <- read.csv("SolidsKey.csv")
sol_dat <- read.csv("May2019_isoC_2plates.csv")

#Merge by sample number
sol_merged <- merge(sol_dat, solidkey, by="Sample..")

#Save merged file
write.csv(sol_merged, file="solid_data_merged.csv")

#Load merged file
sol_merged <- read.csv("solid_data_merged.csv")

#Separate rows for assimilation project
ass_project <- sol_merged[which(is.na(sol_merged$Vessel)),]
sol_merged <- sol_merged[-which(is.na(sol_merged$Vessel)),]

#Calculate standard error 
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

sol_summaryd13 <- summarySE(sol_merged, measurevar="d13C", groupvars=c("Treatment","Sample.Type", "PrePost"))
sol_summaryPercC <- summarySE(sol_merged, measurevar="PercC", groupvars=c("Treatment","Sample.Type", "PrePost"))

#Create a grouped faceted barplot
#Change the levels of pre and post to make the plots appear in the right order
sol_summaryd13$PrePost <- relevel(sol_summaryd13$PrePost, "Pre") ; sol_summaryPercC$PrePost <- relevel(sol_summaryPercC$PrePost, "Pre")

#D13C
tiff(filename="D13C.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(sol_summaryd13, Treatment %in% c("Ambient", "FACE")), aes(fill=Treatment, x=factor(Sample.Type, level=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "ridge-like construction mat. ")), y=d13C))  + 
  #scale_y_reverse() + 
  coord_cartesian(ylim=c(-20,-40)) +
  geom_errorbar(aes(ymin=d13C-se, ymax=d13C+se),
                width=.2, # Width of the error bars
                position=position_dodge(.9)) +
  geom_bar(position="dodge", stat="identity") +
  geom_hline(yintercept = min(filter(sol_summaryd13, Treatment == "Ambient" & Sample.Type == "wood")$d13C), linetype = "dashed") +
  scale_x_discrete(breaks=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "ridge-like construction mat. "),
                   labels=c("WOD","BOD", "GUT", "CON", "FRA", "RID")) +
  facet_grid(.~PrePost, scales="free", space ="free") +
  labs(title=expression("Average"~delta^13~"C in Termite Vessel Substrates"), x="", y=expression(delta^13~"C (\u2030)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#PercC
tiff(filename="PercC.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(sol_summaryPercC, Treatment %in% c("Ambient", "FACE")), aes(fill=Treatment, x=factor(Sample.Type, level=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "ridge-like construction mat. ")), y=PercC))  + 
  coord_cartesian(ylim=c(40, 55)) +
  geom_errorbar(aes(ymin=PercC-se, ymax=PercC+se),
                width=.2, # Width of the error bars
                position=position_dodge(.9)) +
  geom_bar(position="dodge", stat="identity") +
  scale_x_discrete(breaks=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "ridge-like construction mat. "),
                   labels=c("WOD","BOD", "GUT", "CON", "FRA", "RID")) +
  facet_grid(.~PrePost, scales="free", space ="free") +
  labs(title=expression("Average %C in Termite Vessel Substrates"), x="", y="%C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#Make assimilation project graphs
ass_project$PrePost[which(ass_project$PrePost=="")] <- "Post"
#Change the levels of pre and post to make the plots appear in the right order
ass_project$PrePost <- relevel(ass_project$PrePost, "Pre")
ass_summaryd13 <- summarySE(ass_project, measurevar="d13C", groupvars=c("Treatment","Sample.Type", "PrePost"))
ass_summaryPercC <- summarySE(ass_project, measurevar="PercC", groupvars=c("Treatment ","Sample.Type", "PrePost"))

#Create a grouped faceted barplot
#D13C
tiff(filename="ass_D13C.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=ass_summaryd13, aes(fill=Treatment, x=factor(Sample.Type, level=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "plaster-like construction mat. ", "ridge-like construction mat. ")), y=d13C))  + 
  #scale_y_reverse() + 
  coord_cartesian(ylim=c(-20,-40)) +
  geom_errorbar(aes(ymin=d13C-se, ymax=d13C+se),
                width=.2, # Width of the error bars
                position=position_dodge(.9)) +
  geom_bar(position="dodge", stat="identity") +
  geom_hline(yintercept = min(filter(ass_summaryd13, Treatment == "Ambient" & Sample.Type == "wood")$d13C), linetype = "dashed") +
  scale_x_discrete(breaks=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "plaster-like construction mat. ", "ridge-like construction mat. "),
                   labels=c("WOD","BOD", "GUT", "CON", "FRA", "PLA", "RID")) +
  facet_grid(.~PrePost, scales="free", space ="free") +
  labs(title=expression("Average"~delta^13~"C in Assimilation Boxes"), x="", y=expression(delta^13~"C (\u2030)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#PercC
tiff(filename="ass_PercC.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=ass_summaryPercC, aes(fill=Treatment, x=factor(Sample.Type, level=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "plaster-like construction mat. ", "ridge-like construction mat. ")), y=PercC))  + 
  coord_cartesian(ylim=c(40, 65)) +
  geom_errorbar(aes(ymin=PercC-se, ymax=PercC+se),
                width=.2, # Width of the error bars
                position=position_dodge(.9)) +
  geom_bar(position="dodge", stat="identity") +
  scale_x_discrete(breaks=c("wood", "de-gutted bodies", "guts", "construction mat. off  wood",  "frass", "plaster-like construction mat. ", "ridge-like construction mat. "),
                   labels=c("WOD","BOD", "GUT", "CON", "FRA", "PLA", "RID")) +
  facet_grid(.~PrePost, scales="free", space ="free") +
  labs(title=expression("Average %C in Assimilation Boxes"), x="", y="%C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#Assimilation rate line graphs for bodies and guts
ass_bods_summaryd13 <- summarySE(ass_project, measurevar="d13C", groupvars=c("Treatment","Sample.Type", "Ass.day"))
ass_bods_summaryPercC <- summarySE(ass_project, measurevar="PercC", groupvars=c("Treatment ","Sample.Type", "Ass.day"))

#Plot delta13 C for assimilation bodies
#tiff(filename="ass_bodies_D13C.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(ass_bods_summaryd13, Treatment %in% c("Ambient", "FACE") & (Sample.Type %in% c( "de-gutted bodies"))), aes(x=Ass.day, y=d13C, color=Treatment)) +
  geom_line()+
  geom_point() +
  scale_color_discrete(labels=c("FACE","Ambient")) +
  labs(title=expression("Assimilation For De-Gutted Bodies"~delta^13~"C"), x="Time (Days)", y=expression(delta^13~"C (\u2030)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
#dev.off()

#Plot delta13 C for assimilation guts
#tiff(filename="ass_guts_D13C.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(ass_bods_summaryd13, Treatment %in% c("Ambient", "FACE") & (Sample.Type %in% c("guts"))), aes(x=Ass.day, y=d13C, color=Treatment)) +
  geom_line()+
  geom_point() +
  scale_color_discrete(labels=c("FACE","Ambient")) +
  labs(title=expression("Assimilation For Termite Guts"~delta^13~"C"), x="Time (Days)", y=expression(delta^13~"C (\u2030)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
#dev.off()



#Plot percent carbon for assimilation bodies
#Plot percent carbon for assimilation bodies
#tiff(filename="ass_bodies_PercC.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(ass_bods_summaryPercC, Treatment %in% c("Ambient", "FACE") & (Sample.Type %in% c( "de-gutted bodies"))), aes(x=Ass.day, y=PercC, color=Treatment)) +
  geom_line()+
  geom_point() +
  scale_color_discrete(labels=c("FACE","Ambient")) +
  labs(title=expression("Assimilation For De-Gutted Bodies %C"), x="Time (Days)", y="% C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
#dev.off()

#Plot percent carbon for assimilation guts
#tiff(filename="ass_guts_PercC.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(ass_bods_summaryPercC, Treatment %in% c("Ambient", "FACE") & (Sample.Type %in% c("guts"))), aes(x=Ass.day, y=PercC, color=Treatment)) +
  geom_line()+
  geom_point() +
  scale_color_discrete(labels=c("FACE","Ambient")) +
  labs(title=expression("Assimilation For Termite Guts %C"), x="Time (Days)", y="% C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
#dev.off()


#Carbon mass balance data-----
#Import dataset
mass_data <- select(read.csv("FACE_SetUp_data.csv") , -c("Colony"))
#Calculate moles of gas in each container (n = PV/RT)
R_gasconst = 0.08206 #ideal gas constant in L*atm*mol^-1
T_kelvin = 300.15 #temperature of hot room in Kelvin
mutate(mass_data, n_gas = (Pressure_atm*(VesselVol/1000)/(R_gasconst * T_kelvin)), 
                  run_v = paste0("R",Run.,"V",Vessel)) -> mass_data
mass_data <- merge(mass_data, select(gas_merged, -c("Vessel")), by="run_v") 


#Determine fraction of C13 in CO2 and CH4
mutate(mass_data, frac_C13_CO2 = (X13CO2.Mean * mach.dil * syr.dil)/Total_CO2, frac_C13_CH4 = (HR.13CH4.Mean * mach.dil * syr.dil)/Total_CH4) %>%
  mutate(frac_C12_CO2 = 1-frac_C13_CO2, frac_C12_CH4 = 1-frac_C13_CH4) -> mass_data

#Determine mole fraction of CO2 and CH4 in each vessel, accounting for C13 being heavier than C12
mutate(mass_data, mol_frac_CO2 = Total_CO2/1000000, mol_frac_CH4 = Total_CH4/1000000) %>%
  mutate(mol_CO2 = mol_frac_CO2 * n_gas, mol_CH4 = mol_frac_CH4 * n_gas) %>%
  mutate(g_CO2 = (frac_C12_CO2 * mol_CO2 * 44.095) + (frac_C13_CO2 * mol_CO2 * 45.102), #Masses of 1 mole CO2 with C12 and C13 respectively
         g_CH4 = (frac_C12_CH4 * mol_CH4 * 16.04) + (frac_C13_CH4 * mol_CH4 * 17.047)) -> mass_data #Masses of 1 mole CH4 with C12 and C13 respectively

#Plot mass CO2 for termite treatments
tiff(filename="massCO2.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=mass_data[which(mass_data$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.)), y=g_CO2, group=run_v, color=Treatment, pch=Colony)) +
  geom_line()+
  geom_point() +
  scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112)) +
  scale_color_discrete(labels=c("FACE","Ambient")) +
  labs(title=expression("Colony Vessel Mass CO"[2]), x="Time (Hours)", y=expression("Mass CO"[2]~"(g)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#Perform CH4 calculations
mass_data2 <- select(read.csv("FACE_SetUp_data.csv") , -c("Colony"))
#Calculate moles of gas in each container (n = PV/RT)
mutate(mass_data2, n_gas = (Pressure_atm*(VesselVol/1000)/(R_gasconst * T_kelvin)), 
       run_v = paste0("R",Run.,"V",Vessel)) -> mass_data2
mass_data2 <- merge(mass_data2, select(ch4_merged, -c("Vessel")), by="run_v")


#Determine fraction of C13 in CH4
mutate(mass_data2, frac_C13_CH4 = (HR.13CH4.Mean * mach.dil * syr.dil)/Total_CH4) %>%
  mutate(frac_C12_CH4 = 1-frac_C13_CH4) -> mass_data2

#Determine mole fraction of CO2 and CH4 in each vessel, accounting for C13 being heavier than C12
mutate(mass_data2, mol_frac_CH4 = Total_CH4/1000000) %>%
  mutate(mol_CH4 = mol_frac_CH4 * n_gas) %>%
  mutate(g_CH4 = (frac_C12_CH4 * mol_CH4 * 16.04) + (frac_C13_CH4 * mol_CH4 * 17.047)) -> mass_data2 #Masses of 1 mole CH4 with C12 and C13 respectively

#Plot mass CH4 for termite treatments
tiff(filename="massCH4.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=mass_data2[which(mass_data2$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.)), y=g_CH4, group=run_v, color=Treatment, pch=Colony)) +
  geom_line()+
  geom_point() +
  scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112)) +
  scale_color_discrete(labels=c("FACE","Ambient")) +
  labs(title=expression("Colony Vessel Mass CH"[4]), x="Time (Hours)", y=expression("Mass CH"[4]~"(g)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#Plot a mean trendline for each gas by treatment
gas_means <- group_by(mass_data, Treatment, Time..hr.) %>%
              summarise(g_CO2_mean = mean(g_CO2, na.rm=T),
                        g_CO2_sd = sd(g_CO2, na.rm=T)) %>%
              filter(Treatment %in% c("F+T", "L+T"))

gas_means2 <- group_by(mass_data2, Treatment, Time..hr.) %>%
  summarise(g_CH4_mean = mean(g_CH4, na.rm=T),
            g_CH4_sd = sd(g_CH4, na.rm=T))


#Plot mean mass by treatment CO2 for termite treatments
tiff(filename="./Final Figures/Figure 2.tiff", res=300, width=6, height=8, units="in", pointsize=12, compression="lzw", type="cairo")
#tiff(filename="meanbytrt_massCO2.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
Fig2a <- ggplot(data=gas_means[which(gas_means$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.)) + 24, y=g_CO2_mean *1000)) +
          geom_line(aes(color=Treatment)) +
          geom_point(aes(color=Treatment)) +
          stat_smooth(method = "lm", linetype="dashed", lwd=0.5) +
          geom_errorbar(aes(ymin=(g_CO2_mean-g_CO2_sd) *1000, ymax=(g_CO2_mean+g_CO2_sd) * 1000, color=Treatment),
                        width=3, # Width of the error bar crossbars
                        position=position_dodge(.9)) +
          scale_color_discrete(labels=c("FACE","Ambient")) +
          scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112) + 24) +
          labs(title=expression("Mean CO"[2]~"By Treatment"), x="Time (Hours)", y=expression("Mass CO"[2]~"(mg)"), tag = "A") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5)) + 
          annotate("text", x = 70, y = max(gas_means[which(gas_means$Treatment %in% c("F+T", "L+T")),]$g_CO2_mean *1000) * 0.35, label = expression(paste(bold("Adjusted R"^"2")," = 0.978")), hjust = 0) +
          annotate("text", x = 70, y = max(gas_means[which(gas_means$Treatment %in% c("F+T", "L+T")),]$g_CO2_mean *1000) * 0.28, label = expression(paste(bold("Accumulation Rate"), " =")), hjust = 0) +
          annotate("text", x = 90, y = max(gas_means[which(gas_means$Treatment %in% c("F+T", "L+T")),]$g_CO2_mean *1000) * 0.21, label = expression(paste("3.571 mg CO"[2]~"/ hour")), hjust = 0) +
          annotate("text", x = 70, y = max(gas_means[which(gas_means$Treatment %in% c("F+T", "L+T")),]$g_CO2_mean *1000) * 0.14, label = expression(paste(bold("Intercept")," = 0.1218 mg CO"[2])), hjust = 0) 
#dev.off()

#Plot mean mass by treatment CH4 for termite treatments
#tiff(filename="meanbytrt_massCH4.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
Fig2b <- ggplot(data=gas_means2[which(gas_means2$Treatment %in% c("F+T", "L+T")),], aes(x=as.numeric(as.character(Time..hr.)) + 24, y=g_CH4_mean * 1000)) +
          geom_line(aes(color=Treatment)) +
          geom_point(aes(color=Treatment)) +
          stat_smooth(method = "lm", linetype="dashed", lwd=0.5) +
          geom_errorbar(aes(ymin=(g_CH4_mean-g_CH4_sd) *1000 , ymax=(g_CH4_mean+g_CH4_sd) * 1000, color=Treatment),
                        width=3, # Width of the error bar crossbars
                        position=position_dodge(.9)) +
          scale_color_discrete(labels=c("FACE","Ambient")) +
          scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112) + 24) +
          labs(title=expression("Mean CH"[4]~"By Treatment"), x="Time (Hours)", y=expression("Mass CH"[4]~"(mg)"), tag = "B") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5)) + 
          annotate("text", x = 70, y =  max(gas_means2[which(gas_means2$Treatment %in% c("F+T", "L+T")),]$g_CH4_mean *1000) * 0.35, label = expression(paste(bold("Adjusted R"^"2")," = 0.953")), hjust = 0) +
          annotate("text", x = 70, y =  max(gas_means2[which(gas_means2$Treatment %in% c("F+T", "L+T")),]$g_CH4_mean *1000) * 0.28, label = expression(paste(bold("Accumulation Rate"), " =")), hjust = 0) +
          annotate("text", x = 90, y =  max(gas_means2[which(gas_means2$Treatment %in% c("F+T", "L+T")),]$g_CH4_mean *1000) * 0.21, label = expression(paste("0.01598 mg CH"[4]~"/ hour")), hjust = 0) +
          annotate("text", x = 70, y =  max(gas_means2[which(gas_means2$Treatment %in% c("F+T", "L+T")),]$g_CH4_mean *1000) * 0.14, label = expression(paste(bold("Intercept")," = 7.079 x 10"^-4~" mg CH"[4])), hjust = 0) 

Fig2a / Fig2b

dev.off()

#Fit a linear model to mean CO2 and CH4 by treatment to determine flux rates
CO2_lm <- lm(g_CO2_mean ~ as.numeric(as.character(Time..hr.)), data= gas_means[which(gas_means$Treatment %in% c("F+T", "L+T")),])
CH4_lm <- lm(g_CH4_mean ~ as.numeric(as.character(Time..hr.)), data= gas_means2[which(gas_means$Treatment %in% c("F+T", "L+T")),])

summary(CO2_lm)
summary(CH4_lm)

#Calculate %C from wood in the solid samples 
#Gather the data needed
setup_data <- read.csv("FACE_SetUp_data.csv") %>%
  mutate(run_v = paste0("R",Run.,"V",Vessel))    
sol_merged <- mutate(sol_merged, run_v = paste0("R",Run.,"V",Vessel)) %>%
  select(-c("Colony"))


sol_merged <- merge(unique(select(setup_data, c("Colony", "run_v"))), sol_merged, by="run_v")

fracC_data <- data.frame("Colony"=factor(), "Sample.Type"=factor(), "PerWoodC"=numeric())

#Calculate the %C from wood by colony, and sample type
for (i in 2:length(unique(sol_merged$Colony))) {
 for (j in 1:length(unique(droplevels(sol_merged$Sample.Type)))) {
  col = levels(sol_merged$Colony)[i]
  sam = levels(droplevels(sol_merged$Sample.Type))[j]
  wood_vals <- subset(sol_merged, sol_merged$Colony == col & sol_merged$Sample.Type == "wood" & sol_merged$PrePost == "Pre")
  frac_sub <- subset(sol_merged, sol_merged$Colony == col & sol_merged$Sample.Type == sam & sol_merged$PrePost == "Post")
  fracC = (frac_sub[which(frac_sub$Treatment=="FACE"),]$d13C - frac_sub[which(frac_sub$Treatment=="Ambient"),]$d13C) / (wood_vals[which(wood_vals$Treatment=="FACE"),]$d13C - wood_vals[which(wood_vals$Treatment=="Ambient"),]$d13C)
  if (length(fracC) > 0) {
    temp = data.frame("Colony"=col, "Assimilation Day"= NA, "Sample.Type"=sam, "PerWoodC"=fracC*100)
    fracC_data <- rbind(fracC_data, temp)
  }
 }
}

#Calculate %C from wood in gas samples
gas_d13 <- group_by(mass_data, Treatment, Colony, Time..hr.) %>%
  summarise(d13_CO2_mean = mean(Delta.13CO2.Mean, na.rm=T)) %>%
  filter(Treatment %in% c("F+T","L+T"))

gas_d13_2 <- group_by(filter(mass_data2, !is.na(HR.Delta.iCH4.Mean)), Treatment, Colony, Time..hr.) %>%
  summarise(d13_CH4_mean = mean(HR.Delta.iCH4.Mean, na.rm=T))

#Ungroup the factors 
gas_d13 = gas_d13 %>% 
  ungroup() %>%
  mutate_if(is.factor,
            fct_explicit_na,
            na_level = "to_impute")

gas_d13_2 = gas_d13_2 %>% 
  ungroup() %>%
  mutate_if(is.factor,
            fct_explicit_na,
            na_level = "to_impute")

gas_d13_2$Time..hr. = factor(gas_d13_2$Time..hr.)

#Get wood d13 values by colony
gas_wood_vals <- filter(sol_merged, Sample.Type == "wood" & PrePost == "Pre")

gas_fracC_data <- data.frame("Colony"=factor(), "Hour" = factor(), "PerWoodCO2"= numeric())

#Calculate the %C from wood by colony and hour
for (i in 1:length(unique(gas_d13$Time..hr.))) {
  for (j in 1:length(unique(droplevels(gas_d13$Colony)))) {
    hr = levels(gas_d13$Time..hr.)[i]
    col = levels(droplevels(gas_d13$Colony))[j]
    wood_vals <- filter(gas_wood_vals, Colony == col)
    frac_sub <- filter(gas_d13, Time..hr. == hr & Colony == col)
    fracCO2 = (frac_sub[which(frac_sub$Treatment=="F+T"),]$d13_CO2_mean - frac_sub[which(frac_sub$Treatment=="L+T"),]$d13_CO2_mean) / (wood_vals[which(wood_vals$Treatment=="FACE"),]$d13C - wood_vals[which(wood_vals$Treatment=="Ambient"),]$d13C)
    if (length(fracCO2) > 0) {
      temp = data.frame("Colony" = col, "Hour"=hr, "PerWoodCO2"=fracCO2*100)
      gas_fracC_data <- rbind(gas_fracC_data, temp)
    }
    else {
      temp = data.frame("Colony" = col, "Hour"=hr, "PerWoodCO2"=NA)
      gas_fracC_data <- rbind(gas_fracC_data, temp)
    }
  }
}

gas_fracC_data_2 <- data.frame("Colony"=factor(), "Hour" = factor(), "PerWoodCH4"= numeric())

for (i in 1:length(unique(gas_d13_2$Time..hr.))) {
  for (j in 1:length(unique(droplevels(gas_d13_2$Colony)))) {
    hr = levels(gas_d13_2$Time..hr.)[i]
    col = levels(droplevels(gas_d13_2$Colony))[j]
    wood_vals <- filter(gas_wood_vals, Colony == col)
    frac_sub <- filter(gas_d13_2, Time..hr. == hr & Colony == col)
    fracCH4 = (frac_sub[which(frac_sub$Treatment=="F+T"),]$d13_CH4_mean - frac_sub[which(frac_sub$Treatment=="L+T"),]$d13_CH4_mean) / (wood_vals[which(wood_vals$Treatment=="FACE"),]$d13C - wood_vals[which(wood_vals$Treatment=="Ambient"),]$d13C)
    if (length(fracCH4) > 0) {
      temp = data.frame("Colony" = col, "Hour"=hr, "PerWoodCH4"=fracCH4*100)
      gas_fracC_data_2 <- rbind(gas_fracC_data_2, temp)
    }
    else {
      temp = data.frame("Colony" = col, "Hour"=hr, "PerWoodCH4"=NA)
      gas_fracC_data2 <- rbind(gas_fracC_data_2, temp)
    }
  }
}

gas_fracC_data <- merge(gas_fracC_data, gas_fracC_data_2, by=c("Colony", "Hour"), all.x=T)

#Save the percentage wood based carbon data
write.csv(gas_fracC_data, file="gas_PercWoodBasedC.csv")

#Create a dataframe of gas fracC data and mass to determine mass wood-derived C
#Get gas means by colony and treatment
gas_means_col <- group_by(filter(mass_data, Treatment %in% c("F+T", "L+T")), Colony, Time..hr.) %>%
  summarise(mg_CO2_mean = 1000* mean(g_CO2, na.rm=T),
            mg_CO2_sd = 1000* sd(g_CO2, na.rm=T)) 

gas_means2_col <- group_by(mass_data2, Colony, Time..hr.) %>%
  summarise(mg_CH4_mean = 1000* mean(g_CH4, na.rm=T),
            mg_CH4_sd = 1000* sd(g_CH4, na.rm=T))

gas_means_col_merged <- merge(gas_means_col, gas_means2_col, by=c("Colony", "Time..hr."))

#Save the mass of gas data 
write.csv(gas_means_col_merged, file="gas_MassC.csv")

#Average the values by colony 
gas_fracC_data <- group_by(gas_fracC_data, Hour) %>%
  summarise(PerWoodCO2 = mean(PerWoodCO2, na.rm=T), PerWoodCH4 = mean(PerWoodCH4, na.rm=T))

#Rearrange the data for plotting by ggplot2
gas_CO2 <- select(gas_fracC_data, Hour, PerWoodC=PerWoodCO2)
gas_CO2$Gas = rep("CO2", times=length(gas_fracC_data$Hour))
gas_CH4 <- select(gas_fracC_data, Hour, PerWoodC=PerWoodCH4)
gas_CH4$Gas = rep("CH4", times=length(gas_fracC_data$Hour))
gas_fracC_data <- rbind(gas_CO2, gas_CH4)

#Plot the percentage from wood for gas data 
tiff(filename="CO2_CFromWood.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(gas_fracC_data, Gas == "CO2"), aes(x=as.numeric(as.character(Hour)) + 24, y=PerWoodC)) +
  stat_smooth(method = "lm", linetype="dashed", lwd=0.5) +
  geom_line()+
  geom_point() +
  scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112) + 24) +
  labs(title=expression("Percent Carbon in CO"[2]~"From Wood"), x="Time (Hours)", y="% C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

CO2_fracC_lm <-  lm(PerWoodC~ as.numeric(as.character(Hour)), data= filter(gas_fracC_data, Gas == "CO2"))
summary(CO2_fracC_lm)

tiff(filename="CH4_CFromWood.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(gas_fracC_data, Gas == "CH4"), aes(x=as.numeric(as.character(Hour)) + 24, y=PerWoodC)) +
  stat_smooth(method = "lm", linetype="dashed", lwd=0.5) +
  geom_line()+
  geom_point() +
  scale_x_continuous(breaks=c(0, 16, 24, 40, 48, 64, 72, 88, 96, 112) + 24) +
  labs(title=expression("Percent Carbon in CH"[4]~"From Wood"), x="Time (Hours)", y="% C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

CH4_fracC_lm <-  lm(PerWoodC~ as.numeric(as.character(Hour)), data= filter(gas_fracC_data, Gas == "CH4"))
summary(CH4_fracC_lm)


#Calculate %C from wood in assimilation study samples
#Gather the data needed
ass_carb <- group_by(ass_project, Sample.Type, Ass.day, Treatment) %>%
  summarise(mean_d13C = mean(d13C, na.rm=T))
ass_carb$Ass.day <- as.factor(ass_carb$Ass.day)
ass_wood <- filter(ass_carb, Sample.Type == "wood" & Ass.day == "0" & Treatment == "FACE")$mean_d13C  - filter(ass_carb, Sample.Type == "wood" & Ass.day == "0" & Treatment == "Ambient")$mean_d13C

#Ungroup the factors 
ass_carb = ass_carb %>% 
  ungroup() %>%
  mutate_if(is.factor,
     fct_explicit_na,
     na_level = "to_impute")

#Calculate the %C from wood by sample type and day
for (i in 1:length(unique(ass_carb$Ass.day))) {
  for (j in 1:length(unique(droplevels(ass_carb$Sample.Type)))) {
    day = levels(ass_carb$Ass.day)[i]
    sam = levels(droplevels(ass_carb$Sample.Type))[j]
    ass_sub <- subset(ass_carb, ass_carb$Ass.day == day & ass_carb$Sample.Type == sam)
    fracC = (ass_sub[which(ass_sub$Treatment=="FACE"),]$mean_d13C - ass_sub[which(ass_sub$Treatment=="Ambient"),]$mean_d13C) / ass_wood
    if (length(fracC) > 0) {
      temp = data.frame("Colony" = NA, "Assimilation Day"= day, "Sample.Type"=sam, "PerWoodC"=fracC*100)
      fracC_data <- rbind(fracC_data, temp)
    }
  }
}

write.csv(fracC_data, file="PercentCarbon_FromWood.csv")

#Plot the assimilation in bodies and guts
tiff(filename="ass_CFromWood.tiff", res=300, width=6, height=4, units="in", pointsize=12, compression="lzw", type="cairo")
ggplot(data=filter(fracC_data,  is.na(Colony) & Sample.Type %in% c("guts", "de-gutted bodies")), aes(x=as.numeric(as.character(Assimilation.Day)), y=PerWoodC, color=Sample.Type)) +
  stat_smooth(method = "lm", linetype="dashed", lwd=0.5) +
  geom_line()+
  geom_point() +
  scale_colour_discrete(name="Sample Type") +
  scale_x_continuous(breaks=c(6, 12, 18)) +
  labs(title=expression("Percent Carbon From Wood"), x="Time (Days)", y="% C") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()

#Fit a line to each type and obtain the assimilation rate
gut_ass_lm <- lm(PerWoodC ~ as.numeric(Assimilation.Day), data=filter(fracC_data, Sample.Type =="de-gutted bodies"))
bod_ass_lm <- lm(PerWoodC ~ as.numeric(Assimilation.Day), data=filter(fracC_data, Sample.Type =="guts"))

summary(gut_ass_lm)
summary(bod_ass_lm)

#Calculate mass of wood-based carbon in each solid sample type by colony
#Define a helper function to mutate at rows satisfying a condition that we will need later
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
#Start with percent carbon from wood in each sample type by colony
solid_fracC_data <- subset(fracC_data[is.na(fracC_data$Assimilation.Day),], select = -Assimilation.Day)
#Get the mass of carbon in each sample
solid_massC <- subset(sol_merged[!is.na(sol_merged$Colony),], select = c(run_v, Colony, Wt.mg, PercC, Sample.Type, PrePost, Treatment))
#Merge them to add the percent carbon from wood
solid_massC_merged <- merge(solid_massC, solid_fracC_data, by=c("Colony", "Sample.Type")) %>%
  #Adjust the mass of the tissue samples 
  mutate_cond(Sample.Type %in% c("de-gutted bodies", "guts"), Wt.mg = Wt.mg * 75)  %>%
  #Fix percentages
  mutate(PercC= PercC/100, PerWoodC = PerWoodC/100)


#Calculate the frass on wood C mass
#Obtain the post-exposure wood mass from each vessel
post_WoodMass <- subset(subset(solid_massC_merged, select=c(run_v, Wt.mg, PercC, PerWoodC, Sample.Type, PrePost)) %>% filter(Sample.Type == "wood" & PrePost == "Post"), select=-c(Sample.Type, PrePost)) %>%
  mutate(wood.Wt.mg = Wt.mg, wood.PercC = PercC/100, wood.PerWoodC = PerWoodC/100) %>%
  subset(select=-c(Wt.mg, PercC, PerWoodC))
frass_WoodC <- subset(subset(solid_massC_merged, select=c(run_v, PerWoodC, Sample.Type, PrePost)) %>% filter(Sample.Type == "frass" & PrePost == "Post"), select=-c(Sample.Type, PrePost)) %>%
  mutate(frass.PerWoodC = PerWoodC/100) %>%
  subset(select=-c(PerWoodC))

frassWoodMass <- merge(post_WoodMass, frass_WoodC, by="run_v") 

#Save them and send
write.csv(frassWoodMass, file="frassOnWood.csv")
write.csv(solid_massC_merged, file="solidWoodBasedC.csv")

