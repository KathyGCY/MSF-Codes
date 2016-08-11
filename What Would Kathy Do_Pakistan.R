# MSF Protocal for Pakistan Data
# What would Kathy do?


##########################################################################################
# Install Packages
##########################################################################################

install.packages("lattice")
install.packages("pracma")
install.packages("devtools")
install.packages("clusterGeneration")
install.packages("RSNNS")
install.packages("nnet")
install.packages("matcov")
install.packages("plotly")
install.packages("NeuralNetTools")
install.packages("chorddiag")
install.packages("Hmisc")
install.packages("vegan")
install.packages("circlize")
install.packages("circular")
install.packages("base")
install.packages("scatterplot3d")
install.packages("ResourceSelection")
install.packages("car")
install.packages("corrplot")
install.packages("smacof")
install.packages("klaR")
install.packages("MASS")
install.packages("CCA")
install.packages("CCP")
install.packages("GGally")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("MVN")
install.packages("raster")
install.packages("mvnormtest")
install.packages("mclust")
install.packages("vcd")
install.packages("cluster")
install.packages("rpart")
install.packages("Hmisc")
install.packages("foreign")
install.packages("neuralnet")
install.packages("fpc")
install.packages("MTS")
install.packages("plotly")
install.packages("devtools")
install.packages("Rcpp")
install.packages("ggvis")
install.packages("googleVis")
install.packages("rCharts")
install.packages("plotly")
install.packages("dplyr")
install.packages("grDevices")
install.packages("bitops")
install.packages("tidyr")
install.packages("knitr")
install.packages("randomForest")
install.packages("polycor")
install.packages("e1071")
install.packages("nortest")
install.packages("jpeg")
install.packages("nnet")
install.packages("R.utils")
install.packages("mfp")
install.packages("pROC")
install.packages("pROC")


##########################################################################################
# Load Packages
##########################################################################################

library(nnet)
library(devtools)
library(lattice)
library(matcov)
library(plotly)
library(chorddiag)
library(Hmisc)
library(vegan)
library(circlize)
library(circular)
library(base)
library(scatterplot3d)
library(ResourceSelection)
library(NeuralNetTools)
library(car)
library(nnet)
library(NeuralNetTools)
library(RSNNS)
library(corrplot)
library(smacof)
library(klaR)
library(MASS)
library(CCA)
library(CCP)
library(GGally)
library(reshape2)
library(ggplot2)
library(MVN)
library(raster)
library(mvnormtest)
library(mclust)
library(vcd)
library(cluster)
library(rpart)
library(Hmisc)
library(foreign)
library(neuralnet)
library(fpc)
library(MTS)
library(plotly)
library(devtools)
library(Rcpp)
library(ggvis)
library(googleVis)
library(rCharts)
library(plotly)
library(dplyr)
library(grDevices)
library(bitops)
library(tidyr)
library(knitr)
library(randomForest)
library(polycor)
library(e1071)
library(nortest)
library(jpeg)
library(neuralnet)
library(pracma)
library(R.utils)
library(survival)
library(KMsurv)
library(nlme)
library(km.ci)
library(mfp)
library(pROC)
library(pROC)


##########################################################################################
# Import Data
##########################################################################################
lullaby <- read.csv("~/Desktop/MSF/Pakistan/PakistanPeshawarAnti_DATA_LABELS_2016-07-11_1645.csv")


##########################################################################################
# Cleand the data
##########################################################################################


##########################################################################################
# Compute Variables and Define Categoricals
##########################################################################################


##########################################################################################
# Explorative Data Analysis and Visualization
##########################################################################################

# Set Plotting Parameter
par(mfrow=c(1, 1))
par(mar = c(5, 4, 4, 2) + 0.1)

# Pie Chart of Refer
refer = as.factor(lullaby$Referred.from)
levels(refer)
refer.tbl=table(refer)
# refer
# 1. Peshawar    2. Hangu    3. Saada    9. Other 
# 
slices <- c(table(refer)[1], table(refer)[2], table(refer)[3], table(refer)[4])
lbls <- c(names(table(refer))[1], names(table(refer))[2], names(table(refer))[3], names(table(refer))[4])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels = lbls, col = c("#81D8D0", "lightgoldenrod1", "palegreen", "lightpink1"),
    main = "Pie Chart of Refered From")
legend("topright", inset = .05, title = "Refered From", 
       fill =  c("#81D8D0", "lightgoldenrod1", "palegreen", "lightpink1"), 
       c(names(table(refer))[1], names(table(refer))[2], names(table(refer))[3], names(table(refer))[4]), 
       horiz = FALSE)


# Pie Chart of Delivery Type
deliver.type = as.factor(lullaby$Type.of.delivery)
levels(deliver.type)
deliver.tbl = table(deliver.type)
# deliver.type
# 1. Vaginal 2. C-section 
# 131           60 
slices <- c(deliver.tbl[1], deliver.tbl[2])
lbls <- c(names(deliver.tbl)[1], names(deliver.tbl)[2])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels = lbls, col = c("#81D8D0", "lightpink1"),
    main = "Pie Chart of Delivery Type")
legend("topright", inset = .05, title = "Delivery Type", 
       fill =  c("#81D8D0", "lightpink1"), 
       c(names(deliver.tbl)[1], names(deliver.tbl)[2]),
       horiz = FALSE)


# Pie Chart of Presentation
presentation = as.factor(lullaby$Presentation)
levels(presentation)
presentation.tbl = table(presentation)
# presentation
#   1. Cephalic     2. Breech 3. Transverse    9. Unknown 
#           157            22             7             5 
slices <- c(presentation.tbl[1], presentation.tbl[2], presentation.tbl[3], presentation.tbl[4])
lbls <- c(names(presentation.tbl)[1], names(presentation.tbl)[2], names(presentation.tbl)[3],names(presentation.tbl)[4])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels  =  lbls, col = c("#81D8D0", "lightpink1", "lightgoldenrod1", "palegreen"),
    main = "Pie Chart of Presentation")
legend("topright", inset = .05, title = "Presentation", 
       fill =  c("#81D8D0", "lightpink1", "lightgoldenrod1", "palegreen"),
       c(names(presentation.tbl)[1], names(presentation.tbl)[2], names(presentation.tbl)[3],names(presentation.tbl)[4]),
       horiz = FALSE)


# Density of Gravity
par(mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(1, 3))
gravity = lullaby$Gravidity
gravity = na.omit(gravity)
hist(na.omit(gravity), col = "#81D8D0", main = "Histogram of Gravity")
d <- density(na.omit(gravity)) # returns the density data 
plot(d,col = "dodgerblue3", main = "Density of Gravity")
polygon(d, col="#81D8D0", border="dodgerblue3")
boxplot(gravity, notch = TRUE, col = "#81D8D0", main = "Boxplot of  gravity")
abline(h = min(gravity), col = "darkslategray3")
abline(h = max(gravity), col = "darkslategray3")
abline(h = median(gravity), col = "navy")
abline(h = quantile(gravity, c(0.25, 0.75)), col = "dodgerblue3")


# Density of Parity
par(mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(1, 3))
parity = lullaby$Parity
hist(na.omit(parity),col = "#81D8D0",main = "Histogram of Parity")
d <- density(na.omit(parity)) # returns the density data 
plot(d,col = "dodgerblue3",main = "Density of Parity")
polygon(d, col="#81D8D0", border="dodgerblue3")
boxplot(parity,notch = TRUE,col = "#81D8D0" ,main = "Boxplot of  Parity")
abline(h = min(parity), col = "darkslategray3")
abline(h = max(parity), col = "darkslategray3")
abline(h = median(parity), col = "navy")
abline(h = quantile(parity, c(0.25, 0.75)), col = "dodgerblue3")


# Density of Early neonatal death
par(mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(1, 3))
neonatal = as.numeric(lullaby$Early.neonatal.death)
hist(na.omit(neonatal),col = "#81D8D0",main = "Histogram of Neonatal")
d <- density(na.omit(neonatal)) # returns the density data 
plot(d,col = "dodgerblue3",main = "Density of Early.neonatal.death")
polygon(d, col="#81D8D0", border="dodgerblue3")
boxplot(na.omit(neonatal),notch = TRUE,col = "#81D8D0" ,main = "Boxplot of  Neonatal")
abline(h = min(neonatal), col = "darkslategray3")
abline(h = max(neonatal), col = "darkslategray3")
abline(h = median(neonatal), col = "navy")
abline(h = quantile(neonatal, c(0.25, 0.75)), col = "dodgerblue3")

par(mfrow=c(1, 1))
par(mar = c(5, 4, 4, 2) + 0.1)


# Pie Chart of Twin delivery
par(mfrow = c(1, 1))
twin = as.factor(lullaby$Twin.delivery)
levels(twin)
twin.tbl = table(twin)
# twin
#  No Yes  
# 168  23  
slices <- c(twin.tbl[1], twin.tbl[2])
lbls <- c(names(twin.tbl)[1], names(twin.tbl)[2])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices, labels = lbls,col = c("#81D8D0", "lightpink1"),
    main = "Pie Chart of Twin delivery")
legend("topright", inset = .05, title = "Twins", 
       fill =  c("#81D8D0", "lightpink1"), 
       c(names(twin.tbl)[1], names(twin.tbl)[2]),
       horiz = FALSE)


# Pie Chart of Pre-term birth
preterm.birth = as.factor(lullaby$Pre.term.birth)
levels(preterm.birth)
preterm.birth.tbl = table(preterm.birth)
# preterm.birth
#  No Yes 
# 148  43 
slices <- c(preterm.birth.tbl[1], preterm.birth.tbl[2])
lbls <- c(names(preterm.birth.tbl)[1], names(preterm.birth.tbl)[2])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels = lbls, col = c("#81D8D0", "lightpink1"),
    main = "Pie Chart of Pre-term birth")
legend("topright", inset = .05, title = "Pre-term birth", 
       fill =  c("#81D8D0", "lightpink1"), 
       c(names(preterm.birth.tbl)[1], names(preterm.birth.tbl)[2]),
       horiz = FALSE)


# Pie Chart of Mother signs of infection
infect.mother = as.factor(lullaby$Mother.signs.of.infection)
levels(infect.mother)
infect.mother.tbl = table(infect.mother)
names(infect.mother.tbl)[1]= "Unknown"
# infect.mother
#  No Yes Unknown
# 143  45 2 
slices <- c(infect.mother.tbl[1], infect.mother.tbl[2], infect.mother.tbl[3])
lbls <- c(names(infect.mother.tbl)[1], names(infect.mother.tbl)[2], names(infect.mother.tbl)[3])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels = lbls, col = c( "palegreen", "#81D8D0", "lightpink1"),
    main = "Pie Chart of Mother signs of infection")
legend("topright", inset = .05, title = "Mother signs of infection", 
       fill =  c("palegreen", "#81D8D0", "lightpink1"),
       c(names(infect.mother.tbl)[1], names(infect.mother.tbl)[2], names(infect.mother.tbl)[3]),
       horiz = FALSE)


# Pie Chart of PROM above 18h
PROM18 = as.factor(lullaby$PROM.above.18h)
levels(PROM18)
PROM18.tbl = table(PROM18)
PROM18.tbl[4] = PROM18.tbl[4] + PROM18.tbl[1]
PROM18.tbl = PROM18.tbl[-1]

# PROM18
#  No Yes Unknown
# 124  62 2 
slices <- c(PROM18.tbl[1], PROM18.tbl[2], PROM18.tbl[3])
lbls <- c(names(PROM18.tbl)[1], names(PROM18.tbl)[2], names(PROM18.tbl)[3])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels = lbls, col = c("#81D8D0", "lightpink1", "palegreen"),
    main = "Pie Chart of PROM above 18h")
legend("topright", inset = .05, title = "PROM above 18h", 
       fill = c("#81D8D0", "lightpink1", "palegreen"),
       c(names(PROM18.tbl)[1], names(PROM18.tbl)[2], names(PROM18.tbl)[3]),
       horiz = FALSE)


# Pei chart of ATB during labour
atb.lab = as.factor(lullaby$ATB.during.labour)
levels(atb.lab)
atb.lab.tbl = table(atb.lab)
# atb.lab
#  No Yes Unknown
# 128  62 1 
slices <- c(atb.lab.tbl[1], atb.lab.tbl[2], atb.lab.tbl[3])
lbls <- c(names(atb.lab.tbl)[1], names(atb.lab.tbl)[2], names(atb.lab.tbl)[3])
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep = "") # ad % to labels 
pie(slices,labels  =  lbls, col = c("#81D8D0", "lightpink1", "palegreen"),
    main = "Pie Chart of ATB during labour")
legend("topright", inset = .05, title = "ATB during labour", 
       fill = c("#81D8D0", "lightpink1", "palegreen"),
       c(names(atb.lab.tbl)[1], names(atb.lab.tbl)[2], names(atb.lab.tbl)[3]),
       horiz = FALSE)


# Density of Admission Weight 
par(mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(1,3))
weoght.admin = as.numeric(lullaby$Admission.weight)
hist(na.omit(weoght.admin),col = "#81D8D0",main = "Histogram of Admission Weight")
d <- density(na.omit(weoght.admin)) # returns the density data 
plot(d,col = "dodgerblue3",main = "Density of Admission Weight ")
polygon(d, col="#81D8D0", border="dodgerblue3")
boxplot(na.omit(weoght.admin),notch = TRUE,col = "#81D8D0" ,main = "Boxplot of  Neonatal")
abline(h = min(weoght.admin), col = "darkslategray3")
abline(h = max(weoght.admin), col = "darkslategray3")
abline(h = median(weoght.admin), col = "navy")
abline(h = quantile(weoght.admin, c(0.25, 0.75)), col = "dodgerblue3")


# Density of Apgar score
par(mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(1,3))
apgar.score = lullaby$Apgar.score
hist(na.omit(apgar.score),col = "#81D8D0",main = "Histogram of Apgar score")
d <- density(na.omit(apgar.score)) # returns the density data 
plot(d,col = "dodgerblue3",main = "Density of Apgar score")
polygon(d, col="#81D8D0", border="dodgerblue3")
boxplot(apgar.score,notch = TRUE,col = "#81D8D0", main = "Boxplot of  Apgar score")
abline(h = min(apgar.score), col = "darkslategray3")
abline(h = max(apgar.score), col = "darkslategray3")
abline(h = median(apgar.score), col = "navy")
abline(h = quantile(apgar.score, c(0.25, 0.75)), col = "dodgerblue3")


## Grouped Bar Plot of Symphtoms
Asymptomatic                  = lullaby$Main.symptoms..choice.0..Asymptomatic.
Fevers.Hypothermia            = lullaby$Main.symptoms..choice.1..Fevers...Hypothermia.
Hypoglycemia                  = lullaby$Main.symptoms..choice.2..Hypoglycemia.
Hypotonic                     = lullaby$Main.symptoms..choice.3..Hypotonic.
Irritability                  = lullaby$Main.symptoms..choice.4..Irritability.
Poor.Feeding                  = lullaby$Main.symptoms..choice.5..Poor.feeding.
Reduced.conscious.level       = lullaby$Main.symptoms..choice.6..Reduced.conscious.level.
Respiratory.distress          = lullaby$Main.symptoms..choice.7..Respiratory.distress.
Seizures                      = lullaby$Main.symptoms..choice.8..Seizures.
Umbilical.cord_skin.infection = lullaby$Main.symptoms..choice.9..Umbilical.cord...skin.infection.
Other                         = lullaby$Main.symptoms..choice.10..Other.

p <- plot_ly(
  x = c("Possitive", "Negative"),
  y = c(sum(Asymptomatic == "Checked"), sum(Asymptomatic == "Unchecked")), 
  name = "Asymptomatic", 
  type = "bar")
p

p2 <- add_trace(
  p,
  x = c("Possitive", "Negative"),
  y = c(sum(Fevers.Hypothermia == "Checked"), sum(Fevers.Hypothermia == "Unchecked")),
  name = "Fevers.Hypothermia",
  type = "bar")
p2

p3 <- add_trace(
  p2,
  x = c("Possitive", "Negative"),
  y = c(sum(Hypoglycemia == "Checked"), sum(Hypoglycemia == "Unchecked")),
  name = "Hypoglycemia",
  type = "bar")
p3

p4 <- add_trace(
  p3,
  x = c("Possitive", "Negative"),
  y = c(sum(Hypotonic == "Checked"), sum(Hypotonic == "Unchecked")),
  name = "Hypotonic",
  type = "bar")
p4

p5 <- add_trace(
  p4,
  x = c("Possitive", "Negative"),
  y = c(sum(Irritability == "Checked"), sum(Irritability == "Unchecked")),
  name = "Irritability",
  type = "bar")
p5

p6 <- add_trace(
  p5,
  x = c("Possitive", "Negative"),
  y = c(sum(Poor.Feeding == "Checked"), sum(Poor.Feeding == "Unchecked")),
  name = "Poor.Feeding",
  type = "bar")
p6

p7 <- add_trace(
  p6,
  x = c("Possitive", "Negative"),
  y = c(sum(Reduced.conscious.level == "Checked"), sum(Reduced.conscious.level == "Unchecked")),
  name = "Reduced.conscious.level",
  type = "bar")
p7

p8 <- add_trace(
  p7,
  x = c("Possitive", "Negative"),
  y = c(sum(Respiratory.distress == "Checked"), sum(Respiratory.distress == "Unchecked")),
  name = "Respiratory.distress",
  type = "bar")
p8

p9 <- add_trace(
  p8,
  x = c("Possitive", "Negative"),
  y = c(sum(Seizures == "Checked"), sum(Seizures == "Unchecked")),
  name = "Seizures",
  type = "bar")
p9



# Grouped Bar Plot of Diagnosis
Sepsis = lullaby$Diagnosis..choice.1..Sepsis.
Meningitis = lullaby$Diagnosis..choice.2..Meningitis.
NEC = lullaby$Diagnosis..choice.3..NEC.
HIE.Birth.asphyxia = lullaby$Diagnosis..choice.4..HIE.Birth.asphyxia.

p <- plot_ly(
  x = c("Possitive", "Negative"),
  y = c(sum(Sepsis == "Checked"), sum(Sepsis == "Unchecked")),
  name = "Sepsis",
  type = "bar")
p

p2 <- add_trace(
  p,
  x = c("Possitive", "Negative"),
  y = c(sum(Meningitis == "Checked"), sum(Meningitis == "Unchecked")),
  name = "Meningitis",
  type = "bar")
p2

p3 <- add_trace(
  p2,
  x = c("Possitive", "Negative"),
  y = c(sum(NEC == "Checked"), sum(NEC == "Unchecked")),
  name = "NEC",
  type = "bar")
p3

p4 <- add_trace(
  p3,
  x = c("Possitive", "Negative"),
  y = c(sum(HIE.Birth.asphyxia == "Checked"), sum(HIE.Birth.asphyxia == "Unchecked")),
  name = "HIE.Birth.asphyxia",
  type = "bar")
p4

##########################################################################################
# Circle Plot and Heat Plot
##########################################################################################
par(mfrow=c(1,1))

Pakistan.ATB.Trans_12 = matrix(0, nrow = 10, ncol = 10)
Pakistan.ATB.Trans_13 = matrix(0, nrow = 10, ncol = 10)
Pakistan.ATB.Trans_23 = matrix(0, nrow = 10, ncol = 10)

ATB.name <- c("Amikacin", "Ampicillin", "Cefotaxim", "Ceftriaxone","Cloxacillin","Gentamicin","Metronidazole","Vancomicin","Carbapenem","Other")


#####
# Index the ATB-Relavent Columns
index.atb1.commence = which(names(lullaby) == "ATB.of.regimen.1..choice.1..Amikacin.")
index.atb1.fini     = which(names(lullaby) == "ATB.of.regimen.1..choice.98..Other.")
index.atb2.commence = which(names(lullaby) == "ATB.of.regimen.2..choice.1..Amikacin.")
index.atb2.fini     = which(names(lullaby) == "ATB.of.regimen.2..choice.98..Other.")
index.atb3.commence = which(names(lullaby) == "ATB.of.regimen.3....choice.1..Amikacin.")
index.atb3.fini     = which(names(lullaby) == "ATB.of.regimen.3....choice.98..Other.")
atb.Resist.frame=as.data.frame(cbind(lullaby[,index.atb1.commence:index.atb1.fini], lullaby[,index.atb2.commence:index.atb2.fini], lullaby[,index.atb3.commence:index.atb3.fini]))


tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:10)
  {
    for (j in 1:10)
    {
      for (k in 1:length(atb.Resist.frame$ATB.of.regimen.1..choice.1..Amikacin.))
      {
        if ((atb.Resist.frame[k, i] == "Checked")&&(atb.Resist.frame[k, 10+j] == "Checked")) 
        {
          Pakistan.ATB.Trans_12[i, j] = Pakistan.ATB.Trans_12[i, j] + 1
        }
      }
    }
  }
  increase(pb)
}
print("Iteration Complete")
toc()


tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:10)
  {
    for (j in 1:10)
    {
      for (k in 1:length(atb.Resist.frame$ATB.of.regimen.1..choice.1..Amikacin.))
      {
        if ((atb.Resist.frame[k, 10 + i] == "Checked")&&(atb.Resist.frame[k, 20 + j] == "Checked")) 
        {
          Pakistan.ATB.Trans_23[i, j] = Pakistan.ATB.Trans_23[i, j] + 1
        }
      }
    }
  }
  increase(pb)
}
print("Iteration Complete")
toc()


tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:10)
  {
    for (j in 1:10)
    {
      for (k in 1:length(atb.Resist.frame$ATB.of.regimen.1..choice.1..Amikacin.))
      {
        if ((atb.Resist.frame[k, i] == "Checked")&&(atb.Resist.frame[k, 20+j] == "Checked")) 
        {
          Pakistan.ATB.Trans_13[i, j] = Pakistan.ATB.Trans_13[i, j] + 1
        }
      }
    }
  }
  increase(pb)
}
print("Iteration Complete")
toc()


# Plot Phase 1 to Phase 2
groupColors <- c("peachpuff", "#81D8D0", "thistle2","green4", "lawngreen","lightpink1","cornflowerblue","cyan","darkorchid1","violet")
print("Iteration Complete")
m = as.matrix(Pakistan.ATB.Trans_12)
rownames(m) = ATB.name
colnames(m) = ATB.name
chordDiagram(m, grid.col = groupColors, row.col = groupColors)
title("Phase I to Phase II")

# Plot Phase 2 to Phase 3
groupColors <- c("thistle2", "#81D8D0", "lightpink1","green4", "peachpuff","lightpink1","cornflowerblue","cyan","darkorchid1","violet")
m = as.matrix(Pakistan.ATB.Trans_23)
rownames(m) = ATB.name
colnames(m) = ATB.name
chordDiagram(m, grid.col = groupColors, row.col = groupColors)
title("Phase II to Phase III")

# Plot Phase 1 to Phase 3
groupColors <- c("peachpuff", "#81D8D0", "thistle2","green4", "lawngreen","lightpink1","cornflowerblue","cyan","darkorchid1","violet")
m = as.matrix(Pakistan.ATB.Trans_13)
rownames(m) = ATB.name
colnames(m) = ATB.name
chordDiagram(m, grid.col = groupColors, row.col = groupColors)
title("Phase I to Phase III")


############################################################################################
# ATBR
############################################################################################
atbr1.name.start=which(names(lullaby)=="E.Coli.bld.1.Amikacin")
atbr1.name.end=which(names(lullaby)=="Other.bld.1.ESBL.conf")
atbr2.name.start=which(names(lullaby)=="E.Coli.bld.2.Amikacin")
atbr2.name.end=which(names(lullaby)=="Other.bld.2.ESBL.conf")
atbr3.name.start=which(names(lullaby)=="E.Coli.csf.1.Amikacin")
atbr3.name.end=which(names(lullaby)=="Other.csf.1.ESBL.conf")
atbr4.name.start=which(names(lullaby)=="E.Coli.csf.2.Amikacin")
atbr4.name.end=which(names(lullaby)=="Other.csf.2.ESBL.conf")


atbr.frame = as.data.frame(cbind(lullaby[, atbr1.name.start:atbr1.name.end], lullaby[, atbr2.name.start:atbr2.name.end], lullaby[, atbr3.name.start:atbr3.name.end],  lullaby[, atbr4.name.start:atbr4.name.end]))
atbr.frame1 = data.frame(matrix(0, ncol = ncol(atbr.frame), nrow = nrow(atbr.frame))) 
# Compute ATBR into 0-1
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:nrow(atbr.frame))
  {
    for(j in 1:ncol(atbr.frame))
    {
      temp = ifelse(atbr.frame[i, j] == "3. R", 1, 0)
      atbr.frame1[i, j] = temp
    }
  }
  increase(pb)
}
print("Iteration Complete")
toc()

sum(atbr.frame1)


############################################################################################
# Survival
############################################################################################
par(mfrow=c(1,1))

# Determind Duration of Hospitalization
lullaby$Date.of.discharge = factor(lullaby$Date.of.discharge)
lullaby$Date.of.discharge = as.Date(lullaby$Date.of.discharge, format = "%Y-%m-%d")

lullaby$Date.of.LAMA = factor(lullaby$Date.of.LAMA)
lullaby$Date.of.LAMA = as.Date(lullaby$Date.of.LAMA, format = "%Y-%m-%d")

lullaby$Date.of.death = factor(lullaby$Date.of.death)
lullaby$Date.of.death = as.Date(lullaby$Date.of.death, format = "%Y-%m-%d")

lullaby$Date.of.admission = factor(lullaby$Date.of.admission)
lullaby$Date.of.admission = as.Date(lullaby$Date.of.admission, format = "%Y-%m-%d")


dur1=lullaby$Date.of.discharge - lullaby$Date.of.admission
dur2=lullaby$Date.of.LAMA - lullaby$Date.of.admission
dur3=lullaby$Date.of.death - lullaby$Date.of.admission

duration.frame = cbind(dur1, dur2, dur3)
Time = numeric(length = dim(lullaby)[1])

for (i in 1:dim(lullaby)[1])
{
  
  indicator.temp = min(which(!is.na(duration.frame[i,])))
  Time[i] = duration.frame[i, indicator.temp]
}


Patient.ID = lullaby$Patient.ID
lullaby$Baby.outcome = as.factor(lullaby$Baby.outcome)
MORIR = ifelse(lullaby$Baby.outcome == "1. Dead", 1, 0)
Apgar.score = lullaby$Apgar.score
Blood.cultures = as.factor(lullaby$Blood.cultures)

Pakistan_Survival <- as.data.frame(cbind(Patient.ID, MORIR, Time, Apgar.score, Blood.cultures, ATBR, Number.of.CRP, Weight.Loss, Mother.signs.of.infection, PROM.above.18h, ATB.during.labour, Sex))
dim(Pakistan_Survival)
names(Pakistan_Survival)


attach(Pakistan_Survival)
Pakistan_Survival                           = na.omit(Pakistan_Survival)
Pakistan_Survival$MORIR                     = as.factor(ifelse(lullaby$Baby.outcome))
Pakistan_Survival$Baby.outcome              = as.factor(Pakistan_Survival$Baby.outcome)
Pakistan_Survival$Blood.cultures            = as.factor(Pakistan_Survival$Blood.cultures)
Pakistan_Survival$ATBR                      = as.factor(Pakistan_Survival$ATBR)
Pakistan_Survival$Mother.signs.of.infection = as.factor(Pakistan_Survival$Mother.signs.of.infection)
Pakistan_Survival$PROM.above.18h            = as.factor(Pakistan_Survival$PROM.above.18h)
Pakistan_Survival$ATB.during.labour         = as.factor(Pakistan_Survival$ATB.during.labour)
Pakistan_Survival$Sex                       = as.factor(Pakistan_Survival$Sex)
Pakistan_Survival$Resuscitation.at.birth    = as.factor(Pakistan_Survival$Resuscitation.at.birth)
dim(Pakistan_Survival)
# [1] 178  14

# create a Surv object 
# Mortality Study
Pakistan_Survival$MORIR=as.factor(Pakistan_Survival$MORIR)
survobj <- with(Pakistan_Survival, Surv(Time, MORIR == "1"))

# Plot survival distribution of the total sample
# Kaplan-Meier estimator 
# Overall
surv.Pakistan.toy1 = survfit(survobj~1, data = Pakistan_Survival)
summary(surv.Pakistan.toy1)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(surv.Pakistan.toy1, col = "palevioletred3", xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100,
     main = "Survival Distribution (Overall)")
legend("bottomright", title = "Overall", c("Upper", "Value", "Lower"),
       fill = c("lightpink2", "palevioletred3", "lightpink2"),
       lty = c(2, 1, 2))

Pakistan_Survival = na.omit(Pakistan_Survival)
surv.Pakistan.toy = survfit(survobj ~ Apgar.score+Blood.cultures+ATBR+Number.of.CRP+Weight.Loss)
summary(surv.Pakistan.toy)

# Blood Culture Positive
surv.Pakistan.blood = survfit(survobj ~ Blood.cultures, data = Pakistan_Survival)
summary(surv.Pakistan.blood)
plot(surv.Pakistan.blood, xlab="Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("lightpink2","mediumturquoise"),
     main = "Survival Distributions by Blood Culture") 
legend("bottomright", title = "Blood", c("POSITIVE", "NEGATIVE"),
       fill = c("lightpink2", "mediumturquoise"))

# ATBR
surv.Pakistan.ATBR = survfit(survobj ~ ATBR, data = Pakistan_Survival)
summary(surv.Pakistan.ATBR)
plot(surv.Pakistan.ATBR, xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("lightpink2","mediumturquoise"),
     main = "Survival Distributions by ATB Resistance") 
legend("bottomright", title = "ATBR", c("Resist", "NON"),
       fill = c("lightpink2", "mediumturquoise"))

# Num of CRP
Pakistan_Survival$Number.of.CRP = as.factor(Pakistan_Survival$Number.of.CRP)
surv.Pakistan.CRP =  survfit(survobj ~ Number.of.CRP, data = Pakistan_Survival)
summary(surv.Pakistan.CRP)
plot(surv.Pakistan.CRP, xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("cyan", "royalblue3", "pink2", "firebrick4"),
     main = "Survival Distributions by CRP") 
legend("bottomright", title = "CRP", c("0", "1", "2", "3"),
       fill = c("cyan", "royalblue3", "pink2", "firebrick4"))

# Mother.signs.of.infection
surv.Pakistan.motherinfection = survfit(survobj ~ Mother.signs.of.infection, data = Pakistan_Survival)
summary(surv.Pakistan.motherinfection)
plot(surv.Pakistan.motherinfection, xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("lightpink2","mediumturquoise"),
     main = "Survival Distributions by Mother Infection") 
legend("bottomright", title = "Mother Infection", c("INFECTION", "NON"),
       fill = c("lightpink2","mediumturquoise"))

# PROM.above.18h
surv.Pakistan.prom18h = survfit(survobj ~ PROM.above.18h, data = Pakistan_Survival)
summary(surv.Pakistan.prom18h)
plot(surv.Pakistan.prom18h, xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("lightpink2","mediumturquoise"),
     main = "Survival Distributions by PROM > 18h") 
legend("bottomright", title = "Mother Infection", c("OUI", "NON"),
       fill = c("lightpink2","mediumturquoise"))

#ATB.during.labour
surv.Pakistan.labourATB = survfit(survobj ~ ATB.during.labour, data = Pakistan_Survival)
summary(surv.Pakistan.labourATB)
plot(surv.Pakistan.labourATB, xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("lightpink2","mediumturquoise"),
     main = "Survival Distributions by ATB During Labour") 
legend("bottomright", title = "Mother Infection", c("OUI", "NON"),
       fill = c("lightpink2","mediumturquoise"))

#Sex
surv.Pakistan.sex = survfit(survobj ~ Sex, data = Pakistan_Survival)
summary(surv.Pakistan.sex)
plot(surv.Pakistan.sex, xlab = "Survival Time in Days", 
     ylab = "% Surviving", yscale = 100, col = c("mediumturquoise","lightpink2"),
     main = "Survival Distributions by Gender") 
legend("bottomright", title = "Mother Infection", c("M", "F"),
       fill = c("mediumturquoise","lightpink2"))

############################################################################################
# GLM
#Pakistan.glm.050313 <- glm(MORIR ~ Apgar.score + Weight.Loss + Time+ Number.of.CRP+ Blood.cultures + ATBR + Mother.signs.of.infection + PROM.above.18h + ATB.during.labour + Sex + Resuscitation.at.birth, family = "binomial")
Pakistan.glm.050313 <- glm(MORIR ~ Apgar.score + Weight.Loss + Number.of.CRP+ Blood.cultures + ATBR + Mother.signs.of.infection + PROM.above.18h + ATB.during.labour + Sex, data = Pakistan_Survival, family = "binomial")

summary(Pakistan.glm.050313)

# Call:
# glm(formula = MORIR ~ Apgar.score + Weight.Loss + Number.of.CRP + 
#       Blood.cultures + ATBR + Mother.signs.of.infection + PROM.above.18h + 
#       ATB.during.labour + Sex + Resuscitation.at.birth, family = "binomial")
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.8579  -0.4526  -0.3368  -0.2673   2.8125  
# 
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                3.8792739  1.9850024   1.954  0.05067 .  
# Apgar.score               -0.5570064  0.1616727  -3.445  0.00057 ***
# Weight.Loss                0.0020116  0.0008551   2.352  0.01865 *  
# Number.of.CRP              0.2474635  0.4244196   0.583  0.55985    
# Blood.cultures            -1.5754742  1.2367227  -1.274  0.20270    
# ATBR                       0.8021073  0.6877958   1.166  0.24353    
# Mother.signs.of.infection -0.0623218  0.8720595  -0.071  0.94303    
# PROM.above.18h            -0.0505768  0.9564343  -0.053  0.95783    
# ATB.during.labour          0.1699589  1.0541024   0.161  0.87191    
# Sex                        0.0754860  0.5554840   0.136  0.89191    
# Resuscitation.at.birth    -0.8795043  0.7219456  -1.218  0.22313    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 129.68  on 179  degrees of freedom
# Residual deviance: 105.67  on 169  degrees of freedom
# (11 observations deleted due to missingness)
# AIC: 127.67
# 
# Number of Fisher Scoring iterations: 5
# 



sum(Pakistan_Survival$MORIR == "1")
sum(Pakistan_Survival$MORIR == "0")
Cutoff.050313 = sum(Pakistan_Survival$MORIR == "1")/(sum(Pakistan_Survival$MORIR == "1") + sum(Pakistan_Survival$MORIR == "0"))
# [1] 0.8820225
pre <- predict(Pakistan.glm.050313, type='response')

pred.050313  = ifelse(pre > Cutoff.050313, 1, 0)
pred.050313  = ifelse(pre > 0.133, 1, 0)

pred.default = ifelse(pre > 0.5, 1, 0)

table(pred.050313,  Pakistan_Survival$MORIR)
# pred.050313   0   1
#           0 131   8
#           1  26  13

table(pred.default, Pakistan_Survival$MORIR)
# pred.default   0   1
#            0 155  17
#            1   2   4

# ROC Curve
pre <- predict(Pakistan.glm.050313, type='response')
# Probability with Reality
data <- data.frame(prob=pre, obs=Pakistan_Survival$MORIR)
# Sort
data <- data[order(data$prob), ]
n <- nrow(data)
tpr <- fpr <- rep(0, n)
# TPR FPR
for (i in 1:n) 
{
  threshold <- data$prob[i]
  tp <- sum(data$prob > threshold & data$obs == 1)
  fp <- sum(data$prob > threshold & data$obs == 0)
  tn <- sum(data$prob < threshold & data$obs == 0)
  fn <- sum(data$prob < threshold & data$obs == 1)
  tpr[i] <- tp/(tp+fn) # TPR
  fpr[i] <- fp/(tn+fp) # FPR
}
# plot(fpr, tpr, type = 'l', col = "#81D8D0", lwd = 1.5)
# polygon(rbind(fpr, 1), rbind(tpr, 0), col = "grey4")
# title("ROC of GLM")

# plot(fpr,tpr,type = 'l', col = "#81D8D0", lwd = 2)
# polygon(cbind(fpr, 1), cbind(tpr, 0), col = "paleturquoise1")
# title("ROC of GLM")

# Or just use package
library(pROC)
modelroc <- roc(Pakistan_Survival$MORIR, pre)
plot(modelroc, print.auc = TRUE, 
     auc.polygon = TRUE, grid = c(0.1, 0.1),grid.col = c("#81D8D0", "#81D8D0"), 
     max.auc.polygon = TRUE, auc.polygon.col = "paleturquoise1", 
     print.thres = TRUE, print.thres.pch=35, print.thres.cex=par("cex"),
     identity.col = "#81D8D0", identity.lwd=2.5)
title("ROC of GLM", line = 2.75)
legend("topright", title = "AUC", c("AUC = 0.832"),
       fill = c("#81D8D0"))



























































