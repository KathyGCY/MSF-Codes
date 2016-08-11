# MSF Protocal
# What would Kathy do?
# And now we know what happens when a Stage IV OCD codes :)


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
##########################################################################################
# Import Data
##########################################################################################

data.entire <- read.csv("~/Desktop/MSF/Haiti/HaitiDrouillardBurnS_DATA_2016-06-10_2101 copy.csv")


##########################################################################################
# Cleand the data
##########################################################################################

# Compute original data into a redunduncy-eliminated data
# Eliminate Redundent Observations Preserving Information

# NOTES:
# Blood Culture Positive is defined as at least one of the culture is positive

# Let's see what the data looks like
# Check if there's anything wrong like mis-matched variable names
par(mfrow=c(1,1))
names(data.entire)

# We first extract the first obs of all subjects with info and outcome
dim.entire         = dim(data.entire)
nrow.entire        = length(unique(data.entire$bs_pat_id))
ncol.entire        = dim.entire[2]
# 2065 245  
data.entire.name   = names(data.entire)
data.unique        = data.frame(matrix(vector(), nrow.entire, ncol.entire, dimnames = list(c(), data.entire.name)), stringsAsFactors = F)
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for(i in 1:nrow.entire)
  {
    index.temp.DEP   = min(which(data.entire$Organism == i))
    data.unique[i, ] = data.entire[index.temp.DEP, ]
  }
  
  increase(pb)
}
print("Iteration Complete")
toc()

# View(data.unique)

# Now we add a new column indicating ATB Resistance RISN
# ATB Names with RISN is in which columns 
col.start = which(names(data.entire) == "bs_amikacine")
col.end   = which(names(data.entire) == "bs_vanco")

# Iteration
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:nrow.entire)
  {
    # Create temporary data frame with same ID == i
    # Index all patient id                     == i rownumber
    index.temp.1 = which(data.entire$bs_pat_id == i)
    # Extract these rows
    temp.frame.1 = data.entire[index.temp.1, ]
    for (j in col.start:col.end)
    {
      # For all ATBs find the first non N/A value with RISN
      # Index the first non N/A row
      temp.frame.1.index = min(which(temp.frame.1[j] != "N/A"))
      # Find the Colnm   = the certain ATB and row = index
      data.unique[i, j]  = temp.frame.1[temp.frame.1.index, j]
    }
  }
  increase(pb)
}
print("Iteration Complete")
toc()


##########################################################################################
# Compute Variables and Define Categoricals
##########################################################################################

bs_death            = as.factor(data.unique$bs_death)
bs_inj2hosp         = data.unique$bs_hr_bru2hosp
bs_TBSA             = data.unique$bs_surface
bs_age              = data.unique$bs_age
bs_thicken          = as.factor(data.unique$bs_prof)
bs_inhal            = as.factor(data.unique$bs_inhalation)
bs_lieu             = as.factor(data.unique$bs_lieu)
bs_envir            = as.factor(data.unique$bs_environ)
bs_maln             = as.factor(data.unique$bs_malnut)
bs_smoke            = as.factor(data.unique$bs_tabac)
bs_hiv              = as.factor(data.unique$bs_vih)
bs_db               = as.factor(data.unique$bs_diabete)
bs_exam             = as.factor(data.unique$bs_exam)
bs_gender           = as.factor(data.unique$bs_sexe)
bs_etio___1         = as.factor(data.unique$bs_etio___1)
bs_etio___2         = as.factor(data.unique$bs_etio___2)
bs_etio___3         = as.factor(data.unique$bs_etio___3)
bs_etio___4         = as.factor(data.unique$bs_etio___4)
bs_etio___5         = as.factor(data.unique$bs_etio___5)
bs_etio___6         = as.factor(data.unique$bs_etio___6)
bs_etio___7         = as.factor(data.unique$bs_etio___7)
bs_etio___9         = as.factor(data.unique$bs_etio___9)
bs_circonstance___1 = as.factor(data.unique$bs_circonstance___1)
bs_circonstance___2 = as.factor(data.unique$bs_circonstance___2)
bs_circonstance___3 = as.factor(data.unique$bs_circonstance___3)
bs_loc_bru___1      = as.factor(data.unique$bs_loc_bru___1)
bs_loc_bru___2      = as.factor(data.unique$bs_loc_bru___2)
bs_loc_bru___3      = as.factor(data.unique$bs_loc_bru___3)
bs_loc_bru___4      = as.factor(data.unique$bs_loc_bru___4)
bs_loc_bru___5      = as.factor(data.unique$bs_loc_bru___5)
bs_loc_bru___6      = as.factor(data.unique$bs_loc_bru___6)
bs_loc_bru___7      = as.factor(data.unique$bs_loc_bru___7)
bs_loc_bru___8      = as.factor(data.unique$bs_loc_bru___8)
bs_loc_bru___9      = as.factor(data.unique$bs_loc_bru___9)
bs_loc_bru___10     = as.factor(data.unique$bs_loc_bru___10)
bs_loc_bru___11     = as.factor(data.unique$bs_loc_bru___11)
bs_loc_bru___12     = as.factor(data.unique$bs_loc_bru___12)
bs_loc_bru___13     = as.factor(data.unique$bs_loc_bru___13)
bs_loc_bru___14     = as.factor(data.unique$bs_loc_bru___14)
bs_loc_bru___15     = as.factor(data.unique$bs_loc_bru___15)
bs_loc_bru___16     = as.factor(data.unique$bs_loc_bru___16)
bs_loc_bru___17     = as.factor(data.unique$bs_loc_bru___17)
bs_loc_bru___18     = as.factor(data.unique$bs_loc_bru___18)
bs_loc_bru___99     = as.factor(data.unique$bs_loc_bru___99)

bs_loc_risk___1     = as.factor(data.unique$bs_loc_risque___1)
bs_loc_risk___2     = as.factor(data.unique$bs_loc_risque___2)
bs_loc_risk___3     = as.factor(data.unique$bs_loc_risque___3)
bs_loc_risk___4     = as.factor(data.unique$bs_loc_risque___4)
bs_loc_risk___5     = as.factor(data.unique$bs_loc_risque___5)

bs_amikacine        = as.factor(data.unique$bs_amikacine)
bs_amoxyclav        = as.factor(data.unique$bs_amoxyclav)
bs_amoxy            = as.factor(data.unique$bs_amoxy)
bs_ampisulba        = as.factor(data.unique$bs_ampisulba)
bs_ampi             = as.factor(data.unique$bs_ampi)
bs_cefepime         = as.factor(data.unique$bs_cefepime)
bs_cefezidime       = as.factor(data.unique$bs_cefezidime)
bs_cefotaxime       = as.factor(data.unique$bs_cefotaxime)
bs_cefoxit          = as.factor(data.unique$bs_cefoxit)
bs_ceftazid         = as.factor(data.unique$bs_ceftazid)
bs_ceftri           = as.factor(data.unique$bs_ceftri)
bs_cipro            = as.factor(data.unique$bs_cipro)
bs_clinda           = as.factor(data.unique$bs_clinda)
bs_colistine        = as.factor(data.unique$bs_colistine)
bs_cotri            = as.factor(data.unique$bs_cotri)
bs_erythro          = as.factor(data.unique$bs_erythro)
bs_genta            = as.factor(data.unique$bs_genta)
bs_imipen           = as.factor(data.unique$bs_imipen)
bs_meropen          = as.factor(data.unique$bs_meropen)
bs_oxacil           = as.factor(data.unique$bs_oxacil)
bs_peni             = as.factor(data.unique$bs_peni)
bs_pipetazo         = as.factor(data.unique$bs_pipetazo)
bs_tetra            = as.factor(data.unique$bs_tetra)
bs_ticarclav        = as.factor(data.unique$bs_ticarclav)
bs_tobra            = as.factor(data.unique$bs_tobra)
bs_vanco            = as.factor(data.unique$bs_vanco)

# Compute into factor
# With base level 1 = non resist, and level 2 = resist
bs_amikacine.R      = as.factor(ifelse(bs_amikacine   == "1", 1, 0))
bs_amoxyclav.R      = as.factor(ifelse(bs_amoxyclav   == "1", 1, 0))
bs_amoxy.R          = as.factor(ifelse(bs_amoxy       == "1", 1, 0))
bs_ampisulba.R      = as.factor(ifelse(bs_ampisulba   == "1", 1, 0))
bs_ampi.R           = as.factor(ifelse(bs_ampi        == "1", 1, 0))
bs_cefepime.R       = as.factor(ifelse(bs_cefepime    == "1", 1, 0))
bs_cefezidime.R     = as.factor(ifelse(bs_cefezidime  == "1", 1, 0))
bs_cefotaxime.R     = as.factor(ifelse(bs_cefotaxime  == "1", 1, 0))
bs_cefoxit.R        = as.factor(ifelse(bs_cefoxit     == "1", 1, 0))
bs_ceftazid.R       = as.factor(ifelse(bs_ceftazid    == "1", 1, 0))
bs_ceftri.R         = as.factor(ifelse(bs_ceftri      == "1", 1, 0))
bs_cipro.R          = as.factor(ifelse(bs_cipro       == "1", 1, 0))
bs_clinda.R         = as.factor(ifelse(bs_clinda      == "1", 1, 0))
bs_colistine.R      = as.factor(ifelse(bs_colistine   == "1", 1, 0))
bs_cotri.R          = as.factor(ifelse(bs_cotri       == "1", 1, 0))
bs_erythro.R        = as.factor(ifelse(bs_erythro     == "1", 1, 0))
bs_genta.R          = as.factor(ifelse(bs_genta       == "1", 1, 0))
bs_imipen.R         = as.factor(ifelse(bs_imipen      == "1", 1, 0))
bs_meropen.R        = as.factor(ifelse(bs_meropen     == "1", 1, 0))
bs_oxacil.R         = as.factor(ifelse(bs_oxacil      == "1", 1, 0))
bs_peni.R           = as.factor(ifelse(bs_peni        == "1", 1, 0))
bs_pipetazo.R       = as.factor(ifelse(bs_pipetazo    == "1", 1, 0))
bs_tetra.R          = as.factor(ifelse(bs_tetra       == "1", 1, 0))
bs_ticarclav.R      = as.factor(ifelse(bs_ticarclav   == "1", 1, 0))
bs_tobra.R          = as.factor(ifelse(bs_tobra       == "1", 1, 0))
bs_vanco.R          = as.factor(ifelse(bs_vanco       == "1", 1, 0))


##########################################################################################
# Explorative Data Analysis and Visualization
##########################################################################################

# Age
bs_age = data.unique$bs_age
# 5 point summary
summary(bs_age)
# Min. 1st Qu. Median Mean 3rd Qu. Max. 
# 

# Histogramme
hist(bs_age,        col = "cyan")
hist(log10(bs_age), col = "cyan")

# Kernal Density
d              = density(na.omit(bs_age))
plot(d, main   = "Kernel Density of Subject Age")
polygon(d, col = "cyan", border = "dodgerblue4")

# Boxplot (Not displayed in final report)
bs_age   = na.omit(bs_age)
boxplot(bs_age, notch = TRUE, col = "cyan1" , main = "Boxplot of Age")
abline(h = min(bs_age),       col = "Blue")
abline(h = max(bs_age),       col = "Yellow")
abline(h = median(bs_age),    col = "Green")
abline(h = quantile(bs_age, c(0.25, 0.75)), col = "Red")


# Time from injuary 2 hospital
bs_inj2hosp = data.unique$bs_hr_bru2hosp
hist(bs_inj2hosp, col = "cyan")
summary(bs_inj2hosp)

# This part is only for outliers shown of previous histogramme
index.outlier=which.max(bs_inj2hosp)
# [1] 
max(na.omit(bs_inj2hosp))
# 
bs_inj2hosp.r = bs_inj2hosp[-index.outlier]
hist(bs_inj2hosp.r, col = "cyan")
summary(na.omit(bs_inj2hosp.r))
# Min. 1st Qu. Median Mean 3rd Qu. Max. 
# 


#  MSF
bs_MSF           = as.factor(data.unique$bs_struct_MSF)
tb.msf           = table(bs_MSF)
bs_MSF           <- data.frame(labels = c("NON-MSF", "MSF"), 
                               values        = c(tb.msf[1], tb.msf[2]))
plot_ly(bs_MSF, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of MSF Structure")

# Gender
bs_gender        = as.factor(data.unique$bs_sexe)
gender.tbl       = table(bs_gender)
tb.gender        = table(bs_gender)
bs_sex           <- data.frame(labels = c("Male", "Female"), 
                               values        = c(tb.gender[1], tb.gender[2]))
plot_ly(bs_sex, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Gender")

# Outcome
bs_sortie        = as.factor(data.unique$bs_type_sortie)
tb.sortie        = table(bs_sortie)
bs_sortie        <- data.frame(labels = c("Home", "Against Med Advice", "Other Hospital", "Dead"), 
                               values        = c(tb.sortie[1], tb.sortie[2], tb.sortie[3], tb.sortie[4] ))
plot_ly(bs_sortie, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Outcome")

# Exam in a care structure
bs_exam          = as.factor(data.unique$bs_exam)
tb.exam          = table(bs_exam)
bs_exam          <- data.frame(labels = c("Exam", "No Exam"), 
                               values        = c(tb.exam[1], tb.exam[2]))
plot_ly(bs_exam, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Exam")

# HIV
bs_hiv           = as.factor(data.unique$bs_vih)
tb.hiv           = table(bs_hiv)
bs_hiv           <- data.frame(labels = c("Oui", "NON", "Unknown"), 
                               values        = c(tb.hiv[1], tb.hiv[2], tb.hiv[3]))
plot_ly(bs_hiv, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of HIV")

# Diabete
bs_db            = as.factor(data.unique$bs_diabete)
tb.db            = table(bs_db)
bs_db            <- data.frame(labels = c("Diabete", "No Diabete", "Unknown"), 
                               values        = c(tb.db[1], tb.db[2], tb.db[3]))
plot_ly(bs_db, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Diabete")

# Smoking
bs_smoke         = as.factor(data.unique$bs_tabac)
tb.smoke         = table(bs_smoke)
bs_smoke         <- data.frame(labels = c("Smoke", "No Smoke", "Unknown"), 
                               values        = c(tb.smoke[1], tb.smoke[2], tb.smoke[3]))
plot_ly(bs_smoke, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Smoke")

# Malnutrition
bs_maln          = as.factor(data.unique$bs_malnut)
tb.maln          = table(bs_maln)
bs_maln          <- data.frame(labels = c("Malnutrition", "Not Maln", "Unknown"), 
                               values        = c(tb.maln[1], tb.maln[2], tb.maln[3]))
plot_ly(bs_maln, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Malnutrition")

# ETIO
index.start.etio = which(names(data.unique) == "bs_etio___1")
index.end.etio   = which(names(data.unique) == "bs_etio___9")

etio             = as.data.frame (cbind(data.unique[, index.start.etio:index.end.etio]))
num.etio         = colSums (etio, na.rm = TRUE, dims = 1)
percent.etio     = colSums (etio, na.rm = TRUE, dims = 1)*100/nrow(data.unique)
num.etio.sort    = sort.int (num.etio, decreasing = TRUE, index.return = TRUE)

loc              = loc[c(num.loc.sort$ix)]
ETIO.name        = c("Chemical", "Contact hot object", "Electrocution", "Flame", "Gas Explosion", "Hot Liquid", "Raidance", "Others")
ETIO.name        = ETIO.name[c(num.etio.sort$ix)]
frep.Check       = num.etio.sort$x
freq.Uncheck     = nrow(data.unique)-frep.Check
p     <- plot_ly(
  x    = ETIO.name, 
  y    = freq.Uncheck, 
  name = "Unckeck", 
  type = "bar")
p
p2    <- add_trace(
  p, 
  x    = ETIO.name, 
  y    = frep.Check, 
  name = "Check", 
  type = "bar")
p2 %>%
  layout(title = "Etiology")

p     <- plot_ly(
  x    = ETIO.name, 
  y    = frep.Check, 
  name = "Ckeck", 
  type = "bar")
p%>%
  layout(title = "Etiology")

# Circonstance
index.cir.start = which(names(data.unique) == "bs_circonstance___1")
index.cir.end   = which(names(data.unique) == "bs_circonstance___3")

circ            = as.data.frame(cbind(data.unique[, index.cir.start:index.cir.end]))
num.circ        = colSums (circ, na.rm = TRUE, dims = 1)
percent.circ    = colSums (haiti.circ, na.rm = TRUE, dims = 1)*100/nrow(data.unique)
circonstance    = c("Violence", "Accident", "Explosion")
freq.Check      = num.circ
freq.Uncheck    = nrow(data.unique)-freq.Check

p     <- plot_ly(
  x    = circonstance, 
  y    = freq.Check, 
  name = "Ckeck", 
  type = "bar")
p%>%
  layout(title = "Circonstance")

# Accident
index.accident.start = which(names(data.unique) == "bs_accident___1")
index.accident.end   = which(names(data.unique) == "bs_accident___9")

haiti.accident       = as.data.frame(cbind(data.unique[, index.accident.start:index.accident.end]))
num.accident         = colSums (haiti.accident, na.rm = TRUE, dims = 1)
percent.accident     = colSums (haiti.circ, na.rm = TRUE, dims = 1)*100/nrow(data.unique)
accident             = c("Traffic", "Domestic", "Work", "Explosion", "Others")
freq.Check           = num.accident
freq.Uncheck         = nrow(data.unique)-freq.Check

p     <- plot_ly(
  x    = accident, 
  y    = freq.Check, 
  name = "Ckeck", 
  type = "bar")
p%>%
  layout(title = "Accident")

# Environment
bs_envir  = as.factor(data.unique$bs_environ)
tb.envir  = table(bs_envir)
bs_envir  <- data.frame(labels = c("Urban", "Rural"), 
                        values = c(tb.envir[1], tb.envir[2]))
plot_ly(bs_envir, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Environment")

# Particular location
bs_lieu   = as.factor(data.unique$bs_lieu)
tb.lieu   = table(bs_lieu)
bs_lieu   <- data.frame(labels = c("Home", "Public Space", "Work", "Other"), 
                        values = c(tb.lieu[1], tb.lieu[2] , tb.lieu[3] , tb.lieu[4] ))
plot_ly(bs_lieu, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Particular Location")

# Body surface area
bs_TBSA = data.unique$bs_surface
hist(bs_TBSA, col = "cyan")
summary(na.omit(bs_TBSA))
# Min. 1st Qu. Median Mean 3rd Qu. Max. 
# 

# Location of wound
index.loc.start = which(names(data.unique) == "bs_loc_bru___1")
index.loc.end   = which(names(data.unique) == "bs_loc_bru___99")

haiti.loc       = as.data.frame(cbind(data.unique[, index.loc.start:index.loc.end]))
num.loc         = colSums (haiti.loc, na.rm = TRUE, dims = 1)
num.loc.sort    = sort.int(colSums (haiti.loc, na.rm = TRUE, dims = 1), decreasing = TRUE, index.return = TRUE)

loc             = c("Skull", "Face", "Neck", "Right Hand", "Left Hand", "Sup Right Link", "Sup Left Limb", "Anterior Thorax", "Posterior Thorax", "Anterior Abdomen", "Posterior Abdomen", "Buttocks", "Genitals", "Perineum", "Inf Right Limb", "Inf Left Limb", "Right Foot", "Left Foot", "Others")
loc             = loc[c(num.loc.sort$ix)]
freq.Check      = num.loc.sort
freq.Uncheck    = nrow(data.unique)-freq.Check

p     <- plot_ly(
  x    = loc, 
  y    = freq.Check, 
  name = "Ckeck", 
  type = "bar")
p%>%
  layout(title = "Location")

# Location that has Functional risk
index.func.start = which(names(data.unique) == "bs_loc_risque___1")
index.func.end   = which(names(data.unique) == "bs_loc_risque___5")

haiti.locfr      = as.data.frame(cbind(data.unique[, index.func.start:index.func.end]))
num.func         = colSums (haiti.locfr, na.rm = TRUE, dims = 1)

loc.func         = c("Face", "Neck", "Limb", "GenitalsPerineum", "Others")
num.func.sort    = sort.int(num.func, decreasing = TRUE, index.return = TRUE)

loc.func         = loc.func[c(num.func.sort$ix)]

freq.Check       = num.func.sort$x
freq.Uncheck     = nrow(data.unique)-freq.Check

p     <- plot_ly(
  x    = loc.func, 
  y    = freq.Check, 
  name = "Ckeck", 
  type = "bar")
p%>%
  layout(title = "Location Functional Risk")

# Inhalation Injury
bs_inhal = as.factor(data.unique$bs_inhalation)
tb.inhal = table(bs_inhal)
bs_inhal <- data.frame(labels = c("No Inhal Injury", "Inhal Injury"), 
                       values = c(tb.inhal[1] , tb.inhal[2] ))
plot_ly(bs_inhal, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Inhalation Injury")

# Thickeness
bs_thicken = as.factor(data.unique$bs_prof)
tb.thick = table(bs_thicken)
bs_thicken <- data.frame(labels = c("Partial Thickness", "Full-Thickness"), 
                         values = c(tb.thick[1] , tb.thick[2] ))
plot_ly(bs_thicken, labels = labels, values = values, type = "pie") %>%
  layout(title = "Pie Chart Of Prodondeur")


##########################################################################################
# Death Tree
##########################################################################################

# Simplified Tree
tree.death<- rpart(bs_death ~ bs_thicken + bs_inhal + bs_TBSA + bs_lieu + bs_envir + bs_maln + bs_smoke + bs_db + bs_hiv + bs_exam + bs_gender + bs_inj2hosp + bs_age + bs_TBSA, parms = list(split = 'information'))
par(mar = c(.5, .5, .5, .5) , 0.1)
plot(tree.death)
text(tree.death, use.n = TRUE)
title("Death Tree")
par(mar = c(4, 4, 2, 0.5))
summary(tree.death)


##########################################################################################
# Random Forest
##########################################################################################

rf.data=na.omit(cbind(bs_death,bs_thicken,bs_inhal,bs_TBSA,bs_lieu,bs_envir,bs_maln,bs_smoke,bs_db,bs_hiv,bs_exam,bs_gender,bs_inj2hosp,bs_age,bs_TBSA,bs_etio___1,bs_etio___2,bs_etio___3,bs_etio___4,bs_etio___5,bs_etio___6,bs_etio___7,bs_etio___9,bs_circonstance___1,bs_circonstance___2,bs_circonstance___3,bs_loc_bru___1,bs_loc_bru___2,bs_loc_bru___3,bs_loc_bru___4,bs_loc_bru___5,bs_loc_bru___6,bs_loc_bru___7,bs_loc_bru___8,bs_loc_bru___9,bs_loc_bru___10,bs_loc_bru___11,bs_loc_bru___12,bs_loc_bru___13,bs_loc_bru___14,bs_loc_bru___15,bs_loc_bru___16,bs_loc_bru___17,bs_loc_bru___18,bs_loc_bru___99,bs_loc_risk___1,bs_loc_risk___2,bs_loc_risk___3,bs_loc_risk___4,bs_loc_risk___5))

rf.death <- randomForest(bs_death ~ bs_thicken+bs_inhal+bs_TBSA+bs_lieu+bs_envir+bs_maln+bs_smoke+bs_db+bs_hiv+bs_exam+bs_gender+bs_inj2hosp+bs_age+bs_TBSA+bs_etio___1+bs_etio___2+bs_etio___3+bs_etio___4+bs_etio___5+bs_etio___6+bs_etio___7+bs_etio___9+bs_circonstance___1+bs_circonstance___2+bs_circonstance___3+bs_loc_bru___1+bs_loc_bru___2+bs_loc_bru___3+bs_loc_bru___4+bs_loc_bru___5+bs_loc_bru___6+bs_loc_bru___7+bs_loc_bru___8+bs_loc_bru___9+bs_loc_bru___10+bs_loc_bru___11+bs_loc_bru___12+bs_loc_bru___13+bs_loc_bru___14+bs_loc_bru___15+bs_loc_bru___16+bs_loc_bru___17+bs_loc_bru___18+bs_loc_bru___99+bs_loc_risk___1+bs_loc_risk___2+bs_loc_risk___3+bs_loc_risk___4+bs_loc_risk___5,data=rf.data)
out.rf <- predict(rf.death,rf.data[,-1])
print(rf.death)

sort(rf.death$importance,decreasing = TRUE)
n.tree.rf=seq(1:500)
mse.rf=rf.death$mse
plot(rf.death)
p<-plot_ly( x = n.tree.rf, y = mse.rf, mode = "markers",color = mse.rf)
p%>%add_trace(x=seq(0:500),y=rep(0.07,501))%>%
  layout(title = "Random Forest MSE")
outcome.rf.frame=as.data.frame(cbind(out.rf,bs_death))

sum((out.rf>1)&(bs_death==2))
sum(out.rf==1)
summary(rf.death)
importance(rf.death)


##########################################################################################
# ATB Resistance Contengincy Tables and Chi Square Test
##########################################################################################

# Import Data
bs_hemoc_inf     = ifelse(as.factor(data.entire$bs_hemoc_inf)    == "3", 1, 0)
bs_hemoc_inf[is.na(bs_hemoc_inf)] <- 0
bs_compl_infect  = ifelse(as.factor(data.entire$bs_compl_infect) == "1", 1, 0)
bs_compl_infect[is.na(bs_compl_infect)] <- 0
bs_morir         = ifelse(as.factor(data.entire$bs_type_sortie)  == "4", 1, 0)
bs_morir[is.na(bs_morir)] <- 0
bs_pat_id        = data.entire$bs_pat_id
atb.name.start   = which(names(data.entire) == "bs_amikacine")
atb.name.end     = which(names(data.entire) == "bs_vanco")

atbr.frame       = as.data.frame(data.entire[, atb.name.start:atb.name.end])
data.contengincy = as.data.frame(cbind(bs_pat_id, bs_morir, bs_hemoc_inf, bs_compl_infect, atbr.frame))
data.contengincy[is.na(data.contengincy)] <- 0
# Compute ATBR into 0-1
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:nrow(data.contengincy))
  {
    for(j in 5:ncol(data.contengincy))
    {
      data.contengincy[i, j] = ifelse(data.contengincy[i, j] == 1, 1, 0)
    }
  }
  increase(pb)
}
print("Iteration Complete")
toc()

# Initialization
# Creating and Naming Empty Column Vectors
bs_death_blood          = numeric(length(unique(data.contengincy$bs_pat_id)))
bs_bloodculturepositive = numeric(length(unique(data.contengincy$bs_pat_id)))
bs_infectionconfirm     = numeric(length(unique(data.contengincy$bs_pat_id)))
bs_atb_R                = numeric(length(unique(data.contengincy$bs_pat_id)))
# Here is where you could add new variables if interest
# 
# 

# Iteration
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:length(bs_bloodculturepositive))
  {
    # Create temporary data frame with same                      ID == i
    # For example,             i = 1 selects all                "id == 1" rows of data frame
    # Index all patient                                          id == i rownumber
    index.temp.1                 = which(data.contengincy$bs_pat_id == i)
    # Extract these rows
    temp.frame.1                 = data.contengincy[index.temp.1, 1:4]
    temp.resist.frame.1          = data.contengincy[index.temp.1, -c(1:4)]
    # as long as there is one positive blood culture
    # Now we have that temporary data frame of                   id == i
    # Classified as possitive 
    # < = > At least one of the tests is positive
    # < = > Summation of all test results as numerical elements > 0
    temp.death                   = sum(temp.frame.1[, 2])
    temp.blood                   = sum(temp.frame.1[, 3])
    temp.infconfirm              = sum(temp.frame.1[, 4])
    temp.resist                  = sum(temp.resist.frame.1[, -c(1:4)])
    bs_death_blood[i]            = ifelse(temp.death      != 0, 1, 0)
    bs_bloodculturepositive[i]   = ifelse(temp.blood      != 0, 1, 0)
    bs_infectionconfirm[i]       = ifelse(temp.infconfirm != 0, 1, 0)
    bs_atb_R[i]                  = ifelse(temp.resist     != 0, 1, 0)
  }
  increase(pb)
}
print("Iteration Complete")
toc()

bs_blood_id                    = seq(1:length(bs_bloodculturepositive))
data.for.DeathVSInfection      = as.data.frame(cbind(bs_blood_id, bs_bloodculturepositive, bs_complex_inf, bs_death_blood, bs_infectionconfirm, bs_atb_R))
head(data.for.DeathVSInfection)

# Contengincy Table Chi Square Test
tbl.bloodmortality             = table(HaitiBlood_Rfile[, c(2, 4)])
colnames(tbl.bloodmortality)   = c("Survive", "Dead")
rownames(tbl.bloodmortality)   = c("Negative", "Positive")
tbl.bloodmortality
chisq.test(tbl.bloodmortality)

tbl.complexmortality           = table(HaitiBlood_Rfile[, c(3, 4)])
colnames(tbl.complexmortality) = c("Survive", "Dead")
rownames(tbl.complexmortality) = c("NO", "CompleicationInfection")
tbl.complexmortality
chisq.test(tbl.complexmortality)

tbl.confirmmortality           = table(HaitiBlood_Rfile[, c(5, 4)])
colnames(tbl.confirmmortality) = c("Survive", "Dead")
rownames(tbl.confirmmortality) = c("NO", "Confirm")
tbl.confirmmortality
chisq.test(tbl.confirmmortality)


##########################################################################################
# Circle Plot and Heat Plot
##########################################################################################

Haiti.ATB.Trans_12 = matrix(0, nrow = 19, ncol = 19)
Haiti.ATB.Trans_13 = matrix(0, nrow = 19, ncol = 19)
Haiti.ATB.Trans_23 = matrix(0, nrow = 19, ncol = 19)

# First We need to clean the data again
bs_atb1 = data.entire$bs_atb1
bs_atb2 = data.entire$bs_atb2
bs_atb3 = data.entire$bs_atb3
temp.atb.Resist.frame = as.data.frame(cbind(bs_pat_id, bs_atb1, bs_atb2, bs_atb3))
# View(temp.atb.Resist.frame)
atb.Resist.frame = temp.atb.Resist.frame[which(complete.cases(temp.atb.Resist.frame$bs_atb1)), ]
atb.Resist.frame$bs_atb2[is.na(atb.Resist.frame$bs_atb2)] <- 0
atb.Resist.frame$bs_atb3[is.na(atb.Resist.frame$bs_atb3)] <- 0

# Write Markov Chain transit matrix
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:19)
  {
    for (j in 1:19)
    {
      for (k in 1:111)
      {
        if ((atb.Resist.frame$bs_atb1[k] == i)&&(atb.Resist.frame$bs_atb2[k] == j)) 
        {
          Haiti.ATB.Trans_12[i, j] = Haiti.ATB.Trans_12[i, j] + 1
        }
      }
    }
  }
  increase(pb)
  

}
print("Iteration Complete")
toc()

# View(Haiti.ATB.Trans_12)
m1 = Haiti.ATB.Trans_12
rownames(m1) = ATB.name
colnames(m1) = ATB.name
ATB.name <- c("Amikacine", "Amoxycilline", "Azythromycine", "Cefazoline", "Ceftazidime", "Ceftraxone", "Ciprofloxacine", "Clindamycine", "Cloxacilline", "Colistine", "Erythromycine", "FluconazolePO", "FluconazoleIV", "Gentamicin", "Imipenem", "Metronidazole", "PiperacillineTAZO", "Vancomycine", "Others")
dimnames(m1) <- list(have = ATB.name, prefer = ATB.name)
# set the colors
groupColors <- c("firebrick1", "aquamarine ", "gold1", "green4", "burlywood4 ", "darkorange2", "deeppink4 ", "grey1", "bisque3", "aquamarine4", "grey2", "mediumpurple4", "grey3", "lawngreen", "navy", "cornflowerblue", "cyan", "darkorchid1", "lightpink")
chordDiagram(m1, grid.col = groupColors, row.col = groupColors)
title("Phase I to Phase II")

# Phase 2 to Phase3
# Write Markov Chain transit matrix
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:19)
  {
    for (j in 1:19)
    {
      for (k in 1:111)
      {
        if ((atb.Resist.frame$bs_atb2[k] == i)&&(atb.Resist.frame$bs_atb3[k] == j)) 
        {
          Haiti.ATB.Trans_23[i, j] = Haiti.ATB.Trans_23[i, j] + 1
        }
      }
    }
  }
  increase(pb)
  

}
print("Iteration Complete")
toc()

# View(Haiti.ATB.Trans_23)
m2 = Haiti.ATB.Trans_23
rownames(m2) = ATB.name
colnames(m2) = ATB.name
ATB.name     <- c("Amikacine", "Amoxycilline", "Azythromycine", "Cefazoline", "Ceftazidime", "Ceftraxone", "Ciprofloxacine", "Clindamycine", "Cloxacilline", "Colistine", "Erythromycine", "FluconazolePO", "FluconazoleIV", "Gentamicin", "Imipenem", "Metronidazole", "PiperacillineTAZO", "Vancomycine", "Others")
dimnames(m2) <- list(have = ATB.name, prefer = ATB.name)
# set the colors
groupColors  <- c("firebrick1", "aquamarine ", "gold1", "cornflowerblue", "burlywood4 ", "darkorange2", "deeppink4 ", "grey1", "bisque3", "aquamarine4", "grey2", "mediumpurple4", "grey3", "lawngreen", "navy", "lightpink", "aquamarine", "darkorchid1", "lightpink")

chordDiagram(m2, grid.col = groupColors, row.col = groupColors)
title("Phase II to Phase III")

# Phase 1 to Phase3
# Write Markov Chain transit matrix
tic()
pb <- ProgressBar(max=42)
reset(pb)
while (!isDone(pb)) 
{
  for (i in 1:19)
  {
    for (j in 1:19)
    {
      for (k in 1:111)
      {
        if ((atb.Resist.frame$bs_atb1[k] == i)&&(atb.Resist.frame$bs_atb3[k] == j)) 
        {
          Haiti.ATB.Trans_13[i, j] = Haiti.ATB.Trans_13[i, j] + 1}
      }
    }
  }
  increase(pb)
  

}
print("Iteration Complete")
toc()

# View(Haiti.ATB.Trans_13)
m3 = Haiti.ATB.Trans_13
rownames(m3) = ATB.name
colnames(m3) = ATB.name
ATB.name     <- c("Amikacine", "Amoxycilline", "Azythromycine", "Cefazoline", "Ceftazidime", "Ceftraxone", "Ciprofloxacine", "Clindamycine", "Cloxacilline", "Colistine", "Erythromycine", "FluconazolePO", "FluconazoleIV", "Gentamicin", "Imipenem", "Metronidazole", "PiperacillineTAZO", "Vancomycine", "Others")
dimnames(m3) <- list(have = ATB.name, prefer = ATB.name)
# set the colors
groupColors  <- c("firebrick1", "aquamarine ", "gold1", "cornflowerblue", "burlywood4 ", "bisque3", "deeppink4 ", "grey1", "bisque3", "lightpink", "grey2", "mediumpurple4", "grey3", "lawngreen", "navy", "lightpink", "cyan", "plum2", "lightpink")
chordDiagram(m3, grid.col = groupColors, row.col = groupColors)
title("Phase I to Phase III")


# Heat PLot
Haiti.ATB.Heat           = matrix(0, nrow = 19, ncol = 3)
rownames(Haiti.ATB.Heat) = ATB.name
colnames(Haiti.ATB.Heat) = c("Phase_I", "Phase_II", "Phase_III")
for (i in 1:19)
{
  for (j in 1:3)
  {
    Haiti.ATB.Heat[i, j]  = sum(atb.Resist.frame[, j + 1] == i)
  }
}
# # View(Haiti.ATB.Heat)
margin.moi = list(
  l = 150, 
  r = 50, 
  b = 150, 
  t = 100, 
  pad = 4
)
# plot_ly(x = Haiti.ATB.Heat[, 1], y = Haiti.ATB.Heat[, 2], type = "histogram2d")
dimnames(Haiti.ATB.Heat) <- list(have = ATB.name, prefer = c("Phase_I", "Phase_II", "Phase_III"))
# Zoom after display
plot_ly(z = Haiti.ATB.Heat, x = c("Phase_I", "Phase_II", "Phase_III"), y = ATB.name, colorscale = "Hot", type = "heatmap")%>%
  layout(autosize = F, width = 1000, height = 750, margin = margin.moi)%>%
  layout(title = "Heat Plot of ATB")


#########################################################################################
# SVM Prediction of ATB Resistance
#########################################################################################

# Run this part before prediction:
bs_death            = as.factor(data.unique$bs_death)
bs_inj2hosp         = data.unique$bs_hr_bru2hosp
bs_TBSA             = data.unique$bs_surface
bs_age              = data.unique$bs_age
bs_thicken          = as.factor(data.unique$bs_prof)
bs_inhal            = as.factor(data.unique$bs_inhalation)
bs_lieu             = as.factor(data.unique$bs_lieu)
bs_envir            = as.factor(data.unique$bs_environ)
bs_maln             = as.factor(data.unique$bs_malnut)
bs_smoke            = as.factor(data.unique$bs_tabac)
bs_hiv              = as.factor(data.unique$bs_vih)
bs_db               = as.factor(data.unique$bs_diabete)
bs_exam             = as.factor(data.unique$bs_exam)
bs_gender           = as.factor(data.unique$bs_sexe)
bs_etio___1         = as.factor(data.unique$bs_etio___1)
bs_etio___2         = as.factor(data.unique$bs_etio___2)
bs_etio___3         = as.factor(data.unique$bs_etio___3)
bs_etio___4         = as.factor(data.unique$bs_etio___4)
bs_etio___5         = as.factor(data.unique$bs_etio___5)
bs_etio___6         = as.factor(data.unique$bs_etio___6)
bs_etio___7         = as.factor(data.unique$bs_etio___7)
bs_etio___9         = as.factor(data.unique$bs_etio___9)
bs_circonstance___1 = as.factor(data.unique$bs_circonstance___1)
bs_circonstance___2 = as.factor(data.unique$bs_circonstance___2)
bs_circonstance___3 = as.factor(data.unique$bs_circonstance___3)
bs_loc_bru___1      = as.factor(data.unique$bs_loc_bru___1)
bs_loc_bru___2      = as.factor(data.unique$bs_loc_bru___2)
bs_loc_bru___3      = as.factor(data.unique$bs_loc_bru___3)
bs_loc_bru___4      = as.factor(data.unique$bs_loc_bru___4)
bs_loc_bru___5      = as.factor(data.unique$bs_loc_bru___5)
bs_loc_bru___6      = as.factor(data.unique$bs_loc_bru___6)
bs_loc_bru___7      = as.factor(data.unique$bs_loc_bru___7)
bs_loc_bru___8      = as.factor(data.unique$bs_loc_bru___8)
bs_loc_bru___9      = as.factor(data.unique$bs_loc_bru___9)
bs_loc_bru___10     = as.factor(data.unique$bs_loc_bru___10)
bs_loc_bru___11     = as.factor(data.unique$bs_loc_bru___11)
bs_loc_bru___12     = as.factor(data.unique$bs_loc_bru___12)
bs_loc_bru___13     = as.factor(data.unique$bs_loc_bru___13)
bs_loc_bru___14     = as.factor(data.unique$bs_loc_bru___14)
bs_loc_bru___15     = as.factor(data.unique$bs_loc_bru___15)
bs_loc_bru___16     = as.factor(data.unique$bs_loc_bru___16)
bs_loc_bru___17     = as.factor(data.unique$bs_loc_bru___17)
bs_loc_bru___18     = as.factor(data.unique$bs_loc_bru___18)
bs_loc_bru___99     = as.factor(data.unique$bs_loc_bru___99)

bs_loc_risk___1     = as.factor(data.unique$bs_loc_risque___1)
bs_loc_risk___2     = as.factor(data.unique$bs_loc_risque___2)
bs_loc_risk___3     = as.factor(data.unique$bs_loc_risque___3)
bs_loc_risk___4     = as.factor(data.unique$bs_loc_risque___4)
bs_loc_risk___5     = as.factor(data.unique$bs_loc_risque___5)

bs_amikacine        = as.factor(data.unique$bs_amikacine)
bs_amoxyclav        = as.factor(data.unique$bs_amoxyclav)
bs_amoxy            = as.factor(data.unique$bs_amoxy)
bs_ampisulba        = as.factor(data.unique$bs_ampisulba)
bs_ampi             = as.factor(data.unique$bs_ampi)
bs_cefepime         = as.factor(data.unique$bs_cefepime)
bs_cefezidime       = as.factor(data.unique$bs_cefezidime)
bs_cefotaxime       = as.factor(data.unique$bs_cefotaxime)
bs_cefoxit          = as.factor(data.unique$bs_cefoxit)
bs_ceftazid         = as.factor(data.unique$bs_ceftazid)
bs_ceftri           = as.factor(data.unique$bs_ceftri)
bs_cipro            = as.factor(data.unique$bs_cipro)
bs_clinda           = as.factor(data.unique$bs_clinda)
bs_colistine        = as.factor(data.unique$bs_colistine)
bs_cotri            = as.factor(data.unique$bs_cotri)
bs_erythro          = as.factor(data.unique$bs_erythro)
bs_genta            = as.factor(data.unique$bs_genta)
bs_imipen           = as.factor(data.unique$bs_imipen)
bs_meropen          = as.factor(data.unique$bs_meropen)
bs_oxacil           = as.factor(data.unique$bs_oxacil)
bs_peni             = as.factor(data.unique$bs_peni)
bs_pipetazo         = as.factor(data.unique$bs_pipetazo)
bs_tetra            = as.factor(data.unique$bs_tetra)
bs_ticarclav        = as.factor(data.unique$bs_ticarclav)
bs_tobra            = as.factor(data.unique$bs_tobra)
bs_vanco            = as.factor(data.unique$bs_vanco)

bs_inj2hosp=scale(bs_inj2hosp)
bs_TBSA=scale(bs_TBSA)
bs_age=scale(bs_age)
bs_death[is.na(bs_death)] <- 0
bs_atb_R=as.factor(bs_atb_R)
Framework = as.data.frame(cbind(bs_atb_R, bs_death, bs_inj2hosp, bs_TBSA, bs_age, bs_thicken, bs_inhal, bs_lieu, bs_envir, bs_maln, bs_envir, bs_maln, bs_smoke, bs_hiv, bs_db, bs_db, bs_exam, bs_gender, bs_etio___1, bs_etio___2, bs_etio___3, bs_etio___4, bs_etio___5, bs_etio___6, bs_etio___7, bs_etio___9, bs_circonstance___1, bs_circonstance___2, bs_circonstance___3, bs_loc_bru___1, bs_loc_bru___2, bs_loc_bru___3, bs_loc_bru___4, bs_loc_bru___5, bs_loc_bru___6, bs_loc_bru___7, bs_loc_bru___8, bs_loc_bru___9, bs_loc_bru___10, bs_loc_bru___11, bs_loc_bru___12, bs_loc_bru___13, bs_loc_bru___14, bs_loc_bru___15, bs_loc_bru___16, bs_loc_bru___17, bs_loc_bru___18, bs_loc_bru___99, bs_loc_risk___1, bs_loc_risk___2, bs_loc_risk___3, bs_loc_risk___4, bs_loc_risk___5))

# So here's what we can do with this:
# Framework is basically everything in the first sheet 
# CAUTION:
# 
Framework.complete = na.omit(Framework)
pca1.EVAL          = prcomp(Framework.complete)
Rot.mat            = pca1.EVAL$rotation
x                  = pca1.EVAL$x[, 1]
y                  = pca1.EVAL$x[, 2]
z                  = pca1.EVAL$x[, 3]
bs_atb_R.c         = as.factor(Framework.complete$bs_atb_R)
fp                 = as.data.frame(cbind(x, y, z, bs_atb_R.c))
svm.report         <- svm(bs_atb_R.c~x + y + z, type = "C", kernel = "radial", data = fp)

# For the new observation, we extract the info from sheet1 into a vector
# Here is a sample
vec.sample         = rbind(runif(53, -2, 4), runif(53, -2, 4), runif(53, -2, 4), runif(53, -2, 4), runif(53, -2, 4))
vec.project.sample = as.data.frame(vec.sample%*%Rot.mat)
x                  = vec.project.sample[, 1]
y                  = vec.project.sample[, 2]
z                  = vec.project.sample[, 3]
new.data           = as.data.frame(cbind(x, y, z))
# svm.sample.class  <- predict(svm.report, fp[4, 1:3])

svm.sample.class   <- predict(svm.report, new.data[, 1:3])
# 1 2 3 4 5 
#  
# Levels: 0 1
plot(pca1.EVAL$x[, 1], pca1.EVAL$x[, 2], col = c("lawngreen", "lightpink1")[as.factor(fp$bs_atb_R)], pch = c(1, 19)[as.factor(fp$bs_atb_R)])
points(x,           y,                   col = c("turquoise3", "firebrick")[svm.sample.class],       pch = 8)
title("Sample SVM Results")
col.mis = c("lightpink1", "lawngreen", "firebrick", "turquoise3")
pch.mis = c(1, 1, 8, 8)
legend("topleft", inset = .05, title = "Misclassification", 
       c("Observed Positive", "Observed Negative", "Predicted Positive", "Predicted Negative"), col = col.mis, pch = pch.mis, horiz = FALSE)

#########################################################################################
# Neural Network
#########################################################################################
colnames(Framework.complete)[3]="bs_inj2hosp"
colnames(Framework.complete)[4]="bs_TBSA"
colnames(Framework.complete)[5]="bs_age"
n <- names(Framework.complete)
real.atb=Framework.complete$bs_atb_R

# Define Function
f <- as.formula(paste("bs_atb_R ~", paste(n[!n %in% "bs_atb_R"], collapse = " + ")))

# Two Layer neural network
nn2 <- neuralnet(f, data=Framework.complete, hidden=2, err.fct="ce",linear.output=FALSE,likelihood=TRUE)
print(nn2)
# Call: neuralnet(formula = f, data = Framework.complete, hidden = 2,     err.fct = "ce", linear.output = FALSE, likelihood = TRUE)
# 
# 1 repetition was calculated.
# 
# Error         AIC         BIC Reached Threshold Steps
# 1 7.507917539 221.0158351 625.7635168    0.007965057293    21
table(Framework.complete$bs_atb_R,nn2$response)
# 
#     0   1
# 0 361   0
# 1   0  15
plot(nn2, rep = "best")
importance.nn2=garson(nn2)[1]$data
plot(importance.nn2$rel_imp[-c(1:15)])
plot_ly(data = importance.nn2, x = importance.nn2$x_names, y = rel_imp, mode = "markers",color=rel_imp)%>%
  layout(autosize = F, width = 1500, height = 768, margin = margin.moi)%>%
  layout(title="NeuNet_Layer=2_Importance")

# $data
# rel_imp             x_names
# 41 0.005246136041     bs_loc_bru___16
# 48 0.006797279920     bs_loc_risk___4
# 27 0.006863785530      bs_loc_bru___2
# 5  0.008189253511          bs_thicken
# 1  0.009105507966            bs_death
# 3  0.009517456013             bs_TBSA
# 21 0.010942586299         bs_etio___7
# 15 0.011030758977         bs_etio___1
# 47 0.011968542030     bs_loc_risk___3
# 46 0.012025296649     bs_loc_risk___2
# 26 0.012030119593      bs_loc_bru___1
# 18 0.013100120654         bs_etio___4
# 17 0.013304545725         bs_etio___3
# 32 0.014090248878      bs_loc_bru___7
# 42 0.014110317057     bs_loc_bru___17
# 14 0.015554440935           bs_gender
# 16 0.015856871156         bs_etio___2
# 49 0.016965010238     bs_loc_risk___5
# 33 0.017575684262      bs_loc_bru___8
# 28 0.017581252936      bs_loc_bru___3
# 34 0.018865775551      bs_loc_bru___9
# 7  0.018980480518             bs_lieu
# 24 0.019738209242 bs_circonstance___2
# 12 0.019782816951               bs_db
# 40 0.019978004712     bs_loc_bru___15
# 38 0.020189476580     bs_loc_bru___13
# 37 0.020199022175     bs_loc_bru___12
# 4  0.020234056480              bs_age
# 11 0.020353725508              bs_hiv
# 23 0.020666387582 bs_circonstance___1
# 45 0.020770480197     bs_loc_risk___1
# 22 0.021047634712         bs_etio___9
# 19 0.021826537915         bs_etio___5
# 9  0.023643921755             bs_maln
# 43 0.023931564368     bs_loc_bru___18
# 30 0.025023909216      bs_loc_bru___5
# 6  0.025455585711            bs_inhal
# 29 0.026176917687      bs_loc_bru___4
# 31 0.026197629402      bs_loc_bru___6
# 44 0.026472250006     bs_loc_bru___99
# 8  0.027663024321            bs_envir
# 20 0.029371868481         bs_etio___6
# 39 0.029943635786     bs_loc_bru___14
# 25 0.031619568295 bs_circonstance___3
# 10 0.032681659136            bs_smoke
# 2  0.034453324012         bs_inj2hosp
# 13 0.035755025885             bs_exam
# 36 0.037366837569     bs_loc_bru___11
# 35 0.059755455880     bs_loc_bru___10
olden(nn2)
importance1.nn2=olden(nn2)[1]$data
plot(importance1.nn2$importance[-c(1:15)])
plot_ly(data = importance1.nn2, x = importance1.nn2$x_names, y = importance, mode = "markers",color=importance)%>%
  layout(autosize = F, width = 1500, height = 768, margin = margin.moi)%>%
  layout(title="NeuNet_Layer=2_Importance")


# Three Layer neural network
nn3 <- neuralnet(f, data=Framework.complete, hidden=3, err.fct="ce",linear.output=FALSE,likelihood=TRUE)
print(nn3)
# Call: neuralnet(formula = f, data = Framework.complete, hidden = 2,     err.fct = "ce", linear.output = FALSE, likelihood = TRUE)
# 
# 1 repetition was calculated.
# 
# Error         AIC         BIC Reached Threshold Steps
# 1 7.507301782 323.0146036 928.1713316    0.007342241211    20
table(nn3$response,Framework.complete$bs_atb_R)
# 
#     0   1
# 0 361   0
# 1   0  15
plot(nn3, rep = "best")
importance.nn3=garson(nn3)[1]$data
plot(importance.nn3$rel_imp[-c(1:15)])
plot_ly(data = importance.nn3, x = importance.nn3$x_names, y = rel_imp, mode = "markers",color=rel_imp)%>%
  layout(autosize = F, width = 1500, height = 768, margin = margin.moi)%>%
  layout(title="NeuNet_Layer=3_Importance")
olden(nn3)
importance1.nn3=olden(nn3)[1]$data
plot(importance1.nn3$importance[-c(1:15)])
plot_ly(data = importance1.nn3, x = importance1.nn3$x_names, y = importance, mode = "markers",color=importance)%>%
  layout(autosize = F, width = 1500, height = 768, margin = margin.moi)%>%
  layout(title="NeuNet_Layer=3_Importance")

# Method 2
nnn <- nnet(f, data=Framework.complete, size = 4, rang = 0.1,decay = 5e-4, maxit = 1000)
nnntest <- nnet(f, data=Framework.complete,size=8,linout=T)
traindatatest=Framework.complete[,-1]
nnn.predict=predict(nnn, traindatatest, type = "raw")
nnn.predict=nnn.predict[,1]
plot(nnn.predict,Framework.complete$bs_atb_R)
real.atb=Framework.complete$bs_atb_R
nnnfp=as.data.frame(cbind(nnn.predict,real.atb))

plot_ly(data=nnnfp, x = nnn.predict, y = real.atb, mode = "markers",color=nnn.predict)%>%
  layout(title="Neural Networks")

# Method 3
nnnn<-mlp(Framework.complete[,-1],Framework.complete[,1], data=Framework.complete, size = 8,linOut=T)

# PLOT Neural Network
# Write Function plot.nnet
plot.nnet <- function(mod.in,nid=T,all.out=T,all.in=T,bias=T,wts.only=F,rel.rsc=5,circle.cex=5,
                      node.labs=T,var.labs=T,x.lab=NULL,y.lab=NULL,line.stag=NULL,struct=NULL,cex.val=1,
                      alpha.val=1,circle.col='lightblue',pos.col='black',neg.col='grey', max.sp = F, ...){
  require(scales)
  
  # sanity checks
  if('mlp' %in% class(mod.in)) warning('Bias layer not applicable for rsnns object')
  if('numeric' %in% class(mod.in))
  {
    if(is.null(struct)) stop('Three-element vector required for struct')
    if(length(mod.in) != ((struct[1]*struct[2]+struct[2]*struct[3])+(struct[3]+struct[2])))
      stop('Incorrect length of weight matrix for given network structure')
  }
  if('train' %in% class(mod.in))
  {
    if('nnet' %in% class(mod.in$finalModel))
    {
      mod.in<-mod.in$finalModel
      warning('Using best nnet model from train output')
    }
    else stop('Only nnet method can be used with train object')
  }
  
  # gets weights for neural network, output is list
  # if rescaled argument is true, weights are returned but rescaled based on abs value
  nnet.vals<-function(mod.in,nid,rel.rsc,struct.out=struct)
  {
    require(scales)
    require(reshape)
    if('numeric' %in% class(mod.in))
    {
      struct.out<-struct
      wts<-mod.in
    }
    
    # neuralnet package
    if('nn' %in% class(mod.in))
    {
      struct.out<-unlist(lapply(mod.in$weights[[1]],ncol))
      struct.out<-struct.out[-length(struct.out)]
      struct.out<-c(
        length(mod.in$model.list$variables),
        struct.out,
        length(mod.in$model.list$response)
      )    		
      wts<-unlist(mod.in$weights[[1]])   
    }
    
    # nnet package
    if('nnet' %in% class(mod.in))
    {
      struct.out<-mod.in$n
      wts<-mod.in$wts
    }
    
    # RSNNS package
    if('mlp' %in% class(mod.in))
    {
      struct.out<-c(mod.in$nInputs,mod.in$archParams$size,mod.in$nOutputs)
      hid.num<-length(struct.out)-2
      wts<-mod.in$snnsObject$getCompleteWeightMatrix()
      
      # get all input-hidden and hidden-hidden wts
      inps<-wts[grep('Input',row.names(wts)),grep('Hidden_2',colnames(wts)),drop=F]
      inps<-melt(rbind(rep(NA,ncol(inps)),inps))$value
      uni.hids<-paste0('Hidden_',1+seq(1,hid.num))
      for(i in 1:length(uni.hids))
      {
        if(is.na(uni.hids[i+1])) break
        tmp<-wts[grep(uni.hids[i],rownames(wts)),grep(uni.hids[i+1],colnames(wts)),drop=F]
        inps<-c(inps,melt(rbind(rep(NA,ncol(tmp)),tmp))$value)
      }
      
      # get connections from last hidden to output layers
      outs<-wts[grep(paste0('Hidden_',hid.num+1),row.names(wts)),grep('Output',colnames(wts)),drop=F]
      outs<-rbind(rep(NA,ncol(outs)),outs)
      
      # weight vector for all
      wts<-c(inps,melt(outs)$value)
      assign('bias',F,envir=environment(nnet.vals))
    }
    
    if(nid) wts<-rescale(abs(wts),c(1,rel.rsc))
    
    # convert wts to list with appropriate names 
    hid.struct<-struct.out[-c(length(struct.out))]
    row.nms<-NULL
    for(i in 1:length(hid.struct))
    {
      if(is.na(hid.struct[i+1])) break
      row.nms<-c(row.nms,rep(paste('hidden',i,seq(1:hid.struct[i+1])),each=1+hid.struct[i]))
    }
    row.nms<-c(
      row.nms,
      rep(paste('out',seq(1:struct.out[length(struct.out)])),each=1+struct.out[length(struct.out)-1])
    )
    out.ls<-data.frame(wts,row.nms)
    out.ls$row.nms<-factor(row.nms,levels=unique(row.nms),labels=unique(row.nms))
    out.ls<-split(out.ls$wts,f=out.ls$row.nms)
    
    assign('struct',struct.out,envir=environment(nnet.vals))
    out.ls
  }
  
  wts<-nnet.vals(mod.in,nid=F)
  if(wts.only) return(wts)
  
  # circle colors for input, if desired, must be two-vector list, first vector is for input layer
  if(is.list(circle.col))
  {
    circle.col.inp<-circle.col[[1]]
    circle.col<-circle.col[[2]]
  }
  else circle.col.inp<-circle.col
  
  # initiate plotting
  x.range<-c(0,100)
  y.range<-c(0,100)
  
  # these are all proportions from 0-1
  if(is.null(line.stag)) line.stag<-0.011*circle.cex/2
  layer.x<-seq(0.17,0.9,length=length(struct))
  bias.x<-layer.x[-length(layer.x)]+diff(layer.x)/2
  bias.y<-0.95
  circle.cex<-circle.cex
  
  # get variable names from mod.in object
  # change to user input if supplied
  if('numeric' %in% class(mod.in))
  {
    x.names<-paste0(rep('X',struct[1]),seq(1:struct[1]))
    y.names<-paste0(rep('Y',struct[3]),seq(1:struct[3]))
  }
  if('mlp' %in% class(mod.in))
  {
    all.names<-mod.in$snnsObject$getUnitDefinitions()
    x.names<-all.names[grep('Input',all.names$unitName),'unitName']
    y.names<-all.names[grep('Output',all.names$unitName),'unitName']
  }
  if('nn' %in% class(mod.in))
  {
    x.names<-mod.in$model.list$variables
    y.names<-mod.in$model.list$respons
  }
  if('xNames' %in% names(mod.in))
  {
    x.names<-mod.in$xNames
    y.names<-attr(terms(mod.in),'factor')
    y.names<-row.names(y.names)[!row.names(y.names) %in% x.names]
  }
  if(!'xNames' %in% names(mod.in) & 'nnet' %in% class(mod.in))
  {
    if(is.null(mod.in$call$formula))
    {
      x.names<-colnames(eval(mod.in$call$x))
      y.names<-colnames(eval(mod.in$call$y))
    }
    else
    {
      forms<-eval(mod.in$call$formula)
      x.names<-mod.in$coefnames
      facts<-attr(terms(mod.in),'factors')
      y.check<-mod.in$fitted
      if(ncol(y.check)>1) y.names<-colnames(y.check)
      else y.names<-as.character(forms)[2]
    } 
  }
  # change variables names to user sub 
  if(!is.null(x.lab))
  {
    if(length(x.names) != length(x.lab)) stop('x.lab length not equal to number of input variables')
    else x.names<-x.lab
  }
  if(!is.null(y.lab))
  {
    if(length(y.names) != length(y.lab)) stop('y.lab length not equal to number of output variables')
    else y.names<-y.lab
  }
  
  # initiate plot
  plot(x.range,y.range,type='n',axes=F,ylab='',xlab='',...)
  
  # function for getting y locations for input, hidden, output layers
  # input is integer value from 'struct'
  get.ys<-function(lyr, max_space = max.sp)
  {
    if(max_space)
    { 
      spacing <- diff(c(0*diff(y.range),0.9*diff(y.range)))/lyr
    } 
    else 
    {
      spacing<-diff(c(0*diff(y.range),0.9*diff(y.range)))/max(struct)
    }
    
    seq(0.5*(diff(y.range)+spacing*(lyr-1)),0.5*(diff(y.range)-spacing*(lyr-1)),
        length=lyr)
  }
  
  # function for plotting nodes
  # 'layer' specifies which layer, integer from 'struct'
  # 'x.loc' indicates x location for layer, integer from 'layer.x'
  # 'layer.name' is string indicating text to put in node
  layer.points<-function(layer,x.loc,layer.name,cex=cex.val)
  {
    x<-rep(x.loc*diff(x.range),layer)
    y<-get.ys(layer)
    points(x,y,pch=21,cex=circle.cex,col=in.col,bg=bord.col)
    if(node.labs) text(x,y,paste(layer.name,1:layer,sep=''),cex=cex.val)
    if(layer.name=='I' & var.labs) text(x-line.stag*diff(x.range),y,x.names,pos=2,cex=cex.val)      
    if(layer.name=='O' & var.labs) text(x+line.stag*diff(x.range),y,y.names,pos=4,cex=cex.val)
  }
  
  # function for plotting bias points
  # 'bias.x' is vector of values for x locations
  # 'bias.y' is vector for y location
  # 'layer.name' is  string indicating text to put in node
  bias.points<-function(bias.x,bias.y,layer.name,cex,...)
  {
    for(val in 1:length(bias.x))
    {
      points(
        diff(x.range)*bias.x[val],
        bias.y*diff(y.range),
        pch=21,col=in.col,bg=bord.col,cex=circle.cex
      )
      if(node.labs)
        text(
          diff(x.range)*bias.x[val],
          bias.y*diff(y.range),
          paste(layer.name,val,sep=''),
          cex=cex.val
        )
    }
  }
  
  # function creates lines colored by direction and width as proportion of magnitude
  # use 'all.in' argument if you want to plot connection lines for only a single input node
  layer.lines<-function(mod.in,h.layer,layer1=1,layer2=2,out.layer=F,nid,rel.rsc,all.in,pos.col,
                        neg.col,...)
  {
    x0<-rep(layer.x[layer1]*diff(x.range)+line.stag*diff(x.range),struct[layer1])
    x1<-rep(layer.x[layer2]*diff(x.range)-line.stag*diff(x.range),struct[layer1])
    if(out.layer==T)
    {
      y0<-get.ys(struct[layer1])
      y1<-rep(get.ys(struct[layer2])[h.layer],struct[layer1])
      src.str<-paste('out',h.layer)
      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts<-wts[grep(src.str,names(wts))][[1]][-1]
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      wts.rs<-wts.rs[grep(src.str,names(wts.rs))][[1]][-1]
      cols<-rep(pos.col,struct[layer1])
      cols[wts<0]<-neg.col
      if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
      else segments(x0,y0,x1,y1)
    }
    else
    {
      if(is.logical(all.in)) all.in<-h.layer
      else all.in<-which(x.names==all.in)
      y0<-rep(get.ys(struct[layer1])[all.in],struct[2])
      y1<-get.ys(struct[layer2])
      src.str<-paste('hidden',layer1)
      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts<-unlist(lapply(wts[grep(src.str,names(wts))],function(x) x[all.in+1]))
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      wts.rs<-unlist(lapply(wts.rs[grep(src.str,names(wts.rs))],function(x) x[all.in+1]))
      cols<-rep(pos.col,struct[layer2])
      cols[wts<0]<-neg.col
      if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
      else segments(x0,y0,x1,y1)
    }
  }
  bias.lines<-function(bias.x,mod.in,nid,rel.rsc,all.out,pos.col,neg.col,...)
  {
    if(is.logical(all.out)) all.out<-1:struct[length(struct)]
    else all.out<-which(y.names==all.out)
    for(val in 1:length(bias.x))
    {
      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      if(val != length(bias.x)){
        wts<-wts[grep('out',names(wts),invert=T)]
        wts.rs<-wts.rs[grep('out',names(wts.rs),invert=T)]
        sel.val<-grep(val,substr(names(wts.rs),8,8))
        wts<-wts[sel.val]
        wts.rs<-wts.rs[sel.val]
      }
      else
      {
        wts<-wts[grep('out',names(wts))]
        wts.rs<-wts.rs[grep('out',names(wts.rs))]
      }
      cols<-rep(pos.col,length(wts))
      cols[unlist(lapply(wts,function(x) x[1]))<0]<-neg.col
      wts.rs<-unlist(lapply(wts.rs,function(x) x[1]))
      if(nid==F)
      {
        wts.rs<-rep(1,struct[val+1])
        cols<-rep('black',struct[val+1])
      }
      if(val != length(bias.x))
      {
        segments(
          rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
          rep(bias.y*diff(y.range),struct[val+1]),
          rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
          get.ys(struct[val+1]),
          lwd=wts.rs,
          col=cols
        )
      }
      else
      {
        segments(
          rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
          rep(bias.y*diff(y.range),struct[val+1]),
          rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
          get.ys(struct[val+1])[all.out],
          lwd=wts.rs[all.out],
          col=cols[all.out]
        )
      }
    }
  }
  
  # use functions to plot connections between layers
  # bias lines
  if(bias) bias.lines(bias.x,mod.in,nid=nid,rel.rsc=rel.rsc,all.out=all.out,pos.col=alpha(pos.col,alpha.val),
                      neg.col=alpha(neg.col,alpha.val))
  
  # layer lines, makes use of arguments to plot all or for individual layers
  # starts with input-hidden
  # uses 'all.in' argument to plot connection lines for all input nodes or a single node
  if(is.logical(all.in))
  {  
    mapply(
      function(x) layer.lines(mod.in,x,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,
                              all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
      1:struct[1]
    )
  }
  else
  {
    node.in<-which(x.names==all.in)
    layer.lines(mod.in,node.in,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,all.in=all.in,
                pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
  }
  # connections between hidden layers
  lays<-split(c(1,rep(2:(length(struct)-1),each=2),length(struct)),
              f=rep(1:(length(struct)-1),each=2))
  lays<-lays[-c(1,(length(struct)-1))]
  for(lay in lays)
  {
    for(node in 1:struct[lay[1]])
    {
      layer.lines(mod.in,node,layer1=lay[1],layer2=lay[2],nid=nid,rel.rsc=rel.rsc,all.in=T,
                  pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
    }
  }
  # lines for hidden-output
  # uses 'all.out' argument to plot connection lines for all output nodes or a single node
  if(is.logical(all.out))
    mapply(
      function(x) layer.lines(mod.in,x,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
                              all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
      1:struct[length(struct)]
    )
  else
  {
    node.in<-which(y.names==all.out)
    layer.lines(mod.in,node.in,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
                pos.col=pos.col,neg.col=neg.col,all.out=all.out)
  }
  
  # use functions to plot nodes
  for(i in 1:length(struct))
  {
    in.col<-bord.col<-circle.col
    layer.name<-'H'
    if(i==1) { layer.name<-'I'; in.col<-bord.col<-circle.col.inp}
    if(i==length(struct)) layer.name<-'O'
    layer.points(struct[i],layer.x[i],layer.name)
  }
  
  if(bias) bias.points(bias.x,bias.y,'B')
  
}
# Reference
# Dr. Marcus. W. Beck, USEPA, https://github.com/fawda123
# https://beckmw.wordpress.com/


# Plot
plot.nnet(nnnn)
plot.nnet(nn2)
plot.nnet(nn3)
plot.nnet(nnn)

# ROC
plotROC(nn2$response,Framework.complete$bs_atb_R)
title("ROC")
plotROC(nn3$response,Framework.complete$bs_atb_R)
title("ROC")
plotROC(nnn.predict,Framework.complete$bs_atb_R)
title("ROC")
plotROC(nnnn$fitted.values,Framework.complete$bs_atb_R)
title("ROC")



#########################################################################################
#########################################################################################
# Date  = 20160627
#########################################################################################
#########################################################################################


# Okay, try some survival analysis

library(survival)
library(KMsurv)
library(nlme)
library(km.ci)
library(car)
library(mfp)
