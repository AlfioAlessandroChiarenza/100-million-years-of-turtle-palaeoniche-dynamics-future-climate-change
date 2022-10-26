#####################################
## ENSEMBLE SCRIPT  for Turtle ENM ##
#####################################
require(raster) 
require(biomod2)
# load our species data
# Example with Testudinidae data 
DataSpecies <- read.csv(".\\Testudinidae.csv",h=T,sep=",", stringsAsFactors=T) #replace the dot (.) with your preferred directory pathway
head(DataSpecies)

# the name of studied species
myRespName <- 'Testudinidae'

# the presence/absences data for our species 
#
myResp <- as.numeric(DataSpecies[,"Species"])
# the xy coordinates of species data
myRespXY <- DataSpecies[,c("x","y")]

# load the environmental raster layers (for training your ENM)
# Environmental variables imported from your data

myExpl<-stack(list.files(path=".\\layers",pattern='asc',full.names=TRUE )) #replace the dot (.) with your preferred directory pathway


myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
                                     expl.var = myExpl, 
                                     resp.xy = myRespXY, 
                                     resp.name = myRespName,
                                     PA.nb.absences = 1000,
                                     PA.nb.rep = 1,
                                     PA.strategy = 'random',
                                     na.rm = TRUE)

# 2. Defining MAXENT Mododelling options 
#?BIOMOD_ModelingOptions
myBiomodOption <- BIOMOD_ModelingOptions( 
  MAXENT.Phillips = list(path_to_maxent.jar='.\\maxent\\maxent.jar', #replace the dot (.) with your preferred directory pathway (where you keep your maxent files)
                         maximumiterations = 5000, 
                         visible = FALSE, 
                         linear = TRUE,
                         quadratic = TRUE, 
                         product = TRUE, 
                         threshold = TRUE, 
                         hinge = TRUE, 
                         lq2lqptthreshold = 70,
                         l2lqthreshold = 10, 
                         hingethreshold = 15, 
                         beta_threshold = -1, 
                         beta_categorical = -1, 
                         beta_lqp = -1, 
                         beta_hinge = -1, 
                         defaultprevalence = 0.5,   #it was 0.5
                         betamultiplier=2), ### numeric (default 1), multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution. #Check literature on beta multiplier, since a >1 (default) may be recommended for projections (avoid overfitting). It could be from 1 to 20.
  RF = list (ntree=1000))                   #2 is preferred here: see Chiarenza et al. 2020 (https://doi.org/10.1073/pnas.200608711) and Jones et al. 2022 (https://doi.org/10.1038/s41467-022-30793-8) 


# 3. Setting ensemble Modelisation 
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                     models = c('SRE','RF','MAXENT.Phillips'),
                                     models.options = myBiomodOption,
                                     NbRunEval=50,      #set number of replication, here 50 to be consistent with experiments in MaxEnt.java
                                     DataSplit=70,
                                     Prevalence=0.5,
                                     VarImport=3,     #unceratin whether to keep like that
                                     models.eval.meth = c('TSS','ROC','KAPPA'),
                                     SaveObj = TRUE,
                                     do.full.models = TRUE,
                                     rescal.all.models=T,
                                     modeling.id='test')

# files created on hard drive  # adapt to your new path
list.files(".\\ENMmodel") #replace the dot (.) with your preferred directory pathway

# save models evaluation scores and variables importance on hard drive
capture.output(get_evaluations(myBiomodModelOut),
               file=file.path(myRespName, 
                              paste(myRespName,"_formal_models_evaluation.txt", sep="")))

capture.output(get_variables_importance(myBiomodModelOut),
               file=file.path(myRespName, 
                              paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               

# print variable importances                                    
get_variables_importance(myBiomodModelOut)

# 4.1 Projection on TRAINING environemental conditions     #################  WE USED THIS
myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = myExpl,
                                        proj.name = 'tdlae7_sed',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',                               #you can add 'ROC' here to get a further threshold evaluation metric to use
                                        compress = FALSE,
                                        build.clamping.mask = TRUE)

# summary of created object
myBiomodProjection

#Plot projection in other scenario
#plot(myBiomodProjection)           #too heavy for many replications, use instead next lines

test <- get_predictions(myBiomodProjection)
plot(test)                                            #plot all models
plot(test$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(test$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(test$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection

# files created on hard drive
list.files(".Testudinidae_Modern")

################################################### ENSEMBLE
###  6: ensemble modeling

myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by = 'all',
  eval.metric = c('TSS'),                   #you can add 'ROC' here to get a further threshold evaluation metric to use
  eval.metric.quality.threshold = NULL,
  models.eval.meth = c('TSS'),             #you can add 'ROC' here to get a further threshold evaluation metric to use
  prob.mean = TRUE,
  prob.cv = FALSE,
  prob.ci = FALSE,
  prob.ci.alpha = 0.05,
  prob.median = FALSE,
  committee.averaging = FALSE,
  prob.mean.weight = FALSE,
  prob.mean.weight.decay = 'proportional' )


# print summary               

myBiomodEM

# files created on hard drive
list.files(".Testudinidae_Modern")


# get evaluation scores
get_evaluations(myBiomodEM)


##########################################################################
###  7: EnsembleForecasting_current--------->PROJECT TO TRAINING REGION  #

CurrentEnsemble<-BIOMOD_EnsembleForecasting( projection.output = myBiomodProjection,
                                             EM.output = myBiomodEM,
                                             new.env = NULL,
                                             xy.new.env = NULL,
                                             selected.models = 'all',
                                             proj.name = 'tdlae7_sed', #'tdlae7_sed' is the name of the current climatic scenario 
                                             binary.meth = c('TSS'))    #you can add 'ROC' here to get a further threshold evaluation metric to use

# print summary
windows()
plot(CurrentEnsemble)

###################################################
###  8: projection to future conditions RCP 2.6

# tdlai7_sed SCENARIO 
# load environmental variables for the future. 
layers.tdlai7_sed<-stack(list.files(path=".tdlai7_sed\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.tdlai7_sed)
# Select variables to include in the model
names(layers.tdlai7_sed)
#myExpltdlai7_sed<-layers.tdlai7_sed[[c("cmm","drymon","wmm","wetmon")]]
#
myExpltdlai7_sed<-layers.tdlai7_sed[[c("cmm","drymon","wmm","wetmon")]]

# Plot environmental variables on geographical space
plot(myExpltdlai7_sed)

#Project in PROJECTION scenario
myBiomodProjtdlai7_sed <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                            new.env = myExpltdlai7_sed,
                                            proj.name = 'tdlai7_sed',     #'tdlai7_sed' is the name of the future  RCP 2.6 climatic scenario 
                                            selected.models = 'all',
                                            binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                            compress = FALSE,
                                            build.clamping.mask = TRUE)

#plot(myBiomodProjtdlai7_sed)         #don't use with many replications, use lines below

test <- get_predictions(myBiomodProjtdlai7_sed)
plot(test)                                            #plot all models
plot(test$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(test$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(test$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjtdlai7_sed

# files created on hard drive
list.files(".\\Testudinidae\\Mod2F26")

#project ensemble in future scenario

myBiomodProjtdlai7_sedEM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                      EM.output = myBiomodEM,
                                                      new.env = myExpltdlai7_sed,
                                                      xy.new.env = NULL,
                                                      selected.models = 'all',
                                                      proj.name = 'tdlai7_sed',
                                                      binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjtdlai7_sedEM)

# print summary
myBiomodProjtdlai7_sedEM

plot(CurrentEnsemble)  #FOR COMPARISON

## End (tdlai7_sed)

###################################################
###  9: projection to future conditions RCP 4.5

# tdlah7_sed SCENARIO 
# load environmental variables for the future. 
layers.tdlah7_sed<-stack(list.files(path=".\\tdlah7_sed\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.tdlah7_sed)
# Select variables to include in the model
names(layers.tdlah7_sed)
myExpltdlah7_sed<-layers.tdlah7_sed[[c("cmm","drymon","wmm","wetmon")]]
# Plot environmental variables on geographical space
plot(myExpltdlah7_sed)

#Project in PROJECTION scenario
myBiomodProjtdlah7_sed <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                            new.env = myExpltdlah7_sed,
                                            proj.name = 'tdlah7_sed',     #'tdlah7_sed' is the name of the future  RCP 4.5 climatic scenario 
                                            selected.models = 'all',
                                            binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                            compress = FALSE,
                                            build.clamping.mask = TRUE)

#plot(myBiomodProjtdlah7_sed)         #don't use with many replications, use lines below

tdlah7_sed <- get_predictions(myBiomodProjtdlah7_sed)
plot(tdlah7_sed)                                            #plot all models
plot(tdlah7_sed$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(tdlah7_sed$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(tdlah7_sed$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjtdlah7_sed

# files created on hard drive
list.files(".\\Testudinidae\\Mod2F45")

#project ensemble in future scenario

myBiomodProjtdlah7_sedEM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                      EM.output = myBiomodEM,
                                                      new.env = myExpltdlah7_sed,
                                                      xy.new.env = NULL,
                                                      selected.models = 'all',
                                                      proj.name = 'tdlah7_sed',
                                                      binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjtdlah7_sedEM)

# print summary
myBiomodProjtdlah7_sedEM

plot(CurrentEnsemble)  #FOR COMPARISON

## End tdlah7_sed

###################################################
###  9: projection to future conditions RCP 6.0

# tdlag7 SCENARIO 
# load environmental variables for the future. 
layers.tdlag7_sed<-stack(list.files(path=".\\tdlag7_sed\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.tdlag7_sed)
# Select variables to include in the model
names(layers.tdlag7_sed)
myExpltdlag7<-layers.tdlag7_sed[[c("cmm","drymon","wmm","wetmon")]]
# Plot environmental variables on geographical space
plot(myExpltdlag7_sed)

#Project in PROJECTION scenario
myBiomodProjtdlag7 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = myExpltdlag7_sed,
                                        proj.name = 'tdlag7_sed',     #'tdlai7_sed' is the name of the future  RCP 6.0 climatic scenario 
                                        selected.models = 'all',
                                        binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                        compress = FALSE,
                                        build.clamping.mask = TRUE)

#plot(myBiomodProjtdlag7)         #don't use with many replications, use lines below

tdlag7 <- get_predictions(myBiomodProjtdlag7)
plot(tdlag7)                                            #plot all models
plot(tdlag7$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(tdlag7$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(tdlag7$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjtdlag7

# files created on hard drive
list.files(".\\Testudinidae\\Mod2F60")

#project ensemble in future scenario

myBiomodProjtdlag7EM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                  EM.output = myBiomodEM,
                                                  new.env = myExpltdlag7,
                                                  xy.new.env = NULL,
                                                  selected.models = 'all',
                                                  proj.name = 'tdlag7_sed',
                                                  binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjtdlag7EM)

# print summary
myBiomodProjtdlag7EM

plot(CurrentEnsemble)  #FOR COMPARISON

## End tesum

###################################################
###  10: projection to future conditions RCP 8.5

# tdlac7 SCENARIO 
# load environmental variables for the future. 

layers.tdlac7_sed<-stack(list.files(path=".\\tdlac7_sed\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.tdlac7_sed)
# Select variables to include in the model
names(layers.tdlac7_sed)
myExpltdlac7<-layers.tdlac7_sed[[c("cmm","drymon","wmm","wetmon")]]
# Plot environmental variables on geographical space
plot(myExpltdlac7_sed)

#Project in PROJECTION scenario
myBiomodProjtdlac7 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = myExpltdlac7_sed,
                                        proj.name = 'tdlac7_sed',     #'tdlac7_sed' is the name of the future  RCP 8.5 climatic scenario 
                                        selected.models = 'all',
                                        binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                        compress = FALSE,
                                        build.clamping.mask = TRUE)

#plot(myBiomodProjtdlac7)         #don't use with many replications, use lines below

tdlac7 <- get_predictions(myBiomodProjtdlac7)
plot(tdlac7)                                            #plot all models
plot(tdlac7$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(tdlac7$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(tdlac7$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjtdlac7

# files created on hard drive
list.files(".Testudinidae\\Mod2F85")

#project ensemble in future scenario

myBiomodProjtdlac7EM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                  EM.output = myBiomodEM,
                                                  new.env = myExpltdlac7,
                                                  xy.new.env = NULL,
                                                  selected.models = 'all',
                                                  proj.name = 'tdlac7_sed',
                                                  binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjtdlac7EM)

# print summary
myBiomodProjtdlac7EM

plot(CurrentEnsemble)  #FOR COMPARISON

## End tdlac7

###################################################
### 11: projection to BARTONIAN

# teydi1_sed_180 SCENARIO 
# load environmental variables for the future. 

layers.teydi1_sed_180<-stack(list.files(path=".teydi1_sed_180\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.teydi1_sed_180)
# Select variables to include in the model
names(layers.teydi1_sed_180)
myExplteydi1_sed_180<-layers.teydi1_sed_180[[c("cmm","drymon","wmm","wetmon")]]
# Plot environmental variables on geographical space
plot(myExplteydi1_sed_180)

#Project in PROJECTION scenario
myBiomodProjteydi1_sed_180 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                new.env = myExplteydi1_sed_180,
                                                proj.name = 'teydi1_sed_180', #'teydi1_sed_180' is the name of the late Eocene climatic scenario 
                                                selected.models = 'all',
                                                binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                                compress = FALSE,
                                                build.clamping.mask = TRUE)

#plot(myBiomodProjteydi1_sed_180)         #don't use with many replications, use lines below

teydi1_sed_180 <- get_predictions(myBiomodProjteydi1_sed_180)
plot(teydi1_sed_180)                                            #plot all models
plot(teydi1_sed_180$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(teydi1_sed_180$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(teydi1_sed_180$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjteydi1_sed_180

# files created on hard drive
list.files(".Testudinidae\\Mod2BP")

#project ensemble in future scenario

myBiomodProjteydi1_sed_180EM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                          EM.output = myBiomodEM,
                                                          new.env = myExplteydi1_sed_180,
                                                          xy.new.env = NULL,
                                                          selected.models = 'all',
                                                          proj.name = 'teydi1_sed_180',
                                                          binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjteydi1_sed_180EM)

# print summary
myBiomodProjteydi1_sed_180EM

plot(CurrentEnsemble)  #FOR COMPARISON

## End teydi1_sed_180


###################################################
### 12: projection to MAASTRICHTIAN

# teydo1_sed_180 SCENARIO 
# load environmental variables for the future. 

layers.teydo1_sed_180<-stack(list.files(path=".teydo1_sed_180\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.teydo1_sed_180)
# Select variables to include in the model
names(layers.teydo1_sed_180)
myExplteydo1_sed_180<-layers.teydo1_sed_180[[c("cmm","drymon","wmm","wetmon")]]
# Plot environmental variables on geographical space
plot(myExplteydo1_sed_180)

#Project in PROJECTION scenario
myBiomodProjteydo1_sed_180 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                new.env = myExplteydo1_sed_180,
                                                proj.name = 'teydo1_sed_180', #'teydo1_sed_180' is the name of the Maastrichtian climatic scenario
                                                selected.models = 'all',
                                                binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                                compress = FALSE,
                                                build.clamping.mask = TRUE)

#plot(myBiomodProjteydo1_sed_180)         #don't use with many replications, use lines below

teydo1_sed_180 <- get_predictions(myBiomodProjteydo1_sed_180)
plot(teydo1_sed_180)                                            #plot all models
plot(teydo1_sed_180$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(teydo1_sed_180$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(teydo1_sed_180$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjteydo1_sed_180

# files created on hard drive
list.files(".\\Testudinidae\\Mod2Maa")

#project ensemble in future scenario

myBiomodProjteydo1_sed_180EM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                          EM.output = myBiomodEM,
                                                          new.env = myExplteydo1_sed_180,
                                                          xy.new.env = NULL,
                                                          selected.models = 'all',
                                                          proj.name = 'teydo1_sed_180',
                                                          binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjteydo1_sed_180EM)

# print summary
myBiomodProjteydo1_sed_180EM

plot(CurrentEnsemble)  #FOR COMPARISON

## End teydo1_sed_180


###################################################
###  13: projection to TURONIAN

# teyds1_sed_180 SCENARIO 
# load environmental variables for the future. 

layers.teyds1_sed_180<-stack(list.files(path=".\\teyds1_sed_180\\layers",pattern='asc',full.names=TRUE ))

# Descriptive statistics for environmental layers
summary(layers.teyds1_sed_180)
# Select variables to include in the model
names(layers.teyds1_sed_180)
myExplteyds1_sed_180<-layers.teyds1_sed_180[[c("cmm","drymon","wmm","wetmon")]]
# Plot environmental variables on geographical space
plot(myExplteyds1_sed_180)

#Project in PROJECTION scenario
myBiomodProjteyds1_sed_180 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                new.env = myExplteyds1_sed_180,
                                                proj.name = 'teyds1_sed_180', #'teyds1_sed_180' is the name of the Turonian-Coniacian-Santonian climatic scenario
                                                selected.models = 'all',
                                                binary.meth = 'TSS',          #you can add 'ROC' here to get a further threshold evaluation metric to use
                                                compress = FALSE,
                                                build.clamping.mask = TRUE)

#plot(myBiomodProjteyds1_sed_180)         #don't use with many replications, use lines below

teyds1_sed_180 <- get_predictions(myBiomodProjteyds1_sed_180)
plot(teyds1_sed_180)                                            #plot all models
plot(teyds1_sed_180$Testudinidae_PA1_RUN50_SRE)                 #plot a single SRE projection
plot(teyds1_sed_180$Testudinidae_PA1_RUN50_RF)                  #plot a single RF projection
plot(teyds1_sed_180$Testudinidae_PA1_RUN50_MAXENT.Phillips)     #plot a single MAXENT projection


# print summary
myBiomodProjteyds1_sed_180

# files created on hard drive
list.files(".\\Testudinidae\\Mod2TCS")

#project ensemble in future scenario

myBiomodProjteyds1_sed_180EM <-BIOMOD_EnsembleForecasting(projection.output = NULL,
                                                          EM.output = myBiomodEM,
                                                          new.env = myExplteyds1_sed_180,
                                                          xy.new.env = NULL,
                                                          selected.models = 'all',
                                                          proj.name = 'teyds1_sed_180',
                                                          binary.meth = c('TSS'))         #you can add 'ROC' here to get a further threshold evaluation metric to use
windows()
plot(myBiomodProjteyds1_sed_180EM)

# print summary
myBiomodProjteyds1_sed_180EM

plot(CurrentEnsemble)  #FOR COMPARISON

## End teyds1_sed_180
######################################### END of SIMULATIONS