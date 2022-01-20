#######################################################################################
#
#                                                                                     
#   Filename    :	master.R    												  
#                                                                                     
#   Project     :       BiomJ article "Optimal dynamic treatment regime estimation using information extraction from unstructured clinical text"                                                             
#   Authors     :       Nina Zhou, Robert Brook, Ivo Dinov, and Lu Wang                                                                
#   Date        :       11/28/2021
#   Purpose     :       As described in BiomJ article
#																				  
#   R Version   :       4.0.2 (2020-06-22)                                                                
#   RStudio Version    :  	1.3.959
#   Command   :  	R CMD BATCH --no-save master.R master_log.txt  
#
#   Input intermediate results for simulation section :    two_stage_sim.Rdata (objects aa,bb), three_stage_sim.Rdata (objects cc,dd)
#   Input data for application section  :    application_data_day1.csv, application_data_day2.csv
#   Output data files :    master_log.txt (console printouts, includes Table 1), Rplots.pdf (Figure 1);
#                          stage1_dtr.png, stage2_dtr.png (Figure 2 stage 1 and 2 application section), 
#                          temp.html (temporary file to render plots)
#
#   Required R packages :  dplyr, ggplot2, randomForest, missForest, latex2exp, rpart, webshot, htmlwidgets, data.tree, DiagrammeR
#
#   Note  :  requires pandoc to save figures. Installing pandoc at: pandoc.org/installing.html 
#
#
########################################################################################

#Load functions
library(dplyr)
library(ggplot2)
library(randomForest)
library(missForest)
library(latex2exp)
library(rpart)
#setwd() #set working directory to the downloaded Code_and_data folder.
source("../function/Functions.R")

####Simulation Section Code####

#load the text data
text_data <- read.csv("../data/data_simulation/sim_data.csv")

#Check if extracted and real smoking status are matching.
table(text_data$true_smk_status == text_data$smk_extract)#100% matching

#Simulate the weight variable according to extracted info

#summary(text_data$weightlb) 
#42.24lb-228.80lb
#sd(text_data$weightlb, na.rm = TRUE)
#mean = 194.73 sd = 51.36

set.seed(123)
text_data$wt_full <-
  sapply(text_data$weightlb, function(x)
    ifelse(is.na(x),
           rnorm(1, mean = 194.73, sd = 51.36), x))
text_data$wt_full[which(text_data$wt_full <= 0)] = text_data$wt_full[
  which(text_data$wt_full <=0)] + 42.24

# Start simulation for section 3 Simulation.

####################################################
#### simulation 2 Stages 3 Treatments per stage
####################################################

rown = c(
  "Case 1 with text",
  "Case 1 without text",
  "Case 2 with text complete",
  "Case 2 with text imputed",
  "Case 2 without text complete",
  "Case 2 without text imputed"
)

#Load intermediate result aa,bb
load("../data/data_simulation/two_stage_sim.Rdata")

## Avoid running the following code after loading intermediate result
## Can reduce 1:1000 to 1:10 for faster run time
## n = 500

# ptm <- proc.time()
# aa <- sapply(1:1000, function(x)
#   sim_text_2stage(x, 500, text_data))
# proc.time() - ptm #4.68 hours

## n = 1000

# ptm <- proc.time()
# bb <- sapply(1:1000, function(x)
#   sim_text_2stage(x, 1000, text_data))
# proc.time() - ptm #6.14 hours
## Avoid running the code above after loading intermediate result

# Print tables
n500_2stage <- data.frame(matrix(c(
  rowMeans(aa)[1:6],
  apply(aa, 1, sd)[1:6],
  rowMeans(aa)[7:12],
  apply(aa, 1, sd)[7:12]
), ncol = 4))
colnames(n500_2stage) <- c("opt%", "sd", "EY(g^opt)", "sd")
rownames(n500_2stage) <- rown
print("Table 1 - two stage study n = 500")
print(n500_2stage)

n1000_2stage <- data.frame(matrix(c(
  rowMeans(bb)[1:6],
  apply(bb, 1, sd)[1:6],
  rowMeans(bb)[7:12],
  apply(bb, 1, sd)[7:12]
), ncol = 4))
colnames(n1000_2stage) <- c("opt%", "sd", "EY(g^opt)", "sd")
rownames(n1000_2stage) <- rown
print("Table 1 - two stage study n = 1000")
print(n1000_2stage)


#########################################
#Plot Figure 1 two stage - the density of four cases

# Generate true rewards for the test set
kk <- sapply(1:1000, function(x)
  mean_var_calc2stage(x, N2 = 500))

##############
# preparing plot data
EY500 <- aa[c(7, 8, 10, 12), ]
EY1000 <- bb[c(7, 8, 10, 12), ]
EY500 <-
  data.frame(
    Scenario = rep(rep(
      c(
        "Case 1 with IE",
        "Case 2 without IE",
        "Case 2 with IE",
        "Case 1 without IE"
      ),
      each = 1000
    ), 2),
    n = rep(c("n = 500", "n = 1000"), each = 4000),
    EY = c(c(t(EY500)), c(t(EY1000)))
  )
opt_data <-
  data.frame(Scenario = "Optimal",
             n = rep(c("n = 500", "n = 1000"), each = 1000),
             EY = rep(kk[1, ], 4))
EY500 = rbind(EY500, opt_data)
EY500$Scenario <- factor(
  EY500$Scenario,
  levels = c(
    "Case 1 with IE",
    "Case 1 without IE",
    "Case 2 with IE",
    "Case 2 without IE",
    "Optimal"
  )
)
EY500$n <- factor(EY500$n, levels = c("n = 500", "n = 1000"))
EY500$case <-
  c(rep(rep(
    c("Case 1", "Case 2", "Case 2", "Case 1"), each = 1000
  ), 2), rep(c("Case 1", "Case 2"), each = 2000))

A = ggplot(data = EY500, aes(
  x = EY,
  linetype = Scenario,
  color = Scenario
)) +
  geom_line(stat = "density", lwd = 0.8) +
  facet_wrap( ~ n + case) + xlab(expression(paste(widehat(E), "{", Y, "*", 
                                                  (widehat(g)^{opt}), "}"))) +
  ylab("Density") + theme_bw(base_size = 12) + ggtitle("Two-stage Study")

#Plot Figure 1 - two stage
print(A)

########################################################
# simulation 3 Stages 2-3 Treatments per stage
########################################################
rown3 = c(
  "Case 1 with text",
  "Case 1 without text",
  "Case 2 with text imputed",
  "Case 2 without text imputed"
)

#Load intermediate result cc,dd
load("../data/data_simulation/three_stage_sim.Rdata")

## Avoid running the following code after loading intermediate result
## Can reduce 1:1000 to 1:10 for faster run time
## n = 500

# ptm <- proc.time()
# cc <- sapply(1:1000, function(x)
#   sim_text_3stage(x, 500, text_data))
# proc.time() - ptm #3.09 hours
#

## n = 1000

# ptm <- proc.time()
# dd <- sapply(1:1000, function(x)
#   sim_text_3stage(x, 1000, text_data))
# proc.time() - ptm #5.58 hours

## Avoid running the code above after loading intermediate result

# Print tables
n500_3stage <- data.frame(matrix(c(
  rowMeans(cc)[1:4],
  apply(cc, 1, sd)[1:4],
  rowMeans(cc)[5:8],
  apply(cc, 1, sd)[5:8]
),
ncol = 4))
colnames(n500_3stage) <- c("opt%", "sd", "EY(g^opt)", "sd")
rownames(n500_3stage) <- rown3
print("Table 1 - three stage study n = 500")
print(n500_3stage)

n1000_3stage <- data.frame(matrix(c(
  rowMeans(dd)[1:4],
  apply(dd, 1, sd)[1:4],
  rowMeans(dd)[5:8],
  apply(dd, 1, sd)[5:8]
),ncol = 4))
colnames(n1000_3stage) <- c("opt%", "sd", "EY(g^opt)", "sd")
rownames(n1000_3stage) <- rown3
print("Table 1 - three stage study n = 1000")
print(n1000_3stage)

#########################################
#Plot Figure 1 - three stage the density of four cases

# Generate true rewards for the test set
opt <- sapply(1:1000, function(x)
  mean_var_calc3stage(x, N2 = 1000))


##############
# preparing plot data
EY500 <- cc[5:8, ]
EY1000 <- dd[5:8, ]
EY500 <-
  data.frame(
    Scenario = rep(rep(
      c(
        "Case 1 with IE",
        "Case 2 without IE",
        "Case 2 with IE",
        "Case 1 without IE"
      ),
      each = 1000
    ), 2),
    n = rep(c("n = 500", "n = 1000"), each = 4000),
    EY = c(c(t(EY500)), c(t(EY1000)))
  )
opt_data <-
  data.frame(Scenario = "Optimal",
             n = rep(c("n = 500", "n = 1000"), each = 1000),
             EY = rep(opt[1, ], 4))
EY500 = rbind(EY500, opt_data)
EY500$Scenario <- factor(
  EY500$Scenario,
  levels = c(
    "Case 1 with IE",
    "Case 1 without IE",
    "Case 2 with IE",
    "Case 2 without IE",
    "Optimal"
  )
)
EY500$n <- factor(EY500$n, levels = c("n = 500", "n = 1000"))
EY500$case <-
  c(rep(rep(
    c("Case 1", "Case 2", "Case 2", "Case 1"), each = 1000
  ), 2), rep(c("Case 1", "Case 2"), each = 2000))

B = ggplot(data = EY500, aes(
  x = EY,
  linetype = Scenario,
  color = Scenario
)) + geom_line(stat = "density", lwd = 0.8) +
  facet_wrap( ~ n + case) + xlab(expression(paste(widehat(E), "{", Y, "*", 
                                                  (widehat(g) ^ {opt}), "}"))) +
  ylab("Density") + theme_bw(base_size = 12) + ggtitle("Three-stage Study")

# Plot Figure 1 - three stage
print(B)

#####Application Section Code#####

###########Fit the DTR ################
library(webshot)
library(htmlwidgets)
library(data.tree)
library(DiagrammeR)
webshot::install_phantomjs() # install PhatomJS for the first time running the program
source("../function/Functions.R")
source("../function/plot_tree.R")

# load reprocessed application data.
vital_demo_drug_2 <- read.csv("../data/data_hyper_emergency/application_data_day2.csv")
vital_demo_drug_1 <- read.csv("../data/data_hyper_emergency/application_data_day1.csv")

vital_demo_drug_2$HR[which(is.na(vital_demo_drug_2$HR))] = mean(vital_demo_drug_2$HR, na.rm = TRUE)
vital_demo_drug_2$Temp[which(is.na(vital_demo_drug_2$Temp))] = mean(vital_demo_drug_2$Temp, na.rm = TRUE)
vital_demo_drug_2$SpO2[which(is.na(vital_demo_drug_2$SpO2))] = mean(vital_demo_drug_2$SpO2, na.rm = TRUE)
vital_demo_drug_2$trt <- ifelse(vital_demo_drug_2$cal == 1,
                                3,
                                ifelse(
                                  vital_demo_drug_2$beta == 1,
                                  2,
                                  ifelse(vital_demo_drug_2$ace == 1, 1, 4)
                                ))
vital_demo_drug_1$trt <- ifelse(vital_demo_drug_1$cal == 1,
                                3,
                                ifelse(
                                  vital_demo_drug_1$beta == 1,
                                  2,
                                  ifelse(vital_demo_drug_1$ace == 1, 1, 4)
                                ))

vital_demo_drug_2 <-
  vital_demo_drug_2 %>% left_join(vital_demo_drug_1 %>% select(icustay_id, trt), by =
                                    "icustay_id")


Y <- vital_demo_drug_2$sbp -  vital_demo_drug_2$final_sbp
#Characteristics
# HR,Resp,ADMISSION_AGE,wt_avg,GENDER,smoke, kidney, diabetes, copd,neuro
H <-
  cbind(vital_demo_drug_2[, c(
    "dbp",
    "sbp",
    "HR",
    "SpO2",
    "age",
    "wt_avg",
    "GENDER",
    "smoke",
    "kidney",
    "diabetes",
    "copd",
    "hyper"
  )], A1 = vital_demo_drug_2$trt.y)
#1=A,2=B,3=C,4=D
A = vital_demo_drug_2$trt.x

set.seed(58639)
R2surf <-
  randomForest(
    x = cbind(A, H),
    y = Y,
    mtry = 4,
    maxnodes = 20,
    importance = TRUE
  )

#############Continue DTR estimation#############
mu.hat <- cbind(
  predict(R2surf, newdata = cbind(A = rep(1, nrow(
    H
  )), H)),
  predict(R2surf, newdata = cbind(A = rep(2, nrow(
    H
  )), H)),
  predict(R2surf, newdata = cbind(A = rep(3, nrow(
    H
  )), H)),
  predict(R2surf, newdata = cbind(A = rep(4, nrow(
    H
  )), H))
)


pi.hat <- M.propen(A = A, Xs = H)

hyper_tree2 <-
  DTRtree(
    Y = Y,
    A = A,
    H = H,
    pis.hat = pi.hat,
    m.method = "AIPW",
    mus.reg = mu.hat,
    minsplit = 20,
    lambda.pct = 0.3,
    depth = 1
  )

# Plot Figure 2 stage 2
stage2 <- plot(plot.DTR.tree(round(hyper_tree2, 0), H = H))
saveWidget(stage2,"temp.html")
webshot("temp.html","stage2_dtr.png", vwidth = 500, vheight = 500)

rec_trt <- predict.DTR(hyper_tree2, newdata = H)

R2.pred <- predict(R2surf, newdata = cbind(A = rec_trt, H))
R2.true <- predict(R2surf, newdata = cbind(A = A, H))

#Stage 1
#Calculate pseudo outcomes.
Y1 <- vital_demo_drug_2$final_sbp + R2.pred - R2.true

Y1_pseudo = c()
for (i in vital_demo_drug_1$icustay_id) {
  if (i %in% vital_demo_drug_2$icustay_id) {
    Y1_pseudo = c(Y1_pseudo, Y1[which(vital_demo_drug_2$icustay_id == i)])
  } else{
    Y1_pseudo = c(Y1_pseudo, NA)
  }
}

vital_demo_drug_1$Y1_pseudo = Y1_pseudo
vital_demo_drug_1$final_pseudo <-
  ifelse(
    is.na(vital_demo_drug_1$Y1_pseudo),
    vital_demo_drug_1$final_sbp,
    vital_demo_drug_1$Y1_pseudo
  )

vital_demo_drug_1$TempC_Mean[which(is.na(vital_demo_drug_1$TempC_Mean))] = mean(vital_demo_drug_1$TempC_Mean, na.rm = TRUE)
Y <- (vital_demo_drug_1$SysBP_Max - vital_demo_drug_1$final_pseudo)
#Characteristics
# HR,Resp,ADMISSION_AGE,wt_avg,GENDER,smoke, kidney, diabetes, copd, neuro
H <-
  cbind(vital_demo_drug_1[, c(
    "SysBP_Max",
    "HeartRate_Mean",
    "TempC_Mean",
    "SpO2_Mean",
    "age",
    "creatinine_min",
    "wt_avg",
    "GENDER",
    "smoke",
    "kidney",
    "diabetes",
    "copd",
    "hyper"
  )])
#1=A,2=B,3=C,4=D
A = vital_demo_drug_1$trt

set.seed(48934)
R1surf1 <-
  randomForest(
    x = cbind(A, H),
    y = Y,
    ntree = 1000,
    importance = TRUE
  )
R1surf <-
  randomForest(
    x = cbind(A, H),
    y = Y,
    ntree = 500,
    mtry = 4,
    nodesize = 20,
    importance = TRUE
  )
mu.hat <- cbind(
  predict(R1surf, newdata = cbind(A = rep(1, nrow(
    H
  )), H)),
  predict(R1surf, newdata = cbind(A = rep(2, nrow(
    H
  )), H)),
  predict(R1surf, newdata = cbind(A = rep(3, nrow(
    H
  )), H)),
  predict(R1surf, newdata = cbind(A = rep(4, nrow(
    H
  )), H))
)

pi.hat <- M.propen(A = A, Xs = H)

hyper_tree <-
  DTRtree(
    Y = Y,
    A = A,
    H = H,
    pis.hat = pi.hat,
    m.method = "AIPW",
    mus.reg = mu.hat,
    minsplit = 15,
    lambda.pct = 0.1,
    depth = 2
  )

# Plot Figure 2 Stage 1
stage1 <- plot(plot.DTR.tree(round(hyper_tree, 0), H = H))
saveWidget(stage1,"temp.html")
webshot("temp.html","stage1_dtr.png", vwidth = 500, vheight = 500)


rec_trt <- predict.DTR(hyper_tree, newdata = H)

R1.pred <- predict(R1surf, newdata = cbind(A = rec_trt, H))
R1.true <- predict(R1surf, newdata = cbind(A = A, H))
