
library("survival")
library("plyr")
library("survminer")


#####################################Cox####################################

survData <- read.delim("breast.txt" , sep ="\t" , stringsAsFactors=FALSE)

covariates <- c("age", "ethnicity",  "stage", "BRCA1.FPKM", "FPKM.BRCA2.FPKM", "ATM.FPKM" )
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
                        
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = survData)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                         })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
write.table (as.data.frame(res), file="score.txt" , sep="\t", quote=FALSE )



res.cox <- coxph(Surv(time, status) ~ age + stage + ethnicity+ BRCA1.FPKM + BRCA2.FPKM + ATM.FPKM, data =  survData)
summary(res.cox)

################cox for stage##########################


###stage 1 cancer ################
stage1 <- with(survData,
               data.frame(stage = c(0,0),
			              age = rep(mean(age, na.rm = TRUE), 2),
                          BRCA2.FPKM = c(1,2),
						  BRCA1.FPKM = c(1,2),
						  ATM.FPKM = c(1, 2)
                          )
               )
stage1

fit1 <- survfit(res.cox, newdata = stage1)
jpeg(file="stage1_breast_Cox_Plot.jpeg", 500, 500, pointsize=15)
plot(fit1, main="Stage I Breast Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()

##############stage 2 cancer############

stage2 <- with(survData,
               data.frame(stage = c(1,1),
			              age = rep(mean(age, na.rm = TRUE), 2),
                          BRCA2.FPKM = c(1,2),
						  BRCA1.FPKM = c(1,2),
						  ATM.FPKM = c(1, 2)
                          )
               )
stage2


fit2 <- survfit(res.cox, newdata = stage2)

jpeg(file="stage2_breast_Cox_Plot.jpeg.jpeg", 500, 500, pointsize=15)
plot(fit2, main="Stage II Breast Cancer Cox Plott", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()
###############stage 3 cancer################

stage3 <- with(survData,
               data.frame(stage = c(2,2),
			              age = rep(mean(age, na.rm = TRUE), 2),
                          BRCA2.FPKM = c(1,2),
						  BRCA1.FPKM = c(1,2),
						  ATM.FPKM = c(1, 2)
                          )
               )
stage3


fit3 <- survfit(res.cox, newdata = stage3)

jpeg(file="Stage3_breast_Cox_Plot.jpeg.jpeg", 500, 500, pointsize=15)
plot(fit3, main="Stage III Breast Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()

#############stage 4 cancer #################

stage4 <- with(survData,
               data.frame(stage = c(3,3),
			              age = rep(mean(age, na.rm = TRUE), 2),
                          BRCA2.FPKM = c(1,2),
						  BRCA1.FPKM = c(1,2),
						  ATM.FPKM = c(1, 2)
                          )
               )
stage4

fit4 <- survfit(res.cox, newdata = stage4)
jpeg(file="Stage4_breast_Cox_Plot.jpeg.jpeg", 500, 500, pointsize=15)
plot(fit4, main="Stage IV Breast Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()


#######################cox for ethnicity#######################



###white ################
white <- with(survData,
               data.frame(ethnicity = c(0,0),
			              age = rep(mean(age, na.rm = TRUE), 2),
                          BRCA2.FPKM = c(1,2),
						  BRCA1.FPKM = c(1,2),
						  ATM.FPKM = c(1, 2)
                          )
               )
white

fit5 <- survfit(res.cox, newdata = white)
jpeg(file="white_breast_Cox_Plot.jpeg", 500, 500, pointsize=15)
plot(fit5, main="Ethnicity White Breast Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(2:3)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(2:3))
dev.off()

##############black############

black <- with(survData,
               data.frame(ethnicity = c(1,1),
			              age = rep(mean(age, na.rm = TRUE), 2),
                          BRCA2.FPKM = c(1,2),
						  BRCA1.FPKM = c(1,2),
						  ATM.FPKM = c(1, 2)
                          )
               )
black


fit6 <- survfit(res.cox, newdata = black)

jpeg(file="black_breast_Cox_Plot.jpeg.jpeg", 500, 500, pointsize=15)
plot(fit6, main="Ethnicity Black or African American Breast Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(2:3)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(2:3))
dev.off()


################ggsurvplot#########################
######White###########

jpeg(file="ggsurv_White_Breast_Cox_Plot.jpeg", 400, 400, pointsize=25)
ggsurvplot(fit5, data = white, # survfit object with calculated statistics.
  conf.int = FALSE,         # show confidence intervals for 
  ggtheme = theme_test(), # customize plot and risk table with a theme.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "Ethnicity White",
  legend.labs = 
  c("Low FPKM", "High FPKM"),    # change legend labels.
  palette = 
  c( "#2E9FDF", "#E74C3C") # custom color palettes.
)
dev.off()


##########black


jpeg(file="ggsurv_black_Breast_Cox_Plot.jpeg", 400, 400, pointsize=15)
ggsurvplot(fit6, data = black, # survfit object with calculated statistics.
  conf.int = FALSE,         # show confidence intervals for 
  ggtheme = theme_test(), # customize plot and risk table with a theme.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "Ethnicity Black",
  legend.labs = 
  c("Low FPKM", "High FPKM"),    # change legend labels.
  palette = 
  c( "#2E9FDF", "#E74C3C") # custom color palettes.
)
dev.off()

############################stage##########################

######StageI###########

jpeg(file="ggsurv_StageI_breast_Cox_Plot.jpeg")
ggsurvplot(fit1, data = stage1, # survfit object with calculated statistics.
  conf.int = FALSE,         # show confidence intervals for 
  ggtheme = theme_test(), # customize plot and risk table with a theme.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "Stage I Breast Cancer",
  legend.labs = 
  c("Low FPKM", "High FPKM"),    # change legend labels.
  palette = 
  c( "#3E2723", "#4CAF50") # custom color palettes.
)
dev.off()


##########stageII


jpeg(file="ggsurv_stageII_breast_Cox_Plot.jpeg", 500, 500, pointsize=25)
ggsurvplot(fit2, data = stage2, # survfit object with calculated statistics.
  conf.int = FALSE,         # show confidence intervals for 
  ggtheme = theme_test(), # customize plot and risk table with a theme.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "Stage II Breast Cancer",
  legend.labs = 
  c("Low FPKM", "High FPKM"),    # change legend labels.
  palette = 
  c( "#3E2723", "#4CAF50") # custom color palettes.
)
dev.off()

######StageIII###########

jpeg(file="ggsurv_StageIII_breast_Cox_Plot.jpeg", 500, 500, pointsize=25)
ggsurvplot(fit3, data = stage3, # survfit object with calculated statistics.
  conf.int = FALSE,         # show confidence intervals for 
  ggtheme = theme_test(), # customize plot and risk table with a theme.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "Stage IIII Breast Cancer",
  legend.labs = 
  c("Low FPKM", "High FPKM"),    # change legend labels.
  palette = 
  c( "#3E2723", "#4CAF50") # custom color palettes.
)
dev.off()


##########stageIV


jpeg(file="ggsurv_stageIV_breast_Cox_Plot.jpeg", 500, 500, pointsize=25)
ggsurvplot(fit4, data = stage4, # survfit object with calculated statistics.
  conf.int = FALSE,         # show confidence intervals for 
  ggtheme = theme_test(), # customize plot and risk table with a theme.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "Stage IV Breast Cancer",
  legend.labs = 
  c("Low FPKM", "High FPKM"),    # change legend labels.
  palette = 
  c( "#3E2723", "#4CAF50") # custom color palettes.
)
dev.off()






