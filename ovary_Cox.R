
library("survival")
library("plyr")
library("survminer")


#####################################Cox####################################

survData <- read.delim("Ovary.txt" , sep ="\t" , stringsAsFactors=FALSE)
#res.cox <- coxph(Surv(time, status) ~ sex, data = survData)
#res.cox
#summary(res.cox)


covariates <- c("age", "ethnicity", "BRCA1.FPKM", "BRCA2.FPKM", "ATM.FPKM" )
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(time, status)~', x)))
                        
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

res.cox <- coxph(Surv(time, status) ~ age + ethnicity +  BRCA1.FPKM + BRCA2.FPKM + ATM.FPKM, data =  survData)
summary(res.cox)


###visualize Ethnicity##############


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

fit1 <- survfit(res.cox, newdata = white)
jpeg(file="white_ovary_Cox_Plot.jpeg", 500, 500, pointsize=15)
plot(fit1, main="Ethnicity white Ovarian Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(2:3)) 
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


fit2 <- survfit(res.cox, newdata = black)

jpeg(file="black_ovary_Cox_Plot.jpeg.jpeg", 600, 600, pointsize=15)
plot(fit2, main="Ethnicity Black or African American Ovarian Cancer Cox Plot", xlab="Time in Days", ylab="Proportion Surviving", col=c(2:3)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(2:3))
dev.off()



################ggsurvplot#########################
######White###########

jpeg(file="ggsurv_White_ovary_Cox_Plot.jpeg", 500, 500, pointsize=15)
ggsurvplot(fit1, data = white, # survfit object with calculated statistics.
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


jpeg(file="ggsurv_black_ovary_Cox_Plot.jpeg", 500, 500, pointsize=15)
ggsurvplot(fit2, data = black, # survfit object with calculated statistics.
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