
library("survival")
library("plyr")
library("survminer")

###############################################KM################################################
survData <- read.delim("ovary.txt" , sep ="\t" , stringsAsFactors=FALSE)


###BRCA1
fit1 <- survfit(Surv(time, status) ~ BRCA1.FPKM, data = survData)
print(fit1)
jpeg(file="Ovary_FPKM_Plot_BRCA1.jpeg", 500, 500, pointsize=15)
plot(fit1, main="K-M Plot for BRCA1", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()

##BRCA2
fit2 <- survfit(Surv(time, status) ~ BRCA2.FPKM, data = survData)
print(fit2)
jpeg(file="Ovary_FPKM_Plot_BRCA2.jpeg", 500, 500, pointsize=15)
plot(fit2, main="K-M Plot for BRCA2", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()

###ATM
fit3 <- survfit(Surv(time, status) ~ ATM.FPKM, data = survData)
print(fit3)
jpeg(file="Ovary_FPKM_Plot_ATM.jpeg", 500, 500, pointsize=15)
plot(fit3, main="K-M Plot for ATM", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2)) 
legend('topright', c("Low FPKM","High FPKM"), lty=1, col=c(1:2))
dev.off()


#####ethnicity
fit4 <- survfit(Surv(time, status) ~ ethnicity, data = survData)
print(fit4)
jpeg(file="Ovary_ethnicity_Plot.jpeg", 500, 500, pointsize=15)
plot(fit4, main="Ovary_K-M Plot for ethnicity", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:3))
legend('topright', c("white","black","asian"), lty=1, col=c(1:4))
dev.off()


###ggplot#######

###################BRCA1#######################

jpeg(file="BRCA1_Plot_ggplot_ovary.jpeg", 500, 500, pointsize=1)
ggsurvplot(
   fit1,                     # survfit object with calculated statistics.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimaes of survival curves.
   conf.int.style = "step",  # customize style of confidence intervals
   ggtheme = theme_light(), # customize plot and risk table with a theme.
   risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
                            # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "BRCA1 FPKM",
  legend.labs = 
    c("Low", "High"),    # change legend labels.
  palette = 
    c( "#2ECC71", "#E74C3C") # custom color palettes.
)
dev.off()

#################BRCA2################################

jpeg(file="BRCA2_Plot_ggplot_ovary.jpeg", 500, 500, pointsize=1)
ggsurvplot(
   fit2,                     # survfit object with calculated statistics.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimaes of survival curves.
   conf.int.style = "step",  # customize style of confidence intervals
   ggtheme = theme_light(), # customize plot and risk table with a theme.
   risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
                            # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "BRCA2 FPKM",
  legend.labs = 
    c("Low", "High"),    # change legend labels.
  palette = 
    c( "#2ECC71", "#E74C3C") # custom color palettes.
)
dev.off()

########################ATM##################

jpeg(file="ATM_Plot_ggplot_ovary.jpeg", 500, 500, pointsize=1)
ggsurvplot(
   fit3,                     # survfit object with calculated statistics.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimaes of survival curves.
   conf.int.style = "step",  # customize style of confidence intervals
   ggtheme = theme_light(), # customize plot and risk table with a theme.
   risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
                            # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "ATM FPKM",
  legend.labs = 
    c("Low", "High"),    # change legend labels.
  palette = 
    c( "#2ECC71", "#E74C3C") # custom color palettes.
)
dev.off()


###########ethnicity######
jpeg(file="ethnicity_Plot_ggplot_ovary.jpeg", 500, 500, pointsize=1)
ggsurvplot(
   fit4,                     # survfit object with calculated statistics.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimaes of survival curves.
   conf.int.style = "step",  # customize style of confidence intervals
   ggtheme = theme_light(), # customize plot and risk table with a theme.
   risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
                            # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
    legend.title = "Ethnicity",
  legend.labs = 
    c("White", "Black or African American", "Asian"),    # change legend labels.
  palette = 
    c( "#2E9FDF", "#2ECC71", "#E74C3C") # custom color palettes.
)
dev.off()

















