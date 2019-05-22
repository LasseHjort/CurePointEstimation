
#Subset colon cancer patients with available stage variable
Colon2 <- Colon2[Colon2$UICC %in% paste0("UICC stadium ", c("I", "II", "III", "IV")),]

#Create dichotomized stage variable
Colon2$stage <- as.numeric(Colon2$UICC %in% c("UICC stadium III", "UICC stadium IV"))


#Check out the relative survival
# rsfit <- rs.surv(Surv(FU, status)~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
#                  data = Colon2[Colon2$stage == 0,], ratetable = survexp.dk, method = "ederer2")
# rsfit$time <- rsfit$time / ayear
# plot(rsfit)


#Create age groups
Colon2$age_group <- cut(Colon2$age_years, breaks = c(0, 55, 65, 75, 120), 
                        labels = c("-60", "60-70", "70-80", "80-"))

#Create indicators for each strata
stage <- rep(0:1, each = nlevels(Colon2$age_group) * 2)
sex <- rep(rep(c("female", "male"), each = nlevels(Colon2$age_group)), 2)
age <- rep(levels(Colon2$age_group), 2 * 2)

df <- data.frame(stage = stage, sex = sex, age = age)

#Compute the non-parametric relative survival for each strata
L.nonp <- vector("list", length(stage))
for(i in 1:length(L.nonp)){
  cat(i, "\n")
  rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
                   data = Colon2[Colon2$sex == df$sex[i] & 
                                   Colon2$age_group == df$age[i] & 
                                   Colon2$stage == df$stage[i], ], 
                   ratetable = survexp.dk, method = "ederer2") 
  rsfit$time <- rsfit$time / ayear
  D <- data.frame(surv = rsfit$surv, time = rsfit$time, lower = rsfit$lower, upper = rsfit$upper)
  D$stage <- stage[i]
  D$age_group <- age[i]
  L.nonp[[i]] <- D
}

#Plot the relative survival for male patients
which.male <- which(sex == "male")
plot_data.male <- do.call(rbind, L.nonp[which.male])
plot_data.male$stage <- factor(plot_data.male$stage, levels = 0:1, labels = c("I-II", "III-IV"))


ggplot(plot_data.male, aes(x = time, y = surv, colour = stage)) + geom_step() + facet_wrap(~age_group, ncol = 2) + 
  xlab("Follow-up time (years)") + ylab("Relative survival") + theme_bw() + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_colour_discrete(name = "UICC stage") + coord_cartesian(xlim=c(0, 12.5), ylim = c(0, 1)) +   
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = stage, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + ggtitle("Male colon cancer patients")


#Plot the relative survival for female patients
which.female <- which(sex == "female")
plot_data.female <- do.call(rbind, L.nonp[which.female])
plot_data.female$stage <- factor(plot_data.female$stage, levels = 0:1, labels = c("I-II", "III-IV"))


ggplot(plot_data.female, aes(x = time, y = surv, colour = stage)) + geom_step() + facet_wrap(~age_group, ncol = 2) + 
  xlab("Follow-up time (years)") + ylab("Relative survival") + theme_bw() + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_colour_discrete(name = "UICC stage") + coord_cartesian(xlim=c(0, 12.5), ylim = c(0, 1)) +   
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = stage, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + ggtitle("Female colon cancer patients")


#Estimate the Nelson et al. model for each strata
L <- vector("list", length(stage))
for(i in 1:length(L)){
  cat(i, "\n")
  new.data <- Colon2[Colon2$sex == df$sex[i] & 
                       Colon2$age_group == df$age[i] & 
                       Colon2$stage == df$stage[i], ]
  # L[[i]] <- stpm2(Surv(FU_years, status) ~ 1, 
  #                            data = new.data, 
  #                            df = 5, bhazard = new.data$exp_haz) 
  L[[i]] <- GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                             data = new.data, 
                             df = 3, bhazard = "exp_haz", verbose = F) 
}

#Use the predict function to get relative survival probabilities for each strata
pred <- vector("list", length(stage))
time <- seq(0, 13, length.out = 100)
for(i in 1:length(pred)){
  res <- rep(NA, length(time))
  #res[time != 0] <- predict(L[[i]], newdata = data.frame(FU_years = time[time != 0]), keep.attributes = F)
  res[time != 0] <- predict(L[[i]], time = time[time != 0])[[1]]$Estimate
  res[time == 0] <- 1
  pred[[i]] <- data.frame(surv = res, time = time, stage = stage[i], age_group = age[i])
}

#Plot parametric relative surivival for male patients
which.male <- which(sex == "male")
plot_data2.male <- do.call(rbind, pred[which.male])
plot_data2.male$stage <- factor(plot_data2.male$stage, levels = 0:1, labels = c("I-II", "III-IV"))

p <- ggplot(plot_data.male, aes(x = time, y = surv, colour = stage)) + geom_step() + facet_wrap(~age_group, ncol = 2) + 
  xlab("Follow-up time (years)") + ylab("Relative survival") + theme_bw() + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 18), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13)) + 
  scale_colour_discrete(name = "UICC stage") + coord_cartesian(xlim=c(0, 12.5), ylim = c(0, 1)) +   
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = stage, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + ggtitle("Male colon cancer patients") + 
  geom_line(data = plot_data2.male, linetype = "dashed")

pdf(paste0(fig.out, "ColonStratRSMale.pdf"), width = 10, height = 8)
print(p)
dev.off()

#Plot parametric relative surivival for female patients
which.female <- which(sex == "female")
plot_data2.female <- do.call(rbind, pred[which.female])
plot_data2.female$stage <- factor(plot_data2.female$stage, levels = 0:1, labels = c("I-II", "III-IV"))

p <- ggplot(plot_data.female, aes(x = time, y = surv, colour = stage)) + geom_step() + facet_wrap(~age_group, ncol = 2) + 
  xlab("Follow-up time (years)") + ylab("Relative survival") + theme_bw() + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_colour_discrete(name = "UICC stage") + coord_cartesian(xlim=c(0, 12.5), ylim = c(0, 1)) +   
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = stage, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + ggtitle("Female colon cancer patients") + 
  geom_line(data = plot_data2.female, linetype = "dashed")

pdf(paste0(fig.out, "ColonStratRSFeMale.pdf"), width = 10, height = 8)
print(p)
dev.off()


#Calculate conditional probability of cancer-related death for each strata (used in article)
crudeprob <- vector("list", length(stage))
time <- seq(0.0001, 40, length.out = 50)
for(i in 1:length(crudeprob)){
  cat(i, "\n")
  if(df$stage[i] == 0 & (df$age[i] == "80-" | 
                         (df$age[i] == "-60" & df$sex[i] == "female") | 
                         (df$age[i] == "60-70" & df$sex[i] == "female"))){
    link <- "identity"
  } else{
    link <- "loglog"
  }
  crudeprob[[i]] <- calc.Crude(L[[i]], type = "condother", 
                               rmap = list(year = diag_date), 
                               reverse = TRUE, time = time, link = link)[[1]]
  crudeprob[[i]]$time <- time
  crudeprob[[i]]$sex <- sex[i]
  crudeprob[[i]]$stage <- stage[i]
  crudeprob[[i]]$age_group <- age[i]
}

#Create data for plotting
plot_data <- do.call(rbind, crudeprob)
plot_data$stage <- factor(plot_data$stage, levels = 0:1, labels = c("I-II", "III-IV"))
plot_data$sex <- factor(plot_data$sex)
plot_data$stagesex <- plot_data$stage:plot_data$sex

#Create plot for female patients
plot_data2 <- plot_data[plot_data$sex == "female",]
p <- ggplot(plot_data2, aes(x = time, y = Estimate, linetype = stage)) + 
  geom_line(size = 1) + facet_wrap(~ age_group, ncol = 2) + xlab("Years since diagnosis") + 
  ylab("Conditional probability of cancer-related death") + theme_bw() + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 18), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13),
        legend.key.size = unit(2,"line")) + 
  ggtitle("Female colon cancer patients") + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + 
  scale_linetype_discrete(name = "UICC stage")

#Output to file
if(pdf){
  pdf(paste0(fig.out, "ColonStratProbFemale.pdf"), width = 11, height = 9) 
} else {
  tiff(paste0(fig.out, "ColonStratProbFemale.tiff"), res = 300, width = 10 * 300, height = 8 * 300)
}
print(p)
dev.off()

#Create plot for male patients
plot_data2 <- plot_data[plot_data$sex == "male",]
p <- ggplot(plot_data2, aes(x = time, y = Estimate, linetype = stage)) + 
  geom_line(size = 1) + facet_wrap(~ age_group, ncol = 2) + xlab("Years since diagnosis") + 
  ylab("Conditional probability of cancer-related death") + theme_bw() + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 18), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13),
        legend.key.size = unit(2,"line")) + 
  ggtitle("Male colon cancer patients") + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + 
  scale_linetype_discrete(name = "UICC stage")

#Output to file
if(pdf){
  pdf(paste0(fig.out, "ColonStratProbMale.pdf"), width = 11, height = 9) 
} else {
  tiff(paste0(fig.out, "ColonStratProbMale.tiff"), res = 300, width = 10 * 300, height = 8 * 300)
}
print(p)
dev.off()


#Compute the cure point for q = 0.1 using the probability measure
cp1 <- vector("list", length(stage))
for(i in 1:length(stage)){
  cat(i, "\n")
  cp1[[i]] <- calc.Crude.quantile(L[[i]], q = 0.075, rmap = list(year = diag_date), reverse = TRUE)
}

#Compute the cure point for q = 0.05 using the probability measure
cp2 <- vector("list", length(stage))
for(i in 1:length(stage)){
  cat(i, "\n")
  cp2[[i]] <- calc.Crude.quantile(L[[i]], q = 0.05, max.time = 50, rmap = list(year = diag_date), reverse = TRUE)
}

#Compute the cure point for q = 0.025 using the probability measure
cp3 <- vector("list", length(stage))
for(i in 1:length(stage)){
  cat(i, "\n")
  cp3[[i]] <- calc.Crude.quantile(L[[i]], q = 0.025, max.time = 60, rmap = list(year = diag_date), reverse = TRUE)
}

#Create data.frame for plotting the cure points for all strata
thresholds <- c(0.075, 0.05, 0.025)
plot_data <- rbind(do.call(rbind, cp1), do.call(rbind, cp2), do.call(rbind, cp3))
plot_data$threshold <- factor(rep(thresholds, each = length(cp1)))

plot_data$stage <- rep(stage, length(thresholds))
plot_data$stage <- factor(rep(stage, length(thresholds)), levels = 0:1, labels = c("UICC stage I-II", "UICC stage III-IV"))
plot_data$age <- factor(rep(age, length(thresholds)))
plot_data$sex <- factor(rep(sex, length(thresholds)), levels = c("female", "male"), labels = c("Female", "Male"))

plot_data$y <- 4:1
plot_data$y[plot_data$threshold == thresholds[1]] <- plot_data$y[plot_data$threshold == thresholds[1]] + 0.1
plot_data$y[plot_data$threshold == thresholds[3]] <- plot_data$y[plot_data$threshold == thresholds[3]] - 0.1
plot_data$missing <- factor(is.na(plot_data$lower.ci), c(TRUE, FALSE), labels = c("Missing CI", "Available CI"))


#Create plot and output to file
p <- ggplot(plot_data, aes(x = Estimate, y = y, colour = threshold)) + 
  geom_point(size = 2) + facet_grid(stage ~ sex) + 
  ylab("") + xlab("Estimated cure point (years after diagnosis)") + 
  geom_segment(aes(x = lower.ci, xend = upper.ci, y = y, yend = y), linetype = "dashed") + 
  coord_cartesian(xlim=c(0, 30)) + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  scale_y_continuous(labels = c("1" = "80-", "2" = "70-80", "3" = "60-70", "4" = "-60")) + 
  scale_colour_grey() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                              legend.text=element_text(size=15), 
                              axis.title=element_text(size=17),
                              strip.text = element_text(size=15), 
                              axis.text = element_text(size = 14),
                              legend.key.size = unit(2.5,"line"))

if(pdf) {
  pdf(paste0(fig.out, "SubAnalyses.pdf"), width = 11, height = 8.2) 
} else{
  tiff(paste0(fig.out, "SubAnalyses.tiff"), res = 300, width = 9 * 300, height = 6 * 300)
}
print(p)
dev.off()

