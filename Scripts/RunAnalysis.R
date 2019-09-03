##Compute the all-cause survival at 5 years

#DLBCL
sfit <- survfit(Surv(FU_years, status) ~ 1, data = DLBCL)
summary(sfit, 5)

#Colon cancer
sfit <- survfit(Surv(FU_years, status) ~ 1, data = Colon)
summary(sfit, 5)

#AML
sfit <- survfit(Surv(FU_years, status) ~ 1, data = AML)
summary(sfit, 5)


##Fit parametric relative survival models
#DLBCL
#fit_DLBCL <- stpm2(Surv(FU_years, status) ~ 1, data = DLBCL, bhazard = DLBCL$exp_haz, df = 5)
fit_DLBCL <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = DLBCL, bhazard = "exp_haz", df = 3)

#Colon cancer
#fit_Colon <- stpm2(Surv(FU_years, status) ~ 1, data = Colon, bhazard = Colon$exp_haz, df = 5)
fit_Colon <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = Colon, bhazard = "exp_haz", df = 3)

#AML
#fit_AML <- stpm2(Surv(FU_years, status) ~ 1, data = AML, bhazard = AML$exp_haz, df = 5)
fit_AML <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = AML, bhazard = "exp_haz", df = 3)


#Gather parametric models
fits <- list(Colon = list(fit_Colon, data = Colon), 
             AML = list(fit_AML, data = AML), 
             DLBCL = list(fit_DLBCL,data = DLBCL))

#Compute relative survival probs for parametric and non-parametric model
plot_data <- lapply(fits, function(fit){
  rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                   data = fit$data, ratetable = survexp.dk, method = "ederer1")
  #Adjust time to years
  rsfit$time <- rsfit$time / ayear
  
  #Create data.frame with predictions
  D <- data.frame(RS = rsfit$surv, time = rsfit$time, ci.lower = rsfit$lower, ci.upper = rsfit$upper) 
  
  D_para <- predict(fit[[1]], time = D$time, var.type = "n")[[1]]
  names(D_para) <- "Est"
  D_para$time <- D$time
  D_para$model <- "FMC model"
  list(D = D, D_para = D_para)
})


#Create disease variable
diseases <- c("Colon cancer", "Acute myeloid leukemia", "Diffuse large B-cell lymphoma")

#Arrange plotting data and plot relative survival
para_plot_data <- do.call(rbind, lapply(plot_data, function(x) x$D_para))
para_plot_data$disease <- factor(rep(diseases, sapply(plot_data, function(x) nrow(x$D_para))), 
                                 levels = diseases)

npara_plot_data <- do.call(rbind, lapply(plot_data, function(x) x$D))
npara_plot_data$disease <- factor(rep(diseases, sapply(plot_data, function(x) nrow(x$D))), 
                                  levels = diseases)
colnames(npara_plot_data)[1] <- "Est"
npara_plot_data$model <- "Ederer I"

if(pdf){
  pdf(file.path(fig.out, "RSCombined.pdf"), width = 10.8, height = 6)
} else {
  tiff(file.path(fig.out, "RSCombined.tiff"), width = 11 * 300, height = 5 * 300, res = 300)
}

ggplot(data = npara_plot_data, aes(x = time, y = Est, group = model, colour = model)) + geom_step() + 
  facet_grid(.~disease) + geom_step(data = npara_plot_data, aes(x = time, y = ci.lower), linetype = "dashed") + 
  geom_step(data = npara_plot_data, aes(x = time, y = ci.upper), linetype = "dashed") + 
  geom_line(data = para_plot_data, aes(x = time, y = Est), size = 1) + 
  scale_colour_manual(values = c("Ederer I" = "black", 
                                 "FMC model" = "darkgray"), 
                      breaks = c("Ederer I", "FMC model")) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13),
        legend.key.size = unit(2,"line")) + 
  xlab("Time since diagnosis (years)") + 
  ylab("Relative survival") + coord_cartesian(ylim=c(0, 1))
dev.off()


##Create summary table for each disease
data_list <- list(Colon, AML, DLBCL)
#Mean age
mean_age <- sapply(data_list, function(data){
  age_mean <- sprintf("%.0f", mean(data$age_years))
  age_range <- sprintf("%.0f", range(data$age_years))
  paste0(age_mean, "(", age_range[1], "-", age_range[2], ")")
  
})

#Non-parametric relative survival at 5 years
probs_np <- sapply(data_list, function(data){
  rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                   data = data, ratetable = survexp.dk, method = "ederer2")
  s <- summary(rsfit, 5 * ayear)
  pred <- sprintf("%.2f", c(s$surv, s$lower, s$upper))
  paste0(pred[1], "(", pred[2], "-", pred[3], ")")
})

fit_list <- list(fit_Colon, fit_AML, fit_DLBCL)
#Parametric relative survival at 5 years
probs_p <- sapply(fit_list, function(fit){
  #pred <- sprintf("%.2f", predict(fit, newdata = data.frame(FU_years = 5), se.fit = T))
  pred <- sprintf("%.2f", predict(fit, time = 5, var.type = "ci")[[1]])
  paste0(pred[1], "(", pred[2], "-", pred[3], ")")
})

#Parametric cure rate
cure_rate <- sapply(fit_list, function(fit){
  res <- predict(fit, type = "curerate")[[1]]
  pred <- sprintf("%.2f", res)
  paste0(pred[1], "(", pred[2], "-", pred[3], ")")
})

#Crude cure rate
crude_cure_rates <- sapply(fit_list, function(fit){
  res <- calc.Crude(fit, time = 100, type = "disease", rmap = list(year = diag_date))[[1]]
  pred <- sprintf("%.2f", res)
  paste0(pred[1], "(", pred[2], "-", pred[3], ")")
})

#Loss of lifetime at time zero
LL_base <- sapply(fit_list, function(fit){
  LL <- calc.LL(fit, time = 0, rmap = list(year = diag_date))[[1]]
  pred <- sprintf("%.2f", LL)
  paste0(pred[1], "(", pred[2], "-", pred[3], ")")
})

#Create matrix for results
n <- prettyNum(c(nrow(Colon), nrow(AML), nrow(DLBCL)), big.mark=",")
M <- matrix(ncol = 3, nrow = 6)
diseases.short <- c("CC", "AML", "DLBCL")
colnames(M) <- paste0(diseases.short, " (n = ", n, ")")
rownames(M) <- c("Mean age (range)", "5-year RS (Ederer II)", "5-year RS (parametric)", 
                 "Cure proportion", "Probability of dying due to cancer", 
                 "Baseline loss of lifetime (years)")

#Input each of the above results
M[1,] <- mean_age
M[2,] <- probs_np
M[3,] <- probs_p
M[4,] <- cure_rate
M[5,] <- crude_cure_rates
M[6,] <- LL_base

#Print the table to file including legend and caption
print(xtable(M, 
             caption = "Relative survival estimates, cure proportion, probability of dying due to cancer, and baseline loss of 
       lifetime estimates. CC: colon cancer, AML: acute myeloid leukemia, 
             DLBCL: diffuse large B-cell lymphoma, RS: relative survival.", 
             label = "tab:sum", align = "lccc"), file = paste0(tab.out, "Table_sum.tex"), 
      table.placement = "!ht")


##Set time points for evaluating all future comparison measures
time <- seq(1e-5, 15, length.out = 100)


##Calculate cure point for specific epsilon values
filename <- "GeneratedData/cure_points_ci.RData"
if(file.exists(filename)){
  load(filename)
}else{
  LOL_eps <- c(1,2,3) #years
  sc_allLOL <- lapply(LOL_eps, function(q){
    cat(q, "\n")
    lapply(list(fit_Colon, fit_AML, fit_DLBCL), function(fit){
      cuRe::calc.LL.quantile(fit, q = q, max.time = 50, rmap = list(year = diag_date))
    })
  })
  
  crude_eps <- c(0.05, 0.1, 0.15) #probability
  sc_allCrude <- lapply(crude_eps, function(q){
    cat(q, "\n")
    lapply(list(fit_Colon, fit_AML, fit_DLBCL), function(fit){
      cuRe::calc.Crude.quantile(fit, q = q, max.time = 100, rmap = list(year = diag_date), reverse = TRUE)
    })
  })
  #Save to file to avoid excessive running time
  save(LOL_eps, sc_allLOL, crude_eps, sc_allCrude, file = filename)
}


##Calculate loss of lifetime
#Colon cancer
LL_Colon <- calc.LL(fit_Colon, time = time, rmap = list(year = diag_date))

#AML
LL_AML <- calc.LL(fit_AML, time = time, rmap = list(year = diag_date))

#DLBCL
LL_DLBCL <- calc.LL(fit_DLBCL, time = time, rmap = list(year = diag_date))

#Combine results into a single data frame
D <- rbind(LL_Colon[[1]], LL_AML[[1]], LL_DLBCL[[1]])
D$time <- rep(time, 3)
D$disease <- factor(rep(diseases, each = length(time)), levels = diseases)

#Create matrix with cure point information
cpLOL <- sapply(sc_allLOL, function(x) 
  sapply(x, function(y) y$Est))

#Create lines to add to below plot
df_lines1 <- data.frame(x1 = 0, x2 = c(cpLOL), y1 = rep(LOL_eps, each = 3), 
                        y2 = rep(LOL_eps, each = 3), disease = rep(diseases, 3))

df_lines2 <- data.frame(x1 = c(cpLOL), x2 = c(cpLOL), y1 = -2, y2 = rep(LOL_eps, each = 3), 
                        disease = rep(diseases, 3))

lims <- c(min(D$lower.ci), max(D$upper.ci))

#Create plot for the loss of lifetime
p1 <- ggplot(D, aes(x = time, y = Estimate)) + geom_line() + facet_grid(.~disease) + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.3) + 
  theme_bw() + xlab("Years since diagnosis") + ylab("Loss of lifetime (years)") + 
  coord_cartesian(ylim = lims, xlim = c(0, 15)) +  
  geom_segment(data = df_lines1, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "dashed") + 
  geom_segment(data = df_lines2, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "dashed") + 
  theme(axis.title=element_text(size=16),
        strip.text = element_text(size=14), 
        axis.text = element_text(size = 14))


##Probability of cancer related death
#Colon cancer
res_Colon <- calc.Crude(fit_Colon, time = time, type = "condother", rmap = list(year = diag_date), reverse = T)

#AML
res_AML <- calc.Crude(fit_AML, time = time, type = "condother", rmap = list(year = diag_date), reverse = T)

#DLBCL
res_DLBCL <- calc.Crude(fit_DLBCL, time = time, type = "condother", rmap = list(year = diag_date), reverse = T)

#Combine results into a single data frame
D <- rbind(res_Colon[[1]], res_AML[[1]], res_DLBCL[[1]])
D$time <- rep(time, 3)
D$disease <- factor(rep(diseases, each = length(time)), levels = diseases)

#Create matrix with cure point information
cpcrude <- sapply(sc_allCrude, function(x) 
  sapply(x, function(y) y$Est))

#Create lines to add to below plot
df_lines1 <- data.frame(x1 = 0, x2 = c(cpcrude), y1 = rep(crude_eps, each = 3), 
                        y2 = rep(crude_eps, each = 3), disease = rep(diseases, 3))

df_lines2 <- data.frame(x1 = c(cpcrude), x2 = c(cpcrude), y1 = -2, y2 = rep(crude_eps, each = 3), 
                        disease = rep(diseases, 3))

#Create plot for the probability of cancer related death
p2 <- ggplot(D, aes(x = time, y = Estimate)) + geom_line() + facet_grid(.~disease) + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.3) + 
  theme_bw() + xlab("Years since diagnosis") + 
  ylab("Conditional probability of cancer-related death (%)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +  
  geom_segment(data = df_lines1, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "dashed") + 
  geom_segment(data = df_lines2, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "dashed") + 
  theme(axis.title=element_text(size=16),
        strip.text = element_text(size=14), 
        axis.text = element_text(size = 14))


#Plot conditional probability of cancer related death with cure point markers
if(pdf){
  pdf(file.path(fig.out, "ProbDeathLOL.pdf"), width = 11, height = 11) 
} else {
  tiff(file.path(fig.out, "ProbDeathLOL.tiff"), res = 300, width = 11 * 300, height = 10 * 300)
}

grid.arrange(p1, p2, ncol = 1)
dev.off()



##Format the cure points estimated earlier (including confidence intervals)
#Loss of lifetime
res_LOL <- sapply(sc_allLOL, function(x){
  pred <- sprintf("%.2f", x[[1]])
  sc_colon <- paste0(pred[1], "(", pred[3], "-", pred[4], ")")
  pred <- sprintf("%.2f", x[[2]])
  sc_AML <- paste0(pred[1], "(", pred[3], "-", pred[4], ")")
  pred <- sprintf("%.2f", x[[3]])
  sc_DLBCL <- paste0(pred[1], "(", pred[3], "-", pred[4], ")")
  c(sc_colon, sc_AML, sc_DLBCL)
})

#Conditional probability of cancer related death
res_crude <- sapply(sc_allCrude, function(x){
  pred <- sprintf("%.2f", x[[1]])
  sc_colon <- paste0(pred[1], "(", pred[3], "-", pred[4], ")")
  pred <- sprintf("%.2f", x[[2]])
  sc_AML <- paste0(pred[1], "(", pred[3], "-", pred[4], ")")
  pred <- sprintf("%.2f", x[[3]])
  sc_DLBCL <- paste0(pred[1], "(", pred[3], "-", pred[4], ")")
  c(sc_colon, sc_AML, sc_DLBCL)
})

#Combine results in a data frame
M <- rbind(NA, t(res_LOL), NA, t(res_crude))
colnames(M) <- diseases.short
rownames(M) <- c("Loss of lifetime", paste0("$\\epsilon$ = ", LOL_eps), "Probability",
                 paste0("$\\epsilon$ = ", crude_eps))

#Output results to file with caption and legend
print(xtable(M, 
             caption = "Cure point estimates based on the loss of lifetime function and 
             the conditional probability of cancer-related death. 
             The cure point is calculated for varying magins, $\\epsilon$. 
             CC: Colon cancer, AML: Acute myeloid leukemia, DLBCL: Diffuse large B-cell lymphoma.", 
             label = "tab:cure_points", align = "lccc"), file = paste0(tab.out, "Table_curepoints.tex"), 
      sanitize.text.function = identity)



##Calculate cure points for more clinical relevant threshold
filename <- "GeneratedData/cure_points.RData"
if(file.exists(filename)){
  load(filename)
}else{
  LOL_eps2 <- seq(0.5, 4, length.out = 15) #years
  sc_allLOL2 <- lapply(LOL_eps2, function(q){
    cat(q, "\n")
    lapply(list(fit_Colon, fit_AML, fit_DLBCL), function(fit){
      cuRe::calc.LL.quantile(fit, q = q, max.time = 50, rmap = list(year = diag_date))
    })
  })
  
  crude_eps2 <- seq(0.02, 0.16, by = 0.01) #probability
  sc_allCrude2 <- lapply(crude_eps2, function(q){
    cat(q, "\n")
    lapply(list(fit_Colon, fit_AML, fit_DLBCL), function(fit){
      cuRe::calc.Crude.quantile(fit, q = q, max.time = 50, rmap = list(year = diag_date), reverse = TRUE)
    })
  })
  save(LOL_eps2, sc_allLOL2, crude_eps2, sc_allCrude2, file = filename)
}


#Combine loss of lifetime results into a single data frame
plot_data <- lapply(sc_allLOL2, function(x){
  D <- do.call(rbind, x)
  D$disease <- c("CC", "AML", "DLBCL")
  D
})

#Add information to the data frame about threshold
plot_data <- do.call(rbind, plot_data)
plot_data$eps <- rep(LOL_eps2, each = 3)
plot_data$measure <- "Loss of lifetime"

#Comine probability results into a single data frame
plot_data2 <- lapply(sc_allCrude2, function(x){
  D <- do.call(rbind, x)
  D$disease <- c("CC", "AML", "DLBCL")
  D
})

#Add information to the data frame about threshold
plot_data2 <- do.call(rbind, plot_data2)
plot_data2$eps <- rep(crude_eps2, each = 3)
plot_data2$measure <- "Conditional probability of cancer-related death"

#Combine loss of lifetime and probability results and add disease and measure indicators
plot_data_all <- rbind(plot_data, plot_data2)
plot_data_all$disease <- factor(plot_data_all$disease, levels = c("CC", "AML", "DLBCL"))
plot_data_all$measure <- factor(plot_data_all$measure, 
                                levels = c("Loss of lifetime", "Conditional probability of cancer-related death"))

d1 <- max(crude_eps2) - min(crude_eps2)
d2 <- max(LOL_eps2) - min(LOL_eps2)
plot_data_all$wd <- rep(c(d2, d1) / 30, each = nrow(plot_data_all) / 2)
  
#Create the final plot
if(pdf){
  pdf(file.path(fig.out, "StatCureThresholds.pdf"), width = 10.5, height = 7.5) 
} else {
  tiff(file.path(fig.out, "StatCureThresholds.tiff"), width = 10 * 300, height = 7 * 300, res = 300)
}
ggplot(plot_data_all, aes(x = eps, y = Estimate)) + geom_point() + 
  facet_grid(disease ~ measure, scale = "free_x") + 
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, width = wd)) + 
  xlab("Clinical relevant magin") + ylab("Estimated cure point (years after diagnosis)") + 
  coord_cartesian(ylim=c(0, 20)) + 
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank(), 
                     axis.title=element_text(size=17),
                     strip.text = element_text(size=14), 
                     axis.text = element_text(size = 14))
dev.off()

#Display cure point estimates for DLBCL and loss of lifetime with varying MOCR
new.df <- plot_data_all[plot_data_all$measure == "Loss of lifetime",]
new.df <- new.df[new.df$disease == "DLBCL",]
new.df

