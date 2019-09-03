#Do simulations

#Function to compute the expected survival from survexp object. 
exp_function <- function(t, expected){
  s <- summary(expected, t)
  names(s$surv) <- s$time
  a <- s$surv[as.character(t)]
  names(a) <- NULL
  a
}

#Global parameters for gaussian quadrature
gaussxw <- statmod::gauss.quad(100)

#Function for calculating the loss of lifetime based on the true relative survival model
calc.ll <- function(time, eps = 0, sim_data){
  tau <- 60
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_rel_t <- sim_data$rel_surv(time)
  eval_gen_t <- exp_function(time, sim_data$S_exp)
  eval <- rep(NA, length(time))
  for(i in 1:length(eval)){
    points <- scale[i] * gaussxw$nodes  + scale2[i]
    eval_gen <- exp_function(points, sim_data$S_exp)
    eval_rel <- sim_data$rel_surv(points)
    eval[i] <- sum(gaussxw$weights * (eval_gen  - eval_gen * eval_rel / eval_rel_t[i]))
  }
  scale * eval / eval_gen_t - eps
}

#Function for calculating the cure point based on the "true" loss of lifetime function
find_cure_point_LOL <- function(sim_data, eps){
  f <- function(time) calc.ll(time, eps = eps, sim_data = sim_data)
  uniroot(f, lower = 0, upper = 30)$root  
}

#Function for calculating the conditional probability of cancer-related death based 
#on the true relative survival model
calc.prob <- function(time, eps = 0, sim_data){
  tau <- 60
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_rel_t <- sim_data$rel_surv(time)
  eval_gen_t <- exp_function(time, sim_data$S_exp)
  eval <- rep(NA, length(time))
  for(i in 1:length(eval)){
    points <- scale[i] * gaussxw$nodes  + scale2[i]
    eval_gen <- exp_function(points, sim_data$S_exp)
    eval_rel <- sim_data$rel_surv(points)
    eval_haz <- sim_data$haz_fun(points)
    eval[i] <- sum(gaussxw$weights * (eval_gen * eval_rel * eval_haz))
  }
  scale * eval / (eval_gen_t * eval_rel_t) - eps
}

#Function for calculating the cure point based on the "true" conditional probability of cancer related death
find_cure_point_prob <- function(sim_data, eps){
  f <- function(time) calc.prob(time, eps = eps, sim_data = sim_data)
  uniroot(f, lower = 0, upper = 30)$root  
}

#Function for calculating the conditional probability of cure based on the true relative survival model
calc.prob.cure <- function(time, eps, sim_data){
  pi <- sim_data$pars[1]
  surv <- pi + (1 - pi) * exp(-sim_data$pars[3] * time ^ sim_data$pars[2])
  1 - pi / surv - eps
}

#Function for calculating the cure point based on the "true" conditional probability of cure
find_cure_point_prob.cure <- function(sim_data, eps){
  f <- function(time) calc.prob.cure(time, eps = eps, sim_data = sim_data)
  uniroot(f, lower = 0, upper = 30)$root
}


#Simple global function to find equidistant quantiles
getq <- function(k){
  seq(0, 1, length.out = k)
}

#Simulation options
n.sim <- 500
n.obs <- 2000
age <- 60

#Set the number of cores used for the simulations
n.cores <- 48

#Selected margins of clinical relevance
eps_LOL <- 1:3
eps_prob <- c(0.05, 0.1, 0.15)

#Specify simulations
cases_wei_obs <- lapply(cases_wei, function(x) c(x, age = age, n = n.obs))

#List models used in the simulations
models <- c("Nelson et al.", "Andersson et al.", "Jakobsen et al.", "Weibull mixture")


#Create figure with true comparison measures
sim_datas <- list(do.call(sim_surv, cases_wei_obs[[1]]), 
                  do.call(sim_surv, cases_wei_obs[[2]]),
                  do.call(sim_surv, cases_wei_obs[[3]]))
times <- seq(0, 15, length.out = 100)
measures <- c("Conditional probability of cure", 
              "Conditional probability cancer-related death", 
              "Loss of lifetime")
#Calculate true comparison measures for each scenario
generate_results <- lapply(sim_datas, function(sim_data){
  res1 <- calc.prob.cure(time = times, sim_data = sim_data, eps = 0)
  res2 <- calc.prob(time = times, sim_data = sim_data)
  res3 <- calc.ll(time = times, sim_data = sim_data)
  
  data.frame(res = c(res1, res2, res3), time = rep(times, 3), 
             type = rep(measures, each = length(times)))
})
#Combine the values and plot 
res <- do.call(rbind, generate_results)
res$scenario <- paste0("Scenario ", rep(1:3, each = length(times) * 3))

p <- ggplot(data = res, aes(y = res, x = time)) + geom_line() + 
  facet_wrap(scenario ~ type, scales = "free_y") + xlab("Time (years)") + ylab("Comparison measure") + 
  theme_bw()

#Output to file
if(pdf){
  pdf(file.path(fig.out, "TrueComparisonMeasures.pdf"), width = 10, height = 8) 
} else{
  tiff(file.path(fig.out, "TrueComparisonMeasures.tiff"), res = 300, width = 8 * 300, height = 5 * 300)
}
print(p)
dev.off()


#Calculate true loss of lifetime for scenario 2
times <- seq(0, 15, length.out = 100)
res <- calc.ll(time = times, sim_data = do.call(sim_surv, cases_wei_obs[[2]]))
new_res <- data.frame(res = res, time = times, type = "True comparison measure", stringsAsFactors = F)

#Generate simulations and fit parametric model to compute loss of lifetime
sim_paras <- cases_wei_obs[[2]]
sim_data <- do.call(sim_surv, sim_paras)
fit <- fit.cure.model(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz)
#The parameters are preset to make the point of this example more specific
parameters <- c(-0.9306061, -0.1333783, -0.274995)
fit$coefs[[1]] <- parameters[1]
fit$coefs[[2]] <- parameters[2]
fit$coefs[[3]] <- parameters[3]

#Calculate loss of lifetime with the pre-selected model
res.lol <- calc.LL(fit, time = time.points, rmap = list(year = diag_date), var.type = "n")
D <- data.frame(res = res.lol[[1]]$Estimate, time = time.points)
D$type <- "Estimated comparison measure"

#Combine true and estimate loss of lifetime
plot.data <- rbind(new_res, D)

#Plot the loss of lifetime of each model with added lines for the true and estimate cure points
df_lines1 <- data.frame(x1 = -2, x2 = c(10, 3.2), y1 = c(0.2, 4), 
                        y2 = c(0.2, 4))

df_lines2 <- data.frame(x1 = c(8.6, 10, 3.2, 2.85), x2 = c(8.6, 10, 3.2, 2.85), 
                        y1 = -2, y2 = c(0.2, 0.2, 4, 4))

p <- ggplot(data = plot.data, aes(y = res, x = time, colour = type)) + 
  geom_line() + xlab("Time (years)") + ylab("Comparison measure") + theme_bw() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        axis.title=element_text(size=16),
        strip.text = element_text(size=14), 
        axis.text = element_text(size = 14)) + 
  geom_segment(data = df_lines1, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "dashed", colour = "black") + 
  geom_segment(data = df_lines2, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "dashed", colour = "black") + 
  coord_cartesian(ylim = c(0, max(plot.data$res)), xlim = c(0, 15)) + 
  scale_color_manual(values = cbbPalette, limits = unique(plot.data$type))

#Output to file
if(pdf){
  pdf(file.path(fig.out, "EstimationError.pdf"), width = 9, height = 6) 
} else{
  tiff(file.path(fig.out, "EstimationError.pdf"), res = 300, width = 8 * 300, height = 5 * 300)
}
print(p)
dev.off()
  

#Run simulations
filename <- file.path(data.out, "sim_res.RData")
if(file.exists(filename)){
  load(filename)
}else{
  #Set seed
  set.seed(2, "L'Ecuyer")
  L_wei <- vector("list", length(cases_wei_obs))
  for(i in 1:length(cases_wei_obs)){
    cat(i, "\n")
    L <- mclapply(1:n.sim, function(j){
      #cat(j, "\n")
      #Simulate data
      sim_data <- do.call(sim_surv, cases_wei_obs[[i]])

      #Fit Nelson et al. model
      fit1 <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz, df = 5)
      
      #Knots for Andersson et al. model
      d.times <- sim_data$D$FU_years[sim_data$D$status == 1]
      knots <- sort(c(quantile(d.times, probs = sort(c(getq(6), 0.99)))))
      knots <- log(sort(knots))
      inner.knots <- knots[-c(1, length(knots))]
      bdr.knots <- knots[c(1, length(knots))]
      
      #Fit Andersson et al. model
      fit2 <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz,
                    smooth.formula = ~ nsx(log(FU_years), knots = inner.knots, 
                                           Boundary.knots = bdr.knots, cure = T))
      
      #Fit Jakobsen et al. model
      fit3 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = sim_data$D, df = 4, 
                               bhazard = sim_data$D$exp_haz, verbose = F)
      
      #Fit Weibull mixture cure model
      fit4 <- fit.cure.model(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz)
      
      #The following commented code has been used to make various tests
      #Plot relative survival
      # f <- function(t) 0.3 + (1 - 0.3) * exp(- 0.9 * t ^ 0.8)
      # curve(f, from = 0, to = 10, add = t, col = 2, ylim = c(0, 1))
      # plot(fit1, newdata = data.frame(age = 50), add = T, line.col = 1)
      # cure.pred <- predict(fit2, newdata = data.frame(FU_years = 100))
      # g <- function(t) (predict(fit2, newdata = data.frame(FU_years = t)) - cure.pred) / (1 - cure.pred)
      # curve(g, from = 0.0001, to = 10, add = T, col = 1)
      # plot(fit2, newdata = data.frame(age = 50), add = T, line.col = 3)
      # plot(fit3, add = T, col = 4)
      # plot(fit4, add = T, col = 5)
      # 
      # 
      # #Conditional probability of cancer-related death
      # time.points <- seq(0.001, 10, length.out = 50)
      # # 
      # res_true <- calc.prob(time.points, sim_data = sim_data)
      # plot(res_true ~ time.points, ylim = c(-0.1, 1), type = "l")
      # # # 
      # res1 <- calc.Crude(fit1, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = time.points, link = "identity", var.type = "n")
      # plot(res1, col = 2, add = T)
      # # 
      # res2 <- calc.Crude(fit2, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = time.points, link = "identity", var.type = "n")
      # plot(res2, add = T, col = 3)
      # # 
      # res3 <- calc.Crude(fit3, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = time.points, link = "identity", var.type = "n")
      # plot(res3, add = T, col = 4)
      # 
      # res4 <- calc.Crude(fit4, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = time.points, link = "identity", var.type = "n")
      # plot(res4, add = T, col = 5)
      # legend("topright", c("True LOL", models), fill = 1:5)
      # abline(h = 0.05)
      # # 
      # 
      # res2 <- calc.Crude(fit2, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = seq(0.001, 15, length.out = 50), link = "identity")
      # plot(res2, col = 3)
      # 
      # res3 <- calc.Crude(fit3, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = seq(0.001, 15, length.out = 50), link = "identity")
      # plot(res3, add = T, col = 4)
      # 
      # 
      # # 
      # res_true <- calc.ll(time.points, sim_data = sim_data)
      # plot(res_true ~ time.points, ylim = c(-0.1, max(res_true) + 1), type = "l")
      # 
      # res1 <- calc.LL(fit1, time = seq(0.001, 20, length.out = 50),
      #                 rmap = list(year = diag_date), var.type = "n")
      # plot(res1, add = T, col = 2)
      # 
      # res2 <- calc.LL(fit2, time = seq(0.001, 20, length.out = 50),
      #                 rmap = list(year = diag_date), var.type = "n")
      # plot(res2, add = T, col = 3)
      # 
      # res3 <- calc.LL(fit3, time = seq(0.001, 20, length.out = 50),
      #                 rmap = list(year = diag_date), var.type = "n")
      # plot(res3, add = T, col = 4)
      # 
      # res4 <- calc.LL(fit4, time = seq(0.001, 20, length.out = 50),
      #                 rmap = list(year = diag_date), var.type = "n")
      # plot(res4, add = T, col = 5)
      # legend("topright", c("True LOL", models), fill = 1:5)
      
      # res_true <- calc.prob.cure(time.points, sim_data = sim_data, eps = 0)
      # plot(res_true ~ time.points, ylim = c(-0.1, 1), type = "l")
      # 
      # cure.pred <- predict(fit2, newdata = data.frame(FU_years = 100))
      # a <- 1 - cure.pred / predict(fit2, newdata = data.frame(FU_years = time.points))
      # lines(a ~ time.points, col = 2)
      
      #Conditional probability of cancer-related death
      res1 <- lapply(eps_prob, function(eps){
        calc.Crude.quantile(fit1, q = eps, rmap = list(year = diag_date), 
                            reverse = TRUE, max.time = 80)
      })
      
      res2 <- lapply(eps_prob, function(eps){
        calc.Crude.quantile(fit2, q = eps, rmap = list(year = diag_date), 
                            reverse = TRUE, max.time = 80)
      })
      
      res3 <- lapply(eps_prob, function(eps){
        calc.Crude.quantile(fit3, q = eps, rmap = list(year = diag_date), 
                            reverse = TRUE, max.time = 80)
      })
      
      res4 <- lapply(eps_prob, function(eps){
        calc.Crude.quantile(fit4, q = eps, rmap = list(year = diag_date), 
                            reverse = TRUE, max.time = 80)
      })
      
      res_prob <- rbind(do.call(rbind, res1), do.call(rbind, res2), 
                        do.call(rbind, res3), do.call(rbind, res4))
      res_prob$eps <- rep(eps_prob, length(models))
      res_prob$model <- rep(models, each = length(eps_prob))
      
      
      #Loss of lifetime
      res1 <- lapply(eps_LOL, function(eps){
        calc.LL.quantile(fit1, q = eps, rmap = list(year = diag_date), max.time = 80)
      })
      
      res2 <- lapply(eps_LOL, function(eps){
        calc.LL.quantile(fit2, q = eps, rmap = list(year = diag_date), max.time = 80)
      })
      
      res3 <- lapply(eps_LOL, function(eps){
        calc.LL.quantile(fit3, q = eps, rmap = list(year = diag_date), max.time = 80)
      })
      
      res4 <- lapply(eps_LOL, function(eps){
        calc.LL.quantile(fit4, q = eps, rmap = list(year = diag_date), max.time = 80)
      })
      
      res_lol <- rbind(do.call(rbind, res1), do.call(rbind, res2), 
                       do.call(rbind, res3), do.call(rbind, res4))
      res_lol$eps <- rep(eps_LOL, length(models))
      res_lol$model <- rep(models, each = length(eps_LOL))
      
      
      #Conditional probability of cure
      res2 <- lapply(eps_prob, function(eps){
        calc.cure.quantile(fit2, q = eps, max.time = 80)
      })
      
      res3 <- lapply(eps_prob, function(eps){
        calc.cure.quantile(fit3, q = eps, max.time = 80)
      })
      
      res4 <- lapply(eps_prob, function(eps){
        calc.cure.quantile(fit4, q = eps, max.time = 80)
      })
      
      res_probcure <- rbind(do.call(rbind, res2), do.call(rbind, res3), do.call(rbind, res4))
      res_probcure$eps <- rep(eps_prob, length(models) - 1)
      res_probcure$model <- rep(models[-1], each = length(eps_prob))
      
      #Output results as list
      list(cpProbEstimate = res_prob, cpLOLEstimate = res_lol, cpProbCureEstimate = res_probcure)
    }, mc.cores = n.cores)
    
    #Input estimate cure points into list
    L_wei[[i]]$cpLOLEstimate <- lapply(L, function(x) x$cpLOLEstimate)
    L_wei[[i]]$cpProbEstimate <- lapply(L, function(x) x$cpProbEstimate)
    L_wei[[i]]$cpProbCureEstimate <- lapply(L, function(x) x$cpProbCureEstimate)
    
    #Input true cure point into list
    sim_data <- do.call(sim_surv, cases_wei_obs[[i]])
    L_wei[[i]]$cpProb <- sapply(eps_prob, find_cure_point_prob, sim_data = sim_data)
    L_wei[[i]]$cpLOL <- sapply(eps_LOL, find_cure_point_LOL, sim_data = sim_data)
    L_wei[[i]]$cpProbCure <- sapply(eps_prob, find_cure_point_prob.cure, sim_data = sim_data)
    
    cat("Done with ", i, "\n")
    
  }
  save(L_wei, file = filename)
}

#The bias and length of confidence interval are computed in the following

#Format the probability of cure results
D <- lapply(1:length(L_wei), function(i){
  df <- do.call(rbind, L_wei[[i]]$cpProbCureEstimate)
  df$CIlength <- df$upper.ci - df$lower.ci
  df$bias <- abs(df$Est - L_wei[[i]]$cpProbCure)

  a <- data.frame(eps = eps_prob, True = L_wei[[i]]$cpProbCure)
  rownames(a) <- a$eps
  bias <- aggregate(bias ~ eps*model, data = df, FUN = mean)
  CIlength <- aggregate(CIlength ~ eps*model, data = df, FUN = median)
  
  res <- merge(bias, CIlength, by = c("eps", "model"))
  res$True <- a[as.character(res$eps), "True"]
  od <- c("eps", "model", "bias", "CIlength")
  res <- res[,od]
  D1 <- reshape(res, direction = "wide", timevar = "model", idvar = "eps")
  D1 <- cbind(a, D1[,-1])
  D1 <- cbind(D1[, 1:6], "bias.Nelson et al." = NA, "CIlength.Nelson et al." = NA, D1[, 7:8])
  

  rbind(NA, D1)
})

D1 <- do.call(rbind, D)

#Format the conditional probability of cancer-related death results
D <- lapply(1:length(L_wei), function(i){
  df <- do.call(rbind, L_wei[[i]]$cpProbEstimate)
  df$CIlength <- df$upper.ci - df$lower.ci
  df$bias <- abs(df$Est - L_wei[[i]]$cpProb)

  a <- data.frame(eps = eps_prob, True = L_wei[[i]]$cpProb)
  rownames(a) <- a$eps
  bias <- aggregate(bias ~ eps*model, data = df, FUN = mean)
  CIlength <- aggregate(CIlength ~ eps*model, data = df, FUN = median)
  res <- merge(bias, CIlength, by = c("eps", "model"))
  res$True <- a[as.character(res$eps), "True"]
  od <- c("eps", "model", "bias", "CIlength")
  res <- res[,od]
  D1 <- reshape(res, direction = "wide", timevar = "model", idvar = "eps")
  D1 <- cbind(a, D1[,-1])
  rbind(NA, D1)
})

D2 <- do.call(rbind, D)

#Format the loss of lifetim results
D <- lapply(1:length(L_wei), function(i){
  df <- do.call(rbind, L_wei[[i]]$cpLOLEstimate)
  df$CIlength <- df$upper.ci - df$lower.ci
  df$bias <- abs(df$Estimate - L_wei[[i]]$cpLOL)

  a <- data.frame(eps = eps_LOL, True = L_wei[[i]]$cpLOL)
  rownames(a) <- a$eps
  bias <- aggregate(bias ~ eps*model, data = df, FUN = mean)
  CIlength <- aggregate(CIlength ~ eps*model, data = df, FUN = median)
  res <- merge(bias, CIlength, by = c("eps", "model"))
  res$True <- a[as.character(res$eps), "True"]
  od <- c("eps", "model", "bias", "CIlength")
  res <- res[,od]
  D1 <- reshape(res, direction = "wide", timevar = "model", idvar = "eps")
  D1 <- cbind(a, D1[,-1])
  rbind(NA, D1)
})
D3 <- do.call(rbind, D)
names(D3) <- names(D2)

#Check NA values
df <- do.call(rbind, L_wei[[3]]$cpProbCureEstimate)
wh <- which(is.na(df$lower.ci))
df[wh,]
2 / 500 * 100

#Create table with all results
D.all <- rbind(D1, D2, D3)
D.all <- cbind(Scenario = c("Scenario 1", NA, NA, NA, "Scenario 2", NA, NA, NA, "Scenario 3", NA, NA, NA,
                            "Scenario 1", NA, NA, NA, "Scenario 2", NA, NA, NA, "Scenario 3", NA, NA, NA,
                            "Scenario 1", NA, NA, NA, "Scenario 2", NA, NA, NA, "Scenario 3", NA, NA, NA), 
               D.all)

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 12, 12, 24, 24, 36)
addtorow$command <- c("& & &\\multicolumn{2}{|c|}{ARS model} & \\multicolumn{2}{c|}{FMC model} & \\multicolumn{2}{c}{NRS model} & \\multicolumn{2}{|c}{Weibull mixture}\\\\\n", 
                      "& $\\epsilon$& $t_{\\epsilon}$ & Bias & $\\overline{\\text{CIL}}(\\hat t_\\epsilon)$ & Bias & $\\overline{\\text{CIL}}(\\hat t_\\epsilon)$ & Bias & $\\overline{\\text{CIL}}(\\hat t_\\epsilon)$ & Bias & $\\overline{\\text{CIL}}(\\hat t_\\epsilon)$\\\\\n", 
                      "\\hline\n",
                      "\\multicolumn{11}{l}{Conditional probability of cure}\\\\\n", 
                      "\\hline\n",
                      "\\multicolumn{11}{l}{Conditional probability of cancer-related death}\\\\\n",
                      "\\hline\n",
                      "\\multicolumn{11}{l}{Loss of lifetime}\\\\\n", 
                      "\\hline\n")
print(xtable(D.all, align = "clcccccccccc", label = "tab:simres",
             caption = "Bias and CIL of the cure point estimate in simulated data. 
The cure point estimates were based on the probabiltiy of cure, the probability of cancer-related death, and the loss of lifetime function. 
The NRS model was not evaluated for the latter measure since this is not a cure model.
             ARS: relative survival model by Andersson et al.\\ \\cite{Andersson2011}, 
             FMC: flexible mixture cure model by Jakobsen et al.\\ \\cite{Jakobsen20191}, 
             NRS: relative survival model by Nelson et al.\\ \\cite{Nelson2007}, 
             CIL: confidence interval length.", 
             digits = c(1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)),
      add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE,
      sanitize.text.function = identity, hline.after = c(-1),
      file = paste0(tab.out, "SimulationResults.tex"), scalebox = 0.8, table.placement = "!ht")


# #HTML table
# print(xtable(D.all, align = "clcccccccccccccc", label = "tab:simres",
#              caption = "Bias, variance, and coverage of the cure point estimate in simulated data. 
# The cure point estimates were based on the loss of lifetime function, the probability of cancer-related death, 
# and the probabiltiy of cure. The NRS model was not evaluated for the latter measure since this is not a cure model.
#              ARS: relative survival model by Andersson et al.\\ \\cite{Andersson2011}, 
#              FMC: flexible mixture cure model by Jakobsen et al.\\ \\cite{Jakobsen20172}, 
#              NRS: relative survival model by Nelson et al.\\ \\cite{Nelson2007}, 
#              ECP: empirical coverage probability.", 
#              digits = c(1, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1)),
#       add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE,
#       sanitize.text.function = identity, hline.after = c(-1),
#       file = "C:/Users/sw1y/Dropbox/PhDDefence/SimTable.html", scalebox = 0.55, table.placement = "!ht", type = "html")
# 
# 
