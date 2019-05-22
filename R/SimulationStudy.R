#Do simulations

#Function for simulating survival data
sim_surv <- function(pars, age, n){
  #pars is the parameter of the specific relative survival model
  #age is the age used for the general population survival
  #n is the number of samples to simulate
  
  #Generate general population survival. As implemented, age, gender, and year are fixed.
  sim_ages <- rnorm(n, mean = age, sd = 0)
  sim_gender <- factor(rbinom(n, size = 1, prob = 1), levels = 0:1, labels = c("male", "female"))
  dates <- as.Date("1980-01-01") #+ 1:as.numeric(as.Date("1980-01-01") - as.Date("1980-01-01"))
  diag_dates <- sample(dates, n, replace = T)
  D <- data.frame(age = sim_ages * ayear, sex = sim_gender, diag_date = diag_dates)
  S_exp <- survexp(~1, rmap = list(age = age, sex = sex, year = diag_date), 
                   data = D,
                   times = seq(0, (120 - age) * ayear, length.out = 1000),
                   ratetable = survexp.dk, 
                   scale = ayear)
  
  #Function to draw simulations from
  len <- length(S_exp$time)
  sim_fun_gen <- function(x){
    if(x > 1 - S_exp$surv[len]){
      S_exp$time[len]
    }else{
      S_exp$time[findInterval(x, vec = 1 - S_exp$surv) + 1] 
    }
  }
  
  #Set survival of the uncured and relative survival functions
  surv_can_fun <- function(pars, time) exp(-pars[3] * time ^ pars[2]) 
  dens_can_fun <- function(pars, time) exp(-pars[3] * time ^ pars[2]) * pars[3] * pars[2] * time ^ (pars[2] - 1)
  rel_surv <- function(time) pars[1] + (1 - pars[1]) * exp(-pars[3] * time ^ pars[2])
  haz_fun <- function(time) (1 - pars[1]) * dens_can_fun(pars, time) / rel_surv(time) 
  qrel_surv <- function(q){
    res <- rep(Inf, length(q))
    wh <- q < 1 - pars[1]
    res[wh] <- qweibull(q[wh] / (1 - pars[1]), 
                        shape = pars[2], 
                        scale = 1 / pars[3] ^ (1 / pars[2]))
    res
  }
  
  #Simulate uniform variable for general population and disease specific survival
  uni_sim1 <- runif(nrow(D))
  uni_sim2 <- runif(nrow(D))
  
  #Simulate from both distributions
  sim_gen <- sapply(uni_sim1, sim_fun_gen)
  sim_pop <- qrel_surv(uni_sim2)
  D$fu <- pmin(sim_gen, sim_pop)
  
  #Simulate from censoring distribution
  #Set parameters
  max <- 10
  sim_cens <- runif(n = nrow(D), min = 0, max = max)
  
  #Generated follow-up as the minimum of survival time and censoring time
  D$FU <- pmin(D$fu, sim_cens)
  D$status <- as.numeric(D$fu <= sim_cens)
  D$FU[D$FU < 1e-3] <- 1e-3
  # plot(survfit(Surv(FU, status) ~ 1, data = D))
  # curve(rel_surv, from = 0, to = 15, add = T, col = 2)
  #Follow-up in days
  D$FU <- D$FU *  ayear
  # rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
  #                  data = D, ratetable = survexp.dk, method = "ederer2")
  # 
  # plot(rsfit)
  # abline(h = 0.5)
  D$FU_years <- D$FU / ayear
  #Get general population hazard
  D$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                           data = D, ratetable = survexp.dk)
  list(D = D, pars = pars, rel_surv = rel_surv, S_exp = S_exp, haz_fun = haz_fun,
       surv_can_fun = surv_can_fun, dens_can_fun = dens_can_fun)
}


#Plot relative survival scenarios
# cases_wei <- list(list(c(pi = 0.2, shape = 1.2, scale = 1)),
#                   list(c(pi = 0.3, shape = 1.4, scale = 0.5)),
#                   list(c(pi = 0.6, shape = 1, scale = 0.8)))
cases_wei <- list(list(c(pi = 0.2, shape = 1.2, scale = 1)),
                  list(c(pi = 0.3, shape = 1.4, scale = 0.5)),
                  list(c(pi = 0.6, shape = 1.2, scale = 1)))


time.points <- seq(0, 15, length.out = 100)

wei <- function(pars) pars[1] + (1 - pars[1]) * exp(-pars[3] * time.points ^ pars[2])

L <- lapply(cases_wei, function(pars) wei(pars[[1]]))
D_wei <- data.frame(surv = do.call(c, L), time.points = rep(time.points, length(cases_wei)), 
                    Scenario = rep(1:length(cases_wei), each = length(time.points)), 
                    Model = "Weibull")


D <- D_wei
D$Scenario <- factor(D$Scenario)

p <- ggplot(D, aes(x = time.points, y = surv, group = Scenario, linetype = Scenario)) + geom_line() +
  ylab("Relative survival") + xlab("Time") + theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16),
        axis.title=element_text(size=17),
        axis.text = element_text(size = 15),
        legend.key.size = unit(2,"line")) +
  geom_vline(xintercept = 10, linetype = "dashed") + guides(colour = guide_legend(nrow = 1)) + ylim(0,1)


if(pdf){
  pdf(file.path(fig.out, "Cases.pdf"), width = 9, height = 6) 
} else{
  tiff(file.path(fig.out, "Cases.tiff"), res = 300, width = 8 * 300, height = 5 * 300)
}
print(p)
dev.off()


exp_function <- function(t, expected){
  s <- summary(expected, t)
  names(s$surv) <- s$time
  a <- s$surv[as.character(t)]
  names(a) <- NULL
  a
}

gaussxw <- statmod::gauss.quad(100)

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

#Function for calculating the cure point based on the loss of lifetime function
find_cure_point_LOL <- function(sim_data, eps){
  f <- function(time) calc.ll(time, eps = eps, sim_data = sim_data)
  uniroot(f, lower = 0, upper = 30)$root  
}

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

#Function for calculating the cure point based on the conditional probability of cancer related death
find_cure_point_prob <- function(sim_data, eps){
  f <- function(time) calc.prob(time, eps = eps, sim_data = sim_data)
  uniroot(f, lower = 0, upper = 30)$root  
}


calc.prob.cure <- function(time, eps, sim_data){
  pi <- sim_data$pars[1]
  surv <- pi + (1 - pi) * exp(-sim_data$pars[3] * time ^ sim_data$pars[2])
  1 - pi / surv - eps
}

find_cure_point_prob.cure <- function(sim_data, eps){
  f <- function(time) calc.prob.cure(time, eps = eps, sim_data = sim_data)
  uniroot(f, lower = 0, upper = 30)$root
}



getq <- function(k){
  seq(0, 1, length.out = k)
}

#Simulation options
n.sim <- 500
#n.obs <- c(500, 1000, 2000)
n.obs <- 2000
age <- 60

eps_LOL <- 1:3
eps_prob <- c(0.05, 0.1, 0.15)

#Specify simulations
#cases_wei_obs1 <- lapply(cases_wei, function(x) c(x, age = age, n = n.obs[1]))
#cases_wei_obs2 <- lapply(cases_wei, function(x) c(x, age = age, n = n.obs[2]))
#cases_wei_obs3 <- lapply(cases_wei, function(x) c(x, age = age, n = n.obs[3]))
#cases_wei_obs <- c(cases_wei_obs1, cases_wei_obs2, cases_wei_obs3)
cases_wei_obs <- lapply(cases_wei, function(x) c(x, age = age, n = n.obs))
models <- c("Nelson et al.", "Andersson et al.", "Jakobsen et al.", "Weibull mixture")

#Run simulations
filename <- file.path(data.out, "sim_res.RData")
if(file.exists(filename)){
  load(filename)
}else{
  set.seed(2, "L'Ecuyer")
  L_wei <- vector("list", length(cases_wei_obs))
  for(i in 1:length(cases_wei_obs)){
    cat(i, "\n")
    L <- mclapply(1:n.sim, function(j){
      #cat(j, "\n")
      sim_data <- do.call(sim_surv, cases_wei_obs[[i]])
      L <- vector("list", 3)
      n <- 50
      for(i in 1:3){
        cat(i, "\n")
        res <- rep(NA, length.out = n)
        for(j in 1:n){
          sim_data <- do.call(sim_surv, cases_wei_obs[[i]])
          fit3 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = sim_data$D, df = 3, 
                                   bhazard = sim_data$D$exp_haz, verbose = F)
          res[j] <- predict(fit3, type = "curerate")[[1]]$Estimate 
        }
        L[[i]] <- res
      }
      
      L
      
      boxplot(L)
      abline(h = c(0.2, 0.3, 0.6))
      
      mean(L[[1]] - 0.2)
      mean(L[[2]] - 0.3)
      mean(L[[3]] - 0.6)
      
      
      fit <- rs.surv(Surv(FU, status) ~ ratetable(age = age, sex = sex, year = diag_date),
                     ratetable = survexp.dk, data = sim_data$D)
      fit$time <- fit$time / ayear
      plot(fit)
      
      #Fit Jakobsen et al. model
      fit3 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = sim_data$D, df = 3, 
                               bhazard = sim_data$D$exp_haz, verbose = F)
      plot(fit3, add = T)
      
      time <- seq(0.00001, 15, length.out = 100)
      true_rs <- sim_data$rel_surv(time = time)
      plot(true_rs ~ time, type = "l", ylim = c(0, 1), col = 2)
      plot(fit3, add = T, time = time, col = 1)
      
      
      #plot(fit3, time = seq(0.001, 50))
      a <- calc.LL(fit3, time = time, var.type = "n", rmap = list(year = diag_date))
      plot(a, ylim = c(0, 15))
      b <- calc.ll(time = time, sim_data = sim_data)
      lines(time, b, col = 2)
      
      res3 <- lapply(eps_LOL, function(eps){
        calc.LL.quantile(fit3, q = eps, rmap = list(year = diag_date), max.time = 80)
      })
      a <- sapply(res3, function(x) x$Estimate)
      
      b <- sapply(eps_LOL, find_cure_point_LOL, sim_data = sim_data)
      a - b

            
      #Fit Nelson et al. model
      fit1 <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz, df = 4)
      
      #Knots for Andersson et al. model
      d.times <- sim_data$D$FU_years[sim_data$D$status == 1]
      knots <- quantile(d.times, probs = sort(c(getq(5), 0.95)))
      knots[which.max(knots)] <- 10
      #if(knots[length(knots) - 1] < 4){
      #  knots <- c(knots, 6)
      #}
      knots <- log(sort(knots))
      inner.knots <- knots[-c(1, length(knots))]
      bdr.knots <- knots[c(1, length(knots))]
      
      #Fit Andersson et al. model
      fit2 <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz,
                    smooth.formula = ~ nsx(log(FU_years), knots = inner.knots, 
                                           Boundary.knots = bdr.knots, cure = T))
      
      #Fit Jakobsen et al. model
      fit3 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = sim_data$D, df = 3, 
                               bhazard = sim_data$D$exp_haz, verbose = F)
      
      fit4 <- fit.cure.model(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz)
      
      #Plot relative survival
      # plot(fit1, newdata = data.frame(age = 50), add = T, line.col = 2)
      # plot(fit2, newdata = data.frame(age = 50), add = T, line.col = 3)
      # plot(fit3, add = T, col = 4)
      # plot(fit4, add = T, col = 5)
      # 
      # 
      # #Conditional probability of cancer-related death
      # time.points <- seq(0.001, 10, length.out = 50)
      # 
      # res_true <- calc.prob(time.points, sim_data = sim_data)
      # plot(res_true ~ time.points, ylim = c(-0.1, 1), type = "l")
      # 
      # res1 <- calc.Crude(fit1, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = time.points, link = "identity", var.type = "n")
      # plot(res1, col = 2, add = T)
      # 
      # res2 <- calc.Crude(fit2, type = "condother", reverse = T, rmap = list(year = diag_date),
      #                    time = time.points, link = "identity", var.type = "n")
      # plot(res2, add = T, col = 3)
      # 
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
      
      list(cpProbEstimate = res_prob, cpLOLEstimate = res_lol, cpProbCureEstimate = res_probcure)
    }, mc.cores = n.cores)
    
    L_wei[[i]]$cpLOLEstimate <- lapply(L, function(x) x$cpLOLEstimate)
    L_wei[[i]]$cpProbEstimate <- lapply(L, function(x) x$cpProbEstimate)
    L_wei[[i]]$cpProbCureEstimate <- lapply(L, function(x) x$cpProbCureEstimate)
    sim_data <- do.call(sim_surv, cases_wei_obs[[i]])
    L_wei[[i]]$cpProb <- sapply(eps_prob, find_cure_point_prob, sim_data = sim_data)
    L_wei[[i]]$cpLOL <- sapply(eps_LOL, find_cure_point_LOL, sim_data = sim_data)
    L_wei[[i]]$cpProbCure <- sapply(eps_prob, find_cure_point_prob.cure, sim_data = sim_data)
    
    cat("Done with ", i, "\n")
    
  }
  save(L_wei, file = filename)
}


D <- lapply(1:length(L_wei), function(i){
  #Loss of lifetime
  df <- do.call(rbind, L_wei[[i]]$cpLOLEstimate)
  df$bias <- df$Estimate - L_wei[[i]]$cpLOL
  df$var <- df$SE ^ 2
  
  a <- data.frame(eps = eps_LOL, True = L_wei[[i]]$cpLOL)
  b <- aggregate(bias ~ eps*model, data = df, FUN = mean)
  bias <- matrix(b$bias, ncol = 4)
  a <- cbind(a, bias)
  
  for(j in 1:nrow(b)){
    d.new <- df[df$eps == b$eps[j] & df$model == b$model[j],]
    true_cp <- L_wei[[i]]$cpLOL[eps_LOL == b$eps[j]]
    tab <- table(d.new$lower.ci <= true_cp & d.new$upper.ci >= true_cp)
    b$Coverage[j] <- tab["TRUE"] / nrow(d.new) * 100
  }
  
  coverage <- matrix(b$Coverage, ncol = 4)
  
  #df$cidif <- df$upper.ci - df$lower.ci
  conf <- aggregate(var ~ eps*model, data = df, FUN = mean)
  conf <- matrix(conf$var, ncol = 4)
  a <- cbind(a, conf)
  a <- cbind(a, coverage)
  pos <- c(1:3, 7, 11, 4, 8, 12, 5, 9, 13, 6, 10, 14)
  
  #pos <- c(1:3, 7, 4, 8, 5, 9, 6, 10)
  D1 <- a[, pos]
  rbind(NA, D1)
})

D1 <- do.call(rbind, D)


D <- lapply(1:length(L_wei), function(i){
  #Loss of lifetime
  df <- do.call(rbind, L_wei[[i]]$cpProbEstimate)
  df$bias <- df$Est - L_wei[[i]]$cpProb
  df$var <- df$SE ^ 2
  
  a <- data.frame(eps = eps_prob, True = L_wei[[i]]$cpProb)
  b <- aggregate(bias ~ eps*model, data = df, FUN = mean)
  bias <- matrix(b$bias, ncol = 4)
  a <- cbind(a, bias)
  
  for(j in 1:nrow(b)){
    d.new <- df[df$eps == b$eps[j] & df$model == b$model[j],]
    true_cp <- L_wei[[i]]$cpProb[eps_prob == b$eps[j]]
    tab <- table(d.new$lower.ci <= true_cp & d.new$upper.ci >= true_cp)
    b$Coverage[j] <- tab["TRUE"] / nrow(d.new) * 100
  }
  
  coverage <- matrix(b$Coverage, ncol = 4)
  
  #df$cidif <- df$upper.ci - df$lower.ci
  conf <- aggregate(var ~ eps*model, data = df, FUN = mean)
  conf <- matrix(conf$var, ncol = 4)
  a <- cbind(a, conf)
  a <- cbind(a, coverage)
  pos <- c(1:3, 7, 11, 4, 8, 12, 5, 9, 13, 6, 10, 14)
  
  #pos <- c(1:3, 7, 4, 8, 5, 9, 6, 10)
  
  D1 <- a[, pos]
  rbind(NA, D1)
})

D2 <- do.call(rbind, D)


D <- lapply(1:length(L_wei), function(i){
  #Loss of lifetime
  df <- do.call(rbind, L_wei[[i]]$cpProbCureEstimate)
  df$bias <- df$Est - L_wei[[i]]$cpProbCure
  df$var <- df$SE ^ 2
  
  a <- data.frame(eps = eps_prob, True = L_wei[[i]]$cpProbCure)
  b <- aggregate(bias ~ eps*model, data = df, FUN = mean)
  bias <- matrix(b$bias, ncol = length(models) - 1)
  a <- cbind(a, bias)
  
  for(j in 1:nrow(b)){
    d.new <- df[df$eps == b$eps[j] & df$model == b$model[j],]
    true_cp <- L_wei[[i]]$cpProbCure[eps_prob == b$eps[j]]
    tab <- table(d.new$lower.ci <= true_cp & d.new$upper.ci >= true_cp)
    b$Coverage[j] <- tab["TRUE"] / nrow(d.new) * 100
  }
  
  coverage <- matrix(b$Coverage, ncol = length(models) - 1)
  
  #df$cidif <- df$upper.ci - df$lower.ci
  conf <- aggregate(var ~ eps*model, data = df, FUN = mean)
  conf <- matrix(conf$var, ncol = 3)
  a <- cbind(a, conf)
  a <- cbind(a, coverage)
  
  pos1 <- c(1, 2, 3, 6, 9, 4, 7, 10)
  #pos1 <- c(1, 2, 3, 6, 4, 7)
  pos2 <- c(5, 8, 11)
  #pos2 <- c(5, 8)
  
  D1 <- cbind(a[, pos1], matrix(ncol = 3, nrow = nrow(a)), a[, pos2])
  rbind(NA, D1)
})

D3 <- do.call(rbind, D)
names(D3) <- names(D2)

#Check NA values
df <- do.call(rbind, L_wei[[3]]$cpProbCureEstimate)
wh <- which(is.na(df$lower.ci))
df[wh,]
2 / 500 * 100

D.all <- rbind(D1, D2, D3)
D.all <- cbind(Scenario = c("Scenario 1", NA, NA, NA, "Scenario 2", NA, NA, NA, "Scenario 3", NA, NA, NA,
                            "Scenario 1", NA, NA, NA, "Scenario 2", NA, NA, NA, "Scenario 3", NA, NA, NA,
                            "Scenario 1", NA, NA, NA, "Scenario 2", NA, NA, NA, "Scenario 3", NA, NA, NA), 
               D.all)

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 12, 12, 24, 24, 36)
addtorow$command <- c("& & &\\multicolumn{3}{|c|}{ARS model} & \\multicolumn{3}{c|}{FMC model} & \\multicolumn{3}{c}{NRS model} & \\multicolumn{3}{|c}{Weibull mixture}\\\\\n", 
                      "& $\\epsilon$& $t_{\\epsilon}$ & Bias & $\\overline{\\text{Var}}(\\hat t_\\epsilon)$ & ECP $(\\%)$  & Bias & $\\overline{\\text{Var}}(\\hat t_\\epsilon)$ & ECP $(\\%)$ & Bias & $\\overline{\\text{Var}}(\\hat t_\\epsilon)$ & ECP $(\\%)$ & Bias & $\\overline{\\text{Var}}(\\hat t_\\epsilon)$ & ECP $(\\%)$\\\\\n", 
                      "\\hline\n",
                      "\\multicolumn{11}{l}{Loss of lifetime}\\\\\n", 
                      "\\hline\n",
                      "\\multicolumn{11}{l}{Conditional probability of cancer-related death}\\\\\n",
                      "\\hline\n",
                      "\\multicolumn{11}{l}{Conditional probability of cure}\\\\\n", 
                      "\\hline\n")
print(xtable(D.all, align = "clcccccccccccccc", label = "tab:simres",
             caption = "Bias, variance, and coverage of the cure point estimate in simulated data. 
The cure point estimates were based on the loss of lifetime function, the probability of cancer-related death, 
and the probabiltiy of cure. The NRS model was not evaluated for the latter measure since this is not a cure model.
             ARS: relative survival model by Andersson et al.\\ \\cite{Andersson2011}, 
             FMC: flexible mixture cure model by Jakobsen et al.\\ \\cite{Jakobsen20172}, 
             NRS: relative survival model by Nelson et al.\\ \\cite{Nelson2007}, 
             ECP: empirical coverage probability.", 
             digits = c(1, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1)),
      add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE,
      sanitize.text.function = identity, hline.after = c(-1),
      file = paste0(tab.out, "SimulationResults.tex"), scalebox = 0.55, table.placement = "!ht")


#HTML table
print(xtable(D.all, align = "clcccccccccccccc", label = "tab:simres",
             caption = "Bias, variance, and coverage of the cure point estimate in simulated data. 
The cure point estimates were based on the loss of lifetime function, the probability of cancer-related death, 
and the probabiltiy of cure. The NRS model was not evaluated for the latter measure since this is not a cure model.
             ARS: relative survival model by Andersson et al.\\ \\cite{Andersson2011}, 
             FMC: flexible mixture cure model by Jakobsen et al.\\ \\cite{Jakobsen20172}, 
             NRS: relative survival model by Nelson et al.\\ \\cite{Nelson2007}, 
             ECP: empirical coverage probability.", 
             digits = c(1, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1)),
      add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE,
      sanitize.text.function = identity, hline.after = c(-1),
      file = "C:/Users/sw1y/Dropbox/PhDDefence/SimTable.html", scalebox = 0.55, table.placement = "!ht", type = "html")


