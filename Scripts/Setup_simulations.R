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
  
  #Follow-up in days
  D$FU <- D$FU *  ayear
  
  #Follow-up in years
  D$FU_years <- D$FU / ayear
  
  #Get general population hazard
  D$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                           data = D, ratetable = survexp.dk)
  
  #Output in a list
  list(D = D, pars = pars, rel_surv = rel_surv, S_exp = S_exp, haz_fun = haz_fun,
       surv_can_fun = surv_can_fun, dens_can_fun = dens_can_fun)
}


#Plot relative survival scenarios
cases_wei <- list(list(c(pi = 0.2, shape = 1.2, scale = 1)),
                  list(c(pi = 0.3, shape = 0.8, scale = 0.9)),
                  list(c(pi = 0.6, shape = 1.2, scale = 1)))

#Time points used to plot true relative survival
time.points <- seq(0, 15, length.out = 100)

#Plot the true relative survival trajectory
wei <- function(pars) pars[1] + (1 - pars[1]) * exp(-pars[3] * time.points ^ pars[2])
L <- lapply(cases_wei, function(pars) wei(pars[[1]]))
D <- data.frame(surv = do.call(c, L), time.points = rep(time.points, length(cases_wei)), 
                Scenario = rep(1:length(cases_wei), each = length(time.points)), 
                Model = "Weibull")
D$Scenario <- factor(D$Scenario)

p <- ggplot(D, aes(x = time.points, y = surv, group = Scenario, colour = Scenario)) + geom_line(size = 1.1) +
  ylab("Relative survival") + xlab("Time (years)") + theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16),
        axis.title=element_text(size=17),
        axis.text = element_text(size = 15),
        legend.key.size = unit(2,"line")) +
  geom_vline(xintercept = 10, linetype = "dashed") + guides(colour = guide_legend(nrow = 1)) + ylim(0,1) + 
  scale_color_brewer(palette = "Set2")

#Output to file
if(pdf){
  pdf(file.path(fig.out, "Cases.pdf"), width = 9, height = 6) 
} else{
  tiff(file.path(fig.out, "Cases.tiff"), res = 300, width = 8 * 300, height = 5 * 300)
}
print(p)
dev.off()

#Display simulation parameters in table
D <- do.call(rbind, lapply(cases_wei, function(x)x[[1]]))
D <- cbind(Scenario = as.character(1:3), D)

addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- c("Scenario & $\\pi$ & $\\gamma_1$ & $\\gamma_2$\\\\\n")

print(xtable(D, caption = "Parameter values used for simulating survival data.", label = "tab:params", 
             align = "ccccc", digits = 1), include.rownames = F, include.colnames = F, add.to.row = addtorow, 
      santize.text.function = identity, file = file.path(tab.out, "SimParams.tex"))


source("Scripts/SimulationStudy.R")
