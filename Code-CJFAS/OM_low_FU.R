#devtools::install() # install wham and whamMSE from the github repo (https://github.com/lichengxue/MSE_NAA_Project)library(wham)
library(wham)
library(whamMSE)
year_start <- 1993 # starting year in the burn-in period
year_end<- 2022 # end year in the burn-in period
MSE_years  <- 30 # number of years in the feedback loop
# Note: no need to include MSE_years in simulation-estimation 

main.dir = here::here()
setwd(main.dir)
sub.dir <- "Scenario5"
dir.create(file.path(getwd(),sub.dir))

info <- generate_basic_info(n_stocks  = 2, 
                            n_regions = 2, 
                            n_indices = 2, 
                            n_fleets  = 2, 
                            n_seasons = 4, 
                            base.years = year_start:year_end, 
                            n_feedback_years = MSE_years, 
                            life_history = "medium", 
                            n_ages= 12, 
                            Fbar_ages = 12, 
                            recruit_model = 2, 
                            F_info = list(F.year1 = 0.2, Fhist = "F-H-L", Fmax = 2.5, Fmin = 2.5, change_time = 0.5),
                            catch_info = list(catch_cv = 0.1, catch_Neff = 100), 
                            index_info = list(index_cv = 0.1, index_Neff = 100, fracyr_indices = 0.625, q = 0.2), 
                            fracyr_spawn = 0.625, 
                            bias.correct.process = FALSE, 
                            bias.correct.observation = FALSE, 
                            bias.correct.BRPs= FALSE, 
                            mig_type = 0) 

basic_info = info$basic_info # collect basic information
catch_info = info$catch_info # collect fleet catch information
index_info = info$index_info # collect survey information
F_info = info$F # collect fishing information

basic_info <- generate_NAA_where(basic_info = basic_info, move.type = 1) # "bidirectional" movement

move <- generate_move(basic_info = basic_info, 
                      move.type = 1,
                      move.rate = c(0.1, 0), 
                      move.re = "ar1_y",
                      move.sigma = 0.5,
                      prior.sigma = 0.5,
                      move.rho_a = 0.5,
                      move.rho_y = 0.5,
                      use.prior = FALSE)
move$year_re[2,] = "ar1"

n_stocks <- as.integer(basic_info['n_stocks'])
n_regions <- as.integer(basic_info['n_regions'])
n_fleets <- as.integer(basic_info['n_fleets'])
n_indices <- as.integer(basic_info['n_indices'])
n_ages<- as.integer(basic_info['n_ages'])

# Selectivity Configuration
fleet1_pars <- c(5.5, 1)
fleet2_pars <- c(5, 1)
index_pars <- c(2, 1)
sel <- list(model=rep("logistic",n_fleets+n_indices),
            initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))

# M Configuration
M <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))

sigma <- "rec+1"
re_cor <- "iid"
ini.opt<- "equilibrium" # option  <- c("age-specific-fe", "equilibrium")

# Set para. for B-H function
alpha1 <- 7
beta1 <- 1e-4

alpha2 <- 5.5
beta2 <- 1e-4

# Set sigma for NAA
NAA_sig <- 0.2
sigma_vals = array(NAA_sig, dim = c(n_stocks, n_regions, n_ages)) # n_stocks x n_regions x n_ages"
sigma_vals[,,1] = 0.5

# Set initial NAA for each stock
log_N1 <- rep(10, n_stocks) # Create difference between stocks
N1_pars <- generate_ini_N1(log_N1,basic_info,ini.opt)

NAA_re <- list(N1_model=rep(ini.opt,n_stocks),
               sigma=rep(sigma,n_stocks),
               cor=rep(re_cor,n_stocks),
               recruit_model = 3,
               recruit_pars = list(c(alpha1,beta1),c(alpha2,beta2)), # assume same B-H s-r functions for all stocks
               sigma_vals = sigma_vals,
               N1_pars = N1_pars,
               NAA_where = basic_info$NAA_where)

input <- prepare_wham_input(basic_info = basic_info, selectivity = sel, M = M, NAA_re = NAA_re, move = move,
                            catch_info = catch_info, index_info = index_info, F = F_info)

random = input$random # check what processes are random effects
input$random = NULL # so inner optimization won't change simulated RE
om <- fit_wham(input, do.fit = F, do.brps = T, MakeADFun.silent = TRUE)

assess.interval <- 3 # Note: assessment interval is 3 years, given the feedback period is 3 years, there will be only 1 assessment
base.years <- year_start:year_end # Burn-in period
first.year <- head(base.years,1)
terminal.year  <- tail(base.years,1)
assess.years<- seq(terminal.year, tail(om$years,1)-assess.interval,by = assess.interval)

# PAN_NAA EM
library(doParallel)
library(foreach)
detectCores() # check how many cores available
cluster <- makeCluster(50) 
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 1
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(rep(list(fleet1_pars),n_fleets),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model="equilibrium",sigma="rec+1",cor="iid")
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = NULL,
                         age_comp_em = "multinomial", # Note: not multinational!
                         em.opt = list(separate.em = TRUE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval,
                         base_years = base.years,
                         year.use = 20, 
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod1_%03d.RDS",i)))
  
}

stopCluster(cluster)

# FAA_NAA EM
library(doParallel)
library(foreach)
cluster <- makeCluster(50) 
registerDoParallel(cluster)


foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = 1
  n_fleets = n_indices = 2
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model="equilibrium",sigma="rec+1",cor="iid")
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = NULL,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = TRUE, separate.em.type = 2, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod2_%03d.RDS",i)))
  
}

stopCluster(cluster)

# SEP_NAA EM
library(doParallel)
library(foreach)
cluster <- makeCluster(50) 
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 1
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(rep(list(fleet1_pars),n_fleets),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model="equilibrium",sigma="rec+1",cor="iid")
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data,
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = NULL,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = TRUE, separate.em.type = 3, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod3_%03d.RDS",i)))
}

stopCluster(cluster)

# SpD_NAA EM
cluster <- makeCluster(50) 
registerDoParallel(cluster)


foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 2
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model=rep("equilibrium",n_stocks),
                    sigma=rep("rec+1",n_stocks),
                    cor=rep("iid",n_stocks))
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = move,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = FALSE, separate.em.type = NULL, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod4_%03d.RDS",i)))
}

stopCluster(cluster)

# SpE_Fix EM
cluster <- makeCluster(50) 
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 2
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model=rep("equilibrium",n_stocks),
                    sigma=rep("rec+1",n_stocks),
                    cor=rep("iid",n_stocks))
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = move,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = FALSE, separate.em.type = NULL, do.move = TRUE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod5_%03d.RDS",i)))
  
}

stopCluster(cluster)

# SpE_Est EM
cluster <- makeCluster(50)
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 2
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model=rep("equilibrium",n_stocks),
                    sigma=rep("rec+1",n_stocks),
                    cor=rep("iid",n_stocks))
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  move <- generate_move(basic_info = basic_info, 
                        move.type = 1,
                        move.rate = c(0.1, 0), 
                        move.re = "ar1_y",
                        move.sigma = 0.5,
                        prior.sigma = 0.5,
                        move.rho_a = 0.5,
                        move.rho_y = 0.5,
                        use.prior = TRUE)
  
  mod <- loop_through_fn(om = om_with_data,
                         em_info = info,
                         random = random,
                         M_em = M_em,
                         sel_em = sel_em,
                         NAA_re_em = NAA_re_em,
                         move_em = move,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = FALSE, separate.em.type = NULL, do.move = TRUE, est.move = TRUE),
                         assess_years = assess.years,
                         assess_interval = assess.interval,
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod6_%03d.RDS",i)))
  
}

stopCluster(cluster)

# PAN_noNAA EM
library(doParallel)
library(foreach)
detectCores() # check how many cores available
cluster <- makeCluster(50) 
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 1
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(rep(list(fleet1_pars),n_fleets),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model="equilibrium",sigma="rec",cor="iid")
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = NULL,
                         age_comp_em = "multinomial", # Note: not multinational!
                         em.opt = list(separate.em = TRUE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, 
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod7_%03d.RDS",i)))
  
}

stopCluster(cluster)

# FAA_noNAA EM
library(doParallel)
library(foreach)
cluster <- makeCluster(50) 
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = 1
  n_fleets = n_indices = 2
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model="equilibrium",sigma="rec",cor="iid")
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  
  mod <- loop_through_fn(om = om_with_data,
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = NULL,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = TRUE, separate.em.type = 2, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod8_%03d.RDS",i)))
  
}

stopCluster(cluster)

# SEP_noNAA EM
library(doParallel)
library(foreach)
cluster <- makeCluster(50) 
registerDoParallel(cluster)

foreach (i = 1:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 1
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(rep(list(fleet1_pars),n_fleets),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model="equilibrium",sigma="rec",cor="iid")
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = NULL,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = TRUE, separate.em.type = 3, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)
  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod9_%03d.RDS",i)))
}

stopCluster(cluster)

# SpE_noNAA EM
library(doParallel)
library(foreach)
cluster <- makeCluster(50) 
registerDoParallel(cluster)


foreach (i = 51:100) %dopar% {
  
  library(wham)
  library(whamMSE)
  
  om_with_data <- update_om_fn(om, seed = 123+i, random = random)
  
  n_stocks = n_regions = n_fleets = n_indices = 2
  
  sel_em <- list(model=rep("logistic",n_fleets+n_indices),
                 initial_pars=c(list(fleet1_pars, fleet2_pars),rep(list(index_pars),n_indices)))
  
  NAA_re_em <- list(N1_model=rep("equilibrium",n_stocks),
                    sigma=rep("rec",n_stocks),
                    cor=rep("iid",n_stocks))
  
  M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))
  
  mod <- loop_through_fn(om = om_with_data, 
                         em_info = info,
                         random = random,
                         M_em = M_em, 
                         sel_em = sel_em, 
                         NAA_re_em = NAA_re_em, 
                         move_em = move,
                         age_comp_em = "multinomial",
                         em.opt = list(separate.em = FALSE, separate.em.type = NULL, do.move = FALSE, est.move = FALSE),
                         assess_years = assess.years, 
                         assess_interval = assess.interval, 
                         base_years = base.years,
                         year.use = 20, # number of years of data you want to use in the assessment model
                         seed = 123+i,
                         save.sdrep = FALSE,
                         save.last.em = TRUE)  
  saveRDS(mod,file.path(sub.dir,sprintf("Mod10_%03d.RDS",i)))
}

stopCluster(cluster)
