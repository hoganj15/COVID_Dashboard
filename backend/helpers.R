library(deSolve)
library(tidyverse)

countries = read.csv('https://raw.githubusercontent.com/datasets/covid-19/main/data/countries-aggregated.csv') %>% rename(Total = Confirmed)
demographics = read.csv("https://raw.githubusercontent.com/hoganj15/MMA_Assignment_Data/main/COVID/demographics.csv")


######FUNCTION 1 - basic_SIR
#takes a country, and optionally beta and gamma initial values, and L (from Juan's second lecture)
SIR_basic = function(country_chosen, beta = 0.5, gamma = 0.5, L = 10){
  #filter country dataset to include only chosen country
  #add general columns for analysis
  data = countries %>% filter(Country == country_chosen)
  data$new_cases = replace_na((data$Total - lag(data$Total)), 0) #number of new cases over the previous dat
  data = data[min(which(data$new_cases >= 5)):length(data$Country), ] #drop all rows before we see 5+ new cases in a day for the first time
  data$N = rep(as.numeric(demographics %>% filter(Country == country_chosen) %>% select(N)), length(data$Country)) #pulls country population from demographics dataset and creates column from it
  data$day = as.numeric(as.Date(data$Date) - min(as.Date(data$Date))) + 1 #day counter
  #calculate SIR metrics
  data$susceptible = data$N - data$Total
  data$infected = data$Total - data$Recovered - data$Deaths #active cases are all cases minus ones that have recovered or died
  data$removed = data$Recovered + data$Deaths
  
  #create function of differential equations to find rate of infection for first 10 days
  
  Infected = data$infected[1:L]
  Day = data$day[1:L]
  N = as.numeric(demographics %>% filter(Country == country_chosen) %>% select(N))
  
  SIR = function(time, state, parameters){
    par = as.list(c(state, parameters))
    with(par, {
      dS = -beta/N * I * S
      dI = beta/N * I * S - gamma * I
      dR = gamma * I
      list(c(dS, dI, dR))
    })
  }

  #set up initial conditions and fit data to minimize RSS (find curve of best fit)
  
  init = c(S = N-Infected[1], I = Infected[1], R = 0)
  RSS = function(parameters){
    names(parameters) = c("beta", "gamma")
    out = ode(y = init, times = Day, func = SIR, parms = parameters)
    fit = out[ , 3]
    sum((Infected - fit)^2)
  }

  #find beta and gamma which best fir the data
  
  Opt = optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1), control = list(parscale=c(10^-4, 10^-4)))
  Opt_par = setNames(Opt$par, c("beta", "gamma"))

  #calculate R0
  
  R0 = setNames(Opt_par["beta"] / Opt_par["gamma"], "R0")
  data$beta = rep(Opt_par["beta"], length(data$Country))
  data$gamma = rep(Opt_par["gamma"], length(data$Country))
  data$R0 = rep(R0, length(data$Country))
  
  #getting predictions from the ODE
  t = 1:length(data$Country)
  fit = data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
  data$pred_S = fit$S
  data$pred_I = fit$I
  data$pred_R = fit$R
  
  print(paste("The predicted peak is", max(fit$I), "cases occurring on", as.Date(min(data$Date))+which.max(fit$I)))
  #PRR
  a = max(data$infected)
  b = max(data$pred_I)
  PRR = (b - a)/b
  print(paste("Peak Reduction Ratio:", PRR))
  #PDD
  c = which.max(data$infected)
  d = which.max(data$pred_I)
  PDD = c-d
  print(paste("Peak Delay Days:", PDD))
  #CRR
  peak_day = which.max(data$infected) #real peak
  c = data$Total[peak_day] #total cases at the peak infections
  d = data$pred_I[peak_day]+data$pred_R[peak_day] #peak day in worst-case
  CRR = c/d
  print(paste("Case Reduction Ratio:", CRR))
  
  return(data)
}

#########FUNCTION 2 - DEMOGRAPHIC SIR MODEL
#accounts for immigration and emigration
#takes as input the beta and gamma values (calculated above), alpha - the immigration rate, mu - the emigration rate, and number of days to simulate
SIR_demographic = function(Beta, Gamma, Alpha = 1/80, Mu = 1/100, num_days = 1000){
  
  times = seq(0, num_days, by = 1/10)
  parms = c(alpha = Alpha, mu = Mu, N = 1, beta = Beta, gamma = Gamma)
  start = c(S = 0.999, I = 0.001, R = 0)

  sirmod_demographic = function(t, y, parms) {
    # Pull state variables from y vector
    S = y[1]
    I = y[2]
    R = y[3]
    alpha=parms["alpha"]
    beta = parms["beta"]
    mu = parms["mu"]
    gamma = parms["gamma"]
    N = parms["N"]
    ###DEFINE EQUATIONS
    dS = alpha - mu*S - beta * S * I/N
    dI = beta * S * I/N - (mu + gamma) * I
    dR = gamma * I - mu * R
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  }
  
  out = ode(y = start, times = times, func = sirmod_demographic, parms = parms)
  out = as.data.frame(out)
  
  #long-term equilibria
  alpha=parms[1]
  mu=parms[2]
  beta=parms[4]
  gamma=parms[5]
  
  R0=as.numeric(alpha*beta/(mu*(mu+gamma)))
  S_eq=as.numeric((gamma+mu)/beta)
  I_eq=as.numeric((mu/beta)*(R0-1))
  R_eq=as.numeric((gamma/beta)*(R0-1))
  
  out$R0 = rep(R0, length(out$time))
  out$S_eq = rep(S_eq, length(out$time))
  out$I_eq = rep(I_eq, length(out$time))
  out$R_eq = rep(R_eq, length(out$time))
  
  print(paste("R0 in this model:", R0))
  print(paste("Equilibrium S:", S_eq))
  print(paste("Equilibrium I:", I_eq))
  print(paste("Equilibrium R:", R_eq))
  
  return(out)
}

#########FUNCTION 3 - SEIR Model with Vital Dynamics
# 1-k$ of the immigrants arrive infected (other than that it is identical to above)
#takes as input the beta and gamma values, K - percentage of immigrants that arrive infected, alpha - the immigration rate, mu - the emigration rate, and the number of days to simulate
SEIR_vitals = function(Beta, Gamma, K = 0.05, Alpha = 1/80, Mu = 1/100, num_days = 1000){
  times = seq(0, num_days, by = 1/10)
  parms = c(k = K, alpha = Alpha, mu = Mu, N = 1, beta = Beta, gamma = Gamma)
  start = c(S = 0.999, I = 0.001, R = 0)
  
  sirmod_vital = function(t, y, parms) {
    # Pull state variables from y vector
    S = y[1]
    I = y[2]
    R = y[3]
    k = parms["k"]
    alpha=parms["alpha"]
    beta = parms["beta"]
    mu = parms["mu"]
    gamma = parms["gamma"]
    N = parms["N"]
    ###DEFINE EQUATIONS
    dS = alpha * k - mu*S - beta * S * I/N
    dI = alpha * (1-k) + beta * S * I/N - (mu + gamma) * I
    dR = gamma * I - mu * R
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  }
  
  out = ode(y = start, times = times, func = sirmod_vital, parms = parms)
  out = as.data.frame(out)
  
  #long-term equilibria
  alpha=parms[1]
  mu=parms[2]
  beta=parms[4]
  gamma=parms[5]
  
  R0=as.numeric(alpha*beta/(mu*(mu+gamma)))
  S_eq=as.numeric((gamma+mu)/beta)
  I_eq=as.numeric((mu/beta)*(R0-1))
  R_eq=as.numeric((gamma/beta)*(R0-1))
  
  out$R0 = rep(R0, length(out$time))
  out$S_eq = rep(S_eq, length(out$time))
  out$I_eq = rep(I_eq, length(out$time))
  out$R_eq = rep(R_eq, length(out$time))
  
  print(paste("R0 in this model:", R0))
  print(paste("Equilibrium S:", S_eq))
  print(paste("Equilibrium I:", I_eq))
  print(paste("Equilibrium R:", R_eq))
  
  return(out)
}


#########FUNCTION 4 - SIRS MODEL
#takes as input the beta and gamma values, the number of days someone is immune after getting the infection, and the number of days to simulate
SIRS_basic = function(Beta, Gamma, immunity = 200, num_days = 1000){
  times = seq(0, num_days, by = 1/10)
  parms = c(N = 1, beta = Beta , gamma = Gamma, xi= 1/immunity)
  start = c(S = 0.999, I = 0.001, R = 0)
  
  sirsmod_basic = function(t, y, parms) {
    # Pull state variables from y vector
    S = y[1]
    I = y[2]
    R = y[3]
    # Pull parameter values from parms vector
    beta = parms["beta"]
    gamma = parms["gamma"]
    xi = parms["xi"]
    N = parms["N"]
    
    # Define equations
    dS = - beta * S * I/N + xi*R
    dI = beta * S * I/N - ( gamma) * I
    dR = gamma * I - xi*R
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  }
  
  out = ode(y = start, times = times, func = sirsmod_basic, parms = parms)
  out = as.data.frame(out)
  return(out)
}

#########FUNCTION 5 - SIRS MODEL with vaccines
#takes as input the beta and gamma values, the percentage of the population being vaccinated every day, the vaccine effectiveness rate, alpha - the immigration rate, mu - the emigration rate, and the number of days to simulate
SIRS_vaccine = function(Beta, Gamma, vaccination_rate = 0.01, vaccination_eff = 0.95, Alpha = 1/80, Mu = 1/100, num_days = 1000){
  times = seq(0, num_days, by = 1/10)
  parms = c(alpha = Alpha, mu = Mu, N = 1, beta = Beta , gamma = Gamma, rho = vaccination_rate*vaccination_eff)
  start = c(S = 0.999, I = 0.001, R = 0)
  
  sirsmod_vaccine = function(t, y, parms) {
    # Pull state variables from y vector
    S = y[1]
    I = y[2]
    R = y[3]
    # Pull parameter values from parms vector
    alpha=parms["alpha"]
    beta = parms["beta"]
    mu = parms["mu"]
    gamma = parms["gamma"]
    rho = parms["rho"]
    N = parms["N"]
    
    # Define equations
    dS = alpha - mu*S - beta * S * I/N-rho*S
    dI = beta * S * I/N - (mu + gamma) * I
    dR = gamma * I - mu * R+rho*S
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  }
  
  out = ode(y = start, times = times, func = sirsmod_vaccine, parms = parms)
  out = as.data.frame(out)
  return(out)
}


#########FUNCTION 6 - SEASONAL SIR MODEL
#takes as input the gamma value, the transmission rate, the seasonality force, alpha - the immigration rate, mu - the emigration rate, and the number of MONTHS to simulate
SIR_seasonal = function(Gamma, transmission = 20, seasonality = 0.5, Alpha = 1/80, Mu = 1/100, num_months = 200){
  times  = seq(0, num_months, by=1/100)
  params = c(alpha = Alpha, mu = Mu, N = 1, beta0 = transmission, beta1 = seasonality,  gamma = Gamma)
  start = c(S = 0.999,  I = 0.001, R = 0)
  
  sirmod_seasonal = function(t, y, parms){
    S = y[1]
    I = y[2] #function of time
    R = y[3]
    with(as.list(parms),{
      dS =alpha- mu * S  - beta0 * (1+beta1 * cos(2 * pi * t)) * S * I / N
      dI = beta0 * (1 + beta1 * cos(2*pi * t)) *  S * I / N - (mu + gamma) * I
      dR = gamma * I - mu * R
      res = c(dS, dI, dR)
      list(res)
    })
  }
  
  out = as.data.frame(ode(start, times, sirmod_seasonal, params))
  
  alpha=params["alpha"]
  beta0=params["beta0"]
  mu=params["mu"]
  gamma=params["gamma"]
  
  R0=alpha*beta0/(mu*(mu+gamma))
  Istar=  mu*(1-1/R0)*R0/beta0
  as.numeric(R0)
  as.numeric(Istar)
  
  print(paste("Long-term R0:", R0))
  print(paste("Long-term Infection Rate:", Istar))
  
  return(out)
}




#FUNCTIONS FOR PRR, PDD, CRR, as well
#take dataframe as input which has been processed like above
PRR = function(data){
  a = max(data$infected)
  b = max(data$pred_I)
  PRR = (b - a)/b
  return(PRR)
}
PDD = function(data){
  c = which.max(data$infected)
  d = which.max(data$pred_I)
  PDD = c-d
  return(PDD)
}
CRR = function(data){
  peak_day = which.max(data$infected) #real peak
  c = data$Total[peak_day] #total cases at the peak infections
  d = data$pred_I[peak_day]+data$pred_R[peak_day] #peak day in worst-case
  CRR = c/d
  return(CRR)
}


#FUNCTIONS FOR VISUALIZATION
plot_predictions_basic = function(data){
  fig = plot_ly(data, x = ~day)
  fig = fig %>% add_trace(y = ~infected, mode = "lines", name = "Actual")
  fig = fig %>% add_trace(y = ~pred_I, mode = "lines", name = "Predicted")
  fig = fig %>% layout(title = "Basic SIR Predictions", hovermode = 'x unified', xaxis = list(title = "Day"), yaxis = list(title = "Infected"))
  return(fig)
}

plot_predictions_advanced = function(country_chosen, out, model_type = ""){
  data = countries %>% filter(Country == country_chosen)
  data$new_cases = replace_na((data$Total - lag(data$Total)), 0) 
  data = data[min(which(data$new_cases >= 5)):length(data$Country), ] 
  data$N = rep(as.numeric(demographics %>% filter(Country == country_chosen) %>% select(N)), length(data$Country)) 
  data$day = as.numeric(as.Date(data$Date) - min(as.Date(data$Date))) + 1 
  
  data$susceptible = data$N - data$Total
  data$infected = data$Total - data$Recovered - data$Deaths 
  data$removed = data$Recovered + data$Deaths
  
  fig = plot_ly()
  fig = fig %>% add_trace(data = data, x = ~day, y = ~infected, mode = "lines", name = "Actual", hovertemplate = "Day: %{x}<br>Infected: %{y:.3f}")
  fig = fig %>% add_trace(data = out[0:length(data$Country)*10, ], x = ~time, y = ~I*data$N[1], mode = "lines", name = "Predicted", hovertemplate = "Day: %{x}<br>Predicted: %{y:.3f}")
  fig = fig %>% layout(title = paste(model_type, "Model Prediction vs. Actual Data"), xaxis = list(title = "Day"), yaxis = list(title = "Infected"))
  return(fig)
}

plot_predictions_seasonal = function(country_chosen, out){
  data = countries %>% filter(Country == country_chosen)
  data$new_cases = replace_na((data$Total - lag(data$Total)), 0) 
  data = data[min(which(data$new_cases >= 5)):length(data$Country), ] 
  data$N = rep(as.numeric(demographics %>% filter(Country == country_chosen) %>% select(N)), length(data$Country)) 
  data$day = as.numeric(as.Date(data$Date) - min(as.Date(data$Date))) + 1 
  
  data$susceptible = data$N - data$Total
  data$infected = data$Total - data$Recovered - data$Deaths 
  data$removed = data$Recovered + data$Deaths
  
  fig = plot_ly()
  fig = fig %>% add_trace(data = data, x = ~day, y = ~infected, mode = "lines", name = "Actual", hovertemplate = "Day: %{x}<br>Infected: %{y:.3f}")
  fig = fig %>% add_trace(data = out, x = ~time*30, y = ~I*data$N[1], mode = "lines", name = "Predicted", hovertemplate = "Day: %{x}<br>Predicted: %{y:.3f}")
  fig = fig %>% layout(title = "Seasonal SIRS Model Prediction vs. Actual Data", xaxis = list(title = "Day"), yaxis = list(title = "Infected"))
  return(fig)
}


######Function for metrics
plotting_metrics = function(country_chosen){
  data = countries %>% filter(Country == country_chosen)
  data$new_cases = replace_na((data$Total - lag(data$Total)), 0) #number of new cases over the previous dat
  data = data[min(which(data$new_cases >= 5)):length(data$Country), ] #drop all rows before we see 5+ new cases in a day for the first time
  data$N = rep(as.numeric(demographics %>% filter(Country == country_chosen) %>% select(N)), length(data$Country)) #pulls country population from demographics dataset and creates column from it
  data$day = as.numeric(as.Date(data$Date) - min(as.Date(data$Date))) + 1 #day counter
  #calculate SIR metrics
  data$susceptible = data$N - data$Total
  data$infected = data$Total - data$Recovered - data$Deaths #active cases are all cases minus ones that have recovered or died
  data$removed = data$Recovered + data$Deaths
  #PP
  data$point_prevalence = data$infected/data$N
  # Hazard Rate
  Haz_numerator = (data$Total-lag(data$Total))/data$N
  Haz_denominator = 1-(data$Total/data$N)
  data$hazard = Haz_numerator/Haz_denominator
  # Case Fatality
  data$case_fatality = data$Deaths/data$Total
  # Outcome Fatality
  data$outcome_fatality = data$Deaths/(data$Recovered + data$Deaths)
  return(data)
}



