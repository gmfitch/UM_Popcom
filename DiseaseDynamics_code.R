### EXERCISE 1: SIR MODEL ###

SIR_funct = function(N = 1, I0, nu, beta, tEnd, printCode=FALSE)
{
  if(printCode)
  {
    cat(sprintf("SIR_funct(N=%d, I_0=%.2f, beta=%.1f, nu=%d, tEnd=%.1f)\n", N, I_0, beta, nu, tEnd))
  }
  
  S0 = N - I0
  out = sirDeterministic(N=N, S=S0, I=I0, beta=beta, nu=nu, tEnd=tEnd)
  
  return(sirPlot(N, out))
}

##############


sirDeterministic = function(
  # Total population size
  N = 100,
  
  # Initial in susceptible & infected classes
  # All susceptible except 1 if unspecified
  S0 = NULL,
  I0 = NULL,
  
  # End time
  tEnd = 10,
  
  # Duration of infection
  nu = 50,
  
  # Contact rate
  beta = 500
)
{
  require(deSolve)
  
  if(is.null(S0))
  {
    S0 = N - 1
  }
  if(is.null(I0))
  {
    I0 = 1
  }
  
  params = list(N=N, nu=nu, beta=beta)
  initialConditions = c(S=S0, I=I0)
  times = seq(0, tEnd, length.out=1000)
  
  # Define d/dt
  ddtFunc = function(t, states, params)
  {
    with(as.list(c(states, params)),
{
  dSI = beta * S * I
  dIR = I * nu
  
  if(is.nan(dSI) || is.nan(dIR))
  {
    dSI=0
    dIR=0
  }
  
  dSdt = -dSI
  dIdt = dSI - dIR
  
  return(list(c(dSdt, dIdt)))
})
  }

# Simulate ODE
out = as.data.frame(
  ode(
    y = initialConditions,
    times = times,
    func = ddtFunc,
    parms = params,
    method = "ode45"
  )
)

t = out$t
S = out$S
I = out$I
R = N - S - I

return(list(t=t, S=S, I=I, R=R))
}

##########

sirPlot = function(N, out)
{
  ylim = c(0, N)
  plot(out$t, out$S, type='l', lwd=2, lty = 1, col="red", ylim=ylim, xlab='Time (days)', ylab='Relative Population Density', main='S I R Model')
  lines(out$t, out$I, lwd=2, lty = 2, col= "blue")
  lines(out$t, out$R, lwd=2, lty = 3, col= "purple")
  legend("topright", bty="n", c("S","I","R"), col = c('red', 'blue', 'purple'),  lwd = 2, lty = c(1,2,3))
  
  return(invisible(out))
}

#################################
## SIR with demographics

### EXERCISE 2 ###

SIR_dem_funct = function(N = 1, I0, nu, beta, mu, tEnd)
{
  
  S0 = N - I0
  R0 = 0
  out = sir_dem_Deterministic(N=N, S=S0, I=I0, R = R0, beta=beta, nu=nu, mu=mu, tEnd=tEnd)
  
  return(sir_dem_Plot(N, out))
}

##############


sir_dem_Deterministic = function(
  # Total population size
  N = 1,
  
  # Initial in susceptible & infected classes
  # All susceptible except 1 if unspecified
  S0 = NULL,
  I0 = NULL,
  R0 = NULL,
  
  # End time
  tEnd = 10,
  
  # Duration of infection
  nu = 50,
  
  # Contact rate
  beta = 500,
  
  # mu
  mu = 1/(20*365)
  
)
{
  require(deSolve)
  
  if(is.null(S0))
  {
    S0 = N - 1
  }
  if(is.null(I0))
  {
    I0 = 1
  }
  if(is.null(R0))
  {
    R0 = 0
  }
  
  params = list(N=N, nu=nu, beta=beta, mu=mu)
  initialConditions = c(S=S0, I=I0, R=R0)
  times = seq(0, tEnd, length.out=1000)
  
  # Define d/dt
  ddtFunc = function(t, states, params)
  {
    with(as.list(c(states, params)),
{
  
  dSdt = mu*N - beta * S * I / N - mu * S
  dIdt = beta * S * I / N - I * nu - mu * I
  dRdt = I * nu - mu * R
  
  return(list(c(dSdt, dIdt, dRdt)))
})
  }

# Simulate ODE
out = as.data.frame(
  ode(
    y = initialConditions,
    times = times,
    func = ddtFunc,
    parms = params,
    method = "ode45"
  )
)

t = out$t
S = out$S
I = out$I
R = out$R

return(list(t=t, S=S, I=I, R=R))
}

###########

sir_dem_Plot = function(N, out)
{
  ylim = c(0, N)
  plot(out$t, out$S, type='l', lwd=2, lty = 1, col="red", ylim=ylim, xlab='Time (days)', ylab='Relative Population Density', main='S I R Model with births and deaths')
  lines(out$t, out$I, lwd=2, lty = 2, col= "blue")
  lines(out$t, out$R, lwd=2, lty = 3, col= "purple")
  legend("topright", bty="n", c("S","I","R"), col = c('red', 'blue', 'purple'),  lwd = 2, lty = c(1,2,3))
  
  cat("Variables at the end of the simulation: S =", tail(out$S, n=1), ", I =", tail(out$I, n=1), ", R =", round(tail(out$R, n=1)), "\n")

  #return((out))
}

