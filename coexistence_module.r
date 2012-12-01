## Module of functions used for the analysis in Xiao & Fussmann
library(deSolve)

START_VAL = 10 ^ -5 # Global variable, small starting value when needed

AM_eqns = function(t_list, y, p){
  # A-M system with 2 predators and 1 prey
  # Predators can have Holling type 2 or linear functional response
  #   depending on b1 and b2
  # Prey grows logistically
  # t_list - series of time steps
  # y - initial state (starting point)
  # p - list of parameters in sequence of a1, a2, b1, b2, d1, d2
  # For detailed equations and parameterization, see paper
  dR = y[1]*(1-y[1])-p['a1']*y[1]*y[2]/(1+p['b1']*y[1])-p['a2']*y[1]*y[3]/(1+p['b2']*y[1])
  dP1 = y[2]*(p['a1']*y[1]/(1+p['b1']*y[1])-p['d1'])
  dP2 = y[3]*(p['a2']*y[1]/(1+p['b2']*y[1])-p['d2'])
  return(list(c(dR, dP1, dP2)))
}

AM_jacobian = function(t_list, y, p){
  # Jacobian matrix for A-M system
  # Analytical form is given to facilitate the integration process
  # t_list, y, p - see AM_eqns
  matrix(c(1-2*y[1]-p['a1']*y[2]/(1+p['b1']*y[1])^2-p['a2']*y[3]/(1+p['b2']*y[1])^2,
              y[2]*p['a1']/(1+p['b1']*y[1]), y[3]*p['a2']/(1+p['b2']*y[1]),
            -p['a1']*y[1]/(1+p['b1']*y[1]), p['a1']*y[1]/(1+p['b1']*y[1])-p['d1'], 0,
            -p['a2']*y[1]/(1+p['b2']*y[1]), 0, p['a2']*y[1]/(1+p['b2']*y[1])-p['d2']
            ), 3, 3)
}

type_2_growth = function(species, p, R_list){
  # Function to calculate the average growth rate of a consumer with Holling
  #   Type II functional response given the parameters and the resource level
  #   for a time interval
  # species - species identity, can take two values "P1" or "P2"
  # p - see AM_eqns
  # R_list - resource level in consecutive time steps
  # output - average growth rate
  if (species == "P1"){
    return (as.numeric(mean(p['a1'] * R_list / (1 + p['b1'] * R_list)) - p['d1']))
  }
  else {
    return (as.numeric(mean(p['a2'] * R_list / (1 + p['b2'] * R_list)) - p['d2']))
  }
}

coexist_eval = function(eqns, eqns_jac, t, stepsize, y, p){
  # Function to evaluate whether one consumer is able to invade a system with
  #   the other consumer and the resource
  # eqns - differential equations of growth of focal system
  # eqns_jac - Jacobian matrix of focal system
  # t - number of time steps to run
  # p - see AM_eqns
  # y - initial state (starting point). For this function the starting point of
  #   one and only one of the consumers should be zero.
  # stepsize - step size
  # output - long-term growth rate of the invader
  t_list = seq(0, t, stepsize)
  out = as.data.frame(lsoda(y, t_list, eqns, p, jacfunc = eqns_jac, jactype = "fullusr"))
  while (dim(out[complete.cases(out),])[1] != length(t_list)){ # If integration fails
    stepsize = stepsize / 5  # Reduce step size to 1/5
    t_list = seq(0, t, stepsize)
    out = as.data.frame(lsoda(y, t_list, eqns, p, jacfunc = eqns_jac, jactype = "fullusr"))
  }
  # Compare prey level from the last 1/4 steps to R_star
  # Note that this line assumes that the first three differential equations
  #   defines R, P1 and P2, regardless of the system
  R_list = out[(dim(out)[1] - as.integer(dim(out)[1] / 4)) : dim(out)[1], 2]
  if (y[2] == 0){ # If initial system only consists of R and P2
    return (type_2_growth("P1", p, R_list))
  }
  else {return (type_2_growth("P2", p, R_list))}
}

coexist_AM = function(parameters, t = 2000, stepsize = 0.1){
  # Function to determine if two consumers can coexist in the A-M system, assuming
  #   b1 > b2
  # Coexistence is defined as success mutual invasion
  # parameters - vector of parameters in sequence of a1, a2, b1, b2, d1, d2
  # t, stepsize - see coexist_eval
  # output - binary, 1 if coexist, 0 if at least one of the consumers fails to invade
  par_list = as.numeric(parameters)
  a1 = par_list[1]
  a2 = par_list[2]
  b1 = par_list[3]
  b2 = par_list[4]
  d1 = par_list[5]
  d2 = par_list[6]
  out = 0
  if (d1 <= a1 * (b1 - 1) / b1 / (b1 + 1)){  # Check that point is in Region A or C
    if (d2 >= a2 / (a1 / d1 - b1 + b2)){  # Check that point is above curve R1_star = R2_star
      if (d2 <= a2 * (b1 + 1) * d1 / a1 / (b2 + 1)){ # Check that point is below diagonal line
        y1 = c(START_VAL, START_VAL, 0)
        y2 = c(START_VAL, 0, START_VAL)
        if (d2 >= a2 * (b2 - 1) / b2 / (b2 + 1)){ # If point is in region A
          cond2 = 1 # P1 will always be able to invade, according to Theorem 3.5 in Hsu et al. (1978)
        }
        else {cond2 = coexist_eval(AM_eqns, AM_jacobian, t, stepsize, y2, parameters)}

        cond1 = coexist_eval(AM_eqns, AM_jacobian, t, stepsize, y1, parameters)
        if ((cond1 > 0) & (cond2 > 0)){
          out = 1
        }
      }
    }
  }
  return(out)
}

coexist_range = function(func, d1min, d1max, d2min, d2max, par_without_d, out_file, t = 2000){
  # Function to evaluate coexistence given a series of d1 and d2 values
  # The series of both d1 and d2 consists of 100 equally spaced data points
  # func - function to evaluate coexistence for a specific set of parameters
  # d1min, d1max, d2min, d2max - min and max for the range of d1/d2 to be evaluated
  # par_without_d - other parameters in the system of equations, excluding d1 and d2
  # out_file - name of the output file where the results will be recorded.
  #   It has two columns, d1 and d2, and records the (d1, d2) parameterization
  #   that allows coexistence.
  out = data.frame(t(rep(NA, 2)))
  names(out) = c('d1', 'd2')
  out = out[0, ]
  k = 1
  d1_list = seq(d1min, d1max, by = d1max / 100)
  d2_list = seq(d2min, d2max, by = d2max / 100)
  for (i in 1:length(d1_list)){
    for (j in 1:length(d2_list)){
      par_full = c(par_without_d, d1 = d1_list[i], d2 = d2_list[j])
      if (func(par_full, t) == 1){
        out[k, 1] = d1_list[i]
        out[k, 2] = d2_list[j]
        k = k + 1
        write.csv(out, out_file, row.names = F)
      }
    }
  }
}

plot_dynamics = function(eqns, eqns_jac, p, y, t, stepsize, col_list){
  # Function to plot the dynamics of the system
  # eqns, eqns_jac, t, stepsize - see coexist_eval
  # p, y- see AM_eqns
  # col_list - list of colors to be used for each species.
  #   It needs to be of same length as the number of equations.
  t_list = seq(0, t, stepsize)
  out = as.data.frame(lsoda(y, t_list, eqns, p, jacfunc = eqns_jac, jactype = "fullusr"))
  while (dim(out)[1] != length(t_list)){ # If integration fails, reduce step size
    stepsize = stepsize / 5
    t_list = seq(0, t, stepsize)
    out = as.data.frame(lsoda(y, t_list, eqns, p, jacfunc = eqns_jac, jactype = "fullusr"))
  }
  plot(1, 1, xlab = "Timespan", ylab = "Abundance", xlim = c(0, t),
    ylim = c(0, max(out[100:dim(out)[1], 2:dim(out)[2]])), type = "n")
  for (i in 1:(dim(out)[2] - 1)){
    points(t_list, out[, i + 1], col = col_list[i], type = "l", lwd = 1.5)
  }
}
