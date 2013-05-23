## Check that R has stabilized from timestep 7500 to 10000
source('coexistence_module.r')

START_VAL = 10 ^ -5

b_list = seq(1.5, 10, 0.5)
outfile = as.data.frame(matrix(nrow = 0, ncol = 14))
names(outfile) = c('b', 'd', '8000l', '8500l', '9000l', '9500l', '10000l',
  '8000h', '8500h', '9000h', '9500h', '10000h', 'cv_low', 'cv_high')
t = 10000
irow = 1
for (b in b_list){
  d_list = seq(1/(b+1)/100, 1/(b+1), 1/(b+1)/100)
  for (d in d_list){
    if (d < (b-1) / b / (b+1)){
      outfile[irow, 1] = b
      outfile[irow, 2] = d
      stepsize = 0.1
      t_list = seq(0, t, stepsize)
      y = c(START_VAL, START_VAL, 0)
      p = c(a1 = 1, a2 = 1, b1 = b, b2 = 0, d1 = d, d2 = d)
      out = as.data.frame(lsoda(y, t_list, AM_eqns, p, jacfunc = AM_jacobian, jactype = "fullusr"))
      while (dim(out[complete.cases(out),])[1] != length(t_list)){ # If integration fails, reduce step size
        stepsize = stepsize / 5
        t_list = seq(0, t, stepsize)
        out = as.data.frame(lsoda(y, t_list, AM_eqns, p, jacfunc = AM_jacobian, jactype = "fullusr"))
      }
      for (i in 1:5){
        R_list = out[which(t_list >= (7000+500*i) & t_list < (7500+500*i)), 2]
        outfile[irow, i + 2] = min(R_list)
        outfile[irow, i + 7] = max(R_list)
      }
    outfile$cv_low[irow] = sd(as.numeric(outfile[irow, 3:7])) / mean(as.numeric(outfile[irow, 3:7]))
    outfile$cv_high[irow] = sd(as.numeric(outfile[irow, 8:12])) / mean(as.numeric(outfile[irow, 8:12]))
    irow = irow + 1
    }
  write.csv(outfile, 'R_dynamics.csv', row.names = F, quote = F)
  }
}

# Plot dynamics of b-d combinations that has cv larger than 0.05
R_dynamics = read.csv('R_dynamics.csv')
R_dynamics_unstable = R_dynamics[which((abs(R_dynamics$cv_high) > 0.05) | (abs(R_dynamics$cv_low) > 0.05)), ]

pdf_path = 'R_dynamics.pdf'
pdf(file = pdf_path)
par(mfrow = c(3, 2))
for (i in 1:dim(R_dynamics_unstable)[1]){
  b = R_dynamics_unstable$b[i]
  d = R_dynamics_unstable$d[i]
  p = c(a1 = 1, a2 = 1, b1 = 0, b2 = b, d1 = d, d2 = d)
  y = c(START_VAL, 0, START_VAL)
  plot_dynamics(AM_eqns, AM_jacobian, p, y, 10000, 0.1, c('green', 'blue', 'black'))
  title(paste('b = ', b, ', d = ', round(d, digits = 5)))
}
dev.off()


