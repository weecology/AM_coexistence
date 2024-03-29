## Working code to replicate the results in Xiao & Fussmann

library(foreach)
library(doSNOW)
library(deSolve)
library(lattice)
library(gridExtra)
source('coexistence_module.r')

num_range = function(b1, b2){
# Function to evaluate coexistence in the d1-d2 parameter space that allows
#   each consumer to persist with the resource (rectangular region consisting of
#   subregions A, B, C, and D).
# a1 and a2 are assumed to be 1.
# output - a file with the name "NumRange_b1_b2.csv", recording pairs of d1 and d2
#   allowing coexistence for a specific b1-b2 combination.
  par_without_d = c(a1 = 1, a2 = 1, b1 = b1, b2 = b2)
  out_file = paste('NumRange_', b1, '_', b2, '.csv', sep = "")
  coexist_range(coexist_AM, 1/(b1+1)/100, 1/(b1+1), 1/(b2+1)/100, 1/(b2+1), par_without_d, out_file, t = 10000)
}

num_range_a_b = function(a1, a2, b1, b2){
# Modified num_range() which allows a1 and a2 to vary as well as b1 and b2.
  par_without_d = c(a1 = a1, a2 = a2, b1 = b1, b2 = b2)
  out_file = paste('NumRange_', a1, '_', a2, '_', b1, '_', b2, '.csv', sep = "")
  coexist_range(coexist_AM, a1/(b1+1)/100, a1/(b1+1), a2/(b2+1)/100, a2/(b2+1), par_without_d, out_file, t = 10000)
}

growth_rate = function(a, b, d, x){
# Returns the value of the functional response at x, given functional response
#   parameters a and b, and mortality d
  return (a * x / (1 + b * x) - d)
}

compute_A_pot = function(a1, a2, b1, b2){
  # Function to analytically compute the area of the potential coexistence region
  if ((b1 <= b2) | (b1 <= 1)) {A = 0}
  else {
    A = a1*a2*(b1-1)*(b1+b2)/2/(b1-b2)/(b2+1)/b1^2 + a1*a2/(b1-b2)^2*log((b1*b2+2*b1-b2)/b1/(b1+1))
  }
  return(A)
}

plot_pot_rea = function(b1, b2, x_lim, y_lim, main_lab){
  # Function to overlay the potential coexistence region with realized region
  # x_lim, y_lim - upper limit of x and y axes in plot
  # main_lab - subplot label in the title
  dat = read.csv(paste('NumRange_', b1, '_', b2, '.csv', sep = ''), colClasses = rep('numeric', 2))
  plot(dat$d1, dat$d2, xlim = c(0, x_lim), ylim = c(0, y_lim), type = 'p', pch = 20, col = '#BABABA',
    xlab = expression(d[1]), ylab = expression(d[2]),
    main = bquote(paste(.(main_lab), ' ', b[1], '=', .(b1), ', ', b[2], '=', .(b2), sep = '')), cex.lab = 2,
    cex.axis = 2, xaxs = 'i', yaxs = 'i', cex.main = 2)
  curve(1/((1/x)-b1+b2), from = 0, to = (b1-1)/b1/(b1+1), add = T, lty = 'dashed', lwd = 2)
  curve((b1+1)/(b2+1)*x, from = 0, to = (b1-1)/b1/(b1+1), add = T, lty = 'dashed', lwd = 2)
  segments(0, 1/(b2+1), (b1-1)/b1/(b1+1), 1/(b2+1), lwd = 2)
  segments((b1-1)/b1/(b1+1), 0, (b1-1)/b1/(b1+1), 1/(b2+1), lwd = 2)
  if (b2 > 1){
    text((b1-1)/b1/(b1+1)*0.5, (b2-1)/b2/(b2+1)*0.5, "C", cex = 2)
    text((b1-1)/b1/(b1+1)*0.5, (b2-1)/b2/(b2+1)*0.5 + 1/(1+b2)*0.5, "A", cex = 2)
    segments(0, (b2-1)/b2/(b2+1), (b1-1)/b1/(b1+1), (b2-1)/b2/(b2+1), lwd = 2)
  }
  else {
    text((b1-1)/b1/(b1+1)*0.5, 1/(1+b2)*0.5, "A", cex = 2)
  }
}

###################################
############ Analysis #############
###################################

# Create a data frame with all combinations of (b1, b2),
#   with b2 ranging from 0 to 9.5 in intervals of 0.5, and b1 > b2
b_list = as.data.frame(matrix(nrow = 0, ncol = 2))
names(b_list) = c('b1', 'b2')
for (b2 in seq(0, 9.5, 0.5)){
  for (b1 in seq(b2+0.5, 10, 0.5)){
    b_list = rbind(b_list, c(b1, b2))
  }
}

# Create a data frame with all combinations of (a1, a2),
#   with a1 and a2 varying independently from 1 to 5 in intervals of 0.5
a_list = as.data.frame(matrix(nrow = 0, ncol = 2))
names(a_list) = c('a1', 'a2')
for (a1 in seq(1, 10, 0.5)){
  for (a2 in seq(1, 10, 0.5)){
    a_list = rbind(a_list, c(a1, a2))
  }
}

# Evaluate region of coexistence for all (b1, b2) combinations
# This analysis has been caried out in cluster with 8 cores,
#   as it is very time-consuming.
registerDoSNOW(makeCluster(8, type = "SOCK"))
foreach(i=1:dim(b_list)[1]) %dopar% {
  library(deSolve)
  source('coexistence_module.r')
  num_range(b_list[i, 1], b_list[i, 2])
}

# Evaluate region of coexistence for all (b1, b2) combinations
#   with b1 fixed at 10 and b2 fixed at 0 or 3
registerDoSNOW(makeCluster(8, type = "SOCK"))
foreach(i=1:dim(a_list)[1]) %dopar% {
  library(deSolve)
  source('coexistence_module.r')
  num_range_a_b(a_list[i, 1], a_list[i, 2], 10, 0)
  num_range_a_b(a_list[i, 1], a_list[i, 2], 10, 3)
}

##################################
########## Figures ###############
##################################

#### Figure 1 ####
png(filename = 'Fig1.png', width = 800, height = 800)

a1 = 1
a2 = 1
b1 = 10
b2 = 3

plot(1, 1, xlab = '', ylab = '', xlim = c(0, a1 / (1 + b1)), xaxs = 'i', yaxs = 'i',
  ylim = c(0, a2 / (1 + b2)), type = 'n', xaxt = 'n', yaxt = 'n', cex.lab = 1.5)
title(xlab = expression(d[1]), ylab = expression(d[2]), las = 1, cex.lab = 2,
  line = 1)
abline(v = a1 * (b1 - 1) / b1 / (b1 + 1), lwd = 2)
abline(h = a2 * (b2 - 1) / b2 / (b2 + 1), lwd = 2)
curve(a2 * (b1 + 1) / a1 / (b2 + 1) * x, lty = 'dashed', from = 0,
  to = a1 / (1 + b1), lwd = 2, add = T)
curve(a2 / (a1 / x - b1 + b2), lty = 'dashed', from = 0,
  to = a1 / (1 + b1), lwd = 2, add = T)

text(a1 * (b1 - 1) / b1 / (b1 + 1) / 4, a2 * (b2 - 1) / b2 / (b2 + 1) * 0.5,
  'C', cex = 2)
text(a1 * (b1 - 1) / b1 / (b1 + 1) / 4,
  a2 * (b2 - 1) / b2 / (b2 + 1) * 0.5 + a2 / (1 + b2) * 0.5, 'A', cex = 2)
text(a1 * (b1 - 1) / b1 / (b1 + 1) * 0.5 + a1 / (1 + b1) * 0.5,
  a2 * (b2 - 1) / b2 / (b2 + 1) * 0.5, 'D', cex = 2)
text(a1 * (b1 - 1) / b1 / (b1 + 1) * 0.5 + a1 / (1 + b1) * 0.5,
  a2 * (b2 - 1) / b2 / (b2 + 1) * 0.5 + a2 / (1 + b2) * 0.5, 'B', cex = 2)

#shade potential region
x = seq(0, a1 * (b1 - 1) / b1 / (b1 + 1), a1 * (b1 - 1) / b1 / (b1 + 1) / 1000)
y = a2 / (a1 / x - b1 + b2)
x = c(x, a1 * (b1 - 1) / b1 / (b1 + 1), 0)
y = c(y, a2 * (b1 + 1) / a1 / (b2 + 1) * a1 * (b1 - 1) / b1 / (b1 + 1), 0)
polygon(x, y, density = 30, col = 'black', angle = 90)

dev.off()

#### Figure 2 ####
png(filename = 'Fig2.png', width = 800, height = 800)

b_grid = expand.grid(b1 = do.breaks(c(0, 10), 20), b2 = do.breaks(c(0, 10), 20))
for (i in 1:dim(b_grid)[1]){
    b1 = b_grid[i, 1]
    b2 = b_grid[i, 2]
    A_pot = 0  # Area of potential region of coexistence
    A_rea = 0  # Area of realized region of coexistence
    A_ratio = 0 # ratio between A_rea and A_pot
    A_prop = 0  # proportion of A_rea compared to the full rectangular region
    if ((b1 > b2) & (b1 > 1)){
      A_pot = compute_A_pot(1, 1, b1, b2)
      dat = read.csv(paste('NumRange_', b1, '_', b2, '.csv', sep = ''),
        colClasses = rep('numeric', 2))
      num_points = dim(dat)[1]
      A_tot = 1/(b1+1) * 1/(b2+1) # Area of the rectangular region (A+B+C+D)
      A_prop = num_points / 10000
      A_rea = A_prop * A_tot
      if (num_points >= 50){
        A_ratio = A_rea / A_pot
      }
    }
    b_grid$z_pot[i] = A_pot
    b_grid$z_rea[i] = A_rea
    b_grid$z_ratio[i] = A_ratio
    b_grid$z_prop[i] = A_prop
}

z_lab_list = c('A', 'A', 'Ratio', 'Proportion')
title_list = c('(a) Area of Potential Region of Coexistence',
  '(b) Area of Realized Region of Coexistence',
  '(c) Ratio between Potential and Realized Regions',
  '(d) Proportion of Realized Region\n in Full Rectangular Region'
)

plot1 = wireframe(b_grid[, 3] ~ b1 * b2, b_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(b[1]), cex = 1.5),
  ylab = list(expression(b[2]), cex = 1.5), zlab = list(expression(A[potential]), cex = 1.5, rot = 90),
  zlim = c(0, max(b_grid[, 3])), main ='(a) Area of Potential Region of Coexistence', col.regions = '#EDEDED')

plot2 = wireframe(b_grid[, 4] ~ b1 * b2, b_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(b[1]), cex = 1.5),
  ylab = list(expression(b[2]), cex = 1.5), zlab = list(expression(A[realized]), cex = 1.5, rot = 90),
  zlim = c(0, max(b_grid[, 3])), main ='(b) Area of Realized Region of Coexistence', col.regions = '#EDEDED')

plot3 = wireframe(b_grid[, 5] ~ b1 * b2, b_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(b[1]), cex = 1.5),
  ylab = list(expression(b[2]), cex = 1.5), zlab = list('Ratio', cex = 1.5, rot = 90),
  main ='(c) Ratio between Potential and Realized Regions', col.regions = '#EDEDED')

plot4 = wireframe(b_grid[, 6] ~ b1 * b2, b_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(b[1]), cex = 1.5),
  ylab = list(expression(b[2]), cex = 1.5),
  zlab = list('Proportion', cex = 1.5, rot = 90),
  main = '(d) Proportion of Realized Region in Rectangular Region', col.regions = '#EDEDED')

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)

dev.off()

#### Figure 3 ####
png(filename = 'Fig3.png', width = 800, height = 600)

par(mfcol = c(2, 3), oma = c(0, 1, 0, 0), mar = c(6, 5, 4, 2))
b1_list = c(3, 5, 10)
b2_list = c(1, 1.5, 3)
x_lim = max((b1_list - 1)/b1_list/(b1_list + 1))
y_lim = 1 / (min(b2_list) + 1)
main_lab_list = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')
for (i in 1:length(b1_list)){
  b1 = b1_list[i]
  b2 = 1.5
  plot_pot_rea(b1, b2, x_lim, y_lim, main_lab_list[i])
  b1 = 10
  b2 = b2_list[i]
  plot_pot_rea(b1, b2, x_lim, y_lim, main_lab_list[i + 3])
}

dev.off()

#### Figure 4 ####
png(filename = 'Fig4.png', width = 800, height = 800)

a_grid = expand.grid(a1 = do.breaks(c(1, 5), 8), a2 = do.breaks(c(1, 5), 8))
for (i in 1:dim(a_grid)[1]){
    a1 = a_grid[i, 1]
    a2 = a_grid[i, 2]
    A_ratio = 0
    A_pot = compute_A_pot(a1, a2, 10, 0)
    dat = read.csv(paste('NumRange_', a1, '_', a2, '_10_0.csv', sep = ''),
      colClasses = rep('numeric', 2))
    num_points = dim(dat)[1]
    A_tot = a1/(10+1) * a2/(0+1) # Area of the rectangular region (A+B+C+D)
    A_prop = num_points / 10000
    A_rea = A_prop * A_tot
    if (num_points >= 50){
      A_ratio = A_rea / A_pot
    }
    a_grid$z_pot[i] = A_pot
    a_grid$z_rea[i] = A_rea
    a_grid$z_ratio[i] = A_ratio
    a_grid$z_prop[i] = A_prop
}

z_lab_list = c('A', 'A', 'Ratio', 'Proportion')
title_list = c('(a) Area of Potential Region of Coexistence',
  '(b) Area of Realized Region of Coexistence',
  '(c) Ratio between Potential and Realized Regions',
  '(d) Proportion of Realized Region\n in Full Rectangular Region'
)

plot1 = wireframe(a_grid[, 3] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5), zlab = list(expression(A[potential]), cex = 1.5, rot = 90),
  zlim = c(0, max(a_grid[, 3])), main ='(a) Area of Potential Region of Coexistence', col.regions = '#EDEDED')

plot2 = wireframe(a_grid[, 4] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5), zlab = list(expression(A[realized]), cex = 1.5, rot = 90),
  zlim = c(0, max(a_grid[, 3])), main ='(b) Area of Realized Region of Coexistence', col.regions = '#EDEDED')

plot3 = wireframe(a_grid[, 5] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5), zlab = list('Ratio', cex = 1.5, rot = 90),
  zlim = c(max(a_grid[, 5])/2, max(a_grid[, 5])), main ='(c) Ratio between Potential and Realized Regions', col.regions = '#EDEDED')

plot4 = wireframe(a_grid[, 6] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5),
  zlab = list('Proportion', cex = 1.5, rot = 90),
  zlim = c(max(a_grid[, 6])/2, max(a_grid[, 6])), main = '(d) Proportion of Realized Region in Rectangular Region', col.regions = '#EDEDED')

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)

dev.off()

#### Figure 5 ####
png(filename = 'Fig5.png', width = 800, height = 800)

a_grid = expand.grid(a1 = do.breaks(c(1, 5), 8), a2 = do.breaks(c(1, 5), 8))
for (i in 1:dim(a_grid)[1]){
    a1 = a_grid[i, 1]
    a2 = a_grid[i, 2]
    A_ratio = 0 
    A_pot = compute_A_pot(a1, a2, 10, 3)
    dat = read.csv(paste('NumRange_', a1, '_', a2, '_10_3.csv', sep = ''),
      colClasses = rep('numeric', 2))
    num_points = dim(dat)[1]
    A_tot = a1/(10+1) * a2/(3+1) # Area of the rectangular region (A+B+C+D)
    A_prop = num_points / 10000
    A_rea = A_prop * A_tot
    if (num_points >= 50){
      A_ratio = A_rea / A_pot
    }
    a_grid$z_pot[i] = A_pot
    a_grid$z_rea[i] = A_rea
    a_grid$z_ratio[i] = A_ratio
    a_grid$z_prop[i] = A_prop
}

z_lab_list = c('A', 'A', 'Ratio', 'Proportion')
title_list = c('(a) Area of Potential Region of Coexistence',
  '(b) Area of Realized Region of Coexistence',
  '(c) Ratio between Potential and Realized Regions',
  '(d) Proportion of Realized Region\n in Full Rectangular Region'
)

plot1 = wireframe(a_grid[, 3] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5), zlab = list(expression(A[potential]), cex = 1.5, rot = 90),
  zlim = c(0, max(a_grid[, 3])), main ='(a) Area of Potential Region of Coexistence', col.regions = '#EDEDED')

plot2 = wireframe(a_grid[, 4] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5), zlab = list(expression(A[realized]), cex = 1.5, rot = 90),
  zlim = c(0, max(a_grid[, 3])), main ='(b) Area of Realized Region of Coexistence', col.regions = '#EDEDED')

plot3 = wireframe(a_grid[, 5] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5), zlab = list('Ratio', cex = 1.5, rot = 90),
  zlim = c(max(a_grid[, 5])/2, max(a_grid[, 5])), main ='(c) Ratio between Potential and Realized Regions', col.regions = '#EDEDED')

plot4 = wireframe(a_grid[, 6] ~ a1 * a2, a_grid, scales = list(arrows = F, distance = 1.1),
  screen = list(z = -40, x = -70), xlab = list(expression(a[1]), cex = 1.5),
  ylab = list(expression(a[2]), cex = 1.5),
  zlab = list('Proportion', cex = 1.5, rot = 90),
  zlim = c(max(a_grid[, 6])/2, max(a_grid[, 6])), main = '(d) Proportion of Realized Region in Rectangular Region', col.regions = '#EDEDED')

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)

dev.off()
