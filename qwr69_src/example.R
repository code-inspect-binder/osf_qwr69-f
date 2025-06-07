# nGFS for sample_data

# read sample_data.Rdata first
# sample_data is a bootstrap dataset from the empirical dataset used in paper 
# Do not use it for any research purpose

amc <- as.data.frame(sample_data[,c('amc')])
names(amc)[1] <- 'amc' #keep x name
avc <- as.data.frame(sample_data[,c('avc')])
names(avc)[1] <- 'avc' #keep y name
cov <- sample_data[,c('prev_diff','time_needed','age','vac')]

# source nGFS
nGFS(amc, avc, cov) # causal target model
# Unconditional model (no covariates) 
# Estimate Std. Error   t value     Pr(>|t|)
# intercept -0.1216921 0.04836083 -2.516336 1.259196e-02
# amc        0.7096270 0.04706807 15.076609 1.760865e-35
# 
# R-squared and Adjusted R-squared: 
#   [1] 0.5150737
# [1] 0.5128076
# 
# Best z index order:  
#   [1] 2 1 3
# 
# Conditional model (with selected z) 
# Estimate Std. Error   t value     Pr(>|t|)
# intercept   -1.07109122 0.63555349 -1.685289 9.341118e-02
# amc          0.55154383 0.05388593 10.235397 3.191420e-20
# time_needed -0.02084193 0.00791595 -2.632904 9.092125e-03
# prev_diff   -0.33323193 0.09933865 -3.354504 9.425919e-04
# age          0.20099141 0.06520545  3.082433 2.327134e-03
# 
# R-squared and Adjusted R-squared: 
# [1] 0.5739466
# [1] 0.5658698

nGFS(avc, amc, cov) # causal alternative model
# Unconditional model (no covariates) 
# Estimate Std. Error    t value     Pr(>|t|)
# intercept -0.02945615 0.04958752 -0.5940235 5.531238e-01
# avc        0.72583722 0.04814327 15.0766089 1.760865e-35
# 
# R-squared and Adjusted R-squared: 
# [1] 0.5150737
# [1] 0.5128076
# 
# No z (covariates) found
