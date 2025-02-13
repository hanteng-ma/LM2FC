n = 900

###scalar dim

d_scalar = 5

gamma_k = function(k){
  my_vector = matrix(data = c(cos(2*k*pi/2.7),sin(2*k*pi/3.5),cos(2*k*pi/5),sin(2*k*pi/5),sign(sin(k-2.5))), nrow = 5, ncol = 1)
  return(my_vector)
}

X_trajectory = function(t){
  mu_tj = t + sin(t)
  xi_random = c(sample(x = c(rnorm(1,-2,1),rnorm(1,2,1)), size = 1), sample(x = c(rnorm(1,-2,1),rnorm(1,2,1)), size = 1), sample(x = c(rnorm(1,-2,1),rnorm(1,2,1)), size = 1))
  eg_tj = xi_random[1] * ((-1)*sqrt(1/5)*cos(pi*t/5)) + xi_random[2] * (sqrt(1/5)*sin(pi*t/5)) + xi_random[3] * ((-1)*sqrt(1/5)*cos(2*pi*t/5))
  sum_f = mu_tj + eg_tj
  df = list(trajectory = sum_f, cof = xi_random)
  return(df)
}


### KL expansion 求解轨迹的\xi
library(sm)
library(splines)
library(fpca)

file_name_x = 'D:/研究生文件/科研任务/functional_data/sim_data/s5/all_data_0527(5)_'
file_name_y = 'D:/研究生文件/科研任务/functional_data/sim_data/s5/y_label_0527(5)_'

file_name_x_t = 'D:/研究生文件/科研任务/functional_data/sim_data/s5/all_data_test_0527(5)_'
file_name_y_t = 'D:/研究生文件/科研任务/functional_data/sim_data/s5/y_label_test_0527(5)_'

for (iter in 1:100){
  ##########实例生成
  
  my_tj = matrix(, ncol = 3)
  my_xi = matrix(data = NA, nrow = n, ncol = 3)
  for (i in 1:n){
    N = sample(c(5,6,7,8,9,10),1)
    random_point = runif(N, 0, 10)
    sample_tj = X_trajectory(random_point)
    
    my_xi[i,] = sample_tj$cof
    
    data_m = matrix(data = NA, nrow = N, ncol = 3)
    data_m[1:N,1] = rep(i, N)            #标号
    data_m[1:N,2] = sample_tj$trajectory + 0.1 * rnorm(N) #轨迹值 + 随机噪声
    data_m[1:N,3] = random_point         #格点位置
    
    my_tj = rbind(my_tj, data_m)
  }
  my_tj = my_tj[-1,]
  
  ##### 生成标量型变量
  
  
  my_s = matrix(data = NA, nrow = n, ncol = d_scalar)
  for (i in 1:n) {
    my_s[i,1] = runif(1,0.5,0.8) * sign(my_xi[i,1])
    my_s[i,2] = runif(1,0.5,0.8) * sign(my_xi[i,2])
    my_s[i,3] = sample(x = c(rnorm(1,-1,0.25),rnorm(1,1,0.25)), size = 1)
    my_s[i,4] = sample(x = c(rnorm(1,-1,0.25),rnorm(1,1,0.25)), size = 1)
    my_s[i,5] = sample(x = c(rnorm(1,-1,0.25),rnorm(1,1,0.25)), size = 1)
  }
  
  ##### 生成标签
  
  y_label = matrix(data = NA, nrow = n, ncol = 1)#y label
  
  for (i in 1:n) {
    score = c(0,0,0)
    for (k in 1:5) {
      score[k] = cos(2*k*pi/2.7) * my_xi[i,1] + sqrt(2)*sin(2*k*pi/3.5) * my_xi[i,2] +((my_s[i,]) %*% (gamma_k(k)))[1,1] + rnorm(1) ### noise
     
    }
    #print(score)
    y_label[i,1] = which.max(score)
  }
  
  
  ##### test data set 
  
  my_tj_test = matrix(, ncol = 3)
  my_xi_test = matrix(data = NA, nrow = n, ncol = 3)
  for (i in 1:n){
    N = sample(c(5,6,7,8,9,10),1)
    random_point = runif(N, 0, 10)
    sample_tj = X_trajectory(random_point)
    
    my_xi_test[i,] = sample_tj$cof
    
    data_m = matrix(data = NA, nrow = N, ncol = 3)
    data_m[1:N,1] = rep(i, N)            #标号
    data_m[1:N,2] = sample_tj$trajectory + 0.1 * rnorm(N) #轨迹值 + 随机噪声
    data_m[1:N,3] = random_point         #格点位置
    
    my_tj_test = rbind(my_tj_test, data_m)
  }
  my_tj_test = my_tj_test[-1,]
  
  my_s_test = matrix(data = NA, nrow = n, ncol = d_scalar)
  for (i in 1:n) {
    my_s_test[i,1] = runif(1,0.5,0.8) * sign(my_xi_test[i,1])
    my_s_test[i,2] = runif(1,0.5,0.8) * sign(my_xi_test[i,2])
    my_s_test[i,3] = sample(x = c(rnorm(1,-1,0.25),rnorm(1,1,0.25)), size = 1)
    my_s_test[i,4] = sample(x = c(rnorm(1,-1,0.25),rnorm(1,1,0.25)), size = 1)
    my_s_test[i,5] = sample(x = c(rnorm(1,-1,0.25),rnorm(1,1,0.25)), size = 1)
  }
  
  ##### 生成标签
  
  y_label_test = matrix(data = NA, nrow = n, ncol = 1)#y label test
  
  for (i in 1:n) {
    score = c(0,0,0)
    for (k in 1:5) {
      score[k] = cos(2*k*pi/2.7) * my_xi_test[i,1] + sqrt(2)*sin(2*k*pi/3.5) * my_xi_test[i,2] +((my_s_test[i,]) %*% (gamma_k(k)))[1,1] + rnorm(1) ### noise
    }
    y_label_test[i,1] = which.max(score)
  }
  
  ### 对grid归一化
  my_tj_1 = my_tj
  my_tj_1[,3] = my_tj[,3]/10
  
  my_tj_test_1 = my_tj_test
  my_tj_test_1[,3] = my_tj_test[,3]/10
  
  ## candidate models for fitting
  M.set<-c(12,13)
  r.set<-c(4,5,6)
  
  ##parameters for fpca.mle
  ini.method="EM"
  basis.method="bs"
  sl.v=rep(0.5,10)
  max.step=50
  grid.l=seq(0,1,0.01)
  grids=seq(0,1,0.005) #for ini.method = EM
  
  ##fit candidate models by fpca.mle
  result<-fpca.mle(my_tj_1, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
  summary(result)
  
  eigenfest<-result$eigenfunctions
  evalest<-result$eigenvalues
  M<-result$selected_model[1]
  r<-result$selected_model[2]
  grids.new<-result$grid
  sig2est<-result$error_var
  muest<-result$fitted_mean
  
  
  fpcs<-fpca.score(my_tj_1,grids.new,muest,evalest,eigenfest,sig2est,r)
  pred<-fpca.pred(fpcs, muest,eigenfest)
  
  fpcs_test<-fpca.score(my_tj_test_1,grids.new,muest,evalest,eigenfest,sig2est,r)
  pred_test<-fpca.pred(fpcs_test, muest,eigenfest)
  
  N<-length(grids.new)
  
  
  
  all_dim = length(eigenfest[,1]) + d_scalar
  
  all_data = matrix(data = NA, nrow = n, ncol = all_dim)
  all_data[,1:length(eigenfest[,1])] = fpcs
  all_data[,(length(eigenfest[,1])+1):all_dim] = my_s
  
  all_data_test = matrix(data = NA, nrow = n, ncol = all_dim)
  all_data_test[,1:length(eigenfest[,1])] = fpcs_test
  all_data_test[,(length(eigenfest[,1])+1):all_dim] = my_s_test
  
  
  write.csv(x = all_data, file = paste(file_name_x, iter, '.csv',sep = '') )
  write.csv(x = y_label, file = paste(file_name_y, iter, '.csv',sep = '') )
  
  write.csv(x = all_data_test, file = paste(file_name_x_t, iter, '.csv',sep = '') )
  write.csv(x = y_label_test, file = paste(file_name_y_t, iter, '.csv',sep = '') )
  
  print(iter)
}


