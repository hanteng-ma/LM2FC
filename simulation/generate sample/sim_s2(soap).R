n = 400 #300 & 900

###scalar dim

d_scalar = 2

gamma_k = function(k){
  my_vector = matrix(data = c(cos(2*k*pi/2.7),sin(2*k*pi/3.5)), nrow = 2, ncol = 1)
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
source('D:/研究生文件/科研任务/functional_data/SOAP-main/functions_soap.R')
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(locfit))
suppressPackageStartupMessages(library(fdapace))
suppressPackageStartupMessages(library(lsei))
spline.basis=create.bspline.basis(rangeval=c(0,10),nbasis=10,norder=4)
file_name_x = 'D:/研究生文件/科研任务/functional_data/sim_data/s2_soap/all_data_0309(2)_'
file_name_y = 'D:/研究生文件/科研任务/functional_data/sim_data/s2_soap/y_label_0309(2)_'

file_name_x_t = 'D:/研究生文件/科研任务/functional_data/sim_data/s2_soap/all_data_test_0309(2)_'
file_name_y_t = 'D:/研究生文件/科研任务/functional_data/sim_data/s2_soap/y_label_test_0309(2)_'
spline_basis=create.bspline.basis(rangeval=c(0,10),nbasis=10,norder=4)
SOAP = function(observed, timepoints){
  ########################################
  ################PACE###################
  ########################################
  
  res_pace<- FPCA(observed, timepoints,list(dataType='Sparse',error=TRUE, kernel='epan', 
                                            verbose=FALSE,methodBwCov="GCV",methodBwMu="GCV"))
  select_k = SelectK(res_pace, criterion = 'AIC')$K
  pred_pace = predict(res_pace, observed, timepoints,K=select_k )
  error_pace = c()
  
  coef_mat0 = coef(Data2fd(argvals = res_pace$workGrid,y=res_pace$phi,spline.basis))

  ########################################
  ################SOAP###################
  ########################################
  
  observed%>%do.call(c,.)%>%mean
  pc1s = first_FPC(coef_mat0[,1],observed=observed, timepoints=timepoints,minit=6,gamma=1e1,threshold=1e-3)
  previous_beta = list()
  previous_beta[[1]] = pc1s$beta
  pc2s = third_FPC_conditional(coef_mat0[,2], observed=observed, timepoints=timepoints, pc_index=2, gamma=3e4,betalist =previous_beta,threshold=1e-3)
  previous_beta[[2]] = pc2s$beta
  pc3s = third_FPC_conditional(coef_mat0[,3], observed=observed, timepoints=timepoints, pc_index=3, gamma=2e2,betalist =previous_beta,threshold=1e-3)
  previous_beta[[3]] = pc3s$beta
  #pc4s = third_FPC_conditional(coef_mat0[,4], observed=observed, timepoints=timepoints, pc_index=4, gamma=1e4,betalist =previous_beta,threshold=1e-2)
  #previous_beta[[4]] = pc4s$beta
  
  
  #########
  ###AIC###
  #########
  previous_beta0 = previous_beta
  observed2 = observed[which(sapply(observed,length)>1)]
  timepoints2 = timepoints[which(sapply(observed,length)>1)]
  tempy = observed2
  sd_score	 = c()
  beta_samples= list()
  sigma_est = c()
  for (i in 1:length(previous_beta0)){
    print(i)
    res = pred_SOAP_step(previous_beta0[i],	tempy, timepoints2, spline_basis,sigma=0,nminus=0)
    tempy  =res$residuals
    sigma_est = c(sigma_est ,mean((res$residuals%>%do.call(c,.))^2))
    beta_samples[[i]]=as.numeric(res$sfit)
    sd_score = c(sd_score, res$sfit%>%apply(.,2,sd))
  }
  
  N = sapply(observed2,length)%>%sum
  n = length(observed2)
  AIC = N*log(sigma_est ) + N  + 2*n*c(1:length(previous_beta0))
  (k_selet = as.numeric(which.min(AIC)))
  observed2 = observed[which(sapply(observed,length)==5)]
  timepoints2 = timepoints[which(sapply(observed,length)==5)]
  tempy = observed2
  res = pred_SOAP_step(previous_beta0[1:min(4,k_selet)],	tempy, timepoints2, spline_basis,sigma=0,nminus=0)
  (sig = sqrt(mean((res$residuals%>%do.call(c,.))^2)))
  betalist=previous_beta0
  
  ###############################
  ### estimating the individual scores###
  ###############################
  
  scores_est = function(x) {
    library(locfit)
    (yi=observed[[x]])
    (timei=timepoints[[x]])
    xmat = lapply(1:k_selet,function(i){
      pc_fit  = fd(betalist[[i]], spline_basis)
      eval.fd(timei, pc_fit)%>%as.numeric
    })%>%do.call(cbind,.)
    betai = c()
    llike = function(beta10) {
      liklog = dnorm(yi, as.numeric(xmat%*%beta10), sd=sig,log=TRUE)
      sum(liklog)+ sapply(1:k_selet, function(x) {
        log(density.lf(beta_samples[[x]],ev=beta10[x])$y)
      })%>%sum
    }
    beta10 = sapply(beta_samples[1:k_selet],mean)
    betap = lapply(1:k_selet, function(x) {c()})
    for (i in 1:100){
      for (j in 1:k_selet){
        # betap[[j]]=c()
        beta11 = beta10
        # beta11[j] = sample(size=1,x = beta_samples[[j]])
        beta11[j] = rnorm(1,beta10[j],sd= 0.1*sd(beta_samples[[j]]))
        if (runif(1)<exp(llike(beta11)-llike(beta10))){
          beta10=beta11
        } 
        betap[[j]]  =c(betap[[j]],beta10[j])
      }
    }
    sapply(betap, function(x) {x%>%tail(500)%>%mean})%>%return
  }
  
  ## estimated scores 
  score_pred= mclapply(1:length(observed),scores_est)%>%do.call(rbind,.)
  return(score_pred)
}

for (iter in 1:5){
  ##########实例生成
  
  my_tj = matrix(, ncol = 3)
  my_xi = matrix(data = NA, nrow = n, ncol = 3)
  for (i in 1:n){
    N = sample(c(5,6,7,8,9,10),1)
    random_point = sort(runif(N, 0, 10))
    sample_tj = X_trajectory(random_point)
    
    my_xi[i,] = sample_tj$cof
    
    data_m = matrix(data = NA, nrow = N, ncol = 3)
    data_m[1:N,1] = rep(i, N)            #标号
    data_m[1:N,2] = sample_tj$trajectory + 0.1 * rnorm(N) #轨迹值 + 随机噪声
    data_m[1:N,3] = random_point         #格点位置
    
    my_tj = rbind(my_tj, data_m)
  }
  my_tj = my_tj[-1,]
  
  
  timepoints_list = list()
  observed_list = list()
  
  # 获取唯一的轨迹编号
  unique_ids = unique(my_tj[, 1])
  
  # 对每条轨迹进行处理
  for (id in unique_ids) {
    # 找到该轨迹对应的所有行
    traj_data = my_tj[my_tj[, 1] == id, ]
    
    # 提取轨迹值和格点位置
    timepoints = traj_data[, 3]  # 轨迹的格点位置
    observed = traj_data[, 2]    # 轨迹的观测值
    
    # 将结果存储到对应的列表中
    timepoints_list[[id]] = timepoints
    observed_list[[id]] = observed
  }
  
  # 将列表转换为向量，准备后续分析
  timepoints = timepoints_list
  observed = observed_list
  ##### 生成标量型变量

  
  
  my_s = matrix(data = NA, nrow = n, ncol = d_scalar)
  for (i in 1:n) {
    my_s[i,1] = runif(1,0.5,0.8) * sign(my_xi[i,1])
    my_s[i,2] = runif(1,0.5,0.8) * sign(my_xi[i,2])
  }
  
  ##### 生成标签
  
  y_label = matrix(data = NA, nrow = n, ncol = 1)#y label
  
  for (i in 1:n) {
    score = c(0,0,0)
    for (k in 1:5) {
      score[k] = cos(2*k*pi/2.7) * my_xi[i,1] + sqrt(2)*sin(2*k*pi/3.5) * my_xi[i,2] +((my_s[i,]) %*% (gamma_k(k)))[1,1]
    }
    #print(score)
    y_label[i,1] = which.max(score)
  }
  
  
  
  
  
  
  
  fpcs =  SOAP(observed, timepoints)
  
  

  
  
  all_dim = length(fpcs[1,]) + d_scalar
  
  all_data = matrix(data = NA, nrow = n, ncol = all_dim)
  all_data[,1:length(fpcs[1,])] = fpcs
  all_data[,(length(fpcs[1,])+1):all_dim] = my_s
  

  
  
  write.csv(x = all_data[1:300,], file = paste(file_name_x, iter, '.csv',sep = '') )
  write.csv(x = y_label[1:300], file = paste(file_name_y, iter, '.csv',sep = '') )
  
  write.csv(x = all_data[301:400,], file = paste(file_name_x_t, iter, '.csv',sep = '') )
  write.csv(x = y_label[301:400], file = paste(file_name_y_t, iter, '.csv',sep = '') )
  print('iter')
  print(iter)
}
