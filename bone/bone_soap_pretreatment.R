library(loon.data)

source('D:/研究生文件/科研任务/functional_data/SOAP-main/functions_soap.R')
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(locfit))
suppressPackageStartupMessages(library(fdapace))
suppressPackageStartupMessages(library(lsei))
spline.basis=create.bspline.basis(rangeval=c(0,1),nbasis=10,norder=4)
spline_basis=create.bspline.basis(rangeval=c(0,1),nbasis=10,norder=4)
SOAP = function(observed, timepoints){
  ########################################
  ################PACE###################
  ########################################
  spline.basis=create.bspline.basis(rangeval=c(0,1),nbasis=10,norder=4)
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




data("bone_ext")

id_total = max(bone_ext$idnum)
age_l = min(bone_ext$age) - 0.01
age_m = max(bone_ext$age)




my_tj = matrix(, ncol = 3)
my_tj_test = matrix(, ncol = 3)

y_label = c()
y_label_test = c()
my_scalar = matrix(, ncol = 1)
my_scalar_test = matrix(, ncol = 1)

miss_value = 0
miss_value_test = 0
for(i in 1:id_total){
  
  if((i %% 3 == 0) | (i %% 3 == 1) ){
    one_bone = bone_ext[bone_ext[,1]==i,]
    N = length(one_bone[,1])
    
    if(N<=1){
      miss_value = miss_value + 1
      print('miss:')
      print(i)
      next
    }
    
    data_m = matrix(data = NA, nrow = N, ncol = 3)
    data_m[1:N,1] = one_bone$idnum            #标号
    data_m[1:N,2] = one_bone$spnbmd    #轨迹值
    data_m[1:N,3] = (one_bone$age - age_l)/(age_m - age_l) 
    
    my_tj = rbind(my_tj, data_m)
    
    sex = one_bone$sex[1]
    label = one_bone$ethnic[1]
    if(sex == 'male'){
      my_scalar = rbind(my_scalar,matrix(data = 0.00001, 1,1)) ######
    }
    if(sex == 'female'){
      my_scalar = rbind(my_scalar,matrix(data = 0, 1,1))
    }
    
    if(label == 'Asian'){
      y_label = c(y_label,1)
    }
    if(label == 'Black'){
      y_label = c(y_label,2)
    }
    if(label == 'Hispanic'){
      y_label = c(y_label,3)
    }
    if(label == 'White'){
      y_label = c(y_label,4)
    }
  }
  
  if(i %% 3 == 2){
    one_bone = bone_ext[bone_ext[,1]==i,]
    N = length(one_bone[,1])
    
    if(N<=1){
      miss_value_test = miss_value_test + 1
      print('miss:')
      print(i)
      next
    }
    
    data_m = matrix(data = NA, nrow = N, ncol = 3)
    data_m[1:N,1] = one_bone$idnum            #标号
    data_m[1:N,2] = one_bone$spnbmd    #轨迹值
    data_m[1:N,3] = (one_bone$age - age_l)/(age_m - age_l)  
    
    my_tj_test = rbind(my_tj_test, data_m)
    
    sex = one_bone$sex[1]
    label = one_bone$ethnic[1]
    if(sex == 'male'){
      my_scalar_test = rbind(my_scalar_test,matrix(data = 0.00001, 1,1)) #####
    }
    if(sex == 'female'){
      my_scalar_test = rbind(my_scalar_test,matrix(data = 0, 1,1))
    }
    
    if(label == 'Asian'){
      y_label_test = c(y_label_test,1)
    }
    if(label == 'Black'){
      y_label_test = c(y_label_test,2)
    }
    if(label == 'Hispanic'){
      y_label_test = c(y_label_test,3)
    }
    if(label == 'White'){
      y_label_test = c(y_label_test,4)
    }
  }
 
}
my_tj = my_tj[-1,]
my_tj_test = my_tj_test[-1,]
my_scalar = my_scalar[-1,]
my_scalar_test = my_scalar_test[-1,]

##### 绘出样本轨迹
par(mfrow=c(3,3))
for (i in 51:63){
  if((i %% 3 == 0) | (i %% 3 == 1) ){
    id<-i ##for curve i
    t.c<-my_tj[my_tj[,1]==id,3] ##measurement points
    y.c<-my_tj[my_tj[,1]==id,2] ##obs
    #plots
    plot(t.c,y.c,xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id), type = 'p')
    ##points(t.c,y.pred.proj,col=2, pch=2) ##predicted measurements at observed measurement times
  }
  
}




### KL expansion 求解轨迹的\xi
library(sm)
library(splines)
library(fpca)

timepoints_list = list()
observed_list = list()

my_tj = rbind(my_tj, my_tj_test)

# 获取唯一的轨迹编号
unique_ids = unique(my_tj[, 1])

# 对每条轨迹进行处理
id_s = 0
for (id in unique_ids) {
  id_s = id_s + 1
  # 找到该轨迹对应的所有行
  traj_data = my_tj[my_tj[, 1] == id, ]
  
  # 提取轨迹值和格点位置
  timepoints = traj_data[, 3]  # 轨迹的格点位置
  observed = traj_data[, 2]    # 轨迹的观测值
  
  # 将结果存储到对应的列表中
  timepoints_list[[id_s]] = timepoints
  observed_list[[id_s]] = observed
}

# 将列表转换为向量，准备后续分析
timepoints = timepoints_list
observed = observed_list


res_pace<- FPCA(observed, timepoints,list(dataType='Sparse',error=TRUE, 
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
fpcs = score_pred




n = length(y_label)
n_test = length(y_label_test)
d_scalar = 1
d_f = length(fpcs[1,])

all_dim = d_f + d_scalar
all_data = matrix(data = NA, nrow = n, ncol = all_dim)
all_data[,1:d_f] = fpcs[1:n,]
all_data[,(d_f+1):all_dim] = my_scalar




all_data_test = matrix(data = NA, nrow = n_test, ncol = all_dim)
all_data_test[,1:d_f] = fpcs[(n+1):(n+n_test),]
all_data_test[,(d_f+1):all_dim] = my_scalar_test


write.csv(x = all_data, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_data_train_soap_0610.csv')
write.csv(x = y_label, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_y_train_soap_0610.csv')

write.csv(x = all_data_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_data_test_soap_0610.csv')
write.csv(x = y_label_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_y_test_soap_0610.csv')

