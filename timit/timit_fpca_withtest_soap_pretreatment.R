ix_w_train = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/ix_wave_train.txt',header = F)
kcl_w_train = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/kcl_wave_train.txt',header = F)
ao_w_train = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/ao_wave_train.txt',header = F)
ae_w_train = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/ae_wave_train.txt',header = F)


ix_w_test = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/ix_wave_test.txt',header = F)
kcl_w_test = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/kcl_wave_test.txt',header = F)
ao_w_test = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/ao_wave_test.txt',header = F)
ae_w_test = read.csv2(file = 'D:/研究生文件/科研任务/functional_data/py_code/ae_wave_test.txt',header = F)

y_label = rep(1:4,c(length(ix_w_train[,]), length(kcl_w_train[,]), length(ao_w_train[,]), length(ae_w_train[,]) )  )

y_label_test = rep(1:4,c(length(ix_w_test[,]), length(kcl_w_test[,]), length(ao_w_test[,]), length(ae_w_test[,]) )  )

#my_data = rbind(ix_d, dcl_d, ao_d, sh_d)
my_data = rbind(ix_w_train, kcl_w_train, ao_w_train, ae_w_train)
my_data_test = rbind(ix_w_test, kcl_w_test, ao_w_test, ae_w_test)

my_tj = matrix(, ncol = 3)
my_tj_test = matrix(, ncol = 3)

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




n = length(my_data[,])

for(i in 1:length(my_data[,])){
  one_tj = as.double(unlist(strsplit(my_data[i,], split = ',')))
  N = length(one_tj) - 1
  print(c(i,N))
  point_position = (floor(seq(1, N + 1, 1) * (5000/(N + 1))) / 5000)[1:N]
  
  data_m = matrix(data = NA, nrow = N, ncol = 3)
  data_m[1:N,1] = rep(i, N)            #标号
  data_m[1:N,2] = one_tj[-1]    #轨迹值
  data_m[1:N,3] = point_position         #
  data_m = data_m[seq(1, nrow(data_m), by = 100), ]
  my_tj = rbind(my_tj, data_m)
}
my_tj = my_tj[-1,]


n_test = length(my_data_test[,])

for(i in 1:length(my_data_test[,])){
  one_tj = as.double(unlist(strsplit(my_data_test[i,], split = ',')))
  N = length(one_tj) - 1
  print(c(i,N))
  point_position = (floor(seq(1, N + 1, 1) * (5000/(N + 1))) / 5000)[1:N]
  
  data_m = matrix(data = NA, nrow = N, ncol = 3)
  data_m[1:N,1] = rep(i, N) + n            #标号
  data_m[1:N,2] = one_tj[-1]    #轨迹值
  data_m[1:N,3] = point_position         #
  data_m = data_m[seq(1, nrow(data_m), by = 100), ]
  my_tj_test = rbind(my_tj_test, data_m)
}
my_tj_test = my_tj_test[-1,]



#保号开方化

my_tj_sq = my_tj
my_tj_sq[,2] = sign(my_tj[,2]) * sqrt(abs(my_tj[,2]))

my_tj_test_sq = my_tj_test
my_tj_test_sq[,2] = sign(my_tj_test[,2]) * sqrt(abs(my_tj_test[,2]))

##### 绘出开方化后样本轨迹
par(mfrow=c(3,3))
for (i in 81:89){
  id<-i ##for curve i
  t.c<-my_tj_sq[my_tj_sq[,1]==id,3] ##measurement points
  y.c<-my_tj_sq[my_tj_sq[,1]==id,2] ##obs
  #plots
  plot(t.c,y.c,xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id), type = 'l')
  ##points(t.c,y.pred.proj,col=2, pch=2) ##predicted measurements at observed measurement times
}

timepoints_list = list()
observed_list = list()

my_tj = rbind(my_tj_sq, my_tj_test_sq)

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


#fpcs =  SOAP(observed, timepoints)

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
fpcs = score_pred

#####scalar

ix_s_train = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/ix_data_train.txt',header = F)
kcl_s_train = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/kcl_data_train.txt',header = F)
ao_s_train = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/ao_data_train.txt',header = F)
ae_s_train = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/ae_data_train.txt',header = F)

ix_s_test = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/ix_data_test.txt',header = F)
kcl_s_test = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/kcl_data_test.txt',header = F)
ao_s_test = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/ao_data_test.txt',header = F)
ae_s_test = read.table(file = 'D:/研究生文件/科研任务/functional_data/py_code/ae_data_test.txt',header = F)

my_data2 = rbind(ix_s_train, kcl_s_train, ao_s_train, ae_s_train)
my_data2_test = rbind(ix_s_test, kcl_s_test, ao_s_test, ae_s_test)

library(stringr)

### train dataset
my_scalar = matrix(data = NA, n, ncol = 9)
for(i in 1:n){
  my_scalar[i,] = unlist(strsplit(my_data2[i,], split = ','))
}

my_s = matrix(data = NA, nrow = n, ncol = 3)
missing_c = c()

for (i in 1:n) {
  ht = str_sub(my_scalar[i,7],1,-2)
  ht = gsub("'",".",ht)
  ht = as.double(ht)
  print(ht)
  my_s[i,1] = ht/10
  
  if(my_scalar[i,2] == 'M'){
    my_s[i,2] = 0
  }
  if(my_scalar[i,2] == 'F'){
    my_s[i,2] = 0.1
  }
  
  year1 = str_sub(my_scalar[i,6],-2,-1)
  
  if(year1 == '??'){
    missing_c = c(missing_c,i)
    next
  }
  year1 = as.double(year1)
  year2 = str_sub(my_scalar[i,5],-2,-1)
  year2 = as.double(year2)
  age = year2 - year1
  
  my_s[i,3] = age/100
  
}
missing_c

### test dataset

my_scalar_test = matrix(data = NA, n_test, ncol = 9)
for(i in 1:n_test){
  my_scalar_test[i,] = unlist(strsplit(my_data2_test[i,], split = ','))
}

my_s_test = matrix(data = NA, nrow = n_test, ncol = 3)
missing_c_test = c()

for (i in 1:n_test) {
  ht = str_sub(my_scalar_test[i,7],1,-2)
  ht = gsub("'",".",ht)
  ht = as.double(ht)
  print(ht)
  my_s_test[i,1] = ht/10
  
  if(my_scalar_test[i,2] == 'M'){
    my_s_test[i,2] = 0
  }
  if(my_scalar_test[i,2] == 'F'){
    my_s_test[i,2] = 0.1
  }
  
  year1 = str_sub(my_scalar_test[i,6],-2,-1)
  
  if(year1 == '??'){
    missing_c_test = c(missing_c_test,i)
    next
  }
  year1 = as.double(year1)
  year2 = str_sub(my_scalar_test[i,5],-2,-1)
  year2 = as.double(year2)
  age = year2 - year1
  
  my_s_test[i,3] = age/100
  
}
missing_c_test


n = length(my_data[,])
d_scalar = 3
d_f = length(fpcs[1,])
all_dim = d_f + d_scalar
all_data = matrix(data = NA, nrow = n, ncol = all_dim)
all_data[,1:d_f] = fpcs[1:n,]
all_data[,(d_f+1):all_dim] = my_s

y_label = y_label[-missing_c]
all_data = all_data[-missing_c,]

all_dim = d_f + d_scalar
all_data_test = matrix(data = NA, nrow = n_test, ncol = all_dim)
all_data_test[,1:d_f] = fpcs[(n+1):(n+n_test),]
all_data_test[,(d_f+1):all_dim] = my_s_test

y_label_test = y_label_test[-missing_c_test]
all_data_test = all_data_test[-missing_c_test,]

write.csv(x = all_data, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_data_train_soap_0211.csv')
write.csv(x = y_label, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_y_train_soap_0211.csv')

write.csv(x = all_data_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_data_test_soap_0211.csv')
write.csv(x = y_label_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_y_test_soap_0211.csv')
