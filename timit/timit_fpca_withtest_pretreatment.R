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
  data_m[1:N,1] = rep(i, N)            #标号
  data_m[1:N,2] = one_tj[-1]    #轨迹值
  data_m[1:N,3] = point_position         #
  
  my_tj_test = rbind(my_tj_test, data_m)
}
my_tj_test = my_tj_test[-1,]


##### 绘出样本轨迹
par(mfrow=c(3,3))
for (i in 51:59){
  id<-i ##for curve i
  t.c<-my_tj[my_tj[,1]==id,3] ##measurement points
  y.c<-my_tj[my_tj[,1]==id,2] ##obs
  #plots
  plot(t.c,y.c,xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id), type = 'l')
  ##points(t.c,y.pred.proj,col=2, pch=2) ##predicted measurements at observed measurement times
}

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

### KL expansion 求解轨迹的\xi
library(sm)
library(splines)
library(fpca)



## candidate models for fitting
M.set<-c(40)
r.set<-c(5,6)

##parameters for fpca.mle
ini.method="EM"
basis.method="bs"
sl.v=rep(0.5,10)
max.step=50
grid.l=seq(0,1,0.0002)
grids=seq(0,1,0.0002) #for ini.method = EM

##fit candidate models by fpca.mle
result<-fpca.mle(my_tj_sq, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
summary(result)

eigenfest<-result$eigenfunctions
evalest<-result$eigenvalues
M<-result$selected_model[1]
r<-result$selected_model[2]
grids.new<-result$grid
sig2est<-result$error_var
muest<-result$fitted_mean

par(mfrow=c(1,1))
plot(grids.new,muest)

par(mfrow=c(3,3))
for(i in 1:9){
  plot(grids.new,eigenfest[i,],ylim=range(eigenfest),xlab="time",ylab=paste("eigenfunction",i))
  
}
fpcs<-fpca.score(my_tj_sq,grids.new,muest,evalest,eigenfest,sig2est,r)
pred<-fpca.pred(fpcs, muest,eigenfest)

N<-length(grids.new)

par(mfrow=c(3,3))
for (i in 81:89){
  id<-i ##for curve i
  t.c<-my_tj_sq[my_tj_sq[,1]==id,3] ##measurement points
  t.proj<-ceiling(N*t.c) ##measurement points projected on the grid
  y.c<-my_tj_sq[my_tj_sq[,1]==id,2] ##obs
  y.pred.proj<-pred[t.proj,id] ##predicted obs on the measurement points 投影到新格点上的第几个点
  #plots
  plot(t.c*10,y.c,ylim=range(pred[,id]),xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id), type = 'l')
  points(grids.new*10,pred[,id],col=3,type='l')
  ##points(t.c,y.pred.proj,col=2, pch=2) ##predicted measurements at observed measurement times
}

fpcs_test<-fpca.score(my_tj_test_sq,grids.new,muest,evalest,eigenfest,sig2est,r)





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

d_scalar = 3

all_dim = length(eigenfest[,1]) + d_scalar
all_data = matrix(data = NA, nrow = n, ncol = all_dim)
all_data[,1:length(eigenfest[,1])] = fpcs
all_data[,(length(eigenfest[,1])+1):all_dim] = my_s

y_label = y_label[-missing_c]
all_data = all_data[-missing_c,]

all_dim = length(eigenfest[,1]) + d_scalar
all_data_test = matrix(data = NA, nrow = n_test, ncol = all_dim)
all_data_test[,1:length(eigenfest[,1])] = fpcs_test
all_data_test[,(length(eigenfest[,1])+1):all_dim] = my_s_test

y_label_test = y_label_test[-missing_c_test]
all_data_test = all_data_test[-missing_c_test,]

write.csv(x = all_data, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_data_train_0211.csv')
write.csv(x = y_label, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_y_train_0211.csv')

write.csv(x = all_data_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_data_test_0211.csv')
write.csv(x = y_label_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/timit_y_test_0211.csv')
