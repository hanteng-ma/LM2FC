library(loon.data)

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



## candidate models for fitting
M.set<-c(15)
r.set<-c(7,8)

##parameters for fpca.mle
ini.method="EM"
basis.method="bs"
sl.v=rep(0.5,10)
max.step=50
grid.l=seq(0,1,0.01)
grids=seq(0,1,0.001) #for ini.method = EM

##fit candidate models by fpca.mle
result<-fpca.mle(my_tj, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
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

par(mfrow=c(2,2))
for(i in 1:4){
  plot(grids.new,eigenfest[i,],ylim=range(eigenfest),xlab="time",ylab=paste("eigenfunction",i))
  
}
fpcs<-fpca.score(my_tj,grids.new,muest,evalest,eigenfest,sig2est,r)
pred<-fpca.pred(fpcs, muest,eigenfest)

N<-length(grids.new)

par(mfrow=c(3,3))
for (i in 1:13){
  if((i %% 3 == 0) | (i %% 3 == 1) ){
    id<-i ##for curve i
    t.c<-my_tj[my_tj[,1]==id,3] ##measurement points
    t.proj<-ceiling(N*t.c) ##measurement points projected on the grid
    y.c<-my_tj[my_tj[,1]==id,2] ##obs
    y.pred.proj<-pred[t.proj,id] ##predicted obs on the measurement points 投影到新格点上的第几个点
    #plots
    plot(t.c,y.c,xlim = range(0,1), ylim=range(pred[,id]),xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id), type = 'p')
    points(grids.new,pred[,id],col=3,type='l')
    ##points(t.c,y.pred.proj,col=2, pch=2) ##predicted measurements at observed measurement times
  }
}

fpcs_test<-fpca.score(my_tj_test,grids.new,muest,evalest,eigenfest,sig2est,r)  # test sample

n = length(y_label)
n_test = length(y_label_test)
d_scalar = 1

all_dim = length(eigenfest[,1]) + d_scalar
all_data = matrix(data = NA, nrow = n, ncol = all_dim)
all_data[,1:length(eigenfest[,1])] = fpcs
all_data[,(length(eigenfest[,1])+1):all_dim] = my_scalar



all_dim = length(eigenfest[,1]) + d_scalar
all_data_test = matrix(data = NA, nrow = n_test, ncol = all_dim)
all_data_test[,1:length(eigenfest[,1])] = fpcs_test
all_data_test[,(length(eigenfest[,1])+1):all_dim] = my_scalar_test

write.csv(x = all_data, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_data_train_0610.csv')
write.csv(x = y_label, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_y_train_0610.csv')

write.csv(x = all_data_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_data_test_0610.csv')
write.csv(x = y_label_test, file = 'D:/研究生文件/科研任务/functional_data/FSVM_code_py/bone_y_test_0610.csv')

