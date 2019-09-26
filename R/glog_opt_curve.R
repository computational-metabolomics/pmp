attach (testData)
data <- mv_imputation(df=testData$data, method='knn')

out <- glog_transformation (df=data, classes=testData$class, qc_label='QC',
    store_lambda=TRUE)

data_qc <- data[,testData$class=="QC"]

class <- sbcms::sbcdata$class

s <- data.frame(lambda=seq(0,1e10,length.out=100))

k=1
for (lambda in s[,1]) {
    s[k,2]=pmp:::SSE(lambda=lambda,y0=0,y=t(QC$data))
    k=k+1
}

#colnames(s)=c('lambda','SSE')

#G=glog_transform(qc_label='QC',factor_name='class')
#G=method.apply(G,D)
#G$lambda_opt
#G$error_flag
#y=pmp:::SSE(G$lambda_opt,y0 = 0,y=t(QC$data))

#g=ggplot(s,aes(x=lambda,y=SSE))+
#    geom_line()+
#    geom_point(data = data.frame(lambda=G$lambda_opt,SSE=y),mapping=aes(x=lambda,y=SSE),colour='red')
#plot(g)
