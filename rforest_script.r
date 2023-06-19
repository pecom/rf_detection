library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)


colorify = function(df){
    uB = df$u - df$B
    BV = df$B - df$V
    Vr = df$V - df$r
    ri = df$r - df$ip
    iz = df$i - df$zpp
    zY = df$zpp - df$Y
    YJ = df$Y - df$J
    JH = df$J - df$H
    blend = df$blend
    df_new = data.frame(uB, BV, Vr, ri, iz, zY, YJ, JH, df$i,df$blend)
    names(df_new)[9:10] = c("i","blend")
    df_new
}

full_data = FALSE

df_blend = read.csv('./data/rf_data/blend_phot.csv')
if (full_data) {
    df_noblend = read.csv('./data/rf_data/noblend_phot.csv')
} else {
    df_noblend = read.csv('./data/rf_data/noblend_sub_phot.csv')
}
df = rbind(df_blend, df_noblend)

df_vali = read.csv('./data/rf_data/vali_phot.csv')

df.color = colorify(df)
vali.color = colorify(df_vali)

phot.train = df.color[, 1:9]
phot.test = vali.color[, 1:9]

blend.train = df.color[, 10]
blend.test = vali.color[, 10]

train_factors = factor(blend.train)


conf.factors = function(bpred, btrue) {
    # confusion_matrix(obs=btest, pred=bpred, plot=TRUE,
    #              colors = c(low="#eff3ff" , high="#08519c"), unit = "count")
    CrossTable(bpred,btrue)
}
dtree.classify = function(p.test, makeplot=TRUE) {
    out.rpart = rpart(blend.train~.,data=phot.train,
                  minsplit = 1, minbucket = 1, method='class')
    blend.pred = predict(out.rpart,newdata=p.test, type = "class")
    if (makeplot){
        rpart.plot(out.rpart)
    }
    newList <- list("method" = out.rpart, "predict" = blend.pred)
}
dtree.cont = function(p.test, makeplot=TRUE) {
    out.rpart = rpart(blend.train~.,data=phot.train, minsplit = 1, minbucket = 1, method='anova')
    if (makeplot){
        rpart.plot(out.rpart)
    }
    blend.pred = predict(out.rpart,newdata=p.test)
    attr(blend.pred, "names") <- NULL
    newList <- list("method" = out.rpart, "predict" = blend.pred)
}
rf.cont = function(p.test, nt=500) {
    out.rf = randomForest(x=phot.train,y=blend.train,importance=FALSE, ntree=nt)
    blend.pred = predict(out.rf,newdata=p.test)
    attr(blend.pred, "names") <- NULL
    newList <- list("method" = out.rf, "predict" = blend.pred)
}
rf.class = function(p.test, nt=500) {
    out.rf = randomForest(x=phot.train,y=train_factors,importance=TRUE, ntree=10)
    blend.pred = predict(out.rf,newdata=p.test)
    attr(blend.pred, "names") <- NULL
    newList <- list("method" = out.rf, "predict" = blend.pred)
}
score.throw = function(pred, truth) {
    if (sum(pred) == 0) {
        print("Only predicting 0s... there's an issue")
    } else {
        t2 = confusion_matrix(obs=truth, pred=pred)
        newList <- list("blend" = t2[2,2]/(sum(t2[,2])), "total" = sum(t2[2,])/sum(t2))
    }
}

cutoff = .6
ntree = 500

rcont = rf.cont(phot.test, ntree)
blend.predict = rcont$predict

pred.labels = as.integer(blend.predict > .9)
kale = score.throw(pred.labels, blend.test)
cat("Throwing away", kale$total, "% of sample gets", kale$blend, "% of the blends")


cutoffs = seq(0.01, .99, .01)
# datalist = vector("list", length = length(cutoffs))
blendpercs = list()
samplepercs = list()

for (i in seq_along(cutoffs)) {
    pred.labels = as.integer(blend.predict > cutoffs[i])
    dat = score.throw(pred.labels, blend.test)
    # dat$i <- i  # maybe you want to keep track of which iteration produced it?
    # datalist[[i]] <- dat # add it to your list
    blendpercs[[i]] <- dat$blend
    samplepercs[[i]] <- dat$total
}

money.plot = do.call(rbind, Map(data.frame, blend=blendpercs, sample=samplepercs))
if (full_data) {
    save(money.plot, file = "./output/full_split_rf.Rda") 
} else {
    save(money.plot, file = "./output/even_split_rf.Rda") 
}