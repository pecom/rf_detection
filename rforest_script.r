library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)

# Params for the run
even.split = FALSE
include.radius = TRUE
messy.data = FALSE
opt.bands = TRUE

# Functions to change magnitudes to colors
colorify.full = function(df){
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

colorify.optical = function(df){
    uB = df$u - df$B
    BV = df$B - df$V
    Vr = df$V - df$r
    ri = df$r - df$ip
    iz = df$i - df$zpp
    blend = df$blend
    df_new = data.frame(uB, BV, Vr, ri, iz, df$i,df$blend)
    names(df_new)[6:7] = c("i","blend")
    df_new
}
###############################################
# Load Data                                   #
###############################################
if (messy.data) {
    df_blend = read.csv('./data/rf_data/blend_phot_messy.csv')
    df_noblend = read.csv('./data/rf_data/noblend_phot_messy.csv')
    df_vali = read.csv('./data/rf_data/vali_phot_messy.csv')
} else {
    df_blend = read.csv('./data/rf_data/blend_phot.csv')
    df_noblend = read.csv('./data/rf_data/noblend_phot.csv')
    df_vali = read.csv('./data/rf_data/vali_phot.csv')
}

if (even.split) {
    df_noblend = dplyr::sample_n(df_noblend, nrow(df_blend))
}

df = rbind(df_blend, df_noblend)

###############################################
# Massage into final form.                    #
# Change into colors and include flux_radius  #
###############################################
if (opt.bands) {
    df.color = colorify.optical(df)
    vali.color = colorify.optical(df_vali)
} else {
    df.color = colorify.full(df)
    vali.color = colorify.full(df_vali)
}

if (include.radius) {
    df.color = add_column(df.color, FLUX_RADIUS = df$FLUX_RADIUS, .after = "i")
    vali.color = add_column(vali.color, FLUX_RADIUS=df_vali$FLUX_RADIUS, .after="i")
}

df.nparams = ncol(df.color) - 1

phot.train = df.color[, 1:df.nparams]
phot.test = vali.color[, 1:df.nparams]

blend.train = df.color[, df.nparams+1]
blend.test = vali.color[, df.nparams+1]
train_factors = factor(blend.train)

###############################################
# Functions to create decision trees,         #
# random forests, and final plots.            #
###############################################
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
        newList <- list("blend"=0.0, "total"=0.0)
    } else {
        t2 = confusion_matrix(obs=truth, pred=pred)
        newList <- list("blend" = t2[2,2]/(sum(t2[,2])), "total" = sum(t2[2,])/sum(t2))
    }
}

###############################################
# Run random forest!                          #
###############################################
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
cutvals = list()

for (i in seq_along(cutoffs)) {
    pred.labels = as.integer(blend.predict > cutoffs[i])
    dat = score.throw(pred.labels, blend.test)
    # dat$i <- i  # maybe you want to keep track of which iteration produced it?
    # datalist[[i]] <- dat # add it to your list
    blendpercs[[i]] <- dat$blend
    samplepercs[[i]] <- dat$total
    cutvals[[i]] <- i
}

money.plot = do.call(rbind, Map(data.frame, blend=blendpercs, sample=samplepercs, cut=cutvals))
fname = sprintf("./output/rf_split%s_radius%s_messy%s_opt%s.Rda", even.split, include.radius, messy.data, opt.bands)
save(money.plot, file = fname)

tplot = ggplot(data=money.plot, aes(x=sample, y=blend)) + geom_point(color='red', shape=1, size=3)
figname = sprintf("./figs/rf_split%s_radius%s_messy%s_opt%s.png", even.split, include.radius, messy.data, opt.bands)
ggsave(figname, plot=tplot)