# library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)

PSCRATCH = Sys.getenv("PSCRATCH")   
args = commandArgs(trailingOnly=TRUE)
outdir = paste(PSCRATCH, '/output/run', args[1], sep='')

colorify.optical = function(df){
    uB = df$u - df$B
    BV = df$B - df$V
    Vr = df$V - df$r
    ri = df$r - df$ip
    iz = df$i - df$zpp
    blend = df$blend

    df_new = data.frame(uB, BV, Vr, ri, iz, df$i, df$blend)
    names(df_new)[6:7] = c("i", "blend")

#    df_new = data.frame(uB, BV, Vr, ri, iz, df$i, df$FLUX_RADIUS, df$blend)
#    names(df_new)[6:8] = c("i", "FLUX_RADIUS", "blend")
    df_new
}
###############################################
# Load Data                                   #
###############################################

df_blend = read.csv(paste(outdir, '/training_blend_messy.csv', sep=''))
df_noblend = read.csv(paste(outdir, '/training_pure_messy.csv', sep=''))
df_weak = read.csv(paste(outdir, '/vali_weak_messy.csv', sep=''))
df_strong = read.csv(paste(outdir, '/vali_strong_messy.csv', sep=''))

df = rbind(df_blend, df_noblend)

###############################################
# Massage into final form.                    #
# Change into colors and include flux_radius  #
###############################################

df.color = colorify.optical(df)
strong.color = colorify.optical(df_strong)
weak.color = colorify.optical(df_weak)

df.nparams = ncol(df.color) - 1

phot.train = df.color[, 1:df.nparams]
weak.phot = weak.color[, 1:df.nparams]
strong.phot = strong.color[, 1:df.nparams]

blend.train = df.color[, df.nparams+1]
weak.test = weak.color[, df.nparams + 1]
strong.test = strong.color[, df.nparams + 1]

###############################################
# Functions to create and evaluate random     #
# forest.                                     #
###############################################
rf.cont = function(nt=500) {
    out.rf = randomForest(x=phot.train,y=blend.train,importance=TRUE, ntree=nt)
    newList <- list("method" = out.rf)
}

score.throw = function(pred, truth) {
    if (sum(pred) == 0) {
        print("Only predicting 0s... there's an issue")
        newList <- list("blend"=0.0, "total"=0.0)
    } else {
        # t2 = confusion_matrix(obs=truth, pred=pred)
        # newList <- list("blend" = t2[2,2]/(sum(t2[,2])), "total" = sum(t2[2,])/sum(t2))
	t2 = cvms::confusion_matrix(truth, pred)
	newList = list("blend" = t2$Sensitivity, "total" = t2[["Detection Prevalence"]])
    }
}

new.pred = function(rf.obj, p.test) {
    blend.pred = predict(rf.obj,newdata=p.test)
    attr(blend.pred, "names") <- NULL
    newList <- list("predict" = blend.pred)
}

new.moneyplot = function(pred, truth) {
	cutoffs = seq(0.01, .99, .01)
	blendpercs = list()
	samplepercs = list()
	cutvals = list()

	for (i in seq_along(cutoffs)) {
	    pred.labels = as.integer(pred > cutoffs[i])
	    dat = score.throw(pred.labels, truth)
	    # dat$i <- i  # maybe you want to keep track of which iteration produced it?
	    # datalist[[i]] <- dat # add it to your list
	    blendpercs[[i]] <- dat$blend
	    samplepercs[[i]] <- dat$total
	    cutvals[[i]] <- i
	}
	money.plot = do.call(rbind, Map(data.frame, blend=blendpercs, sample=samplepercs, cut=cutvals))
	money.plot
}

###############################################
# Run random forest!                          #
###############################################
ntree = 500
rcont = rf.cont(ntree)

weak.predict = new.pred(rcont$method, weak.phot)
strong.predict = new.pred(rcont$method, strong.phot)

weak.money = new.moneyplot(weak.predict$predict, weak.test)
strong.money = new.moneyplot(strong.predict$predict, strong.test)

save(weak.predict, file=paste(outdir, '/weak_predict_norad.Rda', sep=''))
save(strong.predict, file=paste(outdir, '/strong_predict_norad.Rda', sep=''))

save(weak.money, file=paste(outdir, '/weak_mony_norad.Rda', sep=''))
save(strong.money, file=paste(outdir, '/strong_mony_norad.Rda', sep=''))

save(rcont, file=paste(outdir,"/rfobj_norad.Rda", sep=''))
