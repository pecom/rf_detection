# library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)

PSCRATCH = Sys.getenv("PSCRATCH")   
hdir = Sys.getenv("HOME")
args = commandArgs(trailingOnly=TRUE)
outdir = paste(PSCRATCH, '/output/panic', sep='')
indir = paste(hdir, '/data/rf_data', sep='')
suffix='_norad.Rda'


###############################################
# Functions to create and evaluate random     #
# forest.                                     #
###############################################

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
# load(paste(hdir, '/data/rf_data/testing.Rda', sep=''))
load(paste(hdir, '/data/rf_data/training.Rda', sep=''))

train.phot = training[1:6]
train.lbl = training[,8]

ntree = 500
out.rf = randomForest(x=train.phot,y=train.lbl,importance=TRUE, ntree=ntree, keep.forest=TRUE)
rcont = list("method"=out.rf)

# weak.predict = new.pred(rcont$method, weak.phot)
# strong.predict = new.pred(rcont$method, strong.phot)

# weak.money = new.moneyplot(weak.predict$predict, weak.test)
# strong.money = new.moneyplot(strong.predict$predict, strong.test)

# save(weak.predict, file=paste(outdir, '/weak_predict', suffix, sep=''))
# save(strong.predict, file=paste(outdir, '/strong_predict', suffix, sep=''))

# save(weak.money, file=paste(outdir, '/weak_mony', suffix, sep=''))
# save(strong.money, file=paste(outdir, '/strong_mony', suffix, sep=''))

save(rcont, file=paste(outdir,"/rfobj", suffix, sep=''))
