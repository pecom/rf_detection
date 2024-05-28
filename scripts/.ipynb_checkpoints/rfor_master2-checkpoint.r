# library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)

args = commandArgs(trailingOnly=TRUE)

pdir = Sys.getenv("PSCRATCH")   
hdir = Sys.getenv("HOME")
outdir = paste(pdir, '/output/megaruns', args[1], sep='')
indir = paste(pdir, '/output/run', args[1], sep='')
# suffix=paste('_truestrata_all_fixedbig_i.Rda')
suffix=paste('_9feats', args[1], sep='')
remove.cols = c('FLUX_RADIUS')


###############################################
# Functions to create and evaluate random     #
# forest.                                     #
###############################################

# Functions to change magnitudes to colors
colorify.full = function(df){
    uB = df$u - df$B
    BV = df$B - df$V
    Vr = df$V - df$r
    ri = df$r - df$ip
    iz = df$ip - df$zpp
    zY = df$zpp - df$Y
    YJ = df$Y - df$J
    JH = df$J - df$H
    blend = df$blend
    df_new = data.frame(uB, BV, Vr, ri, iz, zY, YJ, JH, df$ip, df$FLUX_RADIUS, df$blend)
    names(df_new)[9:11] = c("i","FLUX_RADIUS", "blend")
    df_new
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
# Load and format data                        #
###############################################
df_noblend = read.csv(paste(indir, "/training_pure_messy.csv", sep = ""))
df_weak = read.csv(paste(indir, "/training_labels.csv", sep = ""))

df.pure = colorify.full(df_noblend)
df.pure$type = "pure"
df.weak = colorify.full(df_weak)
df.weak$type = df_weak$type

df.master = rbind(df.pure, df.weak)
phot.master = df.master[,1:10]
phot.master = phot.master[,!names(phot.master) %in% remove.cols]


weak.master = as.factor(df.master$type)
strong.master = df.master$type
strong.master[strong.master == 'weak'] = 'pure'
strong.master = as.factor(strong.master)


###############################################
# Run random forest!                          #
###############################################
ntree = 500
out.rf = randomForest(x=phot.master,y=weak.master,ntree=ntree, norm.votes=FALSE, do.trace=TRUE)
rcont = list("method"=out.rf)
save(rcont, file=paste(outdir,"/rfobj", suffix, '_w.Rda', sep=''))

out.rf = randomForest(x=phot.master,y=strong.master,ntree=ntree, norm.votes=FALSE, do.trace=TRUE)
rcont = list("method"=out.rf)
save(rcont, file=paste(outdir,"/rfobj", suffix, '_s.Rda', sep=''))