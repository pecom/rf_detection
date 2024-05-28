# library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)
library(dplyr)
library(pracma)

args = commandArgs(trailingOnly=TRUE)

pdir = Sys.getenv("PSCRATCH")   
hdir = Sys.getenv("HOME")
outdir = paste(pdir, '/output/megaruns/importance', sep='')
indir = paste(pdir, '/output/run1', sep='')
# suffix=paste('_truestrata_all_fixedbig_i.Rda')
remove.cols = c('JH', 'zY', 'YJ', 'uB', 'Vr', 'iz', 'BV', 'ri', 'i')
suffix=paste(paste(remove.cols, collapse=''), '_', sep='')
print(suffix)

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


denom.lin.moneyplot = function(df, blend.filt, sname='score'){
    cutoffs = seq(0.01, .99, .01)
    blendpercs = list()
    samplepercs = list()
    precpercs = list()
    remainblend = list()
    cutvals = list()
    
    blend.denom = sum(blend.filt)
    sample.denom = nrow(df)
    
    for (i in seq_along(cutoffs)) {
        thrown.filt = df[[sname]] > cutoffs[i]
        pred.labels = as.integer(thrown.filt)
        blend.dat = sum(pred.labels[blend.filt])/blend.denom
        sample.dat = sum(pred.labels)/sample.denom
        prec.dat = sum(pred.labels[blend.filt])/sum(pred.labels)
        
        not.thrown = df[!thrown.filt,]
        remainblend[[i]] = sum(not.thrown$blend)/sum(!thrown.filt)


        blendpercs[[i]] <- blend.dat
        samplepercs[[i]] <- sample.dat
        precpercs[[i]] <- prec.dat
        cutvals[[i]] <- i
    }
    money.plot = do.call(rbind, Map(data.frame, recall = blendpercs, cost = samplepercs,
                                    precision=precpercs, cut = cutvals, remain=remainblend))
    money.plot
}


###############################################
# Load and format data                        #
###############################################
df_noblend = read.csv(paste(indir, "/training_pure_messy.csv", sep = ""))
df_noblend = filter(df_noblend, u<30, B<30, V<30, r<30, ip<30, zpp<30, Y<30, J<30,H<30)

df_weak = read.csv(paste(indir, "/training_labels.csv", sep = ""))
df_weak = filter(df_weak, u<30, B<30, V<30, r<30, ip<30, zpp<30, Y<30, J<30,H<30)

df.pure = colorify.full(df_noblend)
df.pure$type = "pure"
df.weak = colorify.full(df_weak)
df.weak$type = df_weak$type

df.master = rbind(df.pure, df.weak)
phot.master = df.master[,1:10]
# phot.master = phot.master[,!names(phot.master) %in% remove.cols]
phot.master = data.frame(FLUX_RADIUS=phot.master$FLUX_RADIUS)
print(str(phot.master))

weak.master = df.master$type
weak.master = as.factor(weak.master)

weak_vali = read.csv(paste(pdir, '/output/run1/vali_weak_som.csv', sep=''))
weak_vali = filter(weak_vali, u<30, B<30, V<30, r<30, ip<30, zpp<30, Y<30, J<30,H<30)
weak.vali = colorify.full(weak_vali)
vali.phot = weak.vali[,1:10]
# vali.phot = vali.phot[,!names(vali.phot) %in% remove.cols]
vali.phot = data.frame(FLUX_RADIUS=vali.phot$FLUX_RADIUS)
print(str(vali.phot))

###############################################
# Run random forest!                          #
###############################################
ntree = 100
out.rf = randomForest(x=phot.master,y=weak.master,ntree=ntree, norm.votes=FALSE, do.trace=TRUE)
rcont = list("method"=out.rf)
save(rcont, file=paste(outdir,"/rfobj", suffix,  '_w.Rda', sep=''))

weak.score = data.frame(predict(out.rf, newdata=vali.phot, type='prob'))
weak.vali$score = 1-weak.score$pure
weak.filt = (weak.vali$blend==1)
weak.plot = denom.lin.moneyplot(weak.vali, weak.filt, 'score')
weak.spline = smooth.spline(weak.plot$cost, weak.plot$recall)
smooth = function(x) predict(weak.spline, x)$y
total = integrate(smooth, 0, 1)$value

finstr = paste("Removing ", suffix, " the area is: ", total, sep='')
print(finstr)
