# library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

pdir = Sys.getenv("PSCRATCH")   
hdir = Sys.getenv("HOME")
outdir = paste(pdir, '/output/multirun/', sep='')
indir = paste(pdir, '/data/subset', sep='')
suffix=paste('_1000pt_10run.Rda')

nruns = 100
nblends = 1000
npure = nblends
# suffix=paste('_truestrata_nofrad_i_iz_yj_zy_bv_jh_ri_ub.Rda')


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

denom.lin.moneyplot = function(df, blend.filt, sname='score'){
    cutoffs = seq(0.01, .99, .01)
    blendpercs = list()
    samplepercs = list()
    cutvals = list()

    blend.denom = sum(blend.filt)
    sample.denom = nrow(df)

    for (i in seq_along(cutoffs)) {
        pred.labels = as.integer(df[[sname]] > cutoffs[i])
        blend.dat = sum(pred.labels[blend.filt])/blend.denom
        sample.dat = sum(pred.labels)/sample.denom

        blendpercs[[i]] <- blend.dat
        samplepercs[[i]] <- sample.dat
        cutvals[[i]] <- i
    }
    money.plot = do.call(rbind, Map(data.frame, blend = blendpercs,
        sample = samplepercs, cut = cutvals))
    money.plot
}

###############################################
# Load and format data                        #
###############################################

df_noblend = read.csv(paste(pdir, "/data/new_run/pure_table_err.csv", sep = ""))
df_weak = read.csv(paste(pdir, "/data/new_run/weak_table_err.csv", sep = ""))
df_strong = read.csv(paste(pdir, "/data/new_run/strong_table_err.csv", sep = ""))

df.pure = colorify.full(df_noblend)
df.pure$type = "pure"
df.weak = colorify.full(df_weak)
df.weak$type = "weak"
df.strong = colorify.full(df_strong)
df.strong$type = "strong"

df.strong$blend = 1
df.full = rbind(df.pure, df.weak, df.strong)

# df.master = read.csv(paste(pdir, '/data/subset/subset_table.csv', sep=''))
# phot.master = phot.master[,!names(phot.master) %in% c("FLUX_RADIUS", "i", "iz", "YJ", "zY", "BV", "JH", "ri", "uB")]

###############################################
# Run random forest!                          #
###############################################
spline.runs = list()
for (i in 1:nruns){
	print("Resampling dataset...")
	pure.ndx = sample(nrow(df.pure), npure, replace=FALSE)
	weak.ndx = sample(nrow(df.weak), nblends%/%2, replace=FALSE)
	strong.ndx = sample(nrow(df.strong), nblends%/%2, replace=FALSE)

	df.train = rbind(df.pure[pure.ndx,], df.weak[weak.ndx,], df.strong[strong.ndx,])
	df.test = rbind(df.pure[-pure.ndx,], df.weak[-weak.ndx,], df.strong[-strong.ndx,])

	test.filt = df.test$type!='pure'

	phot.master = df.train[,1:10]
	type.master = as.factor(df.train$type)

	ntree = 500
	out.rf = randomForest(x=phot.master,y=type.master,ntree=ntree, norm.votes=FALSE, do.trace=FALSE)
	print("Created forest... creating power plot")

	test.scores = data.frame(predict(out.rf, df.test[,1:10], type='prob'))
	test.scores$notpure = 1-test.scores$pure

	single.run = denom.lin.moneyplot(test.scores, test.filt, sname="notpure")

	print("Created plot... smoothing to spline")
	smoothingSpline = smooth.spline(single.run$sample, single.run$blend)
	base.points = predict(smoothingSpline, seq(0.01, 0.99, 0.01))

	base.spline = data.frame(blend = base.points$y,
				 sample = base.points$x,
				 type = i)
	nsplines = length(spline.runs)
	spline.runs[[nsplines+1]] = base.spline
}
mega.spline = do.call("rbind", spline.runs)

save(mega.spline, file=paste(outdir,"mega_spline_", suffix, sep=''))

# save(rcont, file=paste(outdir,"/rfobj", suffix, sep=''))
# save(df.master, file=paste(outdir,"/master_dataframe", suffix, sep=''))
