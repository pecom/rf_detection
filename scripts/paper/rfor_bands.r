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
base.cols = c()
og.bands = c('u', 'B', 'V', 'r', 'i', 'z', 'Y', 'J', 'H')

# bands = c('u')
# bands = c('FLUX_RADIUS')
# remove.cols = c('JH', 'zY', 'YJ', 'uB', 'Vr', 'BV', 'iz', 'FLUX_RADIUS')


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



colorify.funky = function(band.run, col.seq){
	if (length(col.seq)==1){
		phot.master = data.frame(i=band.run[,col.seq[[1]]])
		names(phot.master) = col.seq[[1]]
		print(str(phot.master))
		return(phot.master)
	}

	iter.n = ncol(band.run) - 1
	band.run = band.run[,col.seq]
	band.names = names(band.run)

	temp.data = list()
	temp.names = list()

	for (i in 1:iter.n){
		ndx.1 = i
		ndx.2 = i+1
		col.name = paste(band.names[ndx.1:ndx.2], collapse='')
		temp.names[[i]] = col.name
		temp.data[[i]] = band.run[,ndx.1] - band.run[,ndx.2]
	}
	band.master = data.frame(t(do.call(rbind, temp.data)))
	names(band.master) = temp.names
	band.master
}




###############################################
# Load and format data                        #
###############################################
df_noblend = read.csv(paste(indir, "/training_pure_messy.csv", sep = ""))
df_noblend = filter(df_noblend, u<30, B<30, V<30, r<30, ip<30, zpp<30, Y<30, J<30,H<30)

df_weak = read.csv(paste(indir, "/training_labels.csv", sep = ""))
df_weak = filter(df_weak, u<30, B<30, V<30, r<30, ip<30, zpp<30, Y<30, J<30,H<30)

# df.pure = colorify.full(df_noblend)
df.pure = rename(df_noblend, i=ip, z=zpp)
df.pure$type = "pure"

#df.weak = colorify.full(df_weak)
df.weak = rename(df_weak, i=ip, z=zpp)
df.weak$type = df_weak$type

df.master = rbind(df.pure, df.weak)
phot.master = df.master[,1:10]

weak_vali = read.csv(paste(pdir, '/output/run1/vali_weak_som.csv', sep=''))
weak_vali = filter(weak_vali, u<30, B<30, V<30, r<30, ip<30, zpp<30, Y<30, J<30,H<30)
weak.vali = rename(weak_vali, i=ip, z=zpp)
# weak.vali = colorify.full(weak_vali)
vali.phot = weak.vali[,1:10]
weak.master = df.master$type
weak.master = as.factor(weak.master)

master.recall = list()
master.integral = list()

all.bands = c('u', 'B', 'V', 'r', 'i', 'z', 'Y', 'J', 'H')
bands = c('u', 'B', 'V', 'r', 'z', 'Y', 'J', 'H')
base.cols = c()
bands.importance = list()
for (run.num in 1:7){
	recalls = list()
	integrals = list()
	for (b in bands){
		remove.cols = c(base.cols, b)
		# remove.cols = c('u')
		suffix=paste('bands_', paste(remove.cols, collapse=''), '_', sep='')
		print(paste("Removing bands", paste(remove.cols, collapse='-'), sep=' '))

		phot.run = phot.master[,!names(phot.master) %in% remove.cols]
		band.run = select(phot.run, -FLUX_RADIUS)

		col.seq = og.bands[-match(remove.cols, og.bands)]
		band.master = colorify.funky(band.run, col.seq)

		# band.master[,b] = phot.run[,b]
		if (b!='FLUX_RADIUS'){
			band.master$FLUX_RADIUS = phot.run$FLUX_RADIUS
		}
		band.master$i = phot.run$i
		# phot.master = data.frame(i=phot.master$i)
		print(names(band.master))

		vali.run = vali.phot[,!names(vali.phot) %in% remove.cols]
		vali.run = select(vali.run, -FLUX_RADIUS)

		vali.band = colorify.funky(vali.run, col.seq)
		# vali.band[,b] = vali.phot[,b]
		if (b!='FLUX_RADIUS'){
			vali.band$FLUX_RADIUS = vali.phot$FLUX_RADIUS
		}
		vali.band$i = vali.phot$i
	# 	print(names(vali.band))
		# vali.phot = data.frame(i=vali.phot$i)

		###############################################
		# Run random forest!                          #
		###############################################
		ntree = 100
		file.str = paste(outdir, "/rfobj", suffix, "_w.Rda", sep='')

		load.file=FALSE

		recall.scores = rep(0,5)
		integral.scores = rep(0,5)
		for (i in 1:5){
			if(file.exists(file.str) & load.file){
			#	print("Loading RF...")
				load(file.str)
				out.rf = rcont$method
			} else {
	#			print("Creating RF...")
				out.rf = randomForest(x=band.master,y=weak.master,ntree=ntree, norm.votes=FALSE, do.trace=FALSE)
				rcont = list("method"=out.rf)
			# 	save(rcont, file=paste(outdir,"/rfobj", suffix,  '_w.Rda', sep=''))
			}

			weak.score = data.frame(predict(out.rf, newdata=vali.band, type='prob'))
			weak.vali$score = 1-weak.score$pure
			weak.filt = (weak.vali$blend==1)
			weak.plot = denom.lin.moneyplot(weak.vali, weak.filt, 'score')
			weak.spline = smooth.spline(weak.plot$cost, weak.plot$recall)
			smooth = function(x) predict(weak.spline, x)$y
			total.static = smooth(.1)
			total.int = integrate(smooth, 0, 1)$value

			recall.scores[[i]] = total.static
			integral.scores[[i]] = total.int

			finstr = paste("Static recall ", suffix, " the recall is: ", total.static, sep='')
			# print(finstr)
			finstr = paste("Integrating recall ", suffix, " the recall is: ", total.int, sep='')
			# print(finstr)

			finstr = paste("(", round(total.static*100,2), ", ", round(total.int*100,2), ")", sep='')
			# print(finstr)
		}
		print(paste("Average recall:", mean(recall.scores)))
		print(paste("Average integral:", mean(integral.scores)))
		print("----------------------------------------")

		recalls[b] = mean(recall.scores)
		integrals[b] = mean(integral.scores)
	}

	master.recall[[run.num]] = recalls
	master.integral[[run.num]] = integrals

	max.recall = names(recalls)[which.max(recalls)]
	bands = bands[-match(max.recall, bands)]
	base.cols = c(base.cols, max.recall)
	bands.importance[max.recall] = recalls[which.max(recalls)]
}


print("Master recall")
print(master.recall)
print("----------------------------------------")
print("Bands importance ")
print(bands.importance)
print("----------------------------------------")
print("Base Columns")
print(base.cols)
print("----------------------------------------")





