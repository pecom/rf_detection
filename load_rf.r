library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)

# Params to load
even.split = FALSE
include.radius = TRUE
messy.data = TRUE
opt.bands = TRUE
blend.weak = FALSE

load.data = FALSE

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
    df_new = data.frame(uB, BV, Vr, ri, iz, df$i, df$FLUX_RADIUS, df$blend)
    names(df_new)[6:8] = c("i", "FLUX_RADIUS", "blend")
    df_new
}

colorify.optical.norad = function(df){
    uB = df$u - df$B
    BV = df$B - df$V
    Vr = df$V - df$r
    ri = df$r - df$ip
    iz = df$i - df$zpp
    blend = df$blend
    df_new = data.frame(uB, BV, Vr, ri, iz, df$i, df$blend)
    names(df_new)[6:7] = c("i", "blend")
    df_new
}

if (load.data){
    ###############################################
    # Load Data                                   #
    ###############################################
    if (messy.data) {
        df_blend = read.csv('./data/rf_data/blend_phot_messy.csv')
        df_noblend = read.csv('./data/rf_data/noblend_phot_messy.csv')
        if (blend.weak)  {
            df_vali = read.csv('./data/rf_data/vali_phot_messy.csv')
        } else {
            df_vali = read.csv('./data/rf_data/vali_strong_messy.csv')
        }
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

    # if (include.radius) {
    #     df.color = add_column(df.color, FLUX_RADIUS = df$FLUX_RADIUS, .after = "i")
    #     vali.color = add_column(vali.color, FLUX_RADIUS=df_vali$FLUX_RADIUS, .after="i")
    # }

    df.nparams = ncol(df.color) - 1

    phot.train = df.color[, 1:df.nparams]
    phot.test = vali.color[, 1:df.nparams]

    blend.train = df.color[, df.nparams+1]
    blend.test = vali.color[, df.nparams+1]
    train_factors = factor(blend.train)
}
###############################################
# Functions to create decision trees,         #
# random forests, and final plots.            #
###############################################
vali.load = function(vali) {
    df_vali = read.csv(vali)
    vali.color = colorify.optical(df_vali)
    vali.color = add_column(vali.color, FLUX_RADIUS=df_vali$FLUX_RADIUS, .after="i")
    df.nparams = ncol(vali.color) - 1
    vali = list("phot.test" = vali.color[, 1:df.nparams], "blend.test" = vali.color[,df.nparams+1])
}
    
conf.factors = function(bpred, btrue) {
    # confusion_matrix(obs=btest, pred=bpred, plot=TRUE,
    #              colors = c(low="#eff3ff" , high="#08519c"), unit = "count")
    CrossTable(bpred,btrue)
}

new.pred = function(rf.obj, p.test) {
    blend.pred = predict(rf.obj,newdata=p.test)
    attr(blend.pred, "names") <- NULL
    newList <- list("predict" = blend.pred)
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

fname = sprintf("./output/rf_split%s_radius%s_messy%s_opt%s_weak%s.Rda", even.split, include.radius, messy.data, opt.bands, blend.weak)
rfname = sprintf("./output/rfobj_split%s_radius%s_messy%s_opt%s.Rda", even.split, include.radius, messy.data, opt.bands)
