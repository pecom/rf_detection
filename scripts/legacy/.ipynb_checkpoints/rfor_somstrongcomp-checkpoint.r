


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
outdir = paste(pdir, '/output/run', args[1], sep='')
indir = paste(pdir, '/data/run', args[1], sep='')
# suffix=paste('_truestrata_all_fixedbig_i.Rda')
suffix=paste('_somcomp_strong_label', args[1], '.Rda', sep='')


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

###############################################
# Load and format data                        #
###############################################
# df_noblend = read.csv(paste(outdir, "/training_pure_messy.csv", sep = ""))
# df_weak = read.csv(paste(outdir, "/training_strong_messy.csv", sep = ""))
# # df_strong = read.csv(paste(outdir, "/strong_table_err.csv", sep = ""))

# df.pure = colorify.full(df_noblend)
# df.pure$type = "pure"
# df.weak = colorify.full(df_weak)
# df.weak$type = "strong"

# df.master = rbind(df.pure, df.weak)

df.master = read.csv(paste(outdir, '/iden_strong_rf.csv', sep=''))
df.master$blend[df.master$type=='weak'] = 0

phot.master = df.master[,1:10]
# phot.master = phot.master[,names(phot.master) %in% c("uB", "iz")]
# phot.master = phot.master[,!names(phot.master) %in% c("FLUX_RADIUS", "i", "iz", "YJ", "zY", "BV", "JH", "ri", "uB")]
label.master = as.factor(df.master$blend)
type.master = as.factor(df.master$type)

# sampfac = as.integer(args[1])
# pure.subset = sample_n(df.master[df.master$type=='pure',], 2*sampfac)
# strong.subset = sample_n(df.master[df.master$type=='strong',], sampfac)
# weak.subset = sample_n(df.master[df.master$type=='pure',], sampfac)

###############################################
# Run random forest!                          #
###############################################
ntree = 500
out.rf = randomForest(x=phot.master,y=label.master,ntree=ntree, norm.votes=FALSE, do.trace=TRUE)
rcont = list("method"=out.rf)

save(rcont, file=paste(outdir,"/rfobj", suffix, sep=''))
save(df.master, file=paste(outdir,"/master_dataframe", suffix, sep=''))