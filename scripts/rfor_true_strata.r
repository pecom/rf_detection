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
outdir = paste(pdir, '/output/subset', sep='')
indir = paste(pdir, '/data/subset', sep='')
# suffix=paste('_truestrata_all_', args[1], '.Rda')
suffix=paste('_truestrata_nofrad_i_iz_yj_zy_bv_jh_ri_ub.Rda')


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
# df_noblend = read.csv(paste(pdir, "/data/new_run/pure_table_err.csv", sep = ""))
# df_weak = read.csv(paste(pdir, "/data/new_run/weak_table_err.csv", sep = ""))
# df_strong = read.csv(paste(pdir, "/data/new_run/strong_table_err.csv", sep = ""))
# 
# df.pure = colorify.full(df_noblend)
# df.pure$type = "pure"
# df.weak = colorify.full(df_weak)
# df.weak$type = "weak"
# df.strong = colorify.full(df_strong)
# df.strong$type = "strong"
# 
# df.strong$blend = 1
# df.master = rbind(df.pure, df.weak, df.strong)

df.master = read.csv(paste(pdir, '/data/subset/subset_table.csv', sep=''))

phot.master = df.master[,2:11]
phot.master = phot.master[,!names(phot.master) %in% c("FLUX_RADIUS", "i", "iz", "YJ", "zY", "BV", "JH", "ri", "uB")]
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
out.rf = randomForest(x=phot.master,y=type.master,ntree=ntree, norm.votes=FALSE, do.trace=TRUE)
rcont = list("method"=out.rf)

save(rcont, file=paste(outdir,"/rfobj", suffix, sep=''))
save(df.master, file=paste(outdir,"/master_dataframe", suffix, sep=''))
