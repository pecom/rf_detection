# library(metrica)
library(rpart)
library(rpart.plot)
library(gmodels)
library(randomForest)
library(ggplot2)
library(tibble)
library(cvms)
library(dplyr)
library(kernlab)

args = commandArgs(trailingOnly=TRUE)

pdir = Sys.getenv("PSCRATCH")   
hdir = Sys.getenv("HOME")
outdir = paste(pdir, '/output/megaruns', sep='')
indir = paste(pdir, '/output/run', args[1], sep='')
suffix=paste('_rbf', args[1], '.Rda', sep='')


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
df_noblend = read.csv(paste(indir, "/training_pure_messy.csv", sep = ""))
df_weak = read.csv(paste(indir, "/training_labels.csv", sep = ""))

df.pure = colorify.full(df_noblend)
df.pure$type = "pure"
df.weak = colorify.full(df_weak)
df.weak$type = df_weak$type

df.master = rbind(df.pure, df.weak)
phot.master = df.master[,1:10]

label.master = as.factor(df.master$blend)
type.master = as.factor(df.master$type)

###############################################
# Run random forest!                          #
###############################################
gp = gausspr(x=phot.master, y=type.master, kernel='rbfdot')
gcont = list("method"=gp)

save(gcont, file=paste(outdir,"/gp", suffix, sep=''))
