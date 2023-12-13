hdir = Sys.getenv("HOME")
pdir = Sys.getenv("PSCRATCH")

colorify.full = function(df){
    uB = df$u - df$B
    BV = df$B - df$V
    Vr = df$V - df$r
    ri = df$r - df$ip
    iz = df$ip - df$zpp
    zY = df$zpp - df$Y
    YJ = df$Y - df$J
    JH = df$J - df$H
    i = df$ip
    FLUX_RADIUS = df$FLUX_RADIUS
    blend = df$blend
    df_new = data.frame(uB, BV, Vr, ri, iz, zY, YJ, JH, i, FLUX_RADIUS, blend)
    df_new
}

load.alldata = function(){
    df_noblend = read.csv(paste(pdir, '/data/new_run/pure_table_err.csv', sep=''))
    df_weak = read.csv(paste(pdir, '/data/new_run/weak_table_err.csv', sep=''))
    df_strong = read.csv(paste(pdir, '/data/new_run/strong_table_err.csv', sep=''))

    df.pure = colorify.full(df_noblend)
    df.pure$type = 'pure'
    df.weak = colorify.full(df_weak)
    df.weak$type = 'weak'
    df.strong = colorify.full(df_strong)
    df.strong$type = 'strong'
    df.strong$blend = 1
    df.master = rbind(df.pure, df.weak, df.strong)

    df.master
}
full.fat = load.alldata()

nstrong = 5
nweak = 5
npure = 10

suffix = paste(nstrong, nweak, npure, sep="_") 

strong.subset = (full.fat[full.fat$type == 'strong', ])[sample(sum(full.fat$type == 'strong'), size=nstrong, replace=FALSE),]
weak.subset = (full.fat[full.fat$type == 'weak', ])[sample(sum(full.fat$type == 'weak'), size=nweak, replace=FALSE),]
pure.subset = (full.fat[full.fat$type == 'pure', ])[sample(sum(full.fat$type == 'pure'), size=npure, replace=FALSE),]

df.subset = rbind(strong.subset, weak.subset, pure.subset)
write.csv(df.subset, paste(pdir, "/data/subset/subset_table_", suffix, ".csv", sep = ""))
write.csv(df.subset, paste(pdir, "/data/subset/subset_table.csv", sep = ""))
