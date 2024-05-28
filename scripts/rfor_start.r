library(dplyr)
library(randomForest)
# library(tibble)

# Environment variables used to organize files on Perlmutter
pdir = Sys.getenv("PSCRATCH")   
hdir = Sys.getenv("HOME")
outdir = paste(pdir, '/output/new_run', sep='')


###########################################################
# Option 1: Load in master data frame and split 
# the training-testing manually
###########################################################
# Load in the master dataframe 
df.master = read.csv(paste(pdir,'/output/master.csv', sep=''))

# Separate into pure/strong/weak
pure.df = df.master[df.master$type=='pure',]
weak.df = df.master[df.master$type=='weak',]
strong.df = df.master[df.master$type=='strong',]

# Create training-testing split of each sub dataframe
# Using a 50-50 split for now (try changing this)!
pure.sample = sample(c(TRUE, FALSE), nrow(pure.df), replace=TRUE, prob=c(0.5,0.5))
weak.sample = sample(c(TRUE, FALSE), nrow(weak.df), replace=TRUE, prob=c(0.5,0.5))
strong.sample = sample(c(TRUE, FALSE), nrow(strong.df), replace=TRUE, prob=c(0.5,0.5))

# Apply the split
train.pure = pure.df[pure.sample,]
train.weak = weak.df[weak.sample,]
train.strong = strong.df[strong.sample,]

test.pure = pure.df[!pure.sample,]
test.weak = weak.df[!weak.sample,]
test.strong = strong.df[!strong.sample,]

# Combine the dataframes into one training and testing object
train.df = rbind(train.pure, train.weak, train.strong)
test.df = rbind(test.pure, test.weak, test.strong)

###########################################################
# Option 2: Load in the individual dataframes
# I'll leave this as an exercise for the reader :)
###########################################################


###########################################################
# Split into features + labels and then run RF!
###########################################################

# Let's give the model all the rows for now, maybe try removing
# some rows and see how the performance changes!
train.phot = train.df[,1:10]
test.phot = test.df[,1:10]

# We have some flexibility in what the label is. If it is a numerical
# value the package defaults to a regression tree versus a categorical
# tree if we give factors. Try out different different setups to see
# how the performance changes!
train.type = as.factor(train.df$type)
test.type = as.factor(test.df$type)
# label.master = train.df$blend 

# Create the forest!
ntree = 500
# Setting importance and keep.inbag to FALSE to save on computation
out.rf = randomForest(x=train.phot,y=train.type,importance=FALSE, ntree=ntree, keep.inbag=FALSE, norm.votes=FALSE)
rcont = list("method"=out.rf)

# Save the forest somewhere so you don't have to re-run this every time.
save(rcont, file=paste(outdir,"/rfobj", suffix, sep=''))


###########################################################
# Performance 
# Quantify how we did using the testing set.
###########################################################

# Traditionally we could use the OOB score from RF but to compare
# with SOM/k-NN/other methods we are using a train-test split.
# To compare against those methods we've actually done the training-testing 
# split ahead of time and kept those as static files (instead of what we are
# doing here where we make our own split each time). 

# We quantify our performance with recall-cost graphs. Check the paper
# for the definitions


# Make predictions for each of the test rows and get the probability
# for each label!
test.pred = predict(out.rf, test.phot, type='prob')

# Repackage into dataframe to make the next steps easier
pred.df = data.frame(test.pred)
pred.df$blend = 1 - pred.df$pure # The blend score is 1-pure! 


# An arbitrary threshold that if the blend score is greater than
# that row is labeled a blend.
blend.thresh = .4 
pred.blend = pred.df$blend > blend.thresh
blend.filt = test.type!='pure'


# How many objects do we label as blends?
cost = sum(pred.blend)/nrow(test.df)

# How many blends do we catch?
blend.num = sum(pred.blend[blend.filt])

recall = blend.num/(sum(blend.filt))

# Success!
print(paste("At ", blend.thresh, " threshold we have a cost of ", cost, " and recall of ", recall, sep=''))