#individual
source("./treemix-1.13/src/plotting_funcs.R")
png("tree_allsample.png",width=1200,height=1600,res=100)
plot_tree("tree_allsample")
png("res_allsample.png",width=600,height=600,res=100)
plot_resid("tree_allsample", "sample_order.txt")
#popualtion
source("./treemix-1.13/src/plotting_funcs.R")
png("tree_allfamily.png",width=600,height=600,res=100)
plot_tree("tree_allfamily")
png("res_allfamily.png",width=400,height=400,res=100)
plot_resid("tree_allfamily", "family_order.txt")

print("Done treemix")
