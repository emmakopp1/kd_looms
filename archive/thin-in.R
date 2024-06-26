library(rwty)
library(ape)


phylo = load.trees("/Users/kopp/Documents/kd_looms/data/beast/loom_bcov_1111/loom_bcov_1111.trees",
                   trim = 15)

write.tree(phylo$trees,"/Users/kopp/Documents/kd_looms/data/beast/loom_bcov_1111/loom_bcov_1111_thined.trees")

phylo = load.trees("/Users/kopp/Documents/kd_looms/data/beast/loom_bcov_8421/loom_bcov_8421.trees",
                   trim = 15)

write.tree(phylo$trees,"/Users/kopp/Documents/kd_looms/data/beast/loom_bcov_8421/loom_bcov_8421_thined.trees")
