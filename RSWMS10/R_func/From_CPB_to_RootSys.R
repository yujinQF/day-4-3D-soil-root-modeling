
# set your working directory
# setwd()

#load function
source("./write_RootSys.R")

#import root system
all_roots <- fread("~/../Downloads/root_test.txt")

# Display root system 
plot_RS(all_roots)

# write rootSys file
# With root system above 5000 root segments, it takes a long time
write_RootSys(data = all_roots, out_file = "~/GitHub/GitHub/day-4-3D-soil-root-modeling/RSWMS10/RSWMS/in/RootSys")
