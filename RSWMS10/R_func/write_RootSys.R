
##################################
# R script to convert CPlantOx Rootsysem file .txt into RootSys file for RSWMS
# Adrien Heymans June 2022
#############################


require(tidyverse)
require(data.table)

`%!in%` <- compose(`!`, `%in%`)

trans_rootsys_line <- function(x, raw_file){
  a <- c(unlist(str_split(raw_file$`Time:`[x], " ")),unlist(str_split(raw_file$`Time:`[x+1], " ")))
  a <- as.numeric(a[a != ""])     
  return(a)
}

read_RootSys <- function(path){
  
  raw_file <- read_delim(path, " ", show_col_types = FALSE)
  
  
  po_colnam <- grep("segID#", raw_file$`Time:`)
  po_end <- grep("Total # ", raw_file$`Time:`)[4]
  sta <- po_colnam
  end <- po_end
  cover <- (end-sta)/2
  
  cn <- c(unlist(str_split(raw_file$`Time:`[po_colnam], " ")),unlist(str_split(raw_file$`Time:`[po_colnam+1], " ")))
  cn <- cn[cn != ""]
  cn <- cn[cn != "origination"]
  
  
  data <- matrix(nrow = cover, ncol = 11)
  
  for(i in 1:(cover-1)){
    j <- sta+(i*2)
    data[i,] <- trans_rootsys_line(j, raw_file)
  }
  
  data <- as.data.frame(data)
  colnames(data) <- cn
  
  data <- data %>% filter(!is.na(x))%>%
    mutate(x1 = NA, x2 = x, y1 = NA, y2 = y, z1 = NA, z2 = z)
  
  
  for(i in 2:nrow(data)){
    data$x1[i] <- data$x[data$`segID#` == data$prev[i]]
    data$y1[i] <- data$y[data$`segID#` == data$prev[i]]
    data$z1[i] <- data$z[data$`segID#` == data$prev[i]]
  }
  data <- data %>% filter(!is.na(x1))
  return(data)
}

plot_RS <- function(root_data, type = "time"){
  
  pl <- root_data %>%
    ggplot()+
    geom_segment(aes_string(x = "x1", xend = "x2", y = "z1", yend = "z2", colour = type), size = 1)+
    coord_fixed()+
    theme_dark()+
    viridis::scale_colour_viridis()
  print(pl)
  
}

write_RootSys <- function (data, out_file = "RootSys.1"){
  
  data <- reformat_rootbox(data)
  
  data <- data %>%mutate(length = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2),
                         surface = length*radius*2*pi,
                         mass = pi*radius^2 * length * 0.062)
  
  n_axes <- unique(data$branchID[data$type %!in% c(2,3)])
  age_max = max(data$age)
  
  x1 <- paste0("Time:\n   ",age_max,"\n", "\n",
               "Numbers of seeds\n            1\n",
               "\n",
               "ID, X and Y coordinates of the seeds (one per line)\n     ",
               "1   0    0  \n", "\n",
               "Root DM, shoot DM, leaf area:\n",
               "   ", sum(data$mass), "     0.0000000000000000       0.0000000000000000    \n",
               "\n",
               "Average soil strength and solute concentration experienced by root system:\n",
               "0.29343392611233005   0.0000000000000000\n \n",
               "Total # of axes:\n             ",
               length(n_axes),"\n \n",
               "Total # of branches, including axis(es):\n",
               length(unique(data$branchID)),"\n \n",
               "Total # of segment records:\n",
               nrow(data)+1,"\n \n",
               "segID#    x          y          z      prev or  br#  length   surface  mass\n origination time\n"
  )
  
  tmp <- tibble (node = data$node2ID+1, x = data$x2, y = data$y2, z = data$z2, 
                 prev = data$node1ID+1, or = ifelse(data$type %in% c(1,4,5), 1, 2), 
                 br = data$branchID, length = data$length,
                 surface = data$surface, mass = data$mass)
  tmp <- rbind(tibble(node = 1, x = data$x1[1], y = data$y1[1], z = data$z1[1], 
                      prev = 0, or = 1, br = 1, 
                      length = data$length[1], surface = data$surface[1], mass = data$mass[1]), tmp)
  x2 <- NULL
  for(k in 1:nrow(tmp)){
    if(k == 1){x_tmp <- paste0(paste0(tmp[k,], collapse = "  "), "\n", 0, "\n")}else{
      x_tmp <- paste0(paste0(tmp[k,], collapse = "  "), "\n", data$time[k-1], "\n")
    }
    x2 <- paste0(x2, x_tmp)
  }
  x2 <- paste0(x2,"\n",
               "Total # of growing branch tips:\n        ",
               length(unique(data$branchID)), "\n\n")
  
  x3 <- paste0("tipID#    xg          yg          zg      sg.bhd.tp. ord  br#  tot.br.lgth. axs#\n",
               "overlength  # of estblished points\n",
               "time of establishing (-->)\n")
  
  new <- NULL
  h = 0
  for(i in unique(all_roots$branchID)){
    h = h+1
    
    xend <- all_roots$x2[all_roots$branchID == i]
    xend <- xend[length(xend)]
    yend <- all_roots$y2[all_roots$branchID == i]
    yend <- yend[length(yend)]
    zend <- all_roots$z2[all_roots$branchID == i]
    zend <- zend[length(zend)]
    prev <- all_roots$node1ID[all_roots$branchID == i]
    prev <- prev[length(prev)]
    
    tmp <- tibble(id = h, xend, yend, zend, prev+1, 
                  ord = ifelse(all_roots$type[all_roots$branchID == i][1] %in% c(1,4,5),1,2),br = i,
                  le = sum(all_roots$length[all_roots$branchID == i]), ax = 0)
    x_tmp <- paste0(paste0(tmp[1,], collapse = "  "), "\n", 0,"   ",0, "\n")
    x3 <- paste0(x3, x_tmp)
  }
  
  
  X <- paste0(x1, x2, x3)
  
  write(X, paste0(out_file))
}

reformat_rootbox <- function(all_root){
  
  all_roots <- as.data.table(all_roots)
  all_roots <- all_roots%>%
    transmute(node1ID = node1ID,
              node2ID = node2ID,
              branchID = branchID,
              x1 = x1, y1 = y1, z1 = z1, x2 = x2, y2= y2, z2 = z2,
              radius = radius,
              length = sqrt((all_roots$x2 - all_roots$x1)^2 + (all_roots$y2 - all_roots$y1)^2 + (all_roots$z2 - all_roots$z1)^2),
              time = time,
              type = type,
              age = age)%>%
    arrange(time)
  
  all_roots$node2ID <- 1:nrow(all_roots)
  all_roots$node1ID[all_roots$branchID == 1][1] <- 0
  all_roots$node1ID[all_roots$branchID == 1][-1] <- which(all_roots$branchID == 1)[-length(which(all_roots$branchID == 1))] # tap root ordination
  
  for(i in unique(all_roots$branchID)[-1]){
    all_roots$node1ID[all_roots$branchID == i][-1] <- which(all_roots$branchID == i)[-length(which(all_roots$branchID == i))]
    if(all_roots$type[all_roots$branchID == i][1] %in% c(4,5)){ # connection with the collar
      all_roots$node1ID[all_roots$branchID == i][1] <- 0
    }
    if(all_roots$type[all_roots$branchID == i][1] %in% c(2,3)){ # connection with the parental root
      x1_child <- all_roots$x1[all_roots$branchID == i][1]
      y1_child <- all_roots$y1[all_roots$branchID == i][1]
      z1_child <- all_roots$z1[all_roots$branchID == i][1]
      
      tmp_time <- all_roots$time[all_roots$branchID == i][1]
      
      nearest <- all_roots%>%filter(branchID != i)%>%
        mutate(euc = sqrt((x1-x1_child)^2+ (y1 - y1_child)^2 + (z1 - z1_child)^2))
      nearest <- nearest[nearest$euc == min(nearest$euc), ]
      all_roots$node1ID[all_roots$branchID == i][1] <- nearest$node2ID[1] # oldest segments
    }
  }
  oups <- which(all_roots$node1ID == all_roots$node2ID)
  if(length(oups) > 0){
    for(o in oups){
      self_seg_age <- all_roots$time[all_roots$node2ID == o][1]
      self_seg_id <- all_roots$branchID[all_roots$node2ID == o][1]
      
      nearest <- all_roots%>%filter(branchID == self_seg_id, time < self_seg_age)
      all_roots$node1ID[all_roots$node2ID == o][1] <- nearest$node2ID[nearest$time == max(nearest$time)]
    }
  }
  return(all_roots)
  
}
