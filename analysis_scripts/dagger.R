suppressMessages(library(data.tree))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))

#--------------------------
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(crayon))

### Stuff
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#--------------------------
# utils Fn.
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
​
rej_fn<-function(nodes,lvl,num_rejected,node_pval,ms,ls,l,alpha){
  #pick node effective node and leaf number
  ms_d<-ms[nodes]
  ls_d<-ls[nodes]
  p_vals_d<- node_pval[nodes]
  
  
  ### P-value threshold function
  # r is the considered rank
  crit_func <- function(r,alpha){
    alpha * ls_d * (ms_d + r + num_rejected[as.character(lvl-1)] - 1) / l / ms_d
  }
  
  r <- length(p_vals_d)
  # Determine the appropriate threshold value
  while (sum(p_vals_d <= crit_func(r,alpha)) < r){
    r <- r-1 
  }  
  R <- r 
  tmp_rejected<-nodes[which(p_vals_d <= crit_func(R,alpha))]
  
  return(tmp_rejected)
}
​
detect_inter_cage_cl_fn<-function(feature_Grange,feature_pval_tbl,res_num){
  
  fn_env<-environment()

  feature_bin_n_l<-lapply(1:nrow(feature_pval_tbl),function(x){
    cl_Grange<-   GRanges(seqnames=feature_pval_tbl$chr[x],
                          ranges = IRanges(start=as.numeric(feature_pval_tbl$bins[[x]]),
                                           end=as.numeric(feature_pval_tbl$bins[[x]]) + res_num[feature_pval_tbl$res[x]]-1
                          ))
    return(length(unique(queryHits(findOverlaps(cl_Grange,feature_Grange)))))
    
  })
  feature_pval_tbl<-feature_pval_tbl %>% mutate(feature.bin=unlist(feature_bin_n_l))
  return(feature_pval_tbl)
  
}

produce_depth_tbl_fn<-function(res_cage_set,g_bpt,dagger_roots){
  # Produce CAGE-tree for each resolution
  tmp_g<-induced_subgraph(g_bpt,res_cage_set)
  tmp_g_comp<-components(tmp_g)
  #assign levels/depth to this node set
  ## For each component, depth is the graph distance with the corresponding DAGGER-roots (BPT-leaves)
  comp_set<-which(tmp_g_comp$csize > 1)
  comp_tbl<-do.call(bind_rows,lapply(comp_set,function(x){
    tmp_node<-names(which(tmp_g_comp$membership==x))
    #not working
    tmp_d<-distances(tmp_g,tmp_node[which(tmp_node %in% dagger_roots)],tmp_node,mode = 'in')
    tmp_lvl<-apply(tmp_d,2,function(x)max(x[!(is.infinite(x))])+1)
    tmp_lvl[names(which(apply(tmp_d,2,function(a)any(a==0))))]<-1
    return(tibble(lvl=tmp_lvl,node=names(tmp_lvl),comp=x))
  }))
  gen_lvl<-max(comp_tbl$lvl)
  if(is.infinite(gen_lvl)){gen_lvl<-1}
  isl_set<-which(tmp_g_comp$csize == 1)
  isl_tbl<-do.call(bind_rows,lapply(isl_set,function(x){
    tmp_node<-names(which(tmp_g_comp$membership==x))
    #not working
    return(tibble(lvl=gen_lvl,node=tmp_node,comp=x))
  }))
  node_m_lvl_tbl<-comp_tbl%>%bind_rows(.,isl_tbl)%>%dplyr::rename(m_lvl=lvl)
  return(node_m_lvl_tbl)
  
}

produce_eff_n_and_l_fn<-function(dagger_leaf,node_pval,node_m_lvl_tbl,node_dagger_children,node_dagger_parent){
  
  l<-length(dagger_leaf)
  
  #Compute recursively the effective number of leaves and nodes
  ### Computation for effective leaf/node number is done from DAGGER-leaves to DAGGER-roots
  ms<-rep(0,nrow(node_m_lvl_tbl))
  names(ms)<-node_m_lvl_tbl$node
  ls<-rep(0,nrow(node_m_lvl_tbl))
  names(ls)<-node_m_lvl_tbl$node
  #assign 1 to DAGGER-leaf clusters
  ls[dagger_leaf]<-1
  ms[dagger_leaf]<-1
  #Loop through levels in decreasing DAGGER depth order (starting with DAGGER-leaves)
  for(lvl in rev(sort(unique(node_m_lvl_tbl$m_lvl)))){
    #print(lvl)
    tmp_node<-unlist(node_m_lvl_tbl%>%filter(m_lvl==lvl)%>%dplyr::select(node))
    tmp_leaves_idx<-which(tmp_node %in% dagger_leaf)
    if(length(tmp_leaves_idx)>0){tmp_node<-tmp_node[-tmp_leaves_idx]}
    if(length(tmp_node)<1){next}
    for(p in tmp_node){
      
      ms[p]<-1+sum(ms[node_dagger_children[[p]]]/unlist(lapply(node_dagger_children[[p]],function(x)length(node_dagger_parent[[x]]))))
      
      
      ls[p]<-sum(ls[node_dagger_children[[p]]]/unlist(lapply(node_dagger_children[[p]],function(x)length(node_dagger_parent[[x]]))))
    }
  }
  
  return(list(eff.node=ms,eff.leaf=ls))
}

compute_rejected_test_fn<-function(alpha,node_m_lvl_tbl,node_dagger_parent,node_pval,ms,ls,l,chromo,tmp_res){
  #vector recording the number of rejections at every depth
  num_rejected = rep(0, 1 + max(node_m_lvl_tbl$m_lvl)) 
  names(num_rejected)<-as.character(seq(0, max(node_m_lvl_tbl$m_lvl)))
  #vector recording the actual nodes rejecting the null hypothesis
  rejections<-rep(F,nrow(node_m_lvl_tbl))
  names(rejections)<-node_m_lvl_tbl$node
  #loop through increasing depth levels (starting from DAGGER-roots)
  for (lvl in seq(1, max(node_m_lvl_tbl$m_lvl))){
    #print(lvl)
    nodes_depth_d <- unlist(node_m_lvl_tbl%>%filter(m_lvl==lvl)%>%dplyr::select(node)) 
    
    # Delete the nodes one of whose parents has not been rejected.
    if ( lvl > 1){
      
      nodes_depth_d<-nodes_depth_d[unlist(lapply(nodes_depth_d,function(x){
        all(node_dagger_parent[[x]] %in% names(which(rejections)))
      }))]
      if(any(is.na(node_pval[nodes_depth_d]))){
        #further filter out the nodes for which we don't have p-values
        nodes_depth_d<-names(which(!(is.na(node_pval[nodes_depth_d]))))
        
      }
      
    }
    # Performs the rejection step at depth d.  
    rejected_nodes_depth_d <- rej_fn(nodes_depth_d, lvl, num_rejected,node_pval,ms,ls,l,alpha)
    rejections[rejected_nodes_depth_d] <- T
    num_rejected[as.character(lvl)] <- num_rejected[as.character(lvl-1)] + length(rejected_nodes_depth_d)
    
  }
  return(tibble(chr=chromo,res=tmp_res,node=names(which(rejections)),FDR=alpha,emp.pval=node_pval[names(which(rejections))]))
  
}

mres_DAGGER_fn<-function(chr_pval_tbl,chromo,BHiCect_res_file,alpha_seq,res_num){
  cat(green(chromo), ": DAGGER initialised \n")
  # Build the BHiCect tree
  load(paste0(BHiCect_res_file,chromo,"_spec_res.Rda"))
  cat(green(chromo), ": Build tree \n")
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  #collect leaves and ancestors
  cat(green(chromo), ": Collect children/parent nodes \n")
  tmp_leaves<-chr_bpt$Get('name',filterFun=isLeaf)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  #build the cage-containing sub-tree
  cage_node<-unlist(chr_pval_tbl%>%dplyr::select(cl))
  cage_set<-unique(c(cage_node,unique(unlist(node_ancestor[cage_node])))) 
  Prune(chr_bpt, function(x) x$name %in% cage_set)
  p_node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  
  #rebuild corresponding tree according to DAGGER
  node_dagger_children<-lapply(p_node_ancestor,'[',2)
  #eleminate Root node to "create" DAGGER leaves
  node_dagger_children<-node_dagger_children[-1]
  
  node_dagger_children<-lapply(node_dagger_children,function(x){
    if(x == "Root"){return(NULL)} else{
      return(x)
    }
  })
  #Create DAGGER-parent mapping (immediate BPT children of each node)
  node_dagger_parent<-lapply(cage_set,function(x){
    tmp<-names(which(unlist(node_dagger_children) == x))
    return(unlist(lapply(strsplit(tmp,split="\\."),'[',1)))
  })
  names(node_dagger_parent)<-cage_set
  
  # Build igraph object for CAGE subtree
  g_bpt<-as.igraph.Node(chr_bpt,directed = T,direction = 'climb')
  
  # Loop through each resolution to compute the DAGGER p-value correction
  tmp_res_l<-vector('list',length(unique(chr_pval_tbl$res)))
  names(tmp_res_l)<-unique(chr_pval_tbl$res)
  chr_res_set<-names(sort(res_num[names(tmp_res_l)]))
  for(tmp_res in chr_res_set){
    cat(green(chromo), " DAGGER at ", tmp_res, " \n")
    res_cage_node<-unlist(chr_pval_tbl%>%filter(res==tmp_res)%>%dplyr::select(cl))
    res_cage_set<-unique(c(res_cage_node,grep(paste0(tmp_res),unique(unlist(node_ancestor[res_cage_node])),value=T)))
    #Filter out any cluster without any enriched BPT-children at higher resolution
    if(tmp_res != "5kb"){
      res_cage_set<-res_cage_set[sapply(res_cage_set,function(x){
        any(unlist(lapply(p_node_ancestor[unlist(do.call(bind_rows,tmp_res_l)$node)],function(y){
          x %in% y
        })))
      })]
    }
    if(length(res_cage_set)<2){
      tmp_pval<-unlist(chr_pval_tbl%>%filter(cl%in%res_cage_set)%>%dplyr::select(emp.pval))
      tmp_res_l[[tmp_res]]<-tibble(chr=chromo,res=tmp_res,node=res_cage_set,FDR=NA,emp.pval=tmp_pval)
      next
    }
    ## detect the DAGGER leaves
    ### For each component detect top parent node and label it as leaf
    #### DAGGER leaves are nodes who don't have any BPT-ancestors among the node at the considered resolution! 
    dagger_leaf<-names(which(unlist(lapply(node_dagger_children[res_cage_set],function(x)sum(grepl(tmp_res,x))))<1))
    l<-length(dagger_leaf)
    
    dagger_roots<-res_cage_set[!(res_cage_set %in% unique(unlist(lapply(p_node_ancestor[res_cage_set],'[',-1))))]
    
    node_m_lvl_tbl<-produce_depth_tbl_fn(res_cage_set,g_bpt,dagger_roots)
    
    #Build the node p-value mapping vector
    node_pval<-unlist(chr_pval_tbl%>%filter(cl %in% res_cage_set)%>%dplyr::select(emp.pval))
    names(node_pval)<-unlist(chr_pval_tbl%>%filter(cl %in% res_cage_set)%>%dplyr::select(cl))
    
    eff_quant<-produce_eff_n_and_l_fn(dagger_leaf,node_pval,node_m_lvl_tbl,node_dagger_children,node_dagger_parent)
    
    ms<-eff_quant[[1]]
    ls<-eff_quant[[2]]
    alpha_res_l<-vector('list',length(alpha_seq))
    names(alpha_res_l)<-as.character(alpha_seq)
    
    for ( alpha in alpha_seq){
      
      alpha_res_l[[as.character(alpha)]]<-compute_rejected_test_fn(alpha,node_m_lvl_tbl,node_dagger_parent,node_pval,ms,ls,l,chromo,tmp_res)
    }
    tmp_res_l[[tmp_res]]<-do.call(bind_rows,alpha_res_l)
  }
  
  return(do.call(bind_rows,tmp_res_l))
}

#--------------------------
feature_coord_file<-"./data/features/H1/CAGE/features.bed"
feature_pval_file<-"./results/H1.CAGE.tsv"
BHiCect_res_file<-"./data/rda_clusters/H1/"

out_tsv<-"p22.pvalues.tsv"

# if bed file
feature_coord_tbl<-import.bedGraph(feature_coord_file)

# Read the pvalues table.
feature_pval_tbl<-read_tsv(feature_pval_file) %>% 
  dplyr::rename(chr=chromosome,cl=name,emp.pval=p_value)

#------------------
# chr_set<-unique(feature_pval_tbl$chr)
chr_set <- c("chr22") 

feature_pval_tbl<-do.call(bind_rows,map(chr_set,function(chromo){
  message(chromo)
  load(paste0(BHiCect_res_file,chromo,"_spec_res.Rda"))
  feature_pval_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[cl],
           res=str_split_fixed(cl,"_",4)[,1])
}))
​
dagger_mres_l<-lapply(chr_set,function(chromo){
  chr_feature_coord_tbl<-feature_coord_tbl [seqnames(feature_coord_tbl) == chromo]
  chr_pval_tbl<-feature_pval_tbl %>% filter(chr==chromo)
  cat(green(chromo), ": select clusters with two feature-containing bins \n")
  chr_pval_tbl<-detect_inter_cage_cl_fn(chr_feature_coord_tbl,chr_pval_tbl,res_num) %>% 
    filter(feature.bin>1)
  alpha_seq<-c(0.01)
  
  return(mres_DAGGER_fn(chr_pval_tbl,chromo,BHiCect_res_file,alpha_seq,res_num))
})

dagger_mres_tbl<-do.call(bind_rows,dagger_mres_l)
write_tsv(dagger_mres_tbl,file=out_tsv)
