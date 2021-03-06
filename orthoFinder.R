

orthoFinder<-function (df1, df2, org1, org2, msa, ort = TRUE){
  # Find index of msa having protein ids which both have variants 
  df1_id<-msa$index[msa[[paste0(org1,"_ID")]] %in% unique(df1$Refseq_ID)]
  df2_id<-msa$index[msa[[paste0(org2,"_ID")]] %in% unique(df2$Refseq_ID)]
  common_id<-as.numeric(intersect(df1_id, df2_id))
  df1_df2_ort_list<-vector(mode = "list", length = length(common_id))
  
  j<-0
  pb <- txtProgressBar(min = min(common_id), max = max(common_id), style = 3)
  # Loop over only common alignments
  for (i in common_id){
    # aa positions of variants from selected protein
    df1_aapos<-df1$aapos[df1$Refseq_ID == msa[[paste0(org1,"_ID")]][i]]
    df2_aapos<-df2$aapos[df2$Refseq_ID == msa[[paste0(org2,"_ID")]][i]]
    
    # reference aa of variants
    df1_from<-df1$from[df1$Refseq_ID == msa[[paste0(org1,"_ID")]][i]]
    df2_from<-df2$from[df2$Refseq_ID == msa[[paste0(org2,"_ID")]][i]]
    
    # converted aa
    df1_to<-df1$to[df1$Refseq_ID == msa[[paste0(org1,"_ID")]][i]]
    df2_to<-df2$to[df2$Refseq_ID == msa[[paste0(org2,"_ID")]][i]]
    
    # tables for calculating aa conservation
    df1_conservation<-data.table(from = df1_from, to = df1_to, aapos = df1_aapos)
    df1_conservation<-left_join(df1_conservation, aa_conservation, by = "from")
    df1_conservation<-mutate(df1_conservation, conservation = ifelse(str_detect(df1_conservation$to.y, fixed(df1_conservation$to.x)), 1, 0))
    
    df2_conservation<-data.table(from = df2_from, to = df2_to, aapos = df2_aapos)
    df2_conservation<-left_join(df2_conservation, aa_conservation, by = "from")
    df2_conservation<-mutate(df2_conservation, conservation = ifelse(str_detect(df2_conservation$to.y, fixed(df2_conservation$to.x)), 1, 0))
    
    # get sequence of selected proteins
    df1_seq<-msa[[paste0(org1,"_seq")]][i]
    df2_seq<-msa[[paste0(org2,"_seq")]][i]
    
    # get real (trimmed) sequence of selected proteins
    df1_real_seq<-stri_replace_all_regex(df1_seq, "-", "")
    df2_real_seq<-stri_replace_all_regex(df2_seq, "-", "")
    
    # find index of aa and combine these
    df1_seq_ind<-as.data.table(stri_locate_all_regex(df1_seq, "[A-Z]")[[1]])[[1]]
    df2_seq_ind<-as.data.table(stri_locate_all_regex(df2_seq, "[A-Z]")[[1]])[[1]]
    df1_real_seq_ind<-as.data.table(stri_locate_all_regex(df1_real_seq, "[A-Z]")[[1]])[[1]]
    df2_real_seq_ind<-as.data.table(stri_locate_all_regex(df2_real_seq, "[A-Z]")[[1]])[[1]]
    
    df1_all_seq_ind<-data.table(seq_ind = df1_seq_ind,  aapos = df1_real_seq_ind) %>%
      filter(df1_real_seq_ind %in% df1_aapos)
    
    df2_all_seq_ind<-data.table(seq_ind = df2_seq_ind, aapos = df2_real_seq_ind) %>%
      filter(df2_real_seq_ind %in% df2_aapos)
    
    df1_all_seq_ind<-df1_all_seq_ind[df1_all_seq_ind$seq_ind %in% df2_all_seq_ind$seq_ind]
    df2_all_seq_ind<-df2_all_seq_ind[df2_all_seq_ind$seq_ind %in% df1_all_seq_ind$seq_ind]
    
    # merge by aa positions, get only aa having same reference aa
    if (length(df1_all_seq_ind[[1]]) != 0){
      
      df1_all_seq_ind<-left_join(df1_all_seq_ind, df1_conservation, by = "aapos") %>%
        unique() %>%
        select(!to.y)
      
      df2_all_seq_ind<-left_join(df2_all_seq_ind, df2_conservation, by = "aapos") %>%
        unique() %>%
        select(!to.y)
      
      all_seq_ind<-left_join(df1_all_seq_ind, df2_all_seq_ind, by = "seq_ind") %>%
        filter(from.x == from.y)
      
      # should dismiss conserved aa changes
      if(ort == TRUE){
        all_seq_ind<-filter(all_seq_ind, conservation.x == conservation.y)
      }
      
      # combine both organisms
      j<-j+1
      df1_df2_ort<-setNames(data.table(msa[[paste0(org1,"_ID")]][i],
                                       all_seq_ind$aapos.x,
                                       all_seq_ind$from.x,
                                       all_seq_ind$to.x.x,
                                       msa[[paste0(org2,"_ID")]][i],
                                       all_seq_ind$aapos.y,
                                       all_seq_ind$from.y,
                                       all_seq_ind$to.x.y,
                                       i),
                            c(paste0(org1,"_ID"), paste0(org1,"_aapos"), paste0(org1,"_from"), paste0(org1,"_to"), 
                              paste0(org2,"_ID"), paste0(org2,"_aapos"), paste0(org2,"_from"), paste0(org2,"_to"), "msa_id"))
      
      df1_df2_ort<-df1_df2_ort[df1_df2_ort[[paste0(org1,"_aapos")]] != "",]
      df1_df2_ort_list[[j]]<-df1_df2_ort
    }
    setTxtProgressBar(pb, i)
  }
  df1_df2_ort<-rbindlist(df1_df2_ort_list)
  df1_df2_ort<-unique(df1_df2_ort)
  return(df1_df2_ort)
}
