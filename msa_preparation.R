
# Prepare msa(multiple sequence alignment) file

con<-dbConnect(MySQL(), user='root', password=NULL, dbname='current_project', host='127.0.0.1', port=3306)

msa<-tbl(con, "msa") %>% 
  collect()

# Separate data to columns "Human", "Mouse" and "C.elegans"
msa<-msa %>%
  separate(fasta, c("Human", "Mouse","C.elegans"), sep = "\n>")

# Remove "\n" to get single sequence text
msa<-as.data.frame(apply(msa[,2:4], 2, function(x) gsub("\n", "", x, perl = TRUE)))
msa<-as.data.frame(apply(msa, 2, function(x) gsub(">", "", x, perl = TRUE)))

# If some rows have mixed data (e.g Mouse column having C. elegans sequence data), swap them.
msa1<-msa
msa1$Mouse[grep("elegans", msa1$Mouse)]<-msa1$C.elegans[grep("elegans", msa1$Mouse)]
msa1$C.elegans[grep("musculus", msa$C.elegans)]<-msa$Mouse[grep("musculus", msa$C.elegans)]

# Two of "Human", "Mouse" and "C_elegans" names will be used in place of org1 and org2 in orthoFinder function
msa1<-msa1 %>%
  separate(Human, c("Human_ID","Human_seq"), " \\[Homo sapiens\\]") %>%
  separate(Mouse, c("Mouse_ID","Mouse_seq"), " \\[Mus musculus\\]") %>%
  separate(C.elegans, c("C_elegans_ID","C_elegans_seq"), " \\[Caenorhabditis elegans\\]")

# set rownames as index numbers
msa1$index<-rownames(msa1)

# Find orthologous variants between C. elegans and Human (gnomad) variants
# Note that we will use previously defined names as org1 and org2 organism names
# df1 and df2 are variant data for org1 and org2, respectively.
celegans_gnomad_ort<-orthoFinder(df1 = celegans, df2 = gnomad, 
                                 org1 = "C_elegans", org2 = "Human", 
                                 msa = msa1, ort = TRUE)

