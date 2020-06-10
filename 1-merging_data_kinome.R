####################################
## Loading library              
####################################

library(plyr) 
library(matrixStats)
library (treemap)
library (GO.db)
library (GSEABase)
library (ggplot2)
library (splitstackshape)
library (lattice)
library (gplots)

#######################################################
### Loading UniProt
#######################################################

## Loading UniProt
uniprot <- read.csv ("./uniprot/uniprot-all5.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
names(uniprot)[names(uniprot)=="Entry.name"] <- "ID"
uniprot <- uniprot [c("ID","Gene.names...primary..","Protein.names","Function..CC.","Sequence.similarities","Gene.ontology..GO.","Gene.ontology..biological.process.","Gene.ontology..molecular.function.","Gene.ontology..cellular.component.","Gene.names","Protein.families","Pathway","Keywords","Entry","Length","Mass","Status")]

####################################
## Merging data          
####################################

setwd ("./kinome/")
j_data <- read.csv ("data.csv", quote="", header=TRUE, sep=",", stringsAsFactors = FALSE)
j_table <- read.csv ("table.txt", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)

j_data$Kinase <- toupper(j_data$Kinase)
j_table$RBC.Name <- toupper(j_table$RBC.Name)

j_data_table <- merge (j_data, j_table, by.x="Kinase", by.y="RBC.Name",all.x=TRUE, sort=F)
j_data_table2 <- concat.split.multiple(data = j_data_table, split.cols = "Protein.Accession", seps=";", direction = "long") ## Ignore the warning

id.a <- j_data_table2$Protein.Accession
uniprot_pattern <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
m <- regexpr(uniprot_pattern, id.a)
id.clean <- regmatches(id.a, m)
j_data_table2$Protein.Accession <- id.clean
j_data_table2$Data.mean <- rowMeans(j_data_table2[,c("Data.1","Data.2")])
j_data_table2$Data.SD <- rowSds(as.matrix(j_data_table2[,c("Data.1","Data.2")]))

hist (j_data_table2$Data.mean,n=100,xlim=c(0,50))
boxplot (j_data_table2$Data.mean,ylim=c(0,100))
median(j_data_table2$Data.mean)
  
j_data_table2_selected  <- j_data_table2 [j_data_table2$Data.mean <= 0.5*median(j_data_table2$Data.mean),]
j_data_table2_uniprot <- merge (j_data_table2,uniprot,by.x="Protein.Accession",by.y="Entry",all.x=TRUE,sort=F)
write.csv (j_data_table2_uniprot,"j_data_table2_uniprot.csv")

####################################
## GO enrichment analyis (GO pvalue)
####################################

go_fishers <- function (input, ref) {
  input_keycol <- input$Gene.ontology..biological.process.
  whole_uniprot_keycol <- ref$Gene.ontology..biological.process.
  go_col <- strsplit(input_keycol, split="; ")
  go_col_names <- unique (unlist (go_col))
  go_col <- gsub ("(^.*)(\\[)(.*)(\\])", "\\3", go_col_names)
  go_col_names <- gsub ("( \\[)(.*)(\\])", "", go_col_names)
  total_PD <- nrow (input)
  total_Whole_UniProt <- nrow (ref)
  
  #go_col_slim <- merge (data.frame(GO=go_col), go2goslim, by="GO",all.x=TRUE)
  
  pv.go.total<- NULL
  pv.go<-NULL
  
  for (i in 1:length (go_col)) {
    row <- go_col[i]
    nb_found_in_PD <- length (grep (row, input_keycol))
    id_found_in_PD <- input[grep (row, input_keycol),"ID"]
    nb_found_in_whole_Uniprot <- length (grep (row, whole_uniprot_keycol))
    counts <- matrix (data=c(nb_found_in_PD, (total_PD-nb_found_in_PD), nb_found_in_whole_Uniprot, (total_Whole_UniProt-nb_found_in_whole_Uniprot)), nrow=2)
    pv <- fisher.test(counts)
    pv.go <- data.frame(GO=go_col[i],GO.name=go_col_names[i], NB_REF_TOTAL=total_Whole_UniProt, NB_REF_FOUND=nb_found_in_whole_Uniprot, NB_INPUT_TOTAL=total_PD,NB_INPUT_FOUND=nb_found_in_PD, p.value=pv$p.value, ratio=as.matrix(pv$estimate), ID_INPUT_FOUND=paste0(id_found_in_PD$ID,collapse = "|"),stringsAsFactors = F)
    pv.go.total <- rbind (pv.go.total, pv.go)
    print (paste0("Step:",i,"/",length(go_col)," ",go_col[i]," p-value: ",pv.go$p.value," odds ratio: ", pv.go$ratio," ",go_col_names[i]))
    #print (paste0("Step:",i," ",go_col_names[i],"   ",pv.go))
  }
  
  col_fdr <- p.adjust(pv.go.total$p.value,method="BH")
  pv.go.total <- data.frame(pv.go.total,fdr=col_fdr, stringsAsFactors = F)
  return (pv.go.total)
}



go_fishers_v3 <- function (input, ref, idtrue) {  ## Revised on 8.18.2017, correcting the 2x2 contigency table
  total_PD <- nrow (input)
  total_Whole_UniProt <- nrow (ref)
  input_keycol <- input [,c("ID","Gene.ontology..biological.process.")]
  input_keycol <- concat.split.multiple(input_keycol,split.cols = "Gene.ontology..biological.process.",seps = "; ", direction = "long")
  ref_keycol <- ref [,c("ID","Gene.ontology..biological.process.")]
  ref_keycol <- concat.split.multiple(ref,split.cols = "Gene.ontology..biological.process.",seps = "; ", direction = "long")
  input_go_col <- gsub ("(^.*)(\\[)(.*)(\\])", "\\3", input_keycol$Gene.ontology..biological.process.)
  input_go_col_names <- gsub ("( \\[)(.*)(\\])", "", input_keycol$Gene.ontology..biological.process.)
  ref_go_col <- gsub ("(^.*)(\\[)(.*)(\\])", "\\3", ref_keycol$Gene.ontology..biological.process.)
  ref_go_col_names <- gsub ("( \\[)(.*)(\\])", "", ref_keycol$Gene.ontology..biological.process.)
  go_col <- unique (input_keycol$Gene.ontology..biological.process.)
  go_col_names <- gsub ("( \\[)(.*)(\\])", "", go_col)
  go_col <- gsub ("(^.*)(\\[)(.*)(\\])", "\\3", go_col)
  
  pv.go.total<- NULL
  pv.go<-NULL
  for (i in 1:length (go_col)) {
    row <- go_col[i]
    nb_found_in_PD <- length (grep (row, input_go_col,fixed=T))
    nb_found_in_whole_Uniprot <- length (grep (row, ref_go_col,fixed=T))
    aa <- nb_found_in_PD
    bb <- total_PD
    cc <- nb_found_in_whole_Uniprot
    dd <- total_Whole_UniProt
    counts <- matrix (data=c(aa, (cc-aa), (bb-aa), (dd-cc-(bb-aa))), nrow=2)
    pv <- fisher.test(counts, alternative = "greater")
    
    if (idtrue=="FALSE") {
      pv.go <- data.frame(GO=go_col[i], GO.name=go_col_names[i], NB_REF_TOTAL=total_Whole_UniProt, NB_REF_FOUND=nb_found_in_whole_Uniprot, NB_INPUT_TOTAL=total_PD,NB_INPUT_FOUND=nb_found_in_PD, p.value=pv$p.value, ratio=as.matrix(pv$estimate))
    }
    
    if (idtrue=="TRUE") {
      id_found_in_PD <- input_keycol [grep (row, input_keycol$Gene.ontology..biological.process.),"ID"]
      pv.go <- data.frame(GO=go_col[i], GO.name=go_col_names[i], NB_REF_TOTAL=total_Whole_UniProt, NB_REF_FOUND=nb_found_in_whole_Uniprot, NB_INPUT_TOTAL=total_PD,NB_INPUT_FOUND=nb_found_in_PD, p.value=pv$p.value, ratio=as.matrix(pv$estimate), ID_INPUT_FOUND=paste0(id_found_in_PD$ID,collapse = "|"))
    }
    
    pv.go.total <- rbind (pv.go.total, pv.go)
    # print (paste0("Step:",i,"/",length(go_col)," ",go_col[i]," p-value: ",pv.go$p.value," odds ratio: ", pv.go$ratio," "))
    #print (paste0("Step:",i," ",go_col_names[i],"   ",pv.go))
  }
  col_fdr <- p.adjust(pv.go.total$p.value,method="BH")
  pv.go.total <- data.frame(pv.go.total,fdr=col_fdr)
  return (pv.go.total)
}


input <- j_data_table2_selected
input <- merge (j_data_table2_selected,uniprot,by.x="Protein.Accession",by.y="Entry",all.x=TRUE,sort=F)
input <- input [!is.na(input$ID),]
write.csv (input,"j_data_table2_selected_uniprot.csv")

ref <- uniprot

pv_table_full <- tmp.pv.GO <- go_fishers(input,ref)
write.csv (pv_table_full, "pv_table_full.csv")

####################################
## GO enrichment analyis (Treemap)
####################################
setwd ("..")

fl <- system.file("extdata", "./import/goslim_generic.obo.txt", package="GSEABase")
slim <- getOBOCollection("./import/goslim_generic.obo.txt")

pv_GO_BP <- read.csv ("pv_table_full.csv", quote="\"", header=T, sep=",", stringsAsFactors = FALSE)

###############################################################################
setwd ("./kinome/")

pv_GO_BP_select <- pv_table_full[(pv_table_full$fdr < 0.1 & pv_table_full$ratio >= 1.5), ]
pv_GO_BP_select <- data.frame(pv_GO_BP_select,stringsAsFactors = F)
z_total <- NULL
for (i in 1:nrow(pv_GO_BP_select)) {
  myCollection <- GOCollection (pv_GO_BP_select$GO[i])
  ab <- try (goSlim(myCollection, slim, "BP"),silent=TRUE)
  if (is(ab, "try-error")) {
    z <- "NA"
    z2 <- "NA"
  } else {
    z <- as.character(ab$Term[ab$Count >= 1])
    z2 <- rownames(ab)[ab$Count >= 1]
  }
  z <- data.frame(pv_GO_BP_select[i,],GO.slim.name=z, GO.slim=z2)
  z_total <- rbind (z_total, z)
  print (paste0(i,"/",nrow(pv_GO_BP_select)))
}

z_total <- z_total [!(z_total$GO.slim.name == "biological_process"), ]
z_total <- z_total [!(z_total$GO == z_total$GO.slim),]
z_total_selected <- z_total
z_total_selected_top <- z_total_selected[order(z_total_selected$fdr),]
#z_total_selected_top <- z_total_selected_top [! ((z_total_selected_top$GO.slim == "GO:0044403")), ]
z_total_selected_top <- head(z_total_selected_top,n=30)
z_total_selected_signaling <- z_total_selected[z_total_selected$GO.slim.name == "signal transduction",]
z_total_selected_immune <- z_total_selected[z_total_selected$GO.slim.name == "immune system process",]

write.csv(z_total_selected,"z_total_selected.csv")
pdf ("treemap_selected_Jiang.pdf", width=15, height=10)
z_total_selected$lg_fdr <- log10 (z_total_selected$fdr)
treemap(z_total_selected,index=c("GO.slim.name","GO.name"),vSize="ratio",vColor="lg_fdr",type="value",palette="RdYlBu",mapping=c(-65,-1,0))
dev.off()

pdf ("treemap_selected_Jiang_TOP.pdf", width=8, height=6)
z_total_selected_top$lg_fdr <- log10 (z_total_selected_top$fdr)
treemap(z_total_selected_top,index=c("GO.slim.name","GO.name"),vSize="ratio",vColor="lg_fdr",type="value",palette="RdYlBu",mapping=c(-65,-1,0))
dev.off()

pdf ("treemap_selected_Signaling.pdf", width=15, height=10)
z_total_selected_signaling$lg_fdr <- log10 (z_total_selected_signaling$fdr)
treemap(z_total_selected_signaling,index=c("GO.slim.name","GO.name"),vSize="ratio",vColor="lg_fdr",type="value",palette="RdYlBu",mapping=c(-65,-1,0))
dev.off()

pdf ("treemap_selected_immune.pdf", width=15, height=10)
z_total_selected_immune$lg_fdr <- log10 (z_total_selected_immune$fdr)
treemap(z_total_selected_signaling,index=c("GO.slim.name","GO.name"),vSize="ratio",vColor="lg_fdr",type="value",palette="RdYlBu",mapping=c(-65,-1,0))
dev.off()
