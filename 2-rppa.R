####################################
## Run 1-merging_data_kinome before
## this script
####################################

library (limma)
library (lattice)

####################################
## Loading data (RPPA)
####################################

setwd ("./rppa/")
j_rppa_1_rawlog_1 <- read.csv ("./1/2-normlinear.csv", quote="", header=TRUE, sep=",", stringsAsFactors = FALSE)
j_rppa_1_rawlog_2 <- read.csv ("./2/2-normlinear.csv", quote="", header=TRUE, sep=",", stringsAsFactors = FALSE)
uniprot <- read.csv ("../uniprot/uniprot-all5.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
names(uniprot)[names(uniprot)=="Entry.name"] <- "ID"
#j_rppa_2_normlinear <- read.csv ("2-normlinear.csv", quote="", header=TRUE, sep=",", stringsAsFactors = FALSE)
#j_rppa_3_normlog2 <- read.csv ("3-normlog2.csv", quote="", header=TRUE, sep=",", stringsAsFactors = FALSE)
#j_rppa_4_normlog2_median_centered  <- read.csv ("4-normlog2_median_centered.csv", quote="", header=TRUE, sep=",", stringsAsFactors = FALSE)

## Correcting a gene name:
j_rppa_1_rawlog_1[j_rppa_1_rawlog_1$Antibody.Name.in.Heatmap=="Oct-4-R-C","Gene.Name"] <- "OCT4"
j_rppa_1_rawlog_2[j_rppa_1_rawlog_2$Antibody.Name.in.Heatmap=="Oct-4-R-C","Gene.Name"] <- "OCT4"

j_rppa_input <- merge (j_rppa_1_rawlog_1, j_rppa_1_rawlog_2, by="Antibody.Name.in.Heatmap")
j_rppa_input <- j_rppa_input [!is.na(j_rppa_input$Gene.Name.x),]

j_rppa_input[grep ("GSK3", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "GSK3A"
j_rppa_input[grep ("CD29", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "ITGB1"
j_rppa_input[grep ("H3K9ME2", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "HIST1H3A"
j_rppa_input[grep ("HISTH3", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "HIST1H3A"
j_rppa_input[grep ("OCT4", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "POU5F1"
j_rppa_input[grep ("PDGFR", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "PDGFRA"
j_rppa_input[grep ("PDHK1", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "PDK1"
j_rppa_input[grep ("PIK3BC", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "PIK3CB"
j_rppa_input[grep ("PTGS3", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "COX4I1"
j_rppa_input[grep ("RIP", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "RIPK1"
j_rppa_input[grep ("RPS6K", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "RPS6KA1"
j_rppa_input[grep ("UGT1A", j_rppa_input$Gene.Name.x),"Gene.Name.x"] <- "UGT1A1"
uniprot <- uniprot[!duplicated (uniprot$ID),]
uniprot <- uniprot[!is.na (uniprot$ID),]

j_rppa_input <- merge (j_rppa_input, uniprot, by.x="Gene.Name.x",by.y="Gene.names...primary..",all.x=T, sort = F)

j_rppa_input <- j_rppa_input[!j_rppa_input$ID == "ARF_HUMAN",]
j_rppa_input <- j_rppa_input [!is.na(j_rppa_input$Gene.Name.x),]
#plot (j_rppa_input$control, j_rppa_input$Control.cells,pch=20)
#plot (j_rppa_input$drug.treated.cells, j_rppa_input$drug.30min,pch=20)
#abline (0,1)

j_rppa_input_raw <- j_rppa_input
duplicated (j_rppa_input_raw$Antibody.Name.x)

#j_rppa_input <- j_rppa_input[,c("Control.cells","control","ligand.control.cells","nanoparticle.control.cells","drug.10min","drug.30min","drug.treated.cells","drug.60min")]
row.names(j_rppa_input) <- j_rppa_input_raw[,c("Antibody.Name.x")]

excluded.id <- c("p53","Fibronectin","EMA","Vimentin","CD29","Collagen−VI","PAI-1","Porin","Syk","Bcl2A1","CD171","Cox2","Mcl−1","YB1_pS102","CD44","NDUFB4","PEA−15_pS116","TFRC","PD−L1","CD49b","UGT1A","HER2","TWIST","TIGAR","Src","SDHA","PDGFR−b","Snail","D−a−Tubulin","PR","N−Ras","Bcl−xL","Cyclin−E1","Rab25","LRP6_pS1490","c−Abl","HSP27","PREX1","IGF1R_pY1135_Y1136","Src_pY527","PKC−delta_pS664","Heregulin","PAK4","HER2_pY1248","c−Kit","C−Raf_pS338","ER−a_pS118","c−Met","Elk1_pS383","SCD","MSI2","ER","Sox2","c−Myc","Rictor_pT1135","MIG6","UBAC1","Atg3","Bid","PTEN","Glutamate−D1−2","53BP1","Bad_pS112","Bad_pS112","Granzyme−B","MSH6","Transglutaminase","GAPDH","Rb_pS807_S811","S6_pS240_S244")
excluded.id <- gsub ("−","-",excluded.id)
excluded.id [!excluded.id %in%  row.names(j_rppa_input)]

j_rppa_input <- j_rppa_input[!row.names(j_rppa_input) %in% excluded.id, ]


j_rppa_input$Ctrl <- (j_rppa_input[,c("control")])
j_rppa_input$Ctrl.norm <- j_rppa_input[,c("control")]/(j_rppa_input[,c("control")])
j_rppa_input$Treated <- j_rppa_input[,c("drug.treated.cells")]/j_rppa_input$Control.cells
j_rppa_input$Ligand.control <- j_rppa_input$ligand.control.cells / j_rppa_input$Control.cells
j_rppa_input$Nanoparticle.control <- j_rppa_input$nanoparticle.control.cells / j_rppa_input$Control.cells
#j_rppa_input$Treated.10min <- j_rppa_input$drug.10min / j_rppa_input$Ctrl
#j_rppa_input$Treated.60min <- j_rppa_input$drug.60min / j_rppa_input$Ctrl
write.csv (j_rppa_input,"j_rppa_input.csv")

#j_rppa_hm <- j_rppa_input[,c("Ligand.control","Nanoparticle.control","Treated.10min","Treated.30min","Treated.60min")]
#j_rppa_hm <- j_rppa_input[,c("Nanoparticle.control","Ligand.control","Treated.60min","Treated.30min","Treated.10min","Ctrl.norm")]
j_rppa_hm <- j_rppa_input[,c("Nanoparticle.control","Ligand.control","Treated","Ctrl.norm")]

z <- log2(j_rppa_hm)
#hfit_row <- hclust(dist(z, method = "euclidean"), method="ward.D")
hfit_row <- order (rowMeans(z[,c("Treated.30min","Treated.60min")]))
hfit_row <- order (z[,c("Treated")])

z <- z[,c("Ctrl.norm","Nanoparticle.control","Ligand.control","Treated")]
z <- as.matrix(z[hfit_row,])
j_rppa_hm_z <- j_rppa_input[hfit_row,]
write.csv (j_rppa_hm_z,"j_rppa_hm_z.csv")

zlim1 <- min (z)
zlim2 <- max (z)
pcol <- colorRampPalette (c("blue4","green4","white","white","yellow","red3"))
pcolscale <- 100
brks <- unique (c(seq(-5,-2,length=pcolscale), seq(-2,-0.2,length=pcolscale), seq(-0.2,0.2,length=pcolscale), seq(0.2,2,length=pcolscale), seq(2,5, length=pcolscale)))
zlimscale <- 50
dev.off()
z <- as.matrix (z)
nRow <- nrow (z)
nCol <- ncol (z)
width_f <-  3+0.2*nCol
height_f <- 2+0.11*nRow
pdf ("hm.pdf", width=width_f, height=height_f)
print (contourplot (t(z), xlab='', ylab='', region =TRUE, contour = FALSE, 
                    margin=F, scales=list(tck=0, x=list(rot=45,cex=0.5)), col.regions=pcol, zlim=c(zlim1,zlim2),
                   at=brks, colorkey = list(space = "top")))
dev.off()


z_order <- order (-abs(rowMeans(z)))
z_short <- z [z_order[1:100],]
j_rppa_hm_z_short <- j_rppa_input [z_order[1:100],]
write.csv (j_rppa_hm_z_short,"j_rppa_hm_z_short.csv")

hfit_row <- hclust(dist(z_short, method = "euclidean"), method="ward.D")
z_short <-z_short[hfit_row$order,]
z_short <- z_short[,c("Ctrl.norm","Nanoparticle.control","Ligand.control","Treated")]

z_short <- as.matrix (z_short)
nRow <- nrow (z_short)
nCol <- ncol (z_short)
width_f <-  3+0.2*nCol
height_f <- 2+0.11*nRow
pdf ("hm_main.pdf", width=width_f, height=height_f)
print (contourplot (t(z_short), xlab='', ylab='', region =TRUE, contour = FALSE, 
                    margin=F, scales=list(tck=0, x=list(rot=45,cex=0.5)), col.regions=pcol, zlim=c(zlim1,zlim2),
                    at=brks, colorkey = list(space = "top")))
dev.off()




z_order <- order (-abs(rowMeans(z)))
z_short <- z [z_order[101:nrow(z)],]
j_rppa_hm_z_short2 <- j_rppa_input [z_order[101:nrow(z)],]
write.csv (j_rppa_hm_z_short2, "j_rppa_hm_z_short2.csv")
                                   
hfit_row <- hclust(dist(z_short, method = "euclidean"), method="ward.D")
z_short<-z_short[hfit_row$order,]

z_short <- as.matrix (z_short)
nRow <- nrow (z_short)
nCol <- ncol (z_short)
width_f <-  2+0.11*nRow
height_f <- 2+0.2*nCol
pdf ("hm_other.pdf", width=width_f, height=height_f)
print (contourplot (z_short, xlab='', ylab='', region =TRUE, contour = FALSE, 
                    margin=F, scales=list(tck=0, x=list(rot=45,cex=0.5)), col.regions=pcol, zlim=c(zlim1,zlim2),
                    at=brks, colorkey = list(space = "right")))
dev.off()


j_rppa_kinome <- merge (j_rppa_input_raw, j_data_table2, by.x="Gene.Name.x",by.y="HUGO.symbol",sort=F)
write.csv (j_rppa_kinome,"j_rppa_kinome2.csv")

