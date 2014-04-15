#### Read data files into R ####

bc_A = read.csv("HCT_A.coutb.csv", header = TRUE, sep = "\t")
bc_B = read.csv("HCT_B.coutb.csv", header = TRUE, sep = "\t")
rc_A = read.csv("HCT_A.coutc.csv", header = TRUE, sep = "\t")
rc_B = read.csv("HCT_B.coutc.csv", header = TRUE, sep = "\t")
tc_A = read.csv("HCT_A.coutt.csv", header = TRUE, sep = "\t")
tc_B = read.csv("HCT_B.coutt.csv", header = TRUE, sep = "\t")

bc_C2_cel = read.csv("C2_CEL.coutb.csv", header = TRUE, sep = "\t")
rc_C2_cel = read.csv("C2_CEL.coutc.csv", header = TRUE, sep = "\t")
tc_C2_cel = read.csv("C2_CEL.coutt.csv", header = TRUE, sep = "\t")
bc_C2_ctr = read.csv("C2_CTR.coutb.csv", header = TRUE, sep = "\t")
rc_C2_ctr = read.csv("C2_CTR.coutc.csv", header = TRUE, sep = "\t")
tc_C2_ctr = read.csv("C2_CTR.coutt.csv", header = TRUE, sep = "\t")

samples = read.table("samples.txt", header = FALSE, sep = "\t")
samples_c2_cel = read.table("samples_C2_CEL.txt", header = FALSE, sep = "\t")
samples_c2_ctr = read.table("samples_C2_CTR.txt", header = FALSE, sep = "\t")

#rename the columns to correct cell and ctr numbers

colnames(bc_A)[2:97] <- as.vector(samples$V1)
colnames(bc_B)[2:97] <- as.vector(samples$V1)
colnames(rc_A)[2:97] <- as.vector(samples$V1)
colnames(rc_B)[2:97] <- as.vector(samples$V1)
colnames(tc_A)[2:97] <- as.vector(samples$V1)
colnames(tc_B)[2:97] <- as.vector(samples$V1)

colnames(bc_C2_cel)[2:97] <- as.vector(samples_c2_cel$V1)
colnames(bc_C2_ctr)[2:97] <- as.vector(samples_c2_ctr$V1)
colnames(rc_C2_cel)[2:97] <- as.vector(samples_c2_cel$V1)
colnames(rc_C2_ctr)[2:97] <- as.vector(samples_c2_ctr$V1)
colnames(tc_C2_cel)[2:97] <- as.vector(samples_c2_cel$V1)
colnames(tc_C2_ctr)[2:97] <- as.vector(samples_c2_ctr$V1)

#merge the two libraries for C1, rename rownames and remove geneid column

bc <- merge(bc_A,bc_B,by="GENEID", all=TRUE)
rownames(bc) = bc$GENEID
bc = bc[,-1]
rc <- merge(rc_A,rc_B,by="GENEID",all=TRUE)
rownames(rc) = rc$GENEID
rc = rc[,-1]
tc <- merge(tc_A,tc_B,by="GENEID", all=TRUE)
rownames(tc) = tc$GENEID
tc = tc[,-1]
bc.c2 <- merge(bc_C2_cel,bc_C2_ctr,by="GENEID", all=TRUE)
rownames(bc.c2) = bc.c2$GENEID
bc.c2 = bc.c2[,-1]
rc.c2 <- merge(rc_C2_cel,rc_C2_ctr,by="GENEID", all=TRUE)
rownames(rc.c2) = rc.c2$GENEID
rc.c2 = rc.c2[,-1]
tc.c2 <- merge(tc_C2_cel,tc_C2_ctr,by="GENEID", all=TRUE)
rownames(tc.c2) = tc.c2$GENEID
tc.c2 = tc.c2[,-1]

barplot(as.vector(colSums(objects[[1]]))) # check total reads across samples for 1 experiment

#remove outliers of barcode # 68 (= Dicer control#17)
outliers = grep("DICER_CTR_17",colnames(bc)) #same columns in all dataframes
outliers.c2 = grep("DICER_CTR_34",colnames(bc.c2))

bc = bc[,-outliers]
rc = rc[,-outliers]
tc = tc[,-outliers]

bc.c2 = bc.c2[,-outliers.c2]
rc.c2 = rc.c2[,-outliers.c2]
tc.c2 = tc.c2[,-outliers.c2]

#remove deprecated variables
rm(tc_A,tc_B,bc_A,bc_B,rc_A,rc_B, bc_C2_cel, bc_C2_ctr, rc_C2_cel, rc_C2_ctr, tc_C2_cel, tc_C2_ctr)
rm(samples, samples_c2_cel, samples_c2_ctr, outliers, outliers.c2)

#plot oversequencing of experiments
overseq(rc,bc,"C1")
overseq(rc.c2,bc.c2,"C2")

#Set NA's to 0's
tc = apply(tc,2,function(x){x[is.na(x)]=0;return(x)})
tc.c2 = apply(tc.c2,2,function(x){x[is.na(x)]=0;return(x)})

#make object with all experiments in it. l1=names of experiments (1st row of object) l2=data frames (2nd row)
l1<-list("C1","C2")
l2<-list(tc,tc.c2)
c1.2<-makeobject(l1,l2)
rm(l1,l2)

####filtering and normalization of data####

#plot total reads of samples (outputs barplots, use function combineplots to plot them in 1 plot)
combineplots(1,2,1) # (desired number of rows,columns,1=square plots)


for(i in 1:length(c1.2)){
  totalreads(c1.2,i,2,pdf=FALSE)
} # run totalplots for total reads per sample per condition

dev.off() #return to normal plotting

#remove samples with low total read counts and plot samples leftover
for(i in 1:length(c1.2)){
  c1.2[[i]][[3]]<-filtersamples(c1.2,i,2,10000,70000,pdf=TRUE) # create 3rd row in c1.2, containing filtered data sets
} # run filtsamples to remove any samples with < t total samples

#calculate the conversion factor
ERCC_molecules <- read.csv("ERCC_molecules.csv",sep="\t") # extract molecules spiked in per sample

#calculate conversion factor for c1.2, ith column, jth row in c1.2: total read filtered data frames
for(i in 1:length(c1.2)){
  conversion(c1.2,i,3,pdf=TRUE)
}

#remove ERCCs to avoid clustering artefacts
for(i in 1:length(c1.2)){c1.2[[i]][[4]]<-rmspike(c1.2[[i]][[3]])}

#normalize all genes to median of transcript sums
for(i in 1:length(c1.2)){
  c1.2[[i]][[5]]<-norm.median(c1.2,i,4) #create a 5th row in c1.2, containing normalized data sets
}

#check remaining genes after different cutoffs of median gene expression across samples
for(i in 1:length(c1.2)){
  testcutoff(c1.2,i,5,pdf=TRUE)
}

#cutoff genes that have median expression lower than x
for(i in 1:length(c1.2)){
  c1.2[[i]][[6]]<-filtergenes(c1.2,i,5,5)
}

#remove depracated values
rm(bc,bc.c2,rc,rc.c2,tc,tc.c2)

#so far the reading, normalizing and filtering of data

#the object created here has 5 rows: 1=names of dataset to be used, 
#2=starting datasets, manually constructed by merging datasets of interest and grepping for subsets in it
#3=filtered for total reads
#4=median normalization of dataset
#5=filtered out genes with median lower than arbitrary number
