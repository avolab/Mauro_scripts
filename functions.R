#Functions needed for reading, filtering and normalizing Cel-SEQ data

#remove spike-ins from given data frame
rmspike<-function(x){
  ERCCs<-grep("ERCC-",row.names(x)) # gives vector with row # of spike ins
  data<-x[-ERCCs,] # make new data frame without the specified rows
  return(data) # output new data frame
}

#keep spike-ins from given data frame
keepspike<-function(x){
  ERCCs<-grep("ERCC-",row.names(x)) 
  data<-x[ERCCs,] # make new data frame with only the specified rows
  return(data)
}

# plot oversequencing for combination of datasets,
# x = read count file, y = barcode count file, name = name of data frame
# if pdf= TRUE, will save files in subdirectory "plots" of the current working directory,otherwise will plot to screen
overseq <- function(x,y,name,pdf=FALSE){
  main=paste("oversequencing",name)   # mixes string + name of choice
  xlab=bquote(log[2] ~ "read counts / barcode counts")    # mixes subscript + string
    if(pdf){ # if pdf=TRUE, make pdf
        pdf(paste(getwd(),"/plots/",main,".pdf",sep="")) # make in subdir "/plots" of current working dir, using main for filename
    hist(log2(apply(x,1,sum)/apply(y,1,sum)),breaks=75, col="lightblue", main=main,xlab=xlab)
    dev.off()
  }
    else{hist(log2(apply(x,1,sum)/apply(y,1,sum)),breaks=75, col="lightblue", main=main,xlab=xlab)}
}

# make GENEID the rownames and remove column GENEID
nogeneid<-function(x){
    data <- as.data.frame(x)
    rownames(data) = data[,1]
    data= data[,-1]
    return(data)
}

# construct a data object with all data frames in it, feed it two lists: 
# one with names (of experiments/timepoints etc) and one with data frames
# these will be rearragned to form a 2-dimensional list, where each element contains 1 experiment,
# row 1 = name of the experiment (will be used to lable plots), row 2 = corresponding data frame
# subsequent rows will be used to add normalized/filtered/etc data frames to the object
makeobject<-function(x,y){
  obj <- list()
  if(length(x) != length(y)) stop("lengths of lists do not match")
  
  for(i in 1:length(x)){
    obj[[i]]<-list()
    obj[[i]]<-list(x[i],as.data.frame(y[i]))
  }
  return(obj)
}

#function to combine plots 
combineplots<-function(a,b,square=FALSE){
  par(mfrow=c(a,b)) # will combine plots in grid with a rows and b columns
  if(square){par(pty="s")} # if specified as TRUE, the plots will be made square
}

#plot total number of reads per sample, x = data object, i'th column, j'th row 
#pdf=TRUE to store plot in pdf in /plots directory of current directory
totalreads <- function(x,i,j,pdf=FALSE){
    main=paste("total reads",x[[i]][[1]]) # get name from ith column and first row
  if(pdf){
      pdf(paste(getwd(),"/plots/",main,".pdf",sep=""))
      b<-barplot(colSums(x[[i]][[j]]),xaxt="n",xlab="samples",main=main,col="lightgrey") 
      #no x-axis lables, store barplot values in vector for x-axis tick marks
     axis(1,at=b,labels=c(1:length(x[[i]][[j]]))) # 1=horizontal at=position of tickmarks
    dev.off()
   } 
  else{
    b<-barplot(colSums(x[[i]][[j]]),xaxt="n",xlab="samples",main=main,col="lightgrey") 
    axis(1,at=b,labels=c(1:length(x[[i]][[j]]))) # 1=horizontal at=position of marks
  }
}

#filter all samples that have lower than n total reads
#x=dataframe, i=ith forloop,j=row of object to be used,n=lower transcript bound,m=higher transcript bound)
filtersamples<-function(x,i,j,n,m,pdf=FALSE){
  f <- apply(x[[i]][[j]],2,sum,na.rm=TRUE) > n
  filt<- x[[i]][[j]][,f]
  f<- apply(filt,2,sum,na.rm=TRUE) < m
  filt<-filt[,f]
  main=paste("samples left post filter",x[[i]][[1]])
  xlab=paste(n,"> keep <",m,sep="")
  sub=paste("samples left:",ncol(filt),",removed:",ncol(x[[i]][[j]])-ncol(filt))
  if (pdf){
    pdf(paste(getwd(),"/plots/",main,".pdf",sep=""))
    barplot(colSums(filt), col="green",main=main, xaxt ="n",xlab=xlab, ylab="total reads",sub=sub)
    dev.off()
    return(filt)
  }
  else{
    barplot(colSums(filt), col="green",main=main, xaxt = "n",xlab=xlab, ylab="total reads",sub=sub)
    return(filt)
  }
}

#Create a list of ERCCs found in data file + the number of molecules spiked in (formula from lennart)
conversion <-function(x,i,j,pdf=FALSE){ 
  ERCCs = x[[i]][[j]][grep("ERCC-",row.names(x[[i]][[j]])),] # extract detected ERCCs in samples
  ERCC.list <- merge(ERCC_molecules,ERCCs,by="row.names",all=TRUE) # make table with spike and detected ERCC
  row.names(ERCC.list) <- ERCC.list[,1]
  ERCC.list<- ERCC.list[,-1]
  ERCC.list[is.na(ERCC.list)] <- 0
  conv.list <- lapply(1:ncol(ERCC.list), function(i) { coef(lm(ERCC.list[,1] ~ ERCC.list[,i]))[2]})
  conv.list <- as.vector(unlist(conv.list))
  conv.list <- conv.list[-1]
  conv.list <- 1/conv.list
  main=paste("conversion factor",x[[i]][[1]]) # get name from ith column and first row
  sub=paste("mean conversion factor in %:", mean(conv.list)*100)
  if (pdf){
    pdf(paste(getwd(),"/plots/",main,".pdf",sep=""))  
    b<-barplot(conv.list, col = 'orange', main = main, xlab = "samples", sub=sub,ylab = "conversion factor")
    axis(1,at=b,labels=c(1:length(x[[i]][[j]])))
    dev.off()
  }
  else{
    b<-barplot(conv.list, col = 'orange', main = main, xlab = "samples", sub=sub,ylab = "conversion factor")
    axis(1,at=b,labels=c(1:length(x[[i]][[j]])))
    return(conv.list)
  }
}

#normalization of objects in list of lists to the median
norm.median <-function(x,i,j){
  data<-x[[i]][[j]]
  s <-colSums(data,na.rm=TRUE)
  m <- mean(s)
  norm <-data
  for (l in 1:ncol(data)) {norm[,l] <-data[,l]/s[l]*m
  }
  return(norm)
}

#normalization to median of one dataframe
norm.median.object <-function(x){
  s <-colSums(x,na.rm=TRUE)
  m <- mean(s)
  norm <-x
  for (i in 1:ncol(x)) {norm[,i] <-x[,i]/s[i]*m
  }
  return(norm)
}

#plot number of available transcripts vs cutoffs of median detected transcripts
testcutoff<-function(x,i,j,pdf=FALSE){
  data<-x[[i]][[j]]
  main=paste("genes cutoff test",x[[i]][[1]])
  for(l in 1:20){
    z = apply(data,1,median) > l
    if(l==1){
      rc.cutoff = z
    } else {
       rc.cutoff = cbind(rc.cutoff,z)
    }
  }
  if (pdf){
    pdf(paste(getwd(),"/plots/",main,".pdf",sep="")) 
    plot(apply(rc.cutoff,2,sum),ylab = "number of transcripts",col="black",
    xlab = "cutoff (mean transcript no.)",main=main,type="b",lty=2,pch=19)
    dev.off()
  }
  else{
    plot(apply(rc.cutoff,2,sum),ylab = "number of transcripts",col="black",
    xlab = "cutoff (mean transcript no.)",main=main,type="b",lty=2,pch=19)
  }    
}
    
# filter out genes with median expression lower than X transcripts per cell
# x = data frame, ith column, jth row, n=median # of transcripts
filtergenes <- function(x,i,j,n) {
  data <- x[[i]][[j]]
  f <- apply(data,1,median,na.rm=TRUE) > n # boolean vector for which TRUE: gene has median > n over all samples
  data.filter <- data[f,] # make new data frame for all genes with median > n
  name=x[[i]][[1]][[1]] # assign this to name of dataset in question
  cat("removed", nrow(data)-nrow(data.filter), "genes in",name, sep = " ", fill = TRUE, labels = NULL,
    append = FALSE) # outputs how many genes were removed due to filtering
  return(data.filter)
}

#plot levels of a gene in ETOH and 4OHT samples
#x=dataset,ith column,jth row,n=name of gene
plotgenes<-function(x,i,j,n){
  data<-x[[i]][[j]]
  wt<-data[,grep("HCT_CEL",names(data))]
  ko<-data[,grep("DICER_CEL",names(data))]
  
  name<-grep(n,rownames(data))
  w<-(as.matrix(wt[name,]))
  k<-(as.matrix(ko[name,]))
  par(mfrow=c(1,2))
  boxplot(as.vector(w),main=n,xlab="WT")
  boxplot(as.vector(k),main=n,xlab="KO")
}


# Merge all dataframes of 1 row from object into one
mergeall<-function(x,y){
  new<-merge(x,y, by="row.names", all=TRUE)
  rownames(new)<-new[,1]
  new<-new[,-1]
  return(new)
}

# find coordinates (column numbers) of all samples tagged with experiment tag x in merged dataframe y, ith for loop
#feed this function a vector called conditions, containing names of for example timepoints(or other experimental tags) to be indexed
#will result in list of length(conditions), with each element containing the column numbers of 1 certain condition + experiment(x)
samplecoor<-function(x,y,i){
  a<-grep(x,colnames(y))
  b<-grep(conditions[i],colnames(y))
  c<-a%in%b
  a<-a[c]
  return(a)
}

# correlate samples based on expression of all genes
# x=dataframe, samplist=created in function samplecoor, gives col numbers of right conditions, conditions=list of condition names
# i = ith for loop
cor.samples<-function(x,samplist,conditions,i){
  data<-x[,samplist[[i]]]
  cor.ko.cel <-cor(x[,grep("CEL_O", names(data))])
  cor.wt.cel <-cor(x[,grep("CEL_E", names(data))])
  cor.ko.ctr <-cor(x[,grep("CTR_O", names(data))])
  cor.wt.ctr <-cor(x[,grep("CTR_E", names(data))])
  hist(as.matrix(cor.ko.cel), col="blue", breaks =50, freq=TRUE, main=paste("correlation between cells",conditions[[i]]),xlab="Correlation") 
  hist(as.matrix(cor.wt.cel), add=T, col=rgb(0,1,0,0.8), breaks=50, freq=TRUE)
  hist(as.matrix(cor.ko.ctr), col="blue", breaks =50, freq=TRUE, main=paste("correlation between controls",conditions[[i]]),xlab="Correlation") 
  hist(as.matrix(cor.wt.ctr), add=T, col=rgb(0,1,0,0.8), breaks=50, freq=TRUE)

  legend( "topleft"
  , cex = 1.5, 
  , bty = "n", 
  , legend = c("Wildtype", "Dicer -/-"), 
  , text.col = c("Black", "black"),
  , col = c("Green", "Blue"), 
  , pt.bg = c("black","black")
  , pch = c(20,20)
  )
}


