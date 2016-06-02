#!/usr/bin/env Rscript


################################################################################
##                      
##                             Rscript
##                            
##----------------------- Project: MatrixMaker.R -------------------------------
##
##
## Purpose: This is a RScript that intends to replace the whole CodAnalyze 
##          written in 2007 which job was to build gene models that are 
##          subsequently used by our own software AMIGene. Among improvements, 
##          developements concern the full automation of the process and in 
##          particular, the classification step which in principle requires 
#           human expertise.
##
## Author: Stéphane Cruveiller
## Bug reports: scruveil@genoscope.cns.fr
## First version: 0.1.0 (Aug 5, 2013)
## Current version: 1.0.0 (Apr 15, 2015)
##
##
## Required R libraries
## clusterSim (>=0.42)
## cluster (>=1.14)
## Rscript (>=2.13)
################################################################################ 

## ChangeLog
## 05/08/2013: Initial version
## 06/10/2014: - Bug fixed in the way the discarded codons are handled 
##             - Now the script is wrapped by MMaker.py
## 15/15/2015: - Bug fixed when trying to perform PCA with 0 significany component

## TODO: list of Sid to be tested 
##       - 244 Escherichia coli (3 classes)
##       - 353 Mycobacterium tuberculosis (4 classes)
##       - 26  Borrelia burgdorferii (2 classes)
##       - 1943 Leptospira borgpetersenii (??)


## Get rid of annoying warnings and error messages
suppressPackageStartupMessages(library(getopt))

## Functions for classification analysis

##---------------------------------------------- FUNCTION: AKmeans -------------
##
## Purpose: This function should perform kmeans iteratively (from 2 to 10) to 
##          classify genes into various classes according to their codon usage.
##
## Parameters: x (DFrame) -> A data frame containing the results of the a 
##             factorial correspondence analysis based of codon usage (RSCU, 
##             codon counts)
##             centers (INT) -> the number of centers (i.e. centroids) to test
##             method (STR) -> the distance that should be used to perform 
##             classification (default="euclidean", other values="correlation,
##             manhattan,...")
##             iter.max (INT) -> the maximum number of iterations allowed to 
##             reach convergence
##
## Returns: KM (DFrame) -> A data frame that contains one of the optimal 
##          solutions characterized by intra and inter class variances
##------------------------------------------------------------------------------

## TODO add all the available methods

AKmeans<-function (x, centers, method = "euclidean", iter.max = 200) 
{
	## Load the appropriate packages
	require(amap)
	
	## Convert x to data frame if required
	if (mode(x) == "numeric") 
		x <- data.frame(new.x = x)
	
	## Initiating the process with a first 2 clusters analysis
	KM <- Kmeans(x = x, centers = 2, method = method ,iter.max = iter.max, nstart = 50)
	
	## Looping over ncenters
	for (i in 2:centers) {
		newKM <- Kmeans(x = x, centers = centers, method = method, iter.max = iter.max, nstart = 30)
		cat(paste("Iteration: ",i," -> newSSE: ",sum(newKM$withinss)," vs oldSSE: ",sum(KM$withinss),"\n",sep=""))
		
		## Checking if intraclass distances are smaller than those of the previous analysis 
		if (sum(newKM$withinss) < sum(KM$withinss)) {
			KM <- newKM
		}
	}
	KM$tot.withinss <- sum(KM$withinss)
	xmean <- apply(x, 2, mean)
	centers <- rbind(KM$centers, xmean)
	bss1 <- as.matrix(dist(centers)^2)
	KM$betweenss <- sum(as.vector(bss1[nrow(bss1), ]) * c(KM$size, 0))
	return(KM)
}

##---------------------------------------------- FUNCTION: Evaclust ------------
##
## Purpose: This function should perform kmeans iteratively (from 2 to 9) to classify
##          genes into various classes according to their codon usage. The cluster
##          optimality is based on Davies-Bouldin's index (TODO ref to cite)
##
## Parameters: x (DFrame) -> A data frame containing the results of the a factorial
##             correspondence analysis based of codon usage (RSCU, codon counts)
##             centers (INT) -> the number of centers (i.e. centroids) to test
##             method (STR) -> the distance that should be used to perform 
##             classification (default="euclidean", other values="correlation,manhattan,...")
##             iter.max (INT) -> the maximum number of iterations allowed to reach 
##             convergence
##
## Returns: KM (DFrame) -> A data frame that contains one of the optimal solutions
##          characterized by intra and inter class variances
##------------------------------------------------------------------------------

## TODO add all the available methods

Evaclust<-function (x, centers, method = "euclidean", iter.max = 200) 
{
	## Load the appropriate packages
	suppressPackageStartupMessages(require(clv)) ## required for most internal measures used here
	suppressPackageStartupMessages(require(clusterSim)) ## See values agreement
	
	## Constants definitions
	min_nc<-2
	max_nc<-centers
	
	## Convert x to data frame if required
	if (mode(x) == "numeric") 
		x <- data.frame(new.x = x)
	
	## Create a DFrame that will store BDIs for each iteration
	BDI<-as.data.frame(array(0, c(max_nc-min_nc+1,9)))
	names(BDI)<-c("ncenters","SSE","ASW","BDi","Connectivity","Dunn","CHi","Ci","SDi")
	BDI[,"ncenters"] <- min_nc:max_nc
	
	## Looping over ncenters
	for (nc in min_nc:max_nc) {
		## Performs Kmeans using provided parmaters
		KM<-kmeans(x = x, centers = nc,iter.max = iter.max, nstart = 25)
		## Compute SSE
		BDI[nc-min_nc+1,"SSE"]<- round(sum(KM$withinss),digits=3)
		
		## Obtaining average silhouette width directly from kmeans
		## clustering results: "One way to choose k ‘appropriately’ is to select
		## that value of k for which S(k) is a large as possible." from Rousseeuw
		## Journal of Computational and Applied Mathematics 20(1987): 53-65. 
		BDI[nc-min_nc+1,"ASW"]<- round(index.S(dist(x),KM$cluster),digits=3)
		
		
		## Prepare computation of internal measures
		## 1) Measures that need an object generated by cls.scatt.data 
		## (i.e. Davies-Bouldin, Connectivity, Dunn Index) from clv package
		Imes<-cls.scatt.data(data = x, clust = KM$cluster, dist = method)
		
		## 2) Measures that need an object generated by clv.Scatt (i.e. SD index
		## and SDbw index). Those measures rely on 2 further functions
		## clv.Dis and clv.DensBW
		SDmes<-clv.Scatt(data = x, clust = KM$cluster, dist = method)
		
		## Compute Davies-Bouldin's Index
		## According to Docs found on the web Davies-Bouldin's Index 
		## has a value between 0 and infinity and SHOULD BE MINIMIZED!
		BDI[nc-min_nc+1,"BDi"]<-round(1/index.DB(x,KM$cluster,centrotypes="centroids",p=2,q=1)$DB,digits=3)
		
		## Implementing Connectivity Index
		## According to Rdoc from clValid package connectivity 
		## has a value between 0 and infinity and SHOULD BE MINIMIZED!
		BDI[nc-min_nc+1,"Connectivity"]<-round(1/connectivity(data = x, clust= KM$cluster,neighbour.num=10,dist=method),digits=3)
		
		## Implementing Dunn Index
		## According to Rdoc from clValid package connectivity 
		## has a value between 0 and infinity and SHOULD BE MAXIMIZED!
		BDI[nc-min_nc+1,"Dunn"]<-round(clv.Dunn(Imes,intracls="average", intercls="average"),digits=3)
		
		## Implementing Calinski-Harabasz pseudo F-statistic
		## CHi has a value between 0 and infinity and SHOULD BE MAXIMIZED!
		BDI[nc-min_nc+1,"CHi"]<- round(index.G1(x,KM$cluster,centrotypes="centroids"),digits=3)
		
		## Implementing Cindex from Hubert and Levin (1976)
		## Ci has a value between 0 and infinity and SHOULD BE MINIMIZED!
		BDI[nc-min_nc+1,"Ci"]<- round(1/index.G3(dist(x),KM$cluster),digits=3)
		
		## Implementing BH index from Baker and Hubert (1975)
		## BHi has a value between -1 and 1 and SHOULD BE MAXIMIZED!
		## TODO: problem here can nit find values between -1 and 1
		## BDI[nc-min_nc+1,"BHi"]<- round(index.G2(dist.GDM(x),KM$cluster),digits=3)
		
		## Implementing SD and SDbw
		## Both values SHOULD BE MINIMIZED!
		## See for instance: 
		## http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2848&rep=rep1&type=pdf
		SDdis<-clv.Dis(SDmes$cluster.center)
#		SDdens<-clv.DensBw(data = x, clust = KM$cluster,scatt.obj=SDmes)
		BDI[nc-min_nc+1,"SDi"]<-round(1/clv.SD(SDmes$Scatt,SDdis,alfa = nc),digits=3)
#		BDI[nc-min_nc+1,"SDbw"]<-round(clv.SDbw(SDmes$Scatt,SDdens),digits=3)
		
		## Print for debug purposes
		##print(paste("Indexes for",nc,"clusters: SSE->",BDI[nc-min_nc+1,"SSE"],"Sil->",BDI[nc-min_nc+1,"Sil"],"BDIce->",BDI[nc-min_nc+1,"BDIce"],"BDIav->",BDI[nc-min_nc+1,"BDIav"],"Connectivity->",BDI[nc-min_nc+1,"Connectivity"],"Dunn->",BDI[nc-min_nc+1,"Dunn"],"SD->",BDI[nc-min_nc+1,"SD"],"SDbw->",BDI[nc-min_nc+1,"SDbw"]))
				
	}
	
	## Display results
	##print(paste("min DB for",(min_nc:max_nc)[which.min(BDI[,"BDi"])],"clusters=",min(BDI[,"BDi"])))
	##print(paste("min DB for",(min_nc:max_nc)[which.min(BDI[,"Connectivity"])],"clusters=",min(BDI[,"Connectivity"])))
	##plot(min_nc:max_nc,BDI$SSE, type="p", pch=0, xlab="Number of clusters", ylab="DB")
	##x11()
	##plot(min_nc:max_nc,BDI$BDi, type="p", pch=0,col="red", xlab="Number of clusters", ylab="Validity Index",yaxt="n",ylim=c(0,3))
	##axis(2,at=pretty(range(0,3)))
	##lines(min_nc:max_nc,BDI$BDi,col="red")
	##points(min_nc:max_nc,BDI$Dunn, type="p", pch=0, col="blue")
	##lines(min_nc:max_nc,BDI$Dunn,col="blue")
	
	return(BDI)
}

##---------------------------------------------- FUNCTION: doPCA ---------------
##
## Purpose: This function should perform PCA on codon usage matrix/array 
##          obtained after normalizing the data in some ways. The goal here is
##			to reduce the dimensionality of the data.
##
## Parameters: x (DF) -> A DF/Matrix/Array of codon uasge values.
##             
##
## Returns: x.pca (LST) -> An list that contains coordinates of individuals and
##          the number of significant components (i.e. >5%)
##
##------------------------------------------------------------------------------
doPCA<-function(x){
	
	## Load the appropriate packages
	suppressPackageStartupMessages(require(FactoMineR)) ## required for PCA 
	
	## perform PCA on codon usage and retain all components
	## et the begining
	cat(paste("[INFO] Performing PCA on the whole dataset: ",deparse(substitute(x)),"\n",sep=""))
	x.pca<-PCA(x,ncp=length(colnames(x)),scale.unit=FALSE,graph=FALSE)
	
	##	x.pca
	##	**Results for the Principal Component Analysis (PCA)**
	##			The analysis was performed on 709 individuals, described by 57 variables
	##	*The results are available in the following objects:
	##			
	##			name               description                          
	##	1  "$eig"             "eigenvalues"                        
	##	2  "$var"             "results for the variables"          
	##	3  "$var$coord"       "coord. for the variables"           
	##	4  "$var$cor"         "correlations variables - dimensions"
	##	5  "$var$cos2"        "cos2 for the variables"             
	##	6  "$var$contrib"     "contributions of the variables"     
	##	7  "$ind"             "results for the individuals"        
	##	8  "$ind$coord"       "coord. for the individuals"         
	##	9  "$ind$cos2"        "cos2 for the individuals"           
	##	10 "$ind$contrib"     "contributions of the individuals"   
	##	11 "$call"            "summary statistics"                 
	##	12 "$call$centre"     "mean of the variables"              
	##	13 "$call$ecart.type" "standard error of the variables"    
	##	14 "$call$row.w"      "weights for the individuals"        
	##	15 "$call$col.w"      "weights for the variables"
	
	## Here we need 1 info
	## the number of significant components is deduced from x.pca$eig
	## significant components = higher than 2% of contribution
	sncp <- nrow(x.pca$eig[which(x.pca$eig[,2] >= 5),])
	cat(paste("[INFO] Number of significant components (i.e. > 5%): ",sncp,"\n",sep=""))
	
	## reperform PCA using sncp as number of components if nb of components
	## is strictly positive
	if (sncp < 1) {
		stop("Unable to perform PCA on reduced data!", call.=FALSE)
	}
	else{
		cat(paste("[INFO] Performing PCA on the whole dataset: ",deparse(substitute(x))," with ",sncp," components\n",sep=""))
		x.pca<-PCA(x,ncp = sncp,scale.unit=FALSE,graph=FALSE)
		
		## Return the final object (individual coordinates)
		return(x.pca$ind$coord)
	}
	

}

##---------------------------------------------- FUNCTION: ClustScore ----------
##
## Purpose: This function should compute the score of each proposed solutions 
##          and return the best configuration.
##
## Parameters: x (Dframe) -> A DataFrame/Matrix containing cluster validity 
##             indices.
##
##             
##
## Returns: bestconf (INT) -> The number of clusters for the best partition.
##------------------------------------------------------------------------------
ClustScore<-function(x){
	
	## Display infos
	cat(paste("[INFO] Looking for best partition...\n",sep=""))
	
	## Scores Indices with equal weights i.e. based on ranks
	## so apply rank function by column and then sum by row
	x.ranks<-cbind(x$ncenters,apply(apply(x[, -c(1,2)],2,rank),1,sum))
	
	## Now looking for the best configuration i.e. highest score
	bestconf<-x.ranks[which.max(x.ranks[,2])]
	
	## Display infos
	cat(paste("[INFO] Best partition: ",bestconf," clusters.\n",sep=""))
	
	## Returning the best configuration!
	return(bestconf)
}


##---------------------------------------------- FUNCTION: LoadFile ------------
##
## Purpose: This function should load primary data for the analysis (i.e. codon
##          usage values.
##
## Parameters: x (STR) -> A string corresponding to the path of the file to 
##             load.
##             exclude (STR) -> A string containing the codons to exclude from
##             the analysis
##             
##
## Returns: cua (Matrix) -> A Matrix that contains input data to process.
##------------------------------------------------------------------------------
LoadFile<-function(x,exclude=NULL,normalization="none"){
	
	## Loading table into object
	cat(paste("[INFO] Loading input data (and normalizing if needed!): ",x,"\n",sep=""))
	cua<-read.delim(x,dec=".",header=TRUE,row.names=1)
	# in principle we should get a matrix with nrows and 64 columns!!!
	# If this is not the case, the process stops immediately...
	if (ncol(cua) != 64) {
		## arg call.=FALSE avoids mentioning the function call in the 
		## error stack...
		stop("Unable to load Codon Usage Data properly!", call.=FALSE)
	}
	
	cat(paste("[INFO] ",nrow(cua)," raw objects loaded!\n",sep=""))
	
	## Modifying Data such as to discard unwanted codons
	## cat(paste("[INFO] Excluding codons from the analysis: ",exclude,"\n",sep=""))
	## cua<-cua[,!(names(cua) %in% unlist(strsplit(exclude,',')))]
	
	## convert object to matrix to be able to coerce it to double
	## type needed by functions included in the clv package
	cua.mat<-as.matrix(cua)
	
	## Still here we have a matrix containing integers
	## So coercing it to double
	storage.mode(cua.mat)<-"double"
	if(normalization =="fmax" || normalization =="rscu"){
		
		##Defining Codons list by amino acids
		## Cys, Met, X and Ter codons are automatically discarded
		## Sextets Amino Acids
		Arg<-c("CGA","CGC","CGT","CGG","AGA","AGG")
		Leu<-c("CTA","CTC","CTG","CTT","TTA","TTG")
		Ser<-c("TCA","TCC","TCG","TCT","AGC","AGT")
		
		## Quartets Amino Acids
		Thr<-c("ACA","ACC","ACG","ACT")
		Pro<-c("CCA","CCC","CCG","CCT")
		Ala<-c("GCA","GCC","GCG","GCT")
		Gly<-c("GGA","GGC","GGG","GGT")
		Val<-c("GTA","GTC","GTG","GTT")
		
		## Duets Amino Acids
		Lys<-c("AAA","AAG")
		Asn<-c("AAC","AAT")
		Gln<-c("CAA","CAG")
		His<-c("CAC","CAT")
		Glu<-c("GAA","GAG")
		Asp<-c("GAC","GAT")
		Tyr<-c("TAC","TAT")
		Phe<-c("TTC","TTT")
		
		## Odd Amin Acids
		Ile<-c("ATA","ATC","ATT")
		
		## Create a list of arguments for loop
		## In principle "deparse(substitute(X))" should return "X" 
		## but not evaluated (i.e. STRG format) 
		AAs<- c(deparse(substitute(Arg)),deparse(substitute(Leu)),deparse(substitute(Ser)),deparse(substitute(Thr)),deparse(substitute(Pro)),deparse(substitute(Ala)),deparse(substitute(Gly)),deparse(substitute(Val)),deparse(substitute(Lys)),deparse(substitute(Asn)),deparse(substitute(Gln)),deparse(substitute(His)),deparse(substitute(Glu)),deparse(substitute(Asp)),deparse(substitute(Tyr)),deparse(substitute(Phe)),deparse(substitute(Ile)))
		
		## Create a NULL matrix that will store final results
		## (i.e. normalized counts of codons...
		cua.matnorm<-NULL
		
		## Loop processing args...
		## Looping over AAs list and implicitly over AAs values ;-)
		for( i in AAs){
			argument<-eval(parse(text=i))
			
			if (normalization=="fmax"){
			## This normalization scheme comes from the following paper. In it they the
			## authors reported that this normalization was less sensitive to biases and
			## thus much powerful when used in PCA.  
			## Suzuki H, Saito R, Tomita M. A problem in multivariate analysis of codon
			## usage data and a possible solution. FEBS Lett. 2005 Nov 
			## 21;579(28):6499-504. 
			## Epub 2005 Nov 2. PubMed PMID: 16289058.
				matnorm.tmp<-round(cua.mat[,(select=argument)]/apply(cua.mat[,(select=argument)],1,max),digits=3)
			}
			if(normalization=="rscu"){
				matnorm.tmp<-round(cua.mat[,(select=argument)]/((1/length(argument))*apply(cua.mat[,(select=argument)],1,sum)),digits=3)
			}
			## Building up final Matrix....
			cua.matnorm<-cbind(cua.matnorm,matnorm.tmp)
		}
	}
	
	## Return Input Matrix but...
	## Modifying Data such as to discard unwanted codons
	cat(paste("[INFO] Excluding codons from the analysis: ",exclude,"\n",sep=""))
	if (normalization=="none"){
		## In principle here all values are finite (i.e. INT -ge 0)

		cua.mat<-cua.mat[,!(colnames(cua.mat) %in% unlist(strsplit(exclude,',')))]
		return(cua.mat)
	} else {
		## Here some values may be undefined due to ratio computations of codons
		## that are absent from genes (especially true for short genes!!!). The 
		## following trick will get rid of lines containing at least one 
		## undefined value
		cua.matnorm<-cua.matnorm[,!(colnames(cua.matnorm) %in% unlist(strsplit(exclude,',')))]
		if (nrow(cua.matnorm[complete.cases(cua.matnorm), ]) == 0) {
			stop("Unable to normalize Codon Usage Data properly!",call.=FALSE)
		}
		else{
		cat(paste("[INFO] ",nrow(cua.matnorm[complete.cases(cua.matnorm), ])," objects successfully normalized!\n",sep=""))
		return(cua.matnorm[complete.cases(cua.matnorm), ])
		}
	}
}



##---------------------------------------------- FUNCTION: Get_GOsClass --------
##
## Purpose: This function should return the list of object in each class 
##          determined by the clustering analysis.
##
## Parameters: x (MAT) -> An object of Matrix/DataFrame containing data to cluster.
##			   centers (INT) -> the number of centers to use 
##             
##
## Returns: Nothing (i.e. files are directly written to disk...)
##------------------------------------------------------------------------------
Get_GOsClass<-function(x,nc,iter.max = 200,basename) {

	## throwing a little bit of infos...
	cat(paste("[INFO] Performing definitive partition with ",nc," centers.\n",sep=""))
	
	KMd<-kmeans(x = x, centers = nc,iter.max = iter.max, nstart = 25)
	clusters<-as.data.frame(transform(KMd$cluster))
	for (i in 1:nc) {
		golist<-as.data.frame(rownames(clusters)[which(clusters[,1] == i)])
		names(golist)[1]<-"label"
		if (is.null(basename)){
			filename<-paste("Class_",i,".lst",sep="")
		}
		else {
			filename<-paste(basename,"_Class_",i,".lst",sep="")
		}
		cat(paste("[INFO] Writing items list for class ",i," to ", filename,".\n",sep=""))
		write.table(golist,filename,eol="\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
}




##--------------------------------------------------- MAIN ROUTINE -------------

## Prepare matrix of options to be checked by getopt()
opt <- getopt(matrix(c(
						'cua' , 'c', 1, 'character', 'Codon Usage file (RSCU or Codons counts)',
						'distance','d',1,'character','Distance metric to be used to perform clustering analysis',
						'exclude','e',1,'character','Codons to discard from analysis (comma separated values; Default: "TGC,TGT,AGA,AGG,X,TAA,TAG,TGA"',
						'normed' , 'n', 1, 'character', 'Normalization method; Available options: euclidean,correlation,none; default: none',
						'prefix','p',1,'character','Prefix for files base name',
						'help', 'h', 0, 'logical', 'display help'
				),ncol = 5, byrow = TRUE));

## TODO specify fieldsep?
## TODO specify if headers?

########################## STEP 1 ######################
############# Parsing Command line Arguments ###########
########################################################

## Testing if input exists and how to store it in R
if ( is.null(opt$cua) ){
	stop("[ERROR] One of required file is missing on the command line!")
} else {
	## if no excluded codons list is provided then exclude defaults ones
	if (is.null(opt$exclude)) {
		opt$exclude<-"TGC,TGT,X,TAA,TAG,TGA,TTT"
	}
	
	#Testing validity of normalization method
	if (is.null(opt$normed)) {
		opt$normed<-"none"
	} else {
		NORM<- c("rscu","fmax","none")
		norm <- pmatch(opt$normed, NORM)
		if (is.na(norm)) stop("[ERROR] Invalid normalization method invoked!")
		if (norm == -1) stop("[ERROR] Amiguous normalization method invoked!")
		}
	cua.mat<-LoadFile(opt$cua,exclude=opt$exclude,normalization=opt$normed)
	cat(paste("[INFO] Normalization method: ",opt$normed,"\n",sep=""))
}

## If distance is not specified set it back to default (i.e. euclidean)
if (is.null(opt$distance)){
	opt$distance<-"euclidean"
} else {
	DIST<- c("euclidean","correlation")
	dista <- pmatch(opt$distance, DIST)
	if (is.na(dista)) stop("[ERROR] Invalid distance invoked!")
	if (dista == -1) stop("[ERROR] Amiguous distance invoked!")
}
## Displaying Metric used...
cat(paste("[INFO] Metric used: ",opt$distance,"\n",sep=""))

## First step reducing matrix dimension using PCA
cua.mat.pca<-doPCA(cua.mat)

## Compute and Evaluate quality of clustering results
cua.mat.eva<-Evaclust(cua.mat.pca,6,method=opt$distance)

## Third step computing score for the best partition
nclusts<-ClustScore(cua.mat.eva)

## Fourth step: Generating the list of GOs belonging to the nclusts classes
Get_GOsClass(cua.mat.pca,nclusts,iter.max = 200,basename=opt$prefix)
cat("[INFO] Codon Usage Analysis done.\n")




