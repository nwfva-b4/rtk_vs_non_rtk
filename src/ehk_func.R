#--------------------------------------------------------------------------
# Name:         ehk_func.R
# Author:       Northwest German Forest Research Institute
#               Department Forest Growth
#--------------------------------------------------------------------------

# Function for the Single Tree-Height Imputation using Single Tree Height Curves

# 1. Generate the input table in the required format

# 2. Input variables for height imputation
# 	dm: dg for all trees
# 	ds: dg of trees with height measurements
# 	hs: hg of trees with height measurements
# 	hm: hg all trees (using inverted functions, vgl. Formula 5.2.4.9 BWI-Methodenband 2012)
# 	Calculation sequence for different layers with height data: Plot/Layer/LN/Bagr
# --> If no elevation data is available for the Plot/Layer/LN/Bagr layer, Plot/Layer/LN is used, and so on
#     Global model if no height data is available for the Plot
#_________________________________________________________________________________
# Creates input data in the required format
# id = Unique Plot ID for BWI data, generated from Tnr and Enr
# bagr = Tree species group with coefficeints for EHK (Stand BWI3)
# DBH in cm, Height in m

input_ehk=function(id, bnr, bs, bhd, hoe, nha, bagr){
	
	ein = data.frame(id=id, bnr=bnr, bs=bs, bhd=bhd, hoe=hoe, nha=nha, bagr=bagr)

# Tree species group with coefficeints for EHK  (BWI 2012 Methodenband)
# bagr_ehk: ALN and ALH --> BU for height measurement, i.e., no separate measurements if necessary
# ba_ln: Deciduous tree/conifer. Important for understory or if height measurements are not available for all bagr_ehk entries
	
	z = data.frame(bagr=c("FI", "TA", "DGL", "KI", "LAE", "BU", "EI", "ALH", "ALN"),
							   ba_ln=c("N", "N", "N", "N", "N", "L", "L", "L", "L"),
								 bagr_ehk=c("FI", "TA", "DGL", "KI", "LAE", "BU", "EI", "BU", "BU"),
								 k0=c(0.183,0.079,0.24,0.29,0.074,0.032,0.102,0.122,0.032),
								 k1=c(5.688,3.992,6.033,1.607,3.692,6.04,3.387,5.04,4.42))

	ein=merge(ein,z, by="bagr")

# Define Levels 
# This will be needed later to check which levels have elevation data.
# If elevation data is missing within a level, the next higher level is selected when filling in the gaps
# 1: Tree Species Group/LN/Story Layer/Plot 
# 2.Level: LN/Story Layer/Plot. 
# 3.Level: Story Layer/Plot. 
# 4.Level: Plot
	ein$id4 <-paste(ein$id, sep="")        					 # Plot
	ein$id3 <-paste(ein$id4, ein$bs, sep="_")        # Plot and Story Layer
	ein$id2 <-paste(ein$id3, ein$ba_ln, sep="_")     # Plot, Story Layer and Tree Species Group (LN)                                                       
	ein$id1 <-paste(ein$id2, ein$bagr_ehk, sep="_")  # Plot, Story Layer, Tree Species Group (LN) and Tree Species

	ein=ein[order(ein$id, ein$bnr), 
		c("id", "id4", "id3", "id2", "id1", "bnr", "bs", "bhd", "hoe", "nha", "bagr", "ba_ln", "bagr_ehk", "k0", "k1")]
	return(ein)
}


#_________________________________________________________________________________
# Calculates dm, hm, ds, hs per Level (see above) using tree weighting (nha)
# then individual tree heights via EHK 
# for plots without measured heights, apply the global model h=f(d,bagr,bs)

ehk=function(ein){

  y=ein
	y$g=((pi/4)*y$bhd^2)*y$nha
  
	x1=aggregate(g ~ id1, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id1, data=y, FUN=sum)$nha
	x1$dm1=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id1", all.x=T)
	
	x1=aggregate(g ~ id2, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id2, data=y, FUN=sum)$nha
	x1$dm2=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id2", all.x=T)
	
	x1=aggregate(g ~ id3, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id3, data=y, FUN=sum)$nha
	x1$dm3=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id3", all.x=T)

	x1=aggregate(g ~ id4, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id4, data=y, FUN=sum)$nha
	x1$dm4=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id4", all.x=T)

# Bookmark for printing	
	ein=y

	# dg and hg height measurement trees
	# Keep only the height measurement trees, but merge them with the overall data at the end
  y=y[!is.na(y$hoe) & y$hoe>0,] 
	y$gh=y$g*y$hoe
  
# dg
	x1=aggregate(g ~ id1, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id1, data=y, FUN=sum)$nha
	x1$ds1=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id1", all.x=T)
	
	x1=aggregate(g ~ id2, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id2, data=y, FUN=sum)$nha
	x1$ds2=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id2", all.x=T)
	
	x1=aggregate(g ~ id3, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id3, data=y, FUN=sum)$nha
	x1$ds3=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id3", all.x=T)

	x1=aggregate(g ~ id4, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id4, data=y, FUN=sum)$nha
	x1$ds4=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id4", all.x=T)

# hg
  x1=aggregate(gh ~ id1, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id1, data=y, FUN=sum)$g
  x1$hs1=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id1", all.x=T)
	
  x1=aggregate(gh ~ id2, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id2, data=y, FUN=sum)$g
  x1$hs2=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id2", all.x=T)
	
  x1=aggregate(gh ~ id3, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id3, data=y, FUN=sum)$g
  x1$hs3=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id3", all.x=T)
	
  x1=aggregate(gh ~ id4, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id4, data=y, FUN=sum)$g
  x1$hs4=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id4", all.x=T)

# Keep level with height measurement
	ein$dm=ein$dm4
	ein$dm=ifelse(!is.na(ein$ds3),ein$dm3,ein$dm)
	ein$dm=ifelse(!is.na(ein$ds2),ein$dm2,ein$dm)
	ein$dm=ifelse(!is.na(ein$ds1),ein$dm1,ein$dm)	
	
	ein$ds=ein$ds4
	ein$ds=ifelse(!is.na(ein$ds3),ein$ds3,ein$ds)
	ein$ds=ifelse(!is.na(ein$ds2),ein$ds2,ein$ds)
	ein$ds=ifelse(!is.na(ein$ds1),ein$ds1,ein$ds)

	ein$hs=ein$hs4
	ein$hs=ifelse(!is.na(ein$hs3),ein$hs3,ein$hs)
	ein$hs=ifelse(!is.na(ein$hs2),ein$hs2,ein$hs)
	ein$hs=ifelse(!is.na(ein$hs1),ein$hs1,ein$hs)
	
# Delete columns
	del=c("id1", "id2", "id3", "id4", "g", "dm1", "dm2", "dm3", "dm4", "ds1", "ds2", "ds3", "ds4", "hs1", "hs2", "hs3", "hs4")
		for(i in del)
		{
			ein[,i]=NULL
		}

# Calculate hm using the inverted EHK
	ein$hm=((ein$hs-1.3)/(exp(ein$k0*(1-(ein$dm/ein$ds)))*exp(ein$k0*((1/ein$dm)-(1/ein$ds)))))+1.3

# Height estimation for individual trees using EHK
	ein$hoe_mod = 1.3 + (ein$hm-1.3) * exp(ein$k0*(1-(ein$dm/ein$bhd))) * exp(ein$k1*((1/ein$dm)-(1/ein$bhd)))


# global model for height interpolation when no measurements are available for the plot
	l1=lm(hoe ~ log(bhd) + factor(bagr) + factor(bs), data=y) 
  ein$hoe_mod2=predict(l1, newdata=ein) 
  
  if(any(is.na(ein$hoe_mod))){
  	ein[is.na(ein$hoe_mod),]$hoe_mod=ein[is.na(ein$hoe_mod),]$hoe_mod2
  }
  ein$hoe_mod2=NULL
	
  return(ein)
}



#_________________________________________________________________________________
# Global Model for height imputation




