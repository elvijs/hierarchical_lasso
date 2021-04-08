### Functions for lasso
### the full path
### 131216 160117 270519 220321
###
### Neccesary libraries
library(partitions)
library(multicool)
library(combinat)
library(quadprog)

###
###  Generation of grids
###

## FULL generic grid, meant to be used for models
gridk <- function(p,k=2) ## Recursive grid k^p, all possibilities, all orthants
{
 if (p==1) {  
  res <- c(0:(k-1));   dim(res) <- c(k,1)
 }
 else {
  tmp <- gridk(p-1,k=k); res<-c()
   for(ii in 0:(k-1))  res <- rbind(res, cbind(ii,tmp) )
 }
 return(unname(res))
}

grid3 <- function(p) ## Recursive grid 3^p, all possibilities, all orthants
{
 if (p==1) {  
  res <- c(-1,0,1);   dim(res) <- c(3,1)
 }
 else {
  tmp <- grid3(p-1); res <- rbind(cbind(-1,tmp),cbind(0,tmp),cbind(1,tmp))   
 }
 return(res)
}

grid2 <- function(p) ## Recursive grid 2^p, all possibilities, all orthants
{
 if (p==1) {  
  res <- c(-1,1);   dim(res) <- c(2,1)
 }
 else {
  tmp <- grid2(p-1); res <- rbind(cbind(-1,tmp),cbind(1,tmp))   
 }
 return(res)
}

## imported from work on homology for Magdeburg 251118
ff<-function(k,zo=TRUE) ## generate recursively array size 2^k
{ ## zo=TRUE in levels 01, z0=FALSE in levels -11
  low<-0*zo+(-1)*!zo
  if(k==1) return( cbind(c(low,1)) )
  else return( unname(rbind( cbind(low,ff(k-1,zo)), cbind(1,  ff(k-1,zo)))  )) 
}

## Tests
# gridk(p=3,k=3) ## 27 lattice points {0,1,2}^3
# grid3(p=3) ## 27 lattice points {-1,0,1}^3
# grid2(p=4) ## 16 lattice points {-1,1}^4
# ff(k=4) ## 16 lattice points {0,1}^4

###
###  Functions for monomials <=> integer lattices <=> models
###   A model is given as a matrix with ncol=d and nrow=model size
###

## Give all monomials that divide a *single* DM
## this is a low efficiency version using gridk and removing unneeded terms
lowercorner<-function(DM){
  k<-length(DM); pmax<-max(DM); kandset<-gridk(p=k,k=pmax+1)
  res<-c()
  for(i in 1:nrow(kandset)){
    kual<-kandset[i,]; kmp<-DM>=kual;
    if( prod(kmp)==1  ) res<-rbind(res,unname(kual))
  }
  return((res))
}

## Give all monomials that divide a *list of* DM
## this is a low efficiency version using gridk and removing unneeded terms
## *better* than using the function lowercorner several times
alllowercorner<-function(DM){
  k<-ncol(DM); pmax<-max(c(DM));  kandset<-gridk(p=k,k=pmax+1)
  res<-c()
  for(i in 1:nrow(kandset)){
    kual<-kandset[i,]; 
    kual<-matrix(nrow=nrow(DM),ncol=k,byrow=TRUE,kual)
    kmp<-DM>=kual;
    if( sum(apply(kmp,1,prod))>0  ) res<-rbind(res,unname(kual))
  }
  return((res))
}

## Remove those monomials that are not dividers,
## only keep from LM quadrants given by DM   121119
limpiadivs<-function(LM,DM)
{
  res<-c()
  for(i in 1:nrow(LM))
  {
    testalo<-matrix(ncol=ncol(DM),nrow=nrow(DM),byrow=TRUE,LM[i,])
    testalo<-testalo-DM
    testalo<-apply(testalo<=0,1,prod)==1
    if(sum(testalo)>0) res<-rbind(res,LM[i,])
  }
  return(unname(res))
}

## Generate names (based on rows) of exponent matrix LM
nombres<-function(LM,separa="") t(t(apply(LM,1,paste,sep="",collapse=separa)))

## Generate names for the matrix of coefficients (output of Lasso)
nombracoeffsmatrix<-function(LM,separa="") c("lambda",nombres(LM=LM,separa=separa),"L")

## Remove intercept
clean<-function(LM) LM[nombres(LM)!=c(nombres(matrix(nrow=1,rep(0,ncol(LM)))) ),]

##### Remove duplicated terms, meant for model
klean<-function(LM) LM[!duplicated(nombres(LM=LM,separa="-")),]

## sort the list by degree
sortdeg<-function(LM)  LM[ order(apply(LM,1,sum))  ,]

## trim the model by degree
trimdeg<-function(LM,degree=1) LM[apply(LM,1,sum)<=degree,]

## take away intercept
nointer<-function(LM)
{
  kual<-nombres(matrix(rep(0,ncol(LM)),nrow=1),separa="-") ## to remove
  aliases<-nombres(LM=LM,separa="-")
  res<-LM[aliases!=c(kual),]; if(ncol(LM)==1) res<-matrix(ncol=1,res)
  return(res)
} 


istooclose<-function(res,rest,TOL,p) ## detect moves that are too close
{ 
  resultado<-c(); 
  kand<-res[nrow(res),(1:p)+1]
  for(i in 1:nrow(rest)) resultado<-c(resultado,    sum((kand-rest[i,1:p])^2)<10^(-TOL)   )
  return(resultado)
}


## to compare two vectors
vektorstooclose<-function(v1,v2,TOL)  sum((v1-v2)^2)<10^(-TOL)

## Tests
# lowercorner(DM=c(1,2))  ## model {0,1}x{0,1,2}
# alllowercorner(DM=matrix(ncol=2,byrow=TRUE, c(1,2,2,1))) ## {0,1,2}^2\{2,2} but with duplicated rows
# clean(LM = alllowercorner(DM=matrix(ncol=2,byrow=TRUE, c(1,2,2,1)))) ## without intercept but still with duplicates
# klean(LM = alllowercorner(DM=matrix(ncol=2,byrow=TRUE, c(1,2,2,1)))) ##  without duplicates (8 terms)
# limpiadivs(LM=lowercorner(DM=c(2,1)),DM=lowercorner(DM=c(1,2))) ## {0,1,2}x{0,1} \cap {0,1}x{0,1,2}   =  {0,1}^2
# nombres(LM=lowercorner(DM=c(1,2))) ## list of model term names of {0,1}x{0,1,2}
# sortdeg(LM=lowercorner(DM=c(1,2))) ## order by degree of {0,1}x{0,1,2}
# trimdeg(LM=lowercorner(DM=c(1,2)), degree =2) ## allterms of degree smaller than or equal to two of {0,1}x{0,1,2}
# nointer(LM=lowercorner(DM=c(1,2))) ## The model {0,1}x{0,1,2} without intercept

###
### Generating a multigraded model
###


## Generating terms of a given degree
## k is number of variables, d is degree
graduado<-function(k,d)
{
 t(as.matrix(restrictedparts(n=d,m=k)))->X
 rr<-c()
  for(i in 1:nrow(X)) rr<-rbind(rr,allPerm(initMC(X[i,])))
return(unname(rr))
}

## All model terms up to a given degree
## k is number of variables, d is degree
fullgraduado<-function(k,d)
{
rr<-c(rep(0,k))
for(i in 1:d) rr<-rbind(rr, graduado(k=k,d=i) )
return(unname(rr))
}

## Tests
# graduado(k=4,d=3) ## 20 terms of to degree 3 in 4 variables
# fullgraduado(k=4,d=3) ## 35 terms up to degree 3 in 4 variables

###
### Other functions
###

## Generate a random LH design with option to scale to [0,1]
LHD<-function(d,n,zeroone=FALSE,seed=0)
{ 
  set.seed(seed)
  res<-1:n
  for(i in 1:(d-1)) res<-cbind(res,order(runif(n)))
  if(zeroone) res<-(res-1)/(n-1)
  return(unname(res))
}

## Sequence of linear, with of log values (if loga=TRUE)
## generally used for exploring \lambda in lasso
seqq<-function(minimo=0,delta=0.01,maximo,cuantos,loga=TRUE)
{
  res<-seq(from=minimo,to=maximo,length.out=cuantos); res2<-c()
	if(loga) res2<-10^seq(from=log10(minimo+delta),to=log10(maximo-delta),length.out=cuantos)
  res<-sort(unique(c(res,res2)))
	return(res)
}

## Cumulative sum of list of numbers
cusum<-function(lista) ## given x1,x2,x3 ... retrieves x1,x1+x2,x1+x2+x3,...
{
  res<-lista*0
  for(i in 1:length(lista)) res[i]<-sum(lista[1:i])
  return(res)
}

## Gaps (differences) between elements of list (partial inverse of cusum)
## given x1,x2,x3,... retrieves x2-x1,x3-x2,...
gaps<-function(lista) lista[-1]-lista[-length(lista)]

## Tests
# LHD(d=3,n=4,zeroone = TRUE) ## random LHD in [0,1] with n=4 points in d=3 dimensions
# seqq(minimo=0.1,maximo=5,cuantos=5) ## ten values
# cusum(lista = c(1,4,7)) ## 1 5 12
# gaps(lista = c(1,5,12)) ##  4 7


###
### Manipulation of models
###

## This function embeds LM in a LARGE array of zeroes, in position determined by dondecols
putobject<-function(LM,dondecols=1:ncol(LM),size=ncol(LM))
{
  res<-matrix(ncol=size,nrow=nrow(LM),0); res[,dondecols]<-LM
  return(res)
}

## Turn vector c(1,0,-1) => "+0-" any other than +-1 is turned into zero
## Used with sign() to describe steps of lasso
symbol<-function(chain) 
{
res<-c(); chain<-round(chain,0)
for(i in 1:length(chain)){
  actual<-c("0")
  if(chain[i]==-1) actual<-c("-"); if(chain[i]==1) actual<-c("+");
  res<-c(res,actual)
 }
return(paste(res,collapse=""))
}

## Inverse of symbol(), turning "+0-" => c(1,0,-1)
char2chain<-function(char) 
{
res<-c()
 for(i in strsplit(x=char,split="",useBytes=TRUE)[[1]])
 {
  resc<-0; if(i=="+") resc<-1; if(i=="-") resc<--1;  res<-c(res,resc)
 }
return(res)
}

## Tests
# putobject(LM=lowercorner(DM=c(1,2)),size=4) ## {0,1}x{0,1,2} is embedded in 4 dimensions
# putobject(LM=lowercorner(DM=c(1,2)),dondecols=c(3,1),size=4) ## as above but in columns 3 and 1
# symbol(chain = c(1,-1,0)) ## "+-0"
# char2chain(char=c("00-")) ## 0 0 -1

###
### Functions  lasso
###

## For given Cm, which neighboring orthants increase dimensionality
kands<-function(Cm) ## candidates for movements up
{ ## input is string only
res<-c()
if(prod(Cm!=0)==0)
{
kuales<-(1:length(Cm))[Cm==0]
for(j in kuales){
   Nm<-Cm; Nm[j]<-1;   res<-rbind(res,Nm)
   Nm[j]<--1;   res<-rbind(res,Nm)
    }
}
return(unname(res))
}

## Pseudo inverse of CX^TXC 
## Here C is diagonal of {+-1,0} entries
pseudo<-function(X,C) 
{
XtX<-t(X)%*%X;    BM<-C%*%XtX%*%C  ## big matrix
recover<-( diag(C%*%C)==1  )
## this matrix should be invertible
LM<-BM[recover,recover]
if(sum(recover)==0) LM<-0 else LM<-solve(LM)
RM<-BM*0
RM[recover,recover]<-LM
return(RM)
}

## TRUE if knd is in list of members
membership<-function(knd,members) (prod(knd!=members)==0)

## Manipulation of matrix, remove rows/columns of MMM according to CCC diagonal of {+-1,0} 
## kolumn tells us to do it for single column matrix
## kual=1,2,3  row only, column only, both
## works for non square matrices as well, provided dimensions of C are meaningful
projdown<-function(MMM,CCC,kual=3) 
{
  p<-nrow(CCC)
  if(length(MMM)==p) MMM<-matrix(ncol=1,MMM)
  Indigo<-(1:p)[diag(CCC%*%CCC)==1]
  if(kual==1) res<-MMM[Indigo,]
  if(kual==2) res<-MMM[,Indigo]
  if(kual==3) res<-MMM[Indigo,Indigo] 
  if(length(res)==length(Indigo)) res<-matrix(ncol=1,res)
return(res)
}


## Restore size of element to p rows / both rows and columns
## kual=1,2,3  row only, column only, both
projup<-function(MMM,CCC,kual=3)
{
  p<-nrow(CCC)
  if(length(MMM)<=p) MMM<-matrix(ncol=1,MMM)
  Indigo<-(1:p)[diag(CCC%*%CCC)==1]
  if(kual==1){
    res<-matrix(ncol=1,nrow=p,0);  res[Indigo,]<-MMM
  }
#if(kual==2) res<-MMM[,Indigo] ## uncomment? 291119
  if(kual==3){
     res<-matrix(ncol=p,nrow=p,0); res[Indigo,Indigo]<-MMM
  }
#if(length(res)==pth(Indigo)) res<-matrix(ncol=1,res)
return(res)
}

## Tests
# kands(Cm=c(-1,0,1,0)) ## Four neighboring orthants to (-1,0,1,0)
# pseudo( X=matrix(ncol=3,(1:12)*1), C=diag(c(1,0,1))) ## X lives in 3 factors, but pseudo restricts to x0x
# projdown( MMM =matrix(ncol=3,(1:12)*1), CCC=diag(c(1,0,1)), kual=3) ## only retains a submatrix of size 2x2
# projup(MMM = matrix(nrow=2,1:4), CCC = diag(c(1,0,1,0)),kual=3) ## embeds matrix in a big 4x4


###
### Functions for minimisation
###

## the following needs the R library quadprog
### The minimization of L=1/2||Y-Xb||_2^2 + \lambda ||b||_1
## using the cone approach b=Cu so minimising over cone C^2u>=0
## CM defines cone, XM is design model matrix -minus intercept-
## YM is response vector
## lamm is lambda parameter
mincono<-function(CM,XM,YM,lamm)
{
  Cm<-CM; zero<-prod(diag(Cm)==0)==1
  valor<-0
if(!zero){	
  XtX<-t(XM)%*%XM
  Unovec<-matrix(ncol=1,rep(1,ncol(XM)))
  Dmat<-CM%*%XtX%*%CM
  dvec<-CM%*%t(XM)%*%YM-lamm*CM%*%CM%*%Unovec; 
  Amat<-t(CM%*%CM)
  ## project down elements
  Dmat<-projdown(CCC=Cm,MMM=Dmat,kual=3); 
  dvec<-projdown(CCC=Cm,MMM=dvec,kual=1); 
  Amat<-projdown(CCC=Cm,MMM=Amat,kual=3)
  ## apply quadprog, este es el nucleo de mincono
  AQ<-solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat)
  Uhat<-matrix(ncol=1,AQ$solution)
  ## restore the solution 
  Uhat<-projup(CCC=Cm,MMM=Uhat,kual=1)
  Bhat<-Cm%*%Uhat ## update to true quadrant
  valor<-AQ$value
}
else
{
  Bhat<-matrix(ncol=1,nrow=5,0); 
}
return(  c(diag(Cm),Bhat, valor+1/2*sum(YM^2)) ) ## not Cm!
}

betacono<-function(CM,XM,YM,lamm) ## mincono pero remueve CM del vector de salidas
{
   bb0<-mincono(CM=CM,XM=XM,YM=YM,lamm=lamm)
   bb0<-bb0[-c(1:ncol(XM))]
   return(bb0) ## regresa b^(lamda) and value of function
}

## Tests
# attach(USArrests); mincono(CM=diag(c(1,-1,1)),XM=as.matrix(USArrests[,2:4]),YM=USArrests[,1],lamm=0)
## exactly the same as below. the output includes the cone and criterion  
# lm(USArrests[,1]~as.matrix(USArrests[,2:4])-1)$coefficients
## the following removes cone, only vector and criterion
# attach(USArrests); betacono(CM=diag(c(1,-1,1)),XM=as.matrix(USArrests[,2:4]),YM=USArrests[,1],lamm=0)


## This is exhaustive minimization
## slightly edited 231019 to add GG
lassocono<-function(XM,YM,lmin=0,
     lmax=max(abs(t(XM)%*%YM)),nlambda=10,AA=NULL,lambdaval=NULL,GG=NULL)
{
  if(is.null(lambdaval)) lambdaval<-seq(from=lmin,to=lmax,length.out=nlambda)
  resulta<-c();  nlambda<-length(lambdaval)
  if(is.null(GG))  GG<-grid2(ncol(XM));
  for(i in 1:nlambda){
    BT<-minconoineqquadr(XM=XM,YM=YM,lamm=lambdaval[i],AA=AA,minimize=TRUE,GG=GG)
    resulta<-rbind(resulta,BT)
   }
resulta<-unname(cbind(lambdaval,resulta))
return(resulta)
}

### 231019
vecinos<-function(CC) ## all sign changes, one at a time including CC, for simple search
{
  res<-CC
  for(i in 1:length(CC)){
    tmp<-CC; tmp[i]<--tmp[i]
    res<-rbind(res,tmp)
  }
  return(unique(unname(res))) 
}

## Lasso by quadrants using upsize-downsize movements
## the derivative is estimated numerically so no pseudo, but
## depends on lambda
#lassoqn(XM=DATX,YM=DATY,AA = AA1, TOL=8,eps = 10^-8,buscainicio=TRUE)->Bmo
lassoqn<-function(XM,YM,AA=NULL, TOL=6,eps=10^-10,buscainicio=FALSE)
{
  #  XM<-DATX;YM<-DATY;AA<-AA1; TOL<-6; eps<-10^-10; buscainicio<-TRUE
  rest<-res<-c()
  ### initial step
  Lm<-0;     p<-ncol(XM) ## number of variables
  AA0<-lm(YM~XM-1)
  Cm<-diag(sign(round(AA0$coefficients,TOL))) ### initial Cm
  if(buscainicio){ ## para buscar el punto inicial
    lassocono(XM=XM,YM=YM,lambdaval=c(0),AA=AA)->bb0
    Cm<-diag(bb0[1,(1:ncol(XM))+1])
  }
  Cm[is.na(diag(Cm)),is.na(diag(Cm))]<-0
  Rm<-minconoineq(CM=Cm,XM=XM,YM=YM,AA=AA,lamm=Lm)
  Qm<-Rm[length(Rm)];	Rm<-Rm[-length(Rm)]
  Qm<-sum((YM-XM%*%matrix(ncol=1,Rm))^2)/2
  res<-matrix(nrow=1,c(Lm,Rm,Qm)) ##   (lambda, beta, L)
  visited<-symbol(sign( round(x= Rm,digits=TOL-1)    ))
  ######################
  ### start of main loop
  #while(FALSE){	
  while(!identical(diag(Cm),rep(x=0,times=p))){
    jc<-(1:p)[diag(Cm%*%Cm)!=0] ## candidates to downsize
    jcc<-(1:p)[diag(Cm%*%Cm)==0] ## candidates to upsize
    rest<-c()
    
    ### if there are gaps, first try to go up
    if(length(jc)<p){  ## search candidates to upsize  ## original length(jc)<p
      for(kk in jcc){ ## kk is variable to upsize
        for(kandy in c(-1,1)){ ## candy gives sign
          Cmc<-Cm; Cmc[kk,kk]<-kandy;   ## print(Cmc); ## candidate orthant
          rest<-rbind(rest, newmoves(CM=Cmc,XM=XM,YM=YM,AA=AA,lamm=Lm,eps=eps,TOL=TOL  ))
        }
      } ##end kandy
    } ## end upszing
    
    
    ## then downsizing
    Cm->Cmc;  
    rest<-rbind(rest, newmoves(CM=Cmc,XM=XM,YM=YM,AA=AA,lamm=Lm,eps=eps, TOL=TOL))
    rest<-unname(rest);  
    flagg<-!istooclose(res,rest,TOL=TOL,p=p)
    if(sum(flagg)>=1) rest<-rest[flagg,]
    ## another patch to make sure method does not go down south
    
    Ind<-which.min(rest[,ncol(rest)-2]) ## next move 
    Lm<-rest[Ind,ncol(rest)-2]
    ##V<-c(V,symbol(diag(Cm))) ## abandon current orthant by upsizing/downsizing
    
    ## abandon current orthant by upsizing/downsizing
    
    interest<-rest[Ind,] ## the row to be updated
    ## we only take the values of interest here
    res<-unname(rbind(res,  interest[c( length(interest)-2,1:p, length(interest)-1)]  ))
    
    res[nrow(res),(1:p)+1]<-round(res[nrow(res),(1:p)+1],TOL)
    
    ### may be the bit that's needed, to keep the right lambda
    #  res[nrow(res),1]<-res[nrow(res),1]+res[nrow(res)-1,1]
    ###
    
    ## update
    #Cm	<-diag(sign(rest[Ind,(2*p+1):(3*p)]))
    Cm <- diag(sign(  round(x= interest[1:p], digits= TOL-1)     ) )  
    ## Lm <- interest[length(interest)-2]
    Lm<-res[nrow(res),1]
    visited<-c(visited,symbol(sign( round( x= interest[1:p], digits= TOL-1)     )))
    
  }	
  ### end of main loop
  ####################
  return(res)
} ## output is (lambda, beta, L)



## this function gives a sequence of valid moves
## does not minimize!
## does not remove entries too close to old ones for comfort 
newmoves<-function(CM,XM,YM,AA,lamm,eps,TOL)
{   
  resr<-lambdacc<-c(); p<-ncol(XM) #XXtr
  b0<-minconoineq(CM=CM,XM=XM,YM=YM,AA=AA,lamm=lamm)[-(ncol(XM)+1)] #XXtr
  b0p<-minconoineq(CM=CM,XM=XM,YM=YM,AA=AA,lamm=lamm+eps)[-(ncol(XM)+1)] #XXtr
  d0p<-(b0p-b0)/eps; 
  bcc<-matrix(ncol=ncol(XM),nrow=ncol(XM),0) #XXtr
  lambdacc<--b0/d0p ## candidate lambdas
  for(ii in 1:length(d0p)){  ## first shrink to zero
    if(d0p[ii]!=0){
      newc1<-  matrix(ncol=1,apply(cbind(b0*d0p[ii],-d0p*b0[ii]),1,sum))/d0p[ii]
      #  if(!isvectorAAok(betac=newc,AA=AA) ) 
      newc<-minconoineq(CM=CM,XM=XM,YM=YM,AA=AA,lamm=lamm+lambdacc[ii])[-(ncol(XM)+1)] #XXtr
      if(vektorstooclose(v1=newc1,v2=newc,TOL=TOL)) newc<-newc1
      bcc[,ii]<-newc
    }
  }
  if(!is.null(AA)){
    for(ii in 1:nrow(AA)){ ## then work all the constraints
      c1<-sum(AA[ii,]*d0p); c2<-sum(AA[ii,]*b0)
      if(c1!=0){
        newc1<-  matrix(ncol=1,apply(cbind(b0*c1 ,-d0p*c2),1,sum))/c1
        #  if(!isvectorAAok(betac=newc,AA=AA) ) 
        newc<-minconoineq(CM=CM,XM=XM,YM=YM,AA=AA,lamm=lamm-c2/c1)[-(ncol(XM)+1)] #XXtr
        if(vektorstooclose(v1=newc1,v2=newc,TOL=TOL)) newc<-newc1		
        lambdacc<-c(lambdacc,-c2/c1); bcc<-cbind(bcc,newc)
      }
    }
  }
  if(length(lambdacc)>=1){
    filtroc<-(!is.na(lambdacc)) &  (!is.infinite(lambdacc))  & (lambdacc>lamm*0);## big change & (apply(betacc>=0,2,prod)==1)		
    lambdacc<-lambdacc[filtroc];  ## clear NA, Inf,  <=\lambda ## \beta<0,
    bcc<-bcc[,filtroc]; 
  }
  if(length(bcc)==p) bcc<-matrix(ncol=1,bcc)
  ### update resr
  if(length(lambdacc)>=1){ ## valid lambda
    for(ik in 1:ncol(bcc)){
      Bmc<-matrix(ncol=1,bcc[,ik])
      Qmc<-sum((YM-XM%*%Bmc)^2)/2+(lambdacc[ik]+lamm)*sum(abs(Bmc))
      resr<-rbind(resr,  c(c(Bmc),sum(abs(Bmc)),lambdacc[ik]+lamm,Qmc,-1)  )			
    } ## for columns in betacc
  } 
  
  return(resr)
}

minconor<-function(XM,YM,lamm,AA=NULL,Reggie=1)
{
  XtX<-t(XM)%*%XM
  Unovec<-matrix(ncol=1, rep(1,2*ncol(XM)))
  XtY<-t(XM)%*%YM
  Dmat<-rbind(cbind(XtX,-XtX),cbind(-XtX,XtX+diag( rep(Reggie,ncol(XM))  )))
  dvec<-rbind(XtY,-XtY)-lamm*Unovec
  Amat<-t(rbind( diag(rep(1,2*ncol(XM))),cbind(AA,AA) )) ## crucial inequality
  AQ<-solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat)
  return(c(lamm,AQ$solution,AQ$value+1/2*sum(YM^2))) ## return lamm, betab,L
}

## to contract trajectory from relaxed lasso to TRUE coeffs
contraeresb<-function(resb,extras=TRUE,backcoeffs=FALSE,absolute=FALSE)
{
  if(extras) wrk<-resb[,-c(1,ncol(resb)),drop=FALSE]
  p<-ncol(wrk)/2
  TmpM<-rbind(diag(rep(1,p)),diag(rep(-1,p)))
  if(absolute) TmpM<-abs(TmpM)
  result<-wrk%*%TmpM
  if(!backcoeffs)	result<-unname(cbind(resb[,1],result,resb[,ncol(resb)]))
  return(result)
}


## to contract trajectory from relaxed lasso to ABSOLUTE coeffs
contraeresbabs<-function(resb,extras=TRUE,backcoeffs=FALSE)
{
  if(extras) wrk<-resb[,-c(1,ncol(resb)),drop=FALSE]
  p<-ncol(wrk)/2
  TmpM<-rbind(diag(rep(1,p)),diag(rep(1,p)))
  result<-wrk%*%TmpM
  if(!backcoeffs)	result<-unname(cbind(resb[,1],result,resb[,ncol(resb)]))
  return(result)
}




###################
## minimization constrained by inequalities over orthant
## NOT much meant but can be used with ZEROES in CM
minconoineq<-function(CM,XM,YM,lamm,AA=NULL)
{
  #  XM<-XDt; YM<-YDt; AA<-AA; lamm<-Lm; eps<-10^-10; CM<-Cm
  Cm<-CM; zero<-prod(diag(Cm)==0)==1
  valor<-0
  if(!zero){	
    XtX<-t(XM)%*%XM
    Unovec<-matrix(ncol=1,rep(1,ncol(XM)))
    Dmat<-CM%*%XtX%*%CM
    dvec<-CM%*%t(XM)%*%YM-lamm*CM%*%CM%*%Unovec; 
    Amat<-t(rbind(CM%*%CM,AA)) ## crucial inequality
    ## project down elements  can be used here with all CM
    Dmat<-projdown(CCC=Cm,MMM=Dmat,kual=3); 
    dvec<-projdown(CCC=Cm,MMM=dvec,kual=1); 
    Amat<-projdown(CCC=Cm,MMM=Amat,kual=1);
    ## apply quadprog, este es el nucleo de mincono
    if(sum(diag(Cm^2))>1) 
      AQ<-solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat) else AQ<-solve.QP(Dmat=Dmat,dvec=dvec,Amat=t(Amat))
    Uhat<-matrix(ncol=1,AQ$solution)
    ## restore the solution 
    Uhat<-projup(CCC=Cm,MMM=Uhat,kual=1)
    Bhat<-Cm%*%Uhat ## update to true quadrant
    valor<-AQ$value
  } else
  {
    Bhat<-matrix(ncol=1,nrow=ncol(XM),0);  # 5?
  }
  return(  c(Bhat, valor+1/2*sum(YM^2)) ) ## not Cm!
}


## search over all quadrants
minconoineqquadr<-function(XM,YM,lamm,AA=NULL,minimize=FALSE,GG=grid2(ncol(XM)))
{
  bigmatrix<-c()
  for(i in 1:nrow(GG)){
    BT<-minconoineq(CM=diag(GG[i,]),XM=XM,YM=YM,lamm=lamm,AA=AA)
    bigmatrix<-unname(rbind(bigmatrix,BT))		
  }
  resulta<-cbind(GG,bigmatrix) ## collect results
  if(minimize){
    indica<-which.min(bigmatrix[,ncol(bigmatrix)])
    resulta<-resulta[indica,]
  }
  return(resulta)
}


nextkands<-function(XM,YM,lamm,AA=NULL,Reggie=1,eps=10^-10) ## this gives a collection of lambdas
{
  betasb<-rbind(minconor(XM=XM,YM=YM,lamm=lamm,AA=AA,Reggie=Reggie),
                minconor(XM=XM,YM=YM,lamm=lamm+eps,AA=AA,Reggie=Reggie)  )
  ## first list with contractions
  betas<-contraeresb(resb=betasb,backcoeffs=TRUE)
  dbetas<-(betas[2,]-betas[1,])/eps
  lc<-c(-betas[1,]/dbetas)
  ## now deal with absolute values
  la<-c()
  if(!is.null(AA)){
    betasabs<-contraeresbabs(resb=betasb,backcoeffs=TRUE)
    dbetasabs<-(betasabs[2,]-betasabs[1,])/eps
    la<-c(-betasabs[1,]/dbetasabs)
  }
  lambdak<-round(c(lc,la),3)
  filtroc<-(!is.na(lambdak)) &  (!is.infinite(lambdak))  & (lambdak>lamm*0+eps);
  lambdak<-lambdak[filtroc]
  if(length(lambdak)==0) lambdak<-c()
  if(!is.null(lambdak)) Lm0<-lamm+min(lambdak) else Lm0<-max(abs(t(XM)%*%YM))
  return(Lm0)
}


### 
### TRUE relaxed Lasso, hierarchical
lassonr<-function(XM,YM,AA=NULL,TOL=8,eps=10^-10,Reggie=1)
{
  res<-c(); p<-ncol(XM) ## big result, it has lambda, betab, L
  ### initial step
  Lm<-0 ## initial lambda
  res<-matrix(nrow=1,c(minconor(XM=XM,YM=YM,lamm=Lm,AA=AA,Reggie=Reggie)))
  bc<-contraeresb(res,backcoeffs=TRUE); Cm<-sign( round(bc,TOL)  )
  while(!identical((Cm),rep(x=0,times=p))){
    Lm<-nextkands(XM=XM,YM=YM,lamm=Lm,AA=AA,Reggie=Reggie,eps=eps)
    res<-rbind(res,minconor(XM=XM,YM=YM,lamm=Lm,AA=AA,Reggie=Reggie))
    bc<-contraeresb(res,backcoeffs=TRUE)[nrow(res),];		   
    Cm<-sign( round(bc,TOL)  )
  }			 
  return(res)
}



## to create some lambda values, uniform and log scale
## 231019
espacialo<-function(mmin=0,mmax,nval=10)
{
  #  mmin<-lmin; mmax<-lmax; nval<-10
  res<-seq(from=mmin,to=mmax,length.out=nval)
  res<-c(res,10^seq(from=log10(mmin+(mmax-mmin)/20000),to=log10(mmin+(mmax-mmin)/15),length.out=nval))
  res<-unique(sort(res))
  return(res)
}

## Tests
# vecinos(CC=c(1,1,0,-1)) # three neighboring quadrants
# espacialo(mmax=5) ## 20 values for lambda between 0 and 5


### relaxed Lasso

### criterion for Lasso, relaxed
Lrelax<-function(bv,XM,YM,lamm)
{
## bv is vector stacked b+,b- as column vector
  nn<-length(bv)/2
  BB<-matrix(ncol=1,bv[1:nn]-bv[(nn+1):(2*nn)])
  return(sum((YM-XM%*%BB)^2)/2+lamm*sum(bv))
}

###################################################################################
###################################################################################


### Compute mean squared error using lasso trajectory and
### new data
### Coeffs is assumed to have first column lambda and
###    last column L; extras by default assumes this
### if extras=FALSE then data is only coefficients, by row
validalo<-function(XV,YV,Coeffs=NULL,XM=NULL,YM=NULL,AA=NULL,extras=TRUE)
{
if(is.null(Coeffs)) Coeffs<-lassoqn(XM=XM,YM=YM,AA=AA)
if(extras) Coeffs<-Coeffs[,-c(1,ncol(Coeffs))]
  res<-c()
  for(i in 1:nrow(Coeffs)){
	  err<- sum( ( YV-XV%*%matrix(ncol=1, Coeffs[i,])   )^2  )
    res<-c(res,err) 
  }
#	BM<-matrix(nrow=ncol(Coeffs)-2,ncol=nrow(Coeffs),byrow=FALSE, Coeffs[i,-c(1,ncol(Coeffs))] )
#	BM<-t(Coeffs); BM<-BM[-c(1,ncol(Coeffs)),]
#   err<-matrix(nrow=length(YM),ncol=nrow(Coeffs),byrow=FALSE,YM)-XM%*%BM	
#    apply(err^2,2,sum)/nrow(XV) ## should give the same result
return(unname(res/nrow(XV))) 
}


### Re-estimate using least squares
### extras assumes that Coeffs has extra columns: leading
### column of lambdas and trailing column of L
estimalo<-function(XM,YM,AA=NULL,Coeffs=NULL,extras=TRUE,TOL=6)
{
C0<-Coeffs ## to rescue lambdas and L
if(is.null(Coeffs)) Coeffs<-lassoqn(XM=XM,YM=YM,AA=AA)
if(extras) Coeffs<-Coeffs[,-c(1,ncol(Coeffs))]
  res<-c()
  for(i in 1:nrow(Coeffs)){
	  Cloc<-diag(sign(round(x=Coeffs[i,],digits=TOL)))
    Modlin<-lm(YM~XM%*%Cloc%*%Cloc-1)
		bloc<-c(Modlin$coefficients)
		bloc[is.na(bloc)]<-0
		res<-rbind(res,bloc)
  }
if(extras) res<-cbind(C0[,1],res,C0[,ncol(C0)])
return(unname(res))
}

### Plot
## the Coeffs has the lasso path, and extras clarifies
## if first column is lambda and last is L
## by default plot against s
## simplito gives simple axes labels
graficalo<-function(Coeffs,LM=NULL,extras=TRUE,plotlambda=FALSE,typolinea=1,titulo=NULL,forzalo=FALSE,alfa=0,cexval=0.95,
                        kolor="lightgray",simplito=TRUE)
{
nombrehorizontal<-"|beta|/max|beta|"
if(extras){ 
   lambdav<-Coeffs[,1]; Lv<-Coeffs[,ncol(Coeffs)];
	 Coeffs<-Coeffs[,-c(1,ncol(Coeffs))]	
	 sval<-apply(abs(Coeffs),1,sum); sval<-sval/max(sval)
	 if(plotlambda) {
	     sval<-lambdav; nombrehorizontal<-"lambda"
			}
   } else
	{
	 sval<-apply(abs(Coeffs),1,sum); sval<-sval/max(sval)	
	}
	 hval<-sval
	if(forzalo){
	ordena<-order(sval,decreasing=TRUE); hval<-hval[ordena]; #Coeffs<-Coeffs[ordena,]
	}
if(simplito) plot(range(hval),range(Coeffs),type="n",xlab=nombrehorizontal,ylab="beta",main=titulo) else
  plot(range(hval),range(Coeffs),type="n",
       xlab=expression(  paste("||", beta(lambda), "||",phantom()[{ paste("1")}], "", "","/max", phantom()[{ paste(lambda)}], "", ""       ,"||", beta(lambda), "||",phantom()[{ paste("1")}], "", "")     ),
       ylab=expression(   paste(beta(lambda))   ),main=titulo)
for(i in 1:nrow(Coeffs)) lines(hval[i]*c(1,1),range(Coeffs),col=kolor,lwd=0.5)
for(i in 1:ncol(Coeffs)) lines(hval,Coeffs[,i],col=i,lty=typolinea)
if(!is.null(LM)) text(x=1*(1-plotlambda)-abs(rnorm(nrow(LM),sd=0.1*alfa)),y=Coeffs[1,]+rnorm(nrow(LM),sd=0.1*alfa),labels=nombres(LM),cex=cexval)
}
 

#################################################
#################################################
## functions for hierarchical models and
## restrictions
##

## This function creates the simplified Hasse diagram
## Each row corresponds to a term in LM and numbers in columns are entries that 
## depend on LM, via simple product of unit exponents
hassediagram<-function(LM){
d<-ncol(LM); res<-c()
for(i in 1:nrow(LM))
{
   for(j in 1:d){
        unit<-matrix(ncol=d,0)
        unit[1,j]<-1
        cand<-LM[i,]+unit
        sum((c(nombres(matrix(nrow=1,cand)))==nombres(LM))*(1:nrow(LM)))->result 
        ## if zero entries in output then there is NO nesting
        ## else the index says who is nested #  print(result)
     res<-c(res,result)
      }   
}
return(matrix(ncol=d,res,byrow=TRUE))
}


## This function reweights a constraint
## according to frequency by defaults uses weights
reweight<-function(vektor,weight=TRUE){
if(weight){
  frtable<-table(vektor[vektor!=0])
  ttable<-rbind(as.integer(names(frtable)),unname(frtable))
  if(length(unique(ttable[2,]))>1){
  (valor<-ttable[1,][(ttable[2,])==1])
  (peso<-ttable[2,][(ttable[2,])!=1])
  vektor[vektor==valor]<-vektor[vektor==valor]*peso
  }
}
  return(vektor)
}

## rescale the matrix of restrictions
## this is intended to generate large values
## for positive coefficients
## factor multiplies maximum, unless override
## is TRUE so that it is only factor 
rescale<-function(MAT,factor=1,override=FALSE)
{
  values<-c(MAT); values<-values[values>0]
  maximo<-max(values)
if(override) maximo<-1
  for(i in 1:nrow(MAT)) MAT[i,MAT[i,]>0]<-maximo*factor
  return(MAT)
}

## This function generates the matrix of restrictions
## if type = 1 then return strong hierarchy adding over descendants
## if type = 2 then return weak hierarchy, adding over parents
## otherwise (type = 3) return general strong hierarchy
## by default, weight = TRUE
## if weight=TRUE then all non-zero entries are +-1
##    weight does not matter for type=3
restrict<-function(LM,type=3,weight=TRUE){
 weight<-!weight
 hassemat<-hassediagram(LM); d<-ncol(LM)
 ## this vector determines indexes of parameters
 indexhmat<-matrix(ncol=1,byrow=TRUE,1:((2-1)*nrow(LM)))
 ## the (2-1) introduced to not have relaxation
 truthmat<-apply(hassemat,1,sum)!=0 ## this one determines
 ## where we have non-null hasse branches going down
restrictions2<-restrictions<-c()
for(i in 1:nrow(LM))
{  
   rowsource<-matrix(ncol=(2-1)*nrow(LM),0)  ## the (2-1) introduced to not have relaxation
if(truthmat[i]==TRUE)
{
   rowsource[indexhmat[i,]]<-1 ## this identifies the parent node
   rowsource[ c(indexhmat[ hassemat[i,hassemat[i,]!=0],])]<--1 ## for all edges
   restrictions<-rbind(restrictions,reweight(rowsource,weight=weight))
   rowsource[ c(indexhmat[ hassemat[i,hassemat[i,]!=0],])]<--0 ## clear for all edges
 ## the version for every edge
 for(j in c( hassemat[i,hassemat[i,]!=0]) )
 {
   rowsource[ c(indexhmat[ j,])]<--1
   restrictions2<-rbind(restrictions2,reweight(rowsource,weight=weight))
   rowsource[ c(indexhmat[ j,])]<-0
 }
}
}
descendants<-unique(c(hassemat))
descendants<-descendants[descendants!=0]
restrictions3<-c()
for(aa in descendants){
   rowsource<-matrix(ncol=(2-1)*nrow(LM),0)  ## the (2-1) introduced to not have relaxation
   rowsource[aa]<--1; #  print(aa)
  for(i in 1:d){
#    print(hassemat[,i]==aa);  print((1:nrow(hassemat))[hassemat[,i]==aa])
       kual<-(1:nrow(hassemat))[hassemat[,i]==aa]
     rowsource[kual]<-1
    }
  restrictions3<-rbind(restrictions3,reweight(rowsource,weight=weight))
}
#restrictions ## general restriction as in Bien; adding over descendants
#restrictions2 ## general strict restriction, less sparse
#restrictions3 ## adding over parents sparsest
resultado<-restrictions2
if(type==1) resultado<-restrictions
if(type==2) resultado<-restrictions3
return(resultado)
}


## Tests
# hassediagram(LM=multidegree(d=3,m=2)) ## dependence of higher order monomials on lower order monomials
# restrict(LM=multidegree(d=3,m=2)) ## The set H of constraints with 9 elements because of intercept, also type=3
# restrict(LM=multidegree(d=3,m=2),type=1) ## Set S of constraints with 4 parent elements 
# restrict(LM=multidegree(d=3,m=2),type=2) ## Set W of constraints with 6 son elements 

##
## functions for handling models
## for hierarchical lasso 030913 090414 130619
## 
## formerly in file functions_v2.R
##


## all multilinear terms
## d variables up to degree m
## m<=d if m=d then 2^d terms
multidegree2<-function(d,m)
{
res<- unique(
rbind(rep(0,d),
t((xsimplex(d,m)>0)*1))
)
return(res[order(res%*%rep(1,d)),])
} 

## binary grid of n components size 2^n
bingrid <- function(n)
{
 if (n==1) {  
  res <- c(0,1); dim(res) <- c(2,1)
 }
 else {
  tmp <- bingrid(n-1)
  res <- rbind(cbind(0,tmp),cbind(1,tmp))   
 }
 return(res)
}

## better multidegree using bingrid
multidegree<-function(d,m) ## d variables, degree m
{
 res<-bingrid(d)
 res<-res[res%*%rep(1,d)<=m,]
 return(res[order(res%*%rep(1,d)),])
}

###### design matrix, imported 130813 from old SSM prePC

evalm<-function(monom,desp) ## evaluate monomial at single design point
{ res<-1;
#for(j in 1:length(monom)) res<-res*desp[j]^monom[j]
tmp<-desp^monom
tmp[monom==0]<-1
return(prod(tmp));  
}

evall<-function(LM,desp) ## evaluate list of monomials at single design point
{ res<-c()
for(i in 1:nrow(LM)) res<-rbind(res,evalm(monom=LM[i,],desp=desp))
return(res)
}

desmod<-function(LM,Des) ## SIMPLE design-model matrix
{ res<-c()
for(i in 1:nrow(Des)) res<-rbind(res,c(evall(LM=LM,desp=Des[i,])))
return(res)  
}

buildmodel<-function(LM,Des,standardise=TRUE)
{
  if(standardise) for(i in 1:ncol(Des)) Des[,i]<-(Des[,i]-mean(Des[,i]))/sd(Des[,i])
  DMM<-desmod(LM=LM,Des=Des)
  if(standardise) for(i in 1:ncol(DMM)) ## this below deals with intercept
    if(!identical(LM[i,],rep(0,ncol(LM))))  DMM[,i]<-(DMM[,i]-mean(DMM[,i]))
  return(DMM)
}

#############################################################


### some postprocessing capabilities
## return the polynomial in a more readable form
## lista is the list with variable names
## MD is the list of exponents
givemepoly<-function(lista,MD,ismat=FALSE)
{
  if(!ismat) MD<-matrix(ncol=length(lista),byrow=TRUE,MD)
res<-c()
  for(i in 1:nrow(MD)){
    kual<-MD[i,]!=0; subl<-lista[kual]; subMD<-MD[i,kual]; rr<-c()
  for(kk in 1:length(subl))
     if(subMD[kk]==1) rr<-paste(rr,subl[kk],sep="") else rr<-paste(rr,subl[kk],"^",subMD[kk],sep="")
   res<-c(res,rr)
    }
return(res)
}

