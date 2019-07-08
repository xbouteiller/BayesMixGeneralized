#m14 <- MCMCglmm(Ht_16 ~ Bloc,
#                random =  ~ Alt + Bloc:Alt + Mother + Bloc:Mother,
#                data   = Toul,
#                family = rep("poisson",1),
#                prior = prior.m,
#                nitt = NITT, burnin = BURN, thin = THIN)
library(R2jags)

###################################################################################

ColTrait<-c(12 # 1 - D_14
            ,13 # 2 -  Ht_14
            ,15 # 3 -  Ntot14
            ,17 # 4 -  Mtot14
            ,18 # 5 -  oneacorn 14
            ,19 # 6 -  D_15
            ,20 # 7 -  Ht_15
            ,22 # 8 -  Ntot15
            ,24 # 9 -  Mtot15
            ,25 # 10 -  oneacorn15
            ,26 # 11 -  D_16
            ,27 # 12 -  Ht_16
            ,29 # 13 -  Ntot16
            ,31 # 14 -  Mtot16
            ,32 # 15 -  oneacorn16
            ,37 # 16 -  NtotAll
            ,38 # 17 -  MtotAll
            ,39 # 18 -  oneacornMean
            ,40 # 19 -  Ntot14sans
            ,41 # 20 -  Mtot14sans
            ,42 # 21 -  Ntot15sans
            ,43 # 22 -  Mtot15sans
            ,44 # 23 -  Ntot16sans
            ,45 # 24 -  Mtot16sans
            ,46 # 25 -  NtotAllsans
            ,47 # 26 -  MtotAllsans
            ,49 # 27 -  NtotAllNorm
            ,50 # 28 -  MtotAllNorm
            ,51 # 29 -  GrowthHt
            ,52 # 30 -  GrowthDia
            )


ColTrait<-c(26 # 1 -  D_16
            ,27 # 2 -  Ht_16
            ,39 # 3 -  oneacornMean
            ,51 # 4 -  GrowthHt
            ,52 # 5 -  GrowthDia
)

Mat<-NULL
for (j in 1:length(ColTrait))
{
tab1<- Toul
which(is.na(Toul[ColTrait[j]]))->lesNA
tab1[- lesNA,]->tab2
tab2

##BETWEEN TEMPERATURES (CHAMBERS) MODEL 
##Add comment symbol "#" if you want to use the within chamber model in order to estimate QST
NINDIV=sum(table(tab2$Alt[tab2[ColTrait[j]]>=0])) #NINDIV = Number of individuals
NBLOC=nlevels(tab2$Bloc)
NALT=nlevels(tab2$Alt)
NFAM=nlevels(tab2$Mother)

tab2$Alt<-as.factor(tab2$Alt)

modelstring<-"model{

    
##likelihood
    for(i in 1:NINDIV){
    mod[i]<-mu + a1[Bloc[i]] + A2[Alt[i]] + A3[Mother[i]] + A4[Bloc[i],Alt[i]] + A5[Bloc[i],Mother[i]] 
    Trait[i]~dnorm(mod[i], precision)
    Trait2[i]~dnorm(yMean, TAUTRAIT)
    }
  #Priors sur la précision résiduelle 
    ySigma~dunif(0,1000)
    precision<-1/(ySigma^2)
  
  #Priors sur la moyenne générale
    mu~dnorm(yMean,0.001)
    TAUTRAIT ~ dgamma(0.0001, 0.0001)
    
# Fixed effect
  #Priors sur le facteur Bloc
    for(t in 1:NBLOC) #NBLOC : number of Bloc
    {
      a1[t]~dnorm(0.0, 0.000001)
    }
# Random effects
  # Priors sur le facteur Altitude: Provenance 
    for(p in 1:NALT) # NALT : Numbers of populations
    {
      A2[p]~dnorm(0, tauP)
    }
    #Hyperpriors
    yP ~dunif(0,1000)
    tauP<-1/(yP^2)

    
  # Priors sur le facteur Mother: Family 
    for(q in 1:NFAM) # NFAM : Numbers of family
    {
      A3[q]~dnorm(0, tauF)
    }
    #Hyperpriors
    yF ~dunif(0,1000)
    tauF<-1/(yF^2)

  # Priors sur l'interaction Bloc*Prov 
    for(t in 1:NBLOC) {for(p in 1:NALT){ 
      A4[t,p]~dnorm( 0.0 ,tauINT1 ) ##!!!
      
    }}
    yINT1 ~dunif(0,1000)
    tauINT1<-1/(yINT1^2)

  
  # Priors sur l'interaction Bloc*Family 
    for(t in 1:NBLOC) {for(q in 1:NFAM){ 
      A5[t,q]~dnorm( 0.0 ,tauINT2 ) ##!!!
      
    }}
    yINT2 ~dunif(0,1000)
    tauINT2<-1/(yINT2^2)


##sum to zero constraint (see Kruschke 2015)
for(t in 1:NBLOC){ for(p in 1:NALT) { for(q in 1:NFAM) {

m[t,p,q]<-a1[t] + A2[p] + A3[q] + A4[t,p] + A5[t,q]

}}}

b0<-mu + mean(m[1:NBLOC,1:NALT,1:NFAM])

for(t in 1:NBLOC) {b1[t]<-mean(m[t,1:NALT,1:NFAM])-mean(m[1:NBLOC,1:NALT,1:NFAM])}
for(p in 1:NALT) {B2[p]<-mean(m[1:NBLOC,p,1:NFAM])-mean(m[1:NBLOC,1:NALT,1:NFAM])}
for(q in 1:NFAM) {B3[q]<-mean(m[1:NBLOC,1:NALT,q])-mean(m[1:NBLOC,1:NALT,1:NFAM])}


  ##Estimation of random effect variances
varP<-1/tauP
varF<-1/tauF
varINT1<-1/tauINT1
varINT2<-1/tauINT2
  QST <- varP /(varP+ 8*varF)
resid<-1/precision

  ##HERITABILITE
varTRAIT = 1/ TAUTRAIT
h2<-(4*varF) / varTRAIT
h2manu<-(4*varF) / (varP+ varF+ varINT1+ varINT2+resid)

  
  }"
writeLines(modelstring,con="mixed_model.txt")

###################################################################################
#parametres d'iterations
NECH = 60000
NBURN = 10000
NTHIN = 50


datalist = list(
  Trait = tab2[ColTrait[j]][,1],
  Trait2=tab2[ColTrait[j]][,1],
  Bloc= tab2$Bloc,
  Alt = tab2$Alt,
  Mother = tab2$Mother,
  
  yMean = mean(tab2[ColTrait[j]][,1]),
  
  
  NINDIV=sum(table(tab2$Alt[tab2[ColTrait[j]]>=0])),  #NINDIV = Number of individuals
  NBLOC=nlevels(tab2$Bloc),
  NALT=nlevels(tab2$Alt),
  NFAM=nlevels(tab2$Mother)
  
)


#initial    & chains
NRANGE= nlevels(tab2$Bloc)
NPOP=nlevels(tab2$Alt)
NTREE=nlevels(tab2$Mother)

nchain <- 4


initial<-c(mu = 0, a1 = rep(0,NRANGE) , A2= rep(0,NPOP), A3 = rep(0,NTREE) )


init<-list()
for(k in 1:nchain){init[[k]]<-initial}

# Nombre d'échantillons à extraire par chaîne dans le posterior
Nechantillons <- NECH

# échantillons initiaux à éliminer (burn in) par chaîne
burnin <- NBURN#125000#/100 #50000

# pas d'éclaircissage (thinning)
eclair<-NTHIN #100

# quantités (noeuds) dont nous voulons sauvegarder les échantillons
parameters = c("b0", "b1", "varP", "varF",  "QST", "h2", "h2manu")

# compilation et lancement du calcul  
mix <-jags(model.file="mixed_model.txt", data=datalist, n.chains=nchain, inits=init, parameters.to.save=parameters,
           n.iter=Nechantillons, n.burnin=burnin, n.thin=eclair)


x11()
post<-as.mcmc(mix)
plotQST_alt(post)
segments(x0=mean(loc.Fst.Tot), x1=mean(loc.Fst.Tot),
         y0=-10000, y1=100000,
         col="black",lwd=2)

polygon(x=c(lwr,lwr, upr, upr),
        y=c(0, 100000, 100000, 0),
        col = alpha("orange", 0.3), lty=0)



mat=c(mix$BUGSoutput[11]$mean$QST, mix$BUGSoutput[11]$mean$h2, 
           mix$BUGSoutput[11]$mean$h2manu, mix$BUGSoutput[11]$mean$varF, mix$BUGSoutput[11]$mean$varP)

matInf=c(mix$BUGSoutput[10]$summary[1,3], mix$BUGSoutput[10]$summary[3,3],
             mix$BUGSoutput[10]$summary[4,3],mix$BUGSoutput[10]$summary[5,3], mix$BUGSoutput[10]$summary[6,3])

matSup=c(mix$BUGSoutput[10]$summary[1,7], mix$BUGSoutput[10]$summary[3,7],
         mix$BUGSoutput[10]$summary[4,7],mix$BUGSoutput[10]$summary[5,7], mix$BUGSoutput[10]$summary[6,7])
  
resMat=matrix(c(mat, matInf, matSup), ncol=3)
rownames(resMat)<-c("QST", "h2", "h2manu", "varF", "varP")
colnames(resMat)<-c("mean", "Inf", "Sup")
Mat<-rbind(Mat, resMat)
}
write.table(Mat, file="Toul_QST2_2.csv", sep=";", row.names=T)
              
           
              
              
      plotQST_alt<-function(X, NOM.axeX="", NOM.axeY="", Titre="Qst - Fst", NCHAIN = nchain ){
        which(colnames(X[[1]])=="QST")->pourb0
        list(density(X[[1]][,pourb0], ad=2)$x)->axeX
        list(density(X[[1]][,pourb0], ad=2)$y)->axeY
        min(unlist(axeX))->min.axeX
        max(unlist(axeX))->max.axeX
        min(unlist(axeY))->min.axeY
        max(unlist(axeY))->max.axeY
        
        list()->TOBEPLOTTED
        for(p in 1:NCHAIN){
          TOBEPLOTTED[[p]]<-X[[p]][,pourb0]
        }
        
        hist(unlist(TOBEPLOTTED),# main = colnames(X[[1]])[pourb0],
             xlab= NOM.axeX,
             ylab=NOM.axeY, 
             main= Titre,
             cex.lab = 1.5,
             cex.axis = 2,
             lwd=2,
             xlim=c(0,1),
             breaks=10)
        
        list()->forthequantile
        for(m in 1:NCHAIN){
          forthequantile[[m]]<-X[[m]][,pourb0]
        }
        
        quantile(unlist(forthequantile),  probs = c( 0.025, 0.975))->quantb0
        
        
        
        polygon(x=c(quantb0[1],quantb0[1], quantb0[2], quantb0[2]),
                y=c(0, 1000000, 1000000, 0),
                col = alpha("dodgerblue3", 0.1), lty=0)
        #text(x=quantb0[1],y=max.axeY/20, labels = signif(quantb0[1],3), cex =1.5)
        #text(x=quantb0[2],y=max.axeY/20, labels = signif(quantb0[2],3), cex =1.5)
        
        abline(v =mean(unlist(forthequantile)), lwd=2, lty = 5, col = "black")
        text(x=mean(unlist(forthequantile)),y=max.axeY/7, labels = signif(mean(unlist(forthequantile)),3), cex =1.5)
        
      }
      
      x11()
      post<-as.mcmc(mix)
      plotQST_alt(post)
      segments(x0=mean(loc.Fst.Tot), x1=mean(loc.Fst.Tot),
               y0=-2, y1=10,
               col="black",lwd=2)
      
      polygon(x=c(lwr,lwr, upr, upr),
              y=c(0, 100, 100, 0),
              col = alpha("orange", 0.3), lty=0)
      
###################################################################################



###################################################################################
####                            Modèle Poissson (Ntot)                         #### 
###################################################################################

              
ColTrait<-c(15 # 1 -  Ntot14
            ,22 # 2 -  Ntot15
            ,29 # 3 -  Ntot16
            ,37 # 4 -  NtotAll
            ,40 # 5 -  Ntot14sans
            ,42 # 6 -  Ntot15sans
            ,44 # 7 -  Ntot16sans
            ,46 # 8 -  NtotAllsans
)
              
MatPois<-NULL
for (j in 1:length(ColTrait))
{
  tab1<- Toul
  which(is.na(Toul[ColTrait[j]]))->lesNA
  tab1[- lesNA,]->tab2
  tab2
  
  ##BETWEEN TEMPERATURES (CHAMBERS) MODEL 
  ##Add comment symbol "#" if you want to use the within chamber model in order to estimate QST
  NINDIV=sum(table(tab2$Alt[tab2[ColTrait[j]]>=0])) #NINDIV = Number of individuals
  NBLOC=nlevels(tab2$Bloc)
  NALT=nlevels(tab2$Alt)
  NFAM=nlevels(tab2$Mother)
  
  tab2$Alt<-as.factor(tab2$Alt)

modPoiss<-"model{

#vraisemblance

for(i in 1:NINDIV)  
{
  Trait[i]~dpois(mod[i])
  log(mod[i])<-mu + a1[Bloc[i]] + A2[Alt[i]] + A3[Mother[i]] + A4[Bloc[i],Alt[i]] + A5[Bloc[i],Mother[i]] 
  eps[i]~dnorm(0.0, TauErr)
  Trait2[i]~dnorm(yMean, TAUTRAIT)
}


############################ Priors

  ySigma ~dunif(0,1000)
  TauErr<-1/(ySigma^2)
  mu~dnorm(yMean,0.001)
  TAUTRAIT ~ dgamma(0.0001, 0.0001)

# Fixed effect
  #Priors sur le facteur Bloc
    for(t in 1:NBLOC) #NBLOC : number of Bloc
    {
      a1[t]~dnorm(0.0, 0.000001)
    }
# Random effects
  # Priors sur le facteur Altitude: Provenance 
    for(p in 1:NALT) # NALT : Numbers of populations
    {
      A2[p]~dnorm(0, tauP)
    }
  #Hyperpriors
    yP ~dunif(0,1000)
    tauP<-1/(yP^2)


  # Priors sur le facteur Mother: Family 
    for(q in 1:NFAM) # NFAM : Numbers of family
    {
      A3[q]~dnorm(0, tauF)
    }
  #Hyperpriors
    yF ~dunif(0,1000)
    tauF<-1/(yF^2)

  
  # Priors sur l'interaction Bloc*Prov 
    for(t in 1:NBLOC) {for(p in 1:NALT){ 
    A4[t,p]~dnorm( 0.0 ,tauINT1 ) ##!!!
    
    }}
    yINT1 ~dunif(0,1000)
    tauINT1<-1/(yINT1^2)

  
  # Priors sur l'interaction Bloc*Family 
    for(t in 1:NBLOC) {for(q in 1:NFAM){ 
    A5[t,q]~dnorm( 0.0 ,tauINT2 ) ##!!!
    
    }}
    yINT2 ~dunif(0,1000)
    tauINT2<-1/(yINT2^2)


##sum to zero constraint (see Kruschke 2015)
  for(t in 1:NBLOC){ for(p in 1:NALT) { for(q in 1:NFAM) {
  
  m[t,p,q]<-a1[t] + A2[p] + A3[q] + A4[t,p] + A5[t,q]
  
  }}}

  b0<-mu + mean(m[1:NBLOC,1:NALT,1:NFAM])
  
  for(t in 1:NBLOC) {b1[t]<-mean(m[t,1:NALT,1:NFAM])-mean(m[1:NBLOC,1:NALT,1:NFAM])}
  for(p in 1:NALT) {B2[p]<-mean(m[1:NBLOC,p,1:NFAM])-mean(m[1:NBLOC,1:NALT,1:NFAM])}
  for(q in 1:NFAM) {B3[q]<-mean(m[1:NBLOC,1:NALT,q])-mean(m[1:NBLOC,1:NALT,1:NFAM])}


##Estimation of random effect variances
  varP<-1/tauP
  varF<-1/tauF
  varINT1<-1/tauINT1
  varINT2<-1/tauINT2
  QST <- varP /(varP+ 8*varF)

##HERITABILITE
  varTRAIT = 1/ TAUTRAIT
  h2<-(4*varF) / varTRAIT





}"


writeLines(modPoiss, con="mixed_modelPois.txt")

###################################################################################
#parametres d'iterations
NECH = 60000
NBURN = 10000
NTHIN = 10

datalist = list(
  Trait = tab2[ColTrait[j]][,1],
  Trait2=tab2[ColTrait[j]][,1],
  Bloc= tab2$Bloc,
  Alt = tab2$Alt,
  Mother = tab2$Mother,
  
  yMean = mean(tab2[ColTrait[j]][,1]),
  
  
  NINDIV=sum(table(tab2$Alt[tab2[ColTrait[j]]>=0])),  #NINDIV = Number of individuals
  NBLOC=nlevels(tab2$Bloc),
  NALT=nlevels(tab2$Alt),
  NFAM=nlevels(tab2$Mother)
  
)


#initial    & chains
NRANGE= nlevels(tab2$Bloc)
NPOP=nlevels(tab2$Alt)
NTREE=nlevels(tab2$Mother)

nchain <- 4


initial<-c(mu = 0, a1 = rep(0,NRANGE) , A2= rep(0,NPOP), A3 = rep(0,NTREE) )


init<-list()
for(k in 1:nchain){init[[k]]<-initial}

# Nombre d'échantillons à extraire par chaîne dans le posterior
Nechantillons <- NECH

# échantillons initiaux à éliminer (burn in) par chaîne
burnin <- NBURN#125000#/100 #50000

# pas d'éclaircissage (thinning)
eclair<-NTHIN #100

# quantités (noeuds) dont nous voulons sauvegarder les échantillons
parameters = c("QST", "h2", "h2manu", "varF", "varP")

# compilation et lancement du calcul  
mix <-jags(model.file="mixed_modelPois.txt", data=datalist, n.chains=nchain, inits=init, parameters.to.save=parameters,
           n.iter=Nechantillons, n.burnin=burnin, n.thin=eclair)

mat=c(mix$BUGSoutput[11]$mean$QST, mix$BUGSoutput[11]$mean$h2,
      mix$BUGSoutput[11]$mean$varF,mix$BUGSoutput[11]$mean$varP)

matInf=c(mix$BUGSoutput[10]$summary[1,3],mix$BUGSoutput[10]$summary[3,3],
         mix$BUGSoutput[10]$summary[4,3], mix$BUGSoutput[10]$summary[5,3])

matSup=c(mix$BUGSoutput[10]$summary[1,7],mix$BUGSoutput[10]$summary[3,7],
         mix$BUGSoutput[10]$summary[4,7], mix$BUGSoutput[10]$summary[5,7])

resMat=data.frame(matrix(c(mat, matInf, matSup), ncol=3))
rownames(resMat)<-c("QST", "h2", "varF", "varP")
colnames(resMat)<-c("mean", "Inf", "Sup")
MatPois<-rbind(MatPois, resMat)
}

write.table(MatPois, file="Toul_QST2Pois.csv", sep=";", row.names=T)


###################################################################################
####                      Modèle Négative Binomiale (Ntot)                     #### 
###################################################################################


## Vérification: Le trait suit il une loi binomiale négative
## surpéposition des histograms
###################################################################################
tab1<- Toul
which(is.na(Toul$NtotAllNorm))->lesNA
tab1[- lesNA,]->tab2
tab2

library(MASS)
Param<-fitdistr(round(tab2$NtotAllNorm) ,"negative binomial") # fitter une distrib neg binom sur nos données 

#simule des données
x.negbin<-rnbinom(n = length(tab2$NtotAllNorm), size =Param[1]$estimate[1] ,  mu =Param[1]$estimate[2]) 
# hist(x.poi,main="Poisson distribution")

#comparer les ditrib
hist(tab2$NtotAllNorm,col=rgb(0.8,0.0,0.0,0.5))
hist(x.negbin, col=rgb(0.0,0.0,0.8,0.5), add=T)
###################################################################################




ColTrait<-c(#15 # 3 -  Ntot14
            #17 # 4 -  Mtot14
            #,22 # 8 -  Ntot15
            #,24 # 9 -  Mtot15
            #,29 # 13 -  Ntot16
            #,31 # 14 -  Mtot16
            37 # 16 -  NtotAll
            #,
            #38 # 17 -  MtotAll
            ,49 # 27 -  NtotAllNorm
            #50 # 28 -  MtotAllNorm
)

j=50
MatNegBi<-NULL
for (j in 1:length(ColTrait))
{
  tab1<- Toul
  which(is.na(Toul[ColTrait[j]]))->lesNA
  tab1[- lesNA,]->tab2
  tab2
  
  
  
  ##BETWEEN TEMPERATURES (CHAMBERS) MODEL 
  ##Add comment symbol "#" if you want to use the within chamber model in order to estimate QST
  NINDIV=sum(table(tab2$Alt[tab2[ColTrait[j]]>=0])) #NINDIV = Number of individuals
  NBLOC=nlevels(tab2$Bloc)
  NALT=nlevels(tab2$Alt)
  NFAM=nlevels(tab2$Mother)
  tab2$Alt<-as.factor(tab2$Alt)
  
  
modNegBin<-"model{

#vraisemblance

for(i in 1:NINDIV)  
{
log(mod[i]) <- mu + a1[Bloc[i]] + A2[Alt[i]] + A3[Mother[i]] + A4[Bloc[i],Alt[i]] + A5[Bloc[i],Mother[i]]  + eps[i]
    
	eps[i]~dnorm(0.0, TauErr)
	
    # Transforms mu into p, which is used by the negative binomial distribution
    p[i] <- r/(r + mod[i])
	
	Trait[i] ~ dnegbin(p[i], r)
	Trait2[i] ~ dnegbin(p2, r2)
}


############################ Priors

ySigma ~dunif(0,1000)
TauErr<-1/(ySigma^2)
mu~dnorm(yMean,0.001)
r ~ dunif(0, 50) 


# Priors sur les paramètres généraux

p2 ~ dbeta(1.001,1.001)
r2 ~ dgamma(0.01,0.01)

MeanTRAIT <- r2*(1-p2)/(p2)
varTRAIT <- r2*(1-p2)/(p2*p2)

 
# Fixed effect
#Priors sur le facteur Bloc
for(t in 1:NBLOC) #NBLOC : number of Bloc
{
  a1[t]~dnorm(0.0, 0.000001)
}
# Random effects
# Priors sur le facteur Altitude: Provenance 
for(p in 1:NALT) # NALT : Numbers of populations
{
  A2[p]~dnorm(0, tauP)
}
#Hyperpriors
    yP ~dunif(0,1000)
    tauP<-1/(yP^2)


# Priors sur le facteur Mother: Family 
for(q in 1:NFAM) # NFAM : Numbers of family
{
  A3[q]~dnorm(0, tauF)
}
#Hyperpriors
    yF ~dunif(0,1000)
    tauF<-1/(yF^2)


# Priors sur l'interaction Bloc*Prov 
for(t in 1:NBLOC) {for(p in 1:NALT){ 
A4[t,p]~dnorm( 0.0 ,tauINT1 ) ##!!!

}}
    yINT1 ~dunif(0,1000)
    tauINT1<-1/(yINT1^2)


# Priors sur l'interaction Bloc*Family 
for(t in 1:NBLOC) {for(q in 1:NFAM){ 
A5[t,q]~dnorm( 0.0 ,tauINT2 ) ##!!!

}}
    yINT2 ~dunif(0,1000)
    tauINT2<-1/(yINT2^2)

##sum to zero constraint (see Kruschke 2015)
for(t in 1:NBLOC){ for(p in 1:NALT) { for(q in 1:NFAM) {

m[t,p,q]<-a1[t] + A2[p] + A3[q] + A4[t,p] + A5[t,q]

}}}

b0<- mu + mean(m[1:NBLOC,1:NALT,1:NFAM])

for(t in 1:NBLOC) {b1[t]<-mean(m[t,1:NALT,1:NFAM])-mean(m[1:NBLOC,1:NALT,1:NFAM])}
for(p in 1:NALT) {B2[p]<-mean(m[1:NBLOC,p,1:NFAM])-mean(m[1:NBLOC,1:NALT,1:NFAM])}
for(q in 1:NFAM) {B3[q]<-mean(m[1:NBLOC,1:NALT,q])-mean(m[1:NBLOC,1:NALT,1:NFAM])}

P<-mean(p[1:NINDIV])
MU<-mean(mod[1:NINDIV])

##Estimation of random effect variances
varP<-1/tauP
varF<-1/tauF
varINT1<-1/tauINT1
varINT2<-1/tauINT2
QST <- varP /(varP+ 8*varF)

##HERITABILITE
h2<-(4*varF) / varTRAIT


}"
writeLines(modNegBin, con="mixed_modelNegBi.txt")

###################################################################################
#parametres d'iterations
NECH = 60000
NBURN = 10000
NTHIN = 50


datalist = list(
  Trait = round(tab2[ColTrait[j]][,1]),
  Trait2= round(tab2[ColTrait[j]][,1]),
  Bloc= tab2$Bloc,
  Alt = tab2$Alt,
  Mother = tab2$Mother,
  
  yMean = mean(round(tab2[ColTrait[j]][,1])),
  
  
  NINDIV=sum(table(tab2$Alt[tab2[ColTrait[j]]>=0])),  #NINDIV = Number of individuals
  NBLOC=nlevels(tab2$Bloc),
  NALT=nlevels(tab2$Alt),
  NFAM=nlevels(tab2$Mother)
  
)


#initial    & chains
NRANGE= nlevels(tab2$Bloc)
NPOP=nlevels(tab2$Alt)
NTREE=nlevels(tab2$Mother)

nchain <- 3


initial<-c(mu = 0, a1 = rep(0,NRANGE) , A2= rep(0,NPOP), A3 = rep(0,NTREE) )


init<-list()
for(k in 1:nchain){init[[k]]<-initial}

# Nombre d'échantillons à extraire par chaîne dans le posterior
Nechantillons <- NECH

# échantillons initiaux à éliminer (burn in) par chaîne
burnin <- NBURN#125000#/100 #50000

# pas d'éclaircissage (thinning)
eclair<-NTHIN #100

# quantités (noeuds) dont nous voulons sauvegarder les échantillons
parameters = c("QST", "h2", "varF", "varP")

# compilation et lancement du calcul  
mix <-jags(model.file="mixed_modelNegBi.txt", data=datalist, n.chains=nchain, inits=init, parameters.to.save=parameters,
           n.iter=Nechantillons, n.burnin=burnin, n.thin=eclair)


x11()
post<-as.mcmc(mix)
plotQST_alt(post)
segments(x0=mean(loc.Fst.Tot), x1=mean(loc.Fst.Tot),
         y0=-10000, y1=100000,
         col="black",lwd=2)

polygon(x=c(lwr,lwr, upr, upr),
        y=c(0, 100000, 100000, 0),
        col = alpha("orange", 0.3), lty=0)



mat=c(mix$BUGSoutput[11]$mean$QST, mix$BUGSoutput[11]$mean$h2,
      mix$BUGSoutput[11]$mean$varF,mix$BUGSoutput[11]$mean$varP)

matInf=c(mix$BUGSoutput[10]$summary[1,3],mix$BUGSoutput[10]$summary[3,3],
         mix$BUGSoutput[10]$summary[4,3], mix$BUGSoutput[10]$summary[5,3])

matSup=c(mix$BUGSoutput[10]$summary[1,7],mix$BUGSoutput[10]$summary[3,7],
         mix$BUGSoutput[10]$summary[4,7], mix$BUGSoutput[10]$summary[5,7])

resMat=data.frame(matrix(c(mat, matInf, matSup), ncol=3))
rownames(resMat)<-c("QST", "h2", "varF", "varP")
colnames(resMat)<-c("mean", "Inf", "Sup")
MatNegBi<-rbind(MatNegBi, resMat)

}
write.table(MatNegBi, file="Toul_QST2NegBi_3.csv", sep=";", row.names=T)
