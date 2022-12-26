############################################
############################################
# Regressão Poisson com variabilidade extra
############################################
############################################


# Pacotes:
library(xtable)
library(plotrix)
library(plyr)
library(lattice)
library(HistogramTools)
library(MASS)
library(e1071)

op <- options()


# Diagnóstico e envelope para o modelo Poisson
source("https://www.ime.unicamp.br/~cnaber/diag_pois.r")
source("https://www.ime.unicamp.br/~cnaber/envel_pois.r")

# Diagnóstico e envelope para o modelo binomial negativo
source("https://www.ime.unicamp.br/~cnaber/diag_nbin.r")
source("https://www.ime.unicamp.br/~cnaber/envel_nbin.r")


# Teste CB=M Binomial Negativa
source("https://www.ime.unicamp.br/~cnaber/TestCBMBinNeg.r")



# Definindo variaveis:

mdados <- read.table("https://www.ime.unicamp.br/~cnaber/quine.dat")

etnia <- factor(mdados[,1])
genero <- factor(mdados[,2])
anoe <- factor(mdados[,3])
anoe <- revalue(anoe,c("F0"="8S","F1"="1A","F2" = "2A","F3"="3A"))
desemp <- factor(mdados[,4])
desemp <- revalue(desemp,c("SL"="insuficiente","AL"="suficiente"))
nfaltas <- cbind(mdados[,5])
n <- length(etnia)



####################
####################
# Análise descritiva
####################
####################


datafalt <- data.frame(etnia,genero,anoe,desemp,nfaltas)

cetnia <- ddply(datafalt,.(etnia),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
cgenero <- ddply(datafalt,.(genero),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
canoe <- ddply(datafalt,.(anoe),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
cdesemp <- ddply(datafalt,.(desemp),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))

xtable(cetnia)
xtable(t(cetnia))
xtable(cgenero)
xtable(t(cgenero))
xtable(canoe)
xtable(t(canoe))
xtable(cdesemp)
xtable(t(cdesemp))


X11()
par(mfrow=c(2,2))
boxplot(nfaltas~etnia,cex=1.2,cex.lab=1.2,cex.main=1.2,names=c("aborígene","não-aborígene"))
boxplot(nfaltas~genero,cex=1.2,cex.lab=1.2,cex.main=1.2,names=c("feminino","masculino"))
boxplot(nfaltas~anoe,cex=1.2,cex.lab=1.2,cex.main=1.2,names=c("8 série","1o ano","2o ano","3o ano"))
boxplot(nfaltas~desemp,cex=1.2,cex.lab=1.2,cex.main=1.2,names=c("suficiente","insuficiente"))


# Médias e variâcias para cada grupo

par(mfrow=c(1,1))
ctudo <-  ddply(datafalt,.(etnia,genero,anoe,desemp),summarise,media=mean(nfaltas),vari=var(nfaltas),n=length(nfaltas))
plot(ctudo$media,ctudo$vari,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="média amostral",ylab="variância amostral",pch=19)
abline(0,1,lwd=2)
### OBS > em geral, as variancias possuem valores bem mais altos que as medias, basta ver as escalas
###     > no inicio os pontos estao proximos da reta, mas dps a variancia explode
###     > isso é um indicativo de que a gente pode ter problema ao modelar sobre a dist. poisson


# Gráficos de perfis

par(mfrow=c(1,1))
ez <- qnorm(0.975)

## etnia x gênero
cetgen <- ddply(datafalt,.(etnia,genero),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plotCI(cetgen$media[cetgen$etnia=="A"],uiw=ez*cetgen$dp[cetgen$etnia=="A"]/sqrt(cetgen$n[cetgen$etnia=="A"]),liw=ez*cetgen$dp[cetgen$etnia=="A"]/sqrt(cetgen$n[cetgen$etnia=="A"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="g???nero",ylab="número de faltas",pch=19,col=1,ylim=c(5,30))
lines(cetgen$media[cetgen$etnia=="A"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("feminino","masculino"),cex.axis=1.2)
plotCI(cetgen$media[cetgen$etnia=="N"],uiw=ez*cetgen$dp[cetgen$etnia=="N"]/sqrt(cetgen$n[cetgen$etnia=="N"]),liw=ez*cetgen$dp[cetgen$etnia=="N"]/sqrt(cetgen$n[cetgen$etnia=="N"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cetgen$media[cetgen$etnia=="N"],lwd=2,col=2)
legend(1.5,28,c("aborígene","não-aborígene"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)
### OBS > no sexo feminino nota-se uma diferença em media entre as etnias, a gente observa que os nativos tem um maior numero de faltas em media
###     > no sexo masculino, a gente ja nao nota essa diferença relevante, os valores das medias se aproximam


## etnia x ano
cetano <- ddply(datafalt,.(etnia,anoe),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plotCI(cetano$media[cetano$etnia=="A"],uiw=ez*cetano$dp[cetano$etnia=="A"]/sqrt(cetano$n[cetano$etnia=="A"]),liw=ez*cetano$dp[cetano$etnia=="A"]/sqrt(cetano$n[cetano$etnia=="A"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=1,ylim=c(5,40))
lines(cetano$media[cetano$etnia=="A"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:4,c("8 série","1o ano","2o ano","3o ano"),cex.axis=1.2)
plotCI(cetano$media[cetano$etnia=="N"],uiw=ez*cetano$dp[cetano$etnia=="N"]/sqrt(cetano$n[cetano$etnia=="N"]),liw=ez*cetano$dp[cetano$etnia=="N"]/sqrt(cetano$n[cetano$etnia=="N"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cetano$media[cetano$etnia=="N"],lwd=2,col=2)
legend(1.2,35,c("aborígene","não-aborígene"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)
### OBS > não tem um comportamento muito específico
###     > nota-se uma discrepancia maior nos anos intermediarios, principalmente no segundo ano (aonde teve uma maior diferença entre a media do numero de faltas dos nativo e nao nativos)


## etnia x desemp
cetades <- ddply(datafalt,.(etnia,desemp),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plotCI(cetades$media[cetades$etnia=="A"],uiw=ez*cetades$dp[cetades$etnia=="A"]/sqrt(cetades$n[cetades$etnia=="A"]),liw=ez*cetades$dp[cetades$etnia=="A"]/sqrt(cetades$n[cetades$etnia=="A"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="desempenho",ylab="número de faltas",pch=19,col=1,ylim=c(5,35))
lines(cetades$media[cetades$etnia=="A"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("suficiente","insuficiente"),cex.axis=1.2)
plotCI(cetades$media[cetades$etnia=="N"],uiw=ez*cetades$dp[cetades$etnia=="N"]/sqrt(cetades$n[cetades$etnia=="N"]),liw=ez*cetades$dp[cetades$etnia=="N"]/sqrt(cetades$n[cetades$etnia=="N"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cetades$media[cetades$etnia=="N"],lwd=2,col=2)
legend(1.2,35,c("aborígene","não-aborígene"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)
### OBS > no desempenho "suficiente" nao parece ter uma diferença relevante no numero medio de faltas entre as etnias
###     > no desempenho "insuficiente" nao parece ter uma diferença relevante no numero medio de faltas entre as etnias


## gênero x ano
cgenano <- ddply(datafalt,.(genero,anoe),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plotCI(cgenano$media[cgenano$genero=="F"],uiw=ez*cgenano$dp[cgenano$genero=="F"]/sqrt(cgenano$n[cgenano$genero=="F"]),liw=ez*cgenano$dp[cgenano$genero=="F"]/sqrt(cgenano$n[cgenano$genero=="F"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=1,ylim=c(5,40))
lines(cgenano$media[cgenano$genero=="F"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:4,c("8 série","1o ano","2o ano","3o ano"),cex.axis=1.2)
plotCI(cgenano$media[cgenano$genero=="M"],uiw=ez*cgenano$dp[cgenano$genero=="M"]/sqrt(cgenano$n[cgenano$genero=="M"]),liw=ez*cgenano$dp[cgenano$genero=="M"]/sqrt(cgenano$n[cgenano$genero=="M"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cgenano$media[cgenano$genero=="M"],lwd=2,col=2)
legend(1.2,35,c("feminino","masculino"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)
### OBS > nota-se uma diferença maior no 3 ano
###     > interessante notar que nos dois primeiros anos a media de faltas do sexo feminio pe maior, e nos ultimos dois anos isso inverte
###     > nos tres primeiros anos a diferença da media de faltas entre os sexos não é tao discriminante


# gênero x desemp
cgendes <- ddply(datafalt,.(genero,desemp),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plotCI(cgendes$media[cgendes$genero=="F"],uiw=ez*cgendes$dp[cgendes$genero=="F"]/sqrt(cgendes$n[cgendes$genero=="F"]),liw=ez*cgendes$dp[cgendes$genero=="F"]/sqrt(cgendes$n[cgendes$genero=="F"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="desempenho",ylab="número de faltas",pch=19,col=1,ylim=c(10,30))
lines(cgendes$media[cgendes$genero=="F"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("suficiente","insuficiente"),cex.axis=1.2)
plotCI(cgendes$media[cgendes$genero=="M"],uiw=ez*cgendes$dp[cgendes$genero=="M"]/sqrt(cgendes$n[cgendes$genero=="M"]),liw=ez*cgendes$dp[cgendes$genero=="M"]/sqrt(cgendes$n[cgendes$genero=="M"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cgendes$media[cgendes$genero=="M"],lwd=2,col=2)
legend(1.2,30,c("feminino","masculino"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)
### OBS > no desempenho "suficiente" há uma diferença maior no numero medio de faltas entre os sexos do que no desempenho "insuficiente"
###     > no sexo feminino há um aumento mais relevante na mudança de desempenho do que no sexo masculino que quase chega a ser uma reta


## desempenho x ano
cdesano <- ddply(datafalt,.(desemp,anoe),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plotCI(cdesano$media[cdesano$desemp=="suficiente"],uiw=ez*cdesano$dp[cdesano$desemp=="suficiente"]/sqrt(cdesano$n[cdesano$desemp=="suficiente"]),liw=ez*cdesano$dp[cdesano$desemp=="suficiente"]/sqrt(cdesano$n[cdesano$desemp=="suficiente"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=1,ylim=c(0,35))
lines(cdesano$media[cdesano$desemp=="suficiente"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:4,c("8 série","1o ano","2o ano","3o ano"),cex.axis=1.2)
plotCI(cdesano$media[cdesano$desemp=="insuficiente"],uiw=ez*cdesano$dp[cdesano$desemp=="insuficiente"]/sqrt(cdesano$n[cdesano$desemp=="insuficiente"]),liw=ez*cdesano$dp[cdesano$desemp=="insuficiente"]/sqrt(cdesano$n[cdesano$desemp=="insuficiente"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cdesano$media[cdesano$desemp=="insuficiente"],lwd=2,col=2)
legend(1.2,35,c("suficiente","insuficiente"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)



# Perfis de terceira ordem

## (etnia x gênero) x desempenho
par(mfrow=c(1,2))
cetgendes <- ddply(datafalt,.(etnia,genero,desemp),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))

### Suficiente
plotCI(cetgendes$media[cetgendes$etnia=="A" & cetgendes$desemp=="suficiente"],uiw=ez*cetgendes$dp[cetgendes$etnia=="A"& cetgendes$desemp=="suficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="A"& cetgendes$desemp=="suficiente"]),liw=ez*cetgendes$dp[cetgendes$etnia=="A"& cetgendes$desemp=="suficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="A"& cetgendes$desemp=="suficiente"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="g???nero",ylab="número de faltas",pch=19,col=1,ylim=c(5,38))
title("Desempenho suficiente")
lines(cetgendes$media[cetgendes$etnia=="A"& cetgendes$desemp=="suficiente"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("feminino","masculino"),cex.axis=1.2)
plotCI(cetgendes$media[cetgendes$etnia=="N" & cetgendes$desemp=="suficiente"],uiw=ez*cetgendes$dp[cetgendes$etnia=="N" & cetgendes$desemp=="suficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="N" & cetgendes$desemp=="suficiente"]),liw=ez*cetgendes$dp[cetgendes$etnia=="N" & cetgendes$desemp=="suficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="N" & cetgendes$desemp=="suficiente"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cetgendes$media[cetgendes$etnia=="N" & cetgendes$desemp=="suficiente"],lwd=2,col=2)
legend(1.2,28,c("aborígene","não-aborígene"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)

### Insuficiente
plotCI(cetgendes$media[cetgendes$etnia=="A" & cetgendes$desemp=="insuficiente"],uiw=ez*cetgendes$dp[cetgendes$etnia=="A"& cetgendes$desemp=="insuficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="A"& cetgendes$desemp=="insuficiente"]),liw=ez*cetgendes$dp[cetgendes$etnia=="A"& cetgendes$desemp=="insuficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="A"& cetgendes$desemp=="insuficiente"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="g???nero",ylab="número de faltas",pch=19,col=1,ylim=c(5,38))
title("Desempenho insuficiente")
lines(cetgendes$media[cetgendes$etnia=="A"& cetgendes$desemp=="insuficiente"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("feminino","masculino"),cex.axis=1.2)
plotCI(cetgendes$media[cetgendes$etnia=="N" & cetgendes$desemp=="insuficiente"],uiw=ez*cetgendes$dp[cetgendes$etnia=="N" & cetgendes$desemp=="insuficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="N" & cetgendes$desemp=="insuficiente"]),liw=ez*cetgendes$dp[cetgendes$etnia=="N" & cetgendes$desemp=="insuficiente"]/sqrt(cetgendes$n[cetgendes$etnia=="N" & cetgendes$desemp=="insuficiente"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,pch=17,col=2,add=TRUE)
lines(cetgendes$media[cetgendes$etnia=="N" & cetgendes$desemp=="insuficiente"],lwd=2,col=2)
legend(1.2,32,c("aborígene","não-aborígene"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)

#### OBS > em relação ao desempenho "suficiente", nota-se que ha uma diferença muito pequena na media do numero de faltas no sexo feminino, e uma diferença ate q significativa na media do numero de faltas no sexo masculino
####     > em relação ao desempenho "insuficiente", nota-se que ha uma diferença muito grande na media do numero de faltas no sexo feminino, e uma diferença muito pequena na media do numero de faltas no sexo masculino
####     > esse tipo de analise é util para a gente verificar possiveis interações no modelo, ou seja, entre etnia e sexo aparenta ter uma interação



# Modelo completo

resultM1 <- fit.model <- glm(nfaltas~etnia+genero+anoe+desemp,family=poisson("log"))
summary(resultM1)
xtable(summary(resultM1))

desvioM1 <- deviance(resultM1)
p <- ncol(model.matrix(resultM1))
pvdesvM1 <- 1-pchisq(desvioM1,df=n-p) # percebe-se um p-valor muito baixo, ou seja, se de fato ocorre nessa convergencia para a Qui-quadrado, o ajuste do modelo 
                                      # esta sendo rejeitado, ou seja, tem-se um indicativo de que os dados não tão se adequando ao modelo

diagPoisson(fit.model) # aqui a gente nota o que foi dito sobre a analise de desvio 
                       # percebe-se muitos outliers nos graficos de cima
                       # o histograma ate aparenta estar simetrico em torno do zero, mas apresenta uma cauda bem pesada
                       # nota-se muitos pontos longe da reta em todo o grafico
                       # percebe-se que há uma variabilidade que nap eta sendo captada pelo modelo
                       # observa-se um problema de superdispersao
#abline(0,1)
par(mfrow=c(1,1))
envelPoisson(fit.model,"log") # praticamente todos os pontos fora das bandas

AICM1 <- AIC(fit.model)
BICM1 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(resultM1))$coef
rebetaM1 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM1)
mrebeta <- rbind(rebetaM1)
covbetaM1 <- vcov(fit.model)



# Ajuste dos modelos

## Seleção stepwise

fit.model <- result <- glm(nfaltas~etnia+genero+anoe+desemp,family=poisson("log"))
fit.model0 <- result0 <- glm(nfaltas~1,family=poisson("log"))

### Começando do modelo só com intercpeto

stepAIC(fit.model0,scope=list(upper=result),direction=c("both"))

### Começando do modelo completo

stepAIC(fit.model,scope=list(lower=result0),direction=c("both"))



# Como o modelo Poisson nao se ajustou bem aos dados, vamos implementar o modelo binomial negativo para ver se apresenta melhora

# Binomial Negativa

# modelo com somente os fatores principais

resultBNM1 <- fit.model <- glm.nb(nfaltas~etnia+genero+anoe+desemp,link="log")
summary(resultBNM1)
xtable(summary(resultBNM1))

desvioBNM1 <-deviance(resultBNM1)
p <- ncol(model.matrix(resultBNM1))
pvdesvBNM1 <- 1-pchisq(desvioBNM1,df=n-p) # nota-se um p-valor de 0,047, muito proximo de 0,05, entao ficamos ali em cima do muro pra decidir se rejeita ou nao h0 ; porem, a gente nota sim que teve uma melhora do ajuste significativa

diagnbin(resultBNM1) # nota-se uma melhora nos graficos de cima ; primeiramente em relação a escala, que diminuiu, antes tava de -5 a 5, agr ta de -3 a 3 ; e tbm notamos muito menos pontos fora da banda
                     # histograma aparenta ser mais simetrico em torno do zero, seguindo um dist. Normal
#abline(0,1)
#par(mfrow=c(1,1))
ligacaonbin <- "log"
par(mfrow=c(1,1))
envelnbin(resultBNM1) # muito melhor em relação ao modelo anterior, quase todos os pontos envelopados, ou seja, dentro da banda
                      # essa diferenã mostra como o modelo de poisson restringia a questao da variabilidade quando a gente modelo, pelo mesmo parametro, a media e a variancia

AICBNM1 <- AIC(resultBNM1)
BICBNM1 <- BIC(resultBNM1)

ez < -qnorm(0.975)
rebeta <- (summary(resultBNM1))$coef
rebetaBNM1 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
rphiBNM1 <- cbind(resultBNM1$theta,resultBNM1$SE.theta,resultBNM1$theta-ez*resultBNM1$SE.theta,resultBNM1$theta+ez*resultBNM1$SE.theta,0,0)
xtable(rbind(rebetaBNM1,rphiBNM1))
#xtable(rebetaBNM1)

mrebeta <- rbind(rebetaBNM1)
covbetaBNM1 <- vcov(resultBNM1)



# Teste para retirar os fatores

mC <- rbind(cbind(0,0,1,0,0,0,0),cbind(0,0,0,0,0,0,1))
mM <- rbind(0,0)

testeFCBMBinNeg(resultBNM1,mC,mM)
### OBS > nao rejeitamos a hipotese nula de que os parametros de genero e desempenho são iguais a zero, ou seja, eles nao sao significativos
###     > portanto, ajusta-se um modelo sem esses dois fatores, o resultBNM2, que sera implementadi mais pra frente


# Comparação do modelo binomial 1 com o modelo de Poisson

## Estatísticas de ajuste

AICBIC <- rbind(cbind(AICM1,AICBNM1),cbind(BICM1,BICBNM1),cbind(desvioM1,desvioBNM1),cbind(pvdesvM1,pvdesvBNM1))
colnames(AICBIC) <- c("Poisson","Binomial negativo")
rownames(AICBIC) <- c("AIC","BIC","Desvio","p-valor desvio")
xtable(t(AICBIC))


## Comparação entre gráficos de diagnóstico

#par(mfrow=c(1,2))
diagPoisson(resultM1)
abline(0,1)
diagnbin(resultBNM1)


## Comparação entre os envelopes

par(mfrow=c(1,2))
envelPoisson(resultM1,"log")
title("Modelo de Poisson",cex=1.2)
envelnbin(resultBNM1)
title("Modelo binomial negativo",cex=1.2)



# Modelo reduzido somente com etnia e anoe (modelo binomial 2)

resultBNM2 <- fit.model <- glm.nb(nfaltas~etnia+anoe,link="log")
summary(resultBNM2)
xtable(summary(resultBNM2))

desvioBNM2 <- deviance(resultBNM2)
p <- ncol(model.matrix(resultBNM2))
pvdesvBNM2 <- 1-pchisq(desvioBNM2,df=n-p) # nota-se que o ajuste melhorou em relação ao modelo binomial 1, pois o p-valor aumentou,
                                          # ou seja, menos chance de rejeitar h0, de rejeitar o ajuste

diagnbin(resultBNM2) # nao há tanta diferença
#abline(0,1)
#par(mfrow=c(1,1))
ligacaonbin <- "log"
par(mfrow=c(1,1))
envelnbin(resultBNM2) # muito similiar ao modleo binomial 1

AICBNM2 <- AIC(resultBNM2)
BICBNM2 <- BIC(resultBNM2)

ez <- qnorm(0.975)
rebeta <- (summary(resultBNM2))$coef
rebetaBNM2 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
rphiBNM2 <- cbind(resultBNM2$theta,resultBNM2$SE.theta,resultBNM2$theta-ez*resultBNM2$SE.theta,resultBNM2$theta+ez*resultBNM2$SE.theta,0,0)
xtable(rbind(rebetaBNM2,rphiBNM2))
#xtable(rebetaBNM2)

mrebeta <- rbind(rebetaBNM2)
covbetaBNM2 <- vcov(resultBNM2)



# Notamos que os coeficientes de ano escolar nao foram significativos, mas amos construir um modelo com interação entre eles, antes de tira-los do modelo

# Modelo reduzido com etnia e anoe com interação

resultBNM3 <- fit.model <- glm.nb(nfaltas~etnia+anoe+etnia*anoe,link="log")
mX <- model.matrix(resultBNM3)
summary(resultBNM3)
xtable(summary(resultBNM3))
### OBS > nota-se que a maioria das interações foram significativas, por isso que alguns coeficientes individuais, como o beta2, que antes era significativo
###       deixou de ser significativo, pois as interações puxaram toda essa informação

desvioBNM3 <- deviance(resultBNM3)
p <- ncol(model.matrix(resultBNM3))
pvdesvBNM3 <- 1-pchisq(desvioBNM3,df=n-p) # o ajuste piorou em relação o modelo anterior, o p-valor diminuiu e, se formos rigorosos, ate rejeita a hipotese nula

diagnbin(resultBNM3) # nao teve tanta mudança discrepante em relação ao modelo anterior
#abline(0,1)
#par(mfrow=c(1,1))
ligacaonbin <- "log"
par(mfrow=c(1,1))
envelnbin(resultBNM3) # nao teve muita mudança

AICBNM3 <- AIC(resultBNM3)
BICBNM3 <- BIC(resultBNM3)

ez <- qnorm(0.975)
rebeta <- (summary(resultBNM3))$coef
rebetaBNM3 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
rphiBNM3 <- cbind(resultBNM3$theta,resultBNM3$SE.theta,resultBNM3$theta-ez*resultBNM3$SE.theta,resultBNM3$theta+ez*resultBNM3$SE.theta,0,0)
xtable(rbind(rebetaBNM3,rphiBNM3))
#xtable(rebetaBNM3)

mrebeta <- rbind(rebetaBNM3)
covbetaBNM3 <- vcov(resultBNM3)



# Agora vamos verificar se tem como reduzir mais o modelo, se tem como agrupar esses anos escolares (como o 1 ano e o 2 ano isoladamente, e a sua respectiva interação)

# Teste para redução

mC <- rbind(cbind(0,1,0,0,0,0,0,0),cbind(0,0,1,0,0,0,0,0),cbind(0,0,0,0,1,0,0,0),cbind(0,0,0,0,0,0,0,1))
mM <- rbind(0,0,0,0)

testeFCBMBinNeg(resultBNM3,mC,mM)
### OBS > p-valor alto, não rejeita h0, ou seja, esses fatores nao sao signficativos, ou seja, vamos ajustar um modelo sem eles



# Modelo final 

resultBNM4 <- fit.model <- glm.nb(nfaltas~mX[,4]+mX[,6]+mX[,7],link="log")
mXBNM4 <- model.matrix(resultBNM4)
summary(resultBNM4)
xtable(summary(resultBNM4))

desvioBNM4 <- deviance(resultBNM4)
p <- ncol(model.matrix(resultBNM4))
pvdesvBNM4 <- 1-pchisq(desvioBNM4,df=n-p) # p-valor maior comparando com todos os modelos anteriores, ou seja, esse foi o q apresentou melhor ajuste

diagnbin(resultBNM4) # nao ha nenhuma mudança discrepante, mesmos comentarios dos modelos anteriores
#abline(0,1)
#par(mfrow=c(1,1))
ligacaonbin <- "log"
par(mfrow=c(1,1))
envelnbin(resultBNM4) # nao ha nenhuma mudança discrepante, mesmos comentarios dos modelos anteriores

AICBNM4 <- AIC(resultBNM4)
BICBNM4 <- BIC(resultBNM4)

ez <- qnorm(0.975)
rebeta <- (summary(resultBNM4))$coef
# exp(unique(mXBNM4)%*%rebeta[,1])
rebetaBNM4 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
rphiBNM4 <- cbind(resultBNM4$theta,resultBNM4$SE.theta,resultBNM4$theta-ez*resultBNM4$SE.theta,resultBNM4$theta+ez*resultBNM4$SE.theta,0,0)
xtable(rbind(rebetaBNM4,rphiBNM4))
#xtable(rebetaBNM4)

mrebeta <- rbind(rebetaBNM4)
covbetaBNM4 <- vcov(resultBNM4)


# médias preditas (comparanto somente com os fatores retidos)

predBNM4 <- predict(resultBNM4,type=c("response"),se.fit = TRUE)
mupredBNM4 <- unique(predBNM4$fit)
semupredBNM4 <- unique(predBNM4$se.fit)
liICBNM4 = mupredBNM4-ez*semupredBNM4
lsICBNM4 = mupredBNM4+ez*semupredBNM4



# Como o betas 22 e 23 estão muito proximos, vamos testar se devemos tira-los do modelo

# Redução do modelo

mC <- cbind(0,0,1,-1)
mM <- rbind(0)

testeFCBMBinNeg(resultBNM4,mC,mM)
### OBS > p-valor alto, não rejeita h0, ou seja, esses fatores nao sao signficativos, ou seja, poderiamos reduzir ainda mais esse modelo



# médias preditas (com todos os fatores)

## etnia x ano
par(mfrow=c(1,1))
cetano<-ddply(datafalt,.(etnia,anoe),summarise,media=mean(nfaltas),dp=sqrt(var(nfaltas)),vari=var(nfaltas),ca=skewness(nfaltas),minimo=min(nfaltas),maximo=max(nfaltas),cv=100*((sqrt(var(nfaltas))/mean(nfaltas))),n=length(nfaltas))
plot(cetano$media[cetano$etnia=="A"],axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=1,ylim=c(5,45))
plotCI(c(mupredBNM4[1],mupredBNM4[1],mupredBNM4[2],mupredBNM4[1]),li=c(liICBNM4[1],liICBNM4[1],liICBNM4[2],liICBNM4[1]),ui=c(lsICBNM4[1],lsICBNM4[1],lsICBNM4[2],lsICBNM4[1]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=1,add=TRUE)
lines(c(mupredBNM4[1],mupredBNM4[1],mupredBNM4[2],mupredBNM4[1]),lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:4,c("8 série","1o ano","2o ano","3o ano"),cex.axis=1.2)

lines(cetano$media[cetano$etnia=="N"],type="p",lwd=2,pch=17,cex.lab=1.2,cex.axis=1.2,cex=1.2,col=2)
plotCI(c(mupredBNM4[1],mupredBNM4[3],mupredBNM4[4],mupredBNM4[1]),li=c(liICBNM4[1],liICBNM4[3],liICBNM4[4],liICBNM4[1]),ui=c(lsICBNM4[1],lsICBNM4[3],lsICBNM4[4],lsICBNM4[1]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=17,col=2,add=TRUE)
lines(c(mupredBNM4[1],mupredBNM4[3],mupredBNM4[4],mupredBNM4[1]),lwd=2,col=2)
legend(1.2,35,c("aborígene - observado","não-aborígene - observado","aborígene - predito","não-aborígene - predito"),col=c(1,2,1,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),bty="n",cex=1.2)
### OBS > notamos valores preditos bem proximos dos observados, com exceção da categoria 8 série, aonde o aborigeno observado ja apresenta uma certa dostancia em relação ao aborigeno predito,
###       assim como pro nao-aborigeno, porem com uma distancia menor ; isso tbm acontece na categoria do terceiro ano

## Ano escolar: 8 série
plot(ctudo$media[ctudo$anoe=="8S"],pch=19,cex=1.2,cex.lab=1.2,ylab="média",xlab="grupo",axes=FALSE)
axis(2,cex.axis=1.2)
axis(1,1:8,c("AFS","AFI","AMS","AMI","NFS","NFI","NMS","NMI"),cex.axis=1.2)
plotCI(c(mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1]),
       li=c(liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1]),
       ui=c(lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1]),
       axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=2,add=TRUE)


## Ano escolar: 1 ano
plot(ctudo$media[ctudo$anoe=="1A"],pch=19,cex=1.2,cex.lab=1.2,ylab="média",xlab="grupo",axes=FALSE,ylim=c(0,25))
axis(2,cex.axis=1.2)
axis(1,1:8,c("AFS","AFI","AMS","AMI","NFS","NFI","NMS","NMI"),cex.axis=1.2)
plotCI(c(mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[3],mupredBNM4[3],mupredBNM4[3],mupredBNM4[3]),
       li=c(liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[3],liICBNM4[3],liICBNM4[3],liICBNM4[3]),
       ui=c(lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[3],lsICBNM4[3],lsICBNM4[3],lsICBNM4[3]),
       axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=2,add=TRUE)


## Ano escolar: 2 ano
plot(ctudo$media[ctudo$anoe=="2A"],pch=19,cex=1.2,cex.lab=1.2,ylab="média",xlab="grupo",axes=FALSE,ylim=c(0,45))
axis(2,cex.axis=1.2)
axis(1,1:8,c("AFS","AFI","AMS","AMI","NFS","NFI","NMS","NMI"),cex.axis=1.2)
plotCI(c(mupredBNM4[2],mupredBNM4[2],mupredBNM4[2],mupredBNM4[2],mupredBNM4[4],mupredBNM4[4],mupredBNM4[4],mupredBNM4[4]),
       li=c(liICBNM4[2],liICBNM4[2],liICBNM4[2],liICBNM4[2],liICBNM4[4],liICBNM4[4],liICBNM4[4],liICBNM4[4]),
       ui=c(lsICBNM4[2],lsICBNM4[2],lsICBNM4[2],lsICBNM4[2],lsICBNM4[4],lsICBNM4[4],lsICBNM4[4],lsICBNM4[4]),
       axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=2,add=TRUE)


## Ano escolar: 3 ano
plot(ctudo$media[ctudo$anoe=="3A"],pch=19,cex=1.2,cex.lab=1.2,ylab="média",xlab="grupo",axes=FALSE,ylim=c(10,30))
axis(2,cex.axis=1.2)
axis(1,1:8,c("AFS","AFI","AMS","AMI","NFS","NFI","NMS","NMI"),cex.axis=1.2)
plotCI(c(mupredBNM4[1],mupredBNM4[1],mupredBNM4[1],mupredBNM4[1]),
       li=c(liICBNM4[1],liICBNM4[1],liICBNM4[1],liICBNM4[1]),
       ui=c(lsICBNM4[1],lsICBNM4[1],lsICBNM4[1],lsICBNM4[1]),
       axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="ano",ylab="número de faltas",pch=19,col=2,add=TRUE)