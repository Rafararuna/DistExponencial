########################################
########################################
# Perfis dos clientes
########################################
########################################

# Pacotes:

library(xtable)
library(plotrix)
library(plyr)
library(lattice)
library(HistogramTools)
library(MASS)

op <- options()


# Diagnóstico e envelope para o modelo Poisson
source("https://www.ime.unicamp.br/~cnaber/diag_pois.r")
source("https://www.ime.unicamp.br/~cnaber/envel_pois.r")


# Definindo variaveis

m.dados <- read.table("http://www.ime.unicamp.br/~cnaber/store.dat")
ncli <- m.dados[,1]
ndom <- m.dados[,2]
renda <- m.dados[,3]
idade <- m.dados[,4]
disc <- m.dados[,5]
disl <- m.dados[,6]
n <- length(ncli)

## centralizando nas médias

ndomc  <-  ndom - mean(ndom)
rendac <-  renda - mean(renda)
idadec <- idade - mean(idade)
discc   <- disc - mean(disc)
dislc   <- disl - mean(disl)

ndomcp  <-  ndomc/sd(ndom)
rendacp <-  rendac/sd(renda)
idadecp <- idadec/sd(idade)
disccp   <- discc/sd(disc)
dislcp   <- dislc/sd(disl)


X11()
par(mfrow=c(2,3))
plot(ndom,ncli,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de domicílios",ylab="número de clientes",pch=19) # percebe-se uma tendencia de crescimento, quanto maior o numero de domicilios maior o numero de clientes
plot(renda,ncli,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="renda",ylab="número de clientes",pch=19) # se dermos um corte inicial, nota tambem uma tendencia de crescimento, ou seja, quanto maior a renda maior o numero de clientes
                                                                                                # porem, notamos tbm q, no decorrer no grafico, há obs com renda alta e numero de clietes baixo
plot(idade,ncli,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="idade",ylab="número de clientes",pch=19) # nao nota-se nenhuma tendencia, os dados aperentam ter comportamento aleatorio
plot(disc,ncli,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="distância ao concorrente mais próximo",ylab="número de clientes",pch=19) # nota-se uma tendencia de crescimento, ou seja, quanto maior a distancia ao concorrente proximo maior o numero de clientes
plot(disl,ncli,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="distância à loja",ylab="número de clientes",pch=19) # nota-se uma tendencia de decrescimento, ou seja, quanto maior a distancia à loja menor o numero de clientes


# Medidas-resumo

mcovar <- data.frame(cbind(ndom,renda,idade,disc,disl))
medres <- rbind(rbind(apply(mcovar,2,mean)),rbind(apply(mcovar,2,sd)),100*rbind(apply(mcovar,2,sd))/rbind(apply(mcovar,2,mean)),
                rbind(apply(mcovar,2,quantile,0.5)),rbind(apply(mcovar,2,min)),rbind(apply(mcovar,2,max)))
rownames(medres) <- c("Média","DP","CV(100%)","Mediana","Mínimo","Máximo")
xtable(medres)

### OBS > para que a gente consiga fazer mais comparações entre as variaveis, um alternativa é normalizá-las,
###       ou seja, corrigindo pela media e dividindo pelo desvio-padrao, ai a gente passa a colocar todas as
###       as variaveis com mesma media e desvio-padrao, com media 0 e desvio-padrao 1



# Ajuste do modelo completo

#result<-fit.model<-glm(ncli~ndomc+rendac+idadec+discc+dislc,family=poisson("log"))
result <- fit.model <- glm(ncli~ndomcp+rendacp+idadecp+disccp+dislcp,family=poisson("log"))
summary(result)
xtable(summary(result))

desvioM1 <- deviance(result)
p <- ncol(model.matrix(result))
pvdesvM1 <- 1-pchisq(desvioM1,df=n-p) # p-valor alto, entao podemos assumir que a dist. do modelo é apropriada para a variavel resposta, assumindo a convergencia da função desvio para a dist. Qui-Quadrado
                                      # ou seja, temos um indicativo de o modelo estar bem adequado aos dados

X11()
diagPoisson(fit.model) # comportamento mais dentro do q a gente espera de um modelo bem ajustado
par(mfrow=c(1,1))
envelPoisson(fit.model,"log") # todos pontos estao praticamente dentro da banda

AICM1 <- AIC(fit.model)
BICM1 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaM1 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM1)
mrebeta <- rbind(rebetaM1)
covbetaM1 <- vcov(fit.model)


# Estimativas pontuains e intervalares dos impactos
icov <- exp(rebeta[2:6,1]/c(sd(ndom),sd(renda),sd(idade),sd(disc),sd(disl))) # pegando o exponencial de cada coeficiente e divide pelo desvio padrao
auxcov <-diag(1/c(sd(ndom),sd(renda),sd(idade),sd(disc),sd(disl))) # gradiente (derivando icov)
epicov <- sqrt(diag(auxcov%*%diag(icov)%*%covbetaM1[2:6,2:6]%*%auxcov%*%diag(icov))) # aqui multiplica pela matriz de covariancia e variancia dos coeficientes
liICicov <- icov-ez*epicov # intervalos de confiança das estimativas dos incrementos (METODO DELTA)
lsICicov <- icov+ez*epicov

par(mfrow=c(1,1))
plotCI(icov,li=liICicov,ui=lsICicov,axes=F,xlab="covariável",ylab="impacto no número de clientes",cex=1.2,cex.axis=1.2,pch=19,cex.lab=1.2)
abline(1,0,lwd=2)
axis(2)
axis(1,at=c(1:5),labels=c("ndom","renda","idade","disc","disl"))
### OBS > se ta proximo de um, é um indicativo de q o impacto nao esta sendo relevante dentro da media, caso contrario, ou seja, se esta longe de um, esta sendo relevante

par(mfrow=c(1,3))
plotCI(icov[1],li=liICicov[1],ui=lsICicov[1],axes=F,xlab="covariável",ylab="impacto no número de clientes",cex=1.2,cex.axis=1.2,pch=19,cex.lab=1.2,ylim=c(1,1.001))
abline(1,0,lwd=2)
axis(2)
axis(1,at=c(1),labels=c("ndom"))
plotCI(icov[2],li=liICicov[2],ui=lsICicov[2],axes=F,xlab="covariável",ylab="impacto no número de clientes",cex=1.2,cex.axis=1.2,pch=19,cex.lab=1.2,ylim=c(liICicov[2],1))
abline(1,0,lwd=2)
axis(2)
axis(1,at=c(1),labels=c("renda"))
plotCI(icov[3],li=liICicov[3],ui=lsICicov[3],axes=F,xlab="covariável",ylab="impacto no número de clientes",cex=1.2,cex.axis=1.2,pch=19,cex.lab=1.2)
abline(1,0,lwd=2)
axis(2)
axis(1,at=c(1),labels=c("idade"))
### OBS > ndom: percebe-se um incremento, positivo, acima de 1, ate o limite inferior esta acima de 1
###     > renda: percebe-se um incremendo negativo, abaixo de 1, ate o limite superior ate abaixo de 1
###     > idade: percebe-se um incremendo negativo, abaixo de um 1, ate o limite superior ate abaixo de 1


#plotCI(icov,li=liICicov,ui=lsICicov,axes=F,xlab="covariável",ylab="impacto no número de clientes",cex=1.2,cex.axis=1.2,pch=19,cex.lab=1.2)
#abline(1,0,lwd=2)
#axis(2)
#axis(1,at=c(1:5),labels=c("ndom","renda","idade","disc","disl"))



# Modelo sem a variável sem a idade média

fit.model2 <- glm(ncli~ndomcp+rendacp+disccp+dislcp,family=poisson("log"))

desvioM2 <- deviance(fit.model2)
p <- ncol(model.matrix(fit.model2))
pvdesvM2 <- 1-pchisq(desvioM2,df=n-p) # ainda tem um p-valor alto, ou seja, o modelo esta bem ajustado aos dados

diagPoisson(fit.model2) # nao vemos muita mudança, ainda aparenta estar dentro do esperado para um modelo bem ajustado
par(mfrow=c(1,1))
envelPoisson(fit.model2,"log") # praticamente todos os pontos dentro das bandas

AICM2 <- AIC(fit.model2)
BICM2 <- BIC(fit.model2)


# TRV (para comparar se, de fato, é significativa a introdução da variavel idade no modelo)

lambda <- -2*(as.numeric(logLik(fit.model2))-as.numeric(logLik(fit.model)))
pvalorl <- 1-pchisq(lambda,df=1) # p-valor menor que 0,05, entao a variavel idade é significativa



# Ajuste dos modelos

## Seleção stepwise
fit.model <- result <- glm(ncli~ndomcp+rendacp+idadecp+disccp+dislcp,family=poisson("log"))
fit.model0 <- result0 <-glm(ncli~1,family=poisson("log"))

### Começando do modelo só com intercpeto
stepAIC(fit.model0,scope=list(upper=result),direction=c("both"))

### Começando do modelo completo
stepAIC(fit.model,scope=list(lower=result0),direction=c("both"))

#### OBS > ambos resultaram no modelo com todas as variaveis
####     > como o modelo com a variavel idade apresentou AIC E BIC menores, entao ela deve ficar no modelo 



# Análise preditiva

pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- (pred$fit)
semupred <- (pred$se.fit)
liIC = mupred-ez*semupred
lsIC = mupred+ez*semupred


par(mfrow=c(2,2))
plot(ncli[1:35],cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19,xlab="índice",ylab="número de clientes")
plotCI(mupred[1:35],li=liIC[1:35],ui=lsIC[1:35],pch=17,col=2,add=TRUE)

plot(seq(36,60),ncli[36:60],cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19,xlab="índice",ylim=c(0,40),ylab="número de clientes")
plotCI(seq(36,60),mupred[36:60],li=liIC[36:60],ui=lsIC[36:60],pch=17,col=2,add=TRUE)

plot(seq(61,85),ncli[61:85],cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19,xlab="índice",ylim=c(0,35),ylab="número de clientes")
plotCI(seq(61,85),mupred[61:85],li=liIC[61:85],ui=lsIC[61:85],pch=17,col=2,add=TRUE)

plot(seq(86,110),ncli[86:110],cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19,xlab="índice",ylim=c(0,30),ylab="número de clientes")
plotCI(seq(86,110),mupred[86:110],li=liIC[86:110],ui=lsIC[86:110],pch=17,col=2,add=TRUE)


par(mfrow=c(1,1))
plot(ncli,mupred,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="observado",ylab="predito") # nota-se uma forte associação entre as variaveis (observado e predito)
abline(0,1,lwd=2)
#lines(ncli,liIC,type="p",col=3)
#lines(ncli,lsIC,type="p",col=3)