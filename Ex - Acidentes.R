########################################
########################################
# Número de acidentes
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

# Diagnóstico e envelope para o modelo binomial negativo
source("https://www.ime.unicamp.br/~cnaber/diag_nbin.r")
source("https://www.ime.unicamp.br/~cnaber/envel_nbin.r")



# Número de acidades - 1961

v.x <- c(29,32,20,42,39,25,40,22,40,27,39,8,28,24,21,34,21,21,17,11,18,47,15,15,27,35,36,17,21,37,21,15,20,24,32,25,26,24,15,25,34,21,30)


# Número de acidentes - 1962

v.y <- c(17,17,15,21,26,14,23,12,25,17,16,15,16,7,9,26,15,9,20,11,16,41,12,13,15,25,25,22,13,29,25,12,24,9,17,16,17,16,10,17,22,24,25)


# Total

numac <- c(v.x,v.y) # numero de acidades dos dois anos
n <- length(v.x) # tamanho da amostra - 1961
m <- length(v.y) # tamanho da amostra - 1961
v.n <- c(n,m) 
nt <- n+m # tamanho total da amostra
ano <- c(rep(1961,n),rep(1962,m)) 
anofac <- factor(ano,levels=c("1961","1962")) # fator dos anos


# Medidas descritivas

md.x <- cbind(mean(v.x),var(v.x),sqrt(var(v.x)),100*sqrt(var(v.x))/mean(v.x),min(v.x),as.numeric(quantile(v.x,0.50)),max(v.x))
md.y <- cbind(mean(v.y),var(v.y),sqrt(var(v.y)),100*sqrt(var(v.y))/mean(v.y),min(v.y),as.numeric(quantile(v.y,0.50)),max(v.y))
xtable(rbind(md.x,md.y)) # media, variancia, dp, cv, minimo, mediana, maximo

m.X <- rbind(matrix(0,n,1),matrix(1,m,1))
boxplot(rbind(cbind(v.x),cbind(v.y))~m.X,cex=1.3,cex.lab=1.3,names=c("1961","1962"))

par(mfrow=c(2,1))
a1<-hist(v.x,xlab="número de acidentes (1961)",ylab="frequência relativa",main="",cex=1.2,cex.lab=1.2,cex.axis=1.2,probability=TRUE)
a2<-hist(v.y,xlab="número de acidentes (1962)",ylab="frequência relativa",main="",cex=1.2,cex.lab=1.2,cex.axis=1.2,probability=TRUE)
a1$density <- a1$counts/sum(a1$counts)*100
plot(a1,freq=F,main="",ylab="percentual",xlab="número de acidentes (1961)",cex=1.3,cex.axis=1.3,cex.lab=1.3)
a2$density <- a2$counts/sum(a2$counts)*100
plot(a2,freq=F,main="",ylab="percentual",xlab="número de acidentes (1962)",cex=1.3,cex.axis=1.3,cex.lab=1.3)



# Inferência

resultM1 <- fit.model <- result <- glm(numac~anofac,family=poisson(link="log"))
summary(result)
xtable(summary(result))

desvioM1 <- deviance(result) # desvio do modelo
p <- ncol(model.matrix(result)) # numero de parametros
pvdesvM1 <- 1-pchisq(desvioM1,df=nt-p) # estatística qui- quadrado, assumindo a convergencia da função desvio para a dist. Qui-quadrado
                                       # percebe-se um p-valor muito baixo, ou seja, se de fato ocorre nessa convergencia para a Qui-quadrado, o ajuste do modelo 
                                       # esta sendo rejeitado, ou seja, tem-se um indicativo de que os dados não tão se adequando ao modelo

diagPoisson(fit.model) # percebe-se alguns pontos discrepantes
                       # pelo histograma, parece ter uma dist. simetrica em torno do zero
                       # acerca do ultimo grafica, nota-se que as obs foram agrupadas, e o ideal é que elas estivessem proximo da diagonal
abline(0,1)
par(mfrow=c(1,1))
envelPoisson(fit.model,"log") # percebe-se que realmente o modelo n se ajustou bem aos dados, pois quase todos os pontos estao fora das bandas, ou seja, nao tem um ajuste adequado do modelo

AICM1 <- AIC(fit.model)
BICM1 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaM1 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM1)
mrebeta <- rbind(rebetaM1)
covbetaM1 <- vcov(fit.model)

pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredM1 <- unique(pred$fit)
semupredM1 <- unique(pred$se.fit)
liICM1 = mupredM1-ez*semupredM1
lsICM1 = mupredM1+ez*semupredM1

resultpred <- cbind(c(mupredM1),c(semupredM1),c(liICM1),c(lsICM1))
xtable(resultpred)



# Análise preditiva (gerando a resposta a partir da media estimada)

ypred <- rpois(nt,lambda=mupredM1)

par(mfrow=c(1,2))
boxplot(rbind(cbind(v.x),cbind(v.y))~m.X,cex=1.3,cex.lab=1.3,names=c("1961","1962"),ylim=c(0,60))
boxplot(rbind(cbind(ypred[ano==1961]),cbind(ypred[ano==1962]))~m.X,cex=1.3,cex.lab=1.3,names=c("1961","1962"),add=TRUE,border=2)
legend(1.5,60,c("observado","predito"),col=c(1,2),lty=c(1,1),bty="n",cex=1.3)
### OBS > nota-se que os boxplots preditos é muito mais homogeneo do que os observados, ou seja, apresentam uma variabilidade menor, em ambos os anos
###     > nota-se uma variabilidade predito bem abaixo do observado 
###     > entao, provavelmente, o problema da falta de ajuste desse modelo, é a superdispersão


par(mfrow=c(2,1))
a1<-hist(v.x,xlab="número de acidentes (1961)",ylab="frequência relativa", main="Sobre os dados observados",cex=1.2,cex.lab=1.2,cex.axis=1.2,probability=TRUE)
a2<-hist(v.y,xlab="número de acidentes (1962)",ylab="frequência relativa",main="",cex=1.2,cex.lab=1.2,cex.axis=1.2,probability=TRUE)
a1pred<-hist(ypred[ano==1961],xlab="número de acidentes (1961)",ylab="frequência relativa",main="Sobre os dados preditos",cex=1.2,cex.lab=1.2,cex.axis=1.2,probability=TRUE)
a2pred<-hist(ypred[ano==1962],xlab="número de acidentes (1962)",ylab="frequência relativa",main="",cex=1.2,cex.lab=1.2,cex.axis=1.2,probability=TRUE)

PlotRelativeFrequency (a1,main="",ylab="percentual",xlab="número de acidentes (1961)",cex=1.3,cex.axis=1.3,cex.lab=1.3,ylim=c(0,0.40))
PlotRelativeFrequency (a1pred,main="",ylab="percentual",xlab="número de acidentes (1961)",cex=1.3,cex.axis=1.3,cex.lab=1.3,add=TRUE,border=2)
legend(30,0.50,c("observado","predito"),col=c(1,2),lty=c(1,1),bty="n",cex=1.3)
PlotRelativeFrequency (a2,main="",ylab="percentual",xlab="número de acidentes (1962)",cex=1.3,cex.axis=1.3,cex.lab=1.3,ylim=c(0,0.50))
PlotRelativeFrequency (a2pred,main="",ylab="percentual",xlab="número de acidentes (1962)",cex=1.3,cex.axis=1.3,cex.lab=1.3,add=TRUE,border=2)
legend(25,0.60,c("observado","predito"),col=c(1,2),lty=c(1,1),bty="n",cex=1.3)
### OBS > as valores preditos ficam mais proximos dos valores observados no de 1962 do que no ano de 1961
###     > nota-se tbm quem os valores observados estao mais disperos, mais distribuidos ao longo do histograma, enquando que os valores preditos estao mais concentrados em um regiao   

### OBS > observamos que o modelo nao se ajustou bem aos dados, o que pode ser causado pela superdispersão
###     > a superdispersão por sua vez pode provocar erros-padrao subestimados, resultando em intervalos de confiança menores do que eles seriam



# Outras funções de ligação (envelopes)

fit.model <- glm(numac~anofac,family=poisson(link="identity"))

par(mfrow=c(1,1))
envelPoisson(fit.model,"identity") # tambem nao esta bom, quase todos os pontos estao fora da banda

fit.model <- glm(numac~anofac,family=poisson(link="sqrt"))
par(mfrow=c(1,1))
envelPoisson(fit.model,"sqrt") # tambem nao esta bom, quase todos os pontos estao fora da banda

### OBS > ou seja, o problema não é a função de ligação, e sim a superdispersão



# Para tentar corrigir essa superdispersao, vamos implementar o modelo binomaial negativo

# Modelo Binomial negativo

resultBN <- fit.model <- glm.nb(numac~anofac,link="log")
summary(resultBN)
xtable(summary(resultBN))

desvioM2 <- deviance(resultBN)
p <- ncol(model.matrix(resultBN))
pvdesvM2 <- 1-pchisq(desvioM2,df=nt-p) # aqui temos um p-valor alto, ou seja, não rejeita h0, ou seja, nao rejeita o ajuste,
                                       # assim, ja notamos um avanço

diagnbin(fit.model) # ja notamos uma melhora em relaçao aos dois grafcos de cima, que antes tinham varios pontos do lado de fora, e agr tem-se muito poucos, apenas 4
                    # em relação aos demais graficos n teve tanta diferença
#abline(0,1)
par(mfrow=c(1,1))
ligacaonbin <- "log"
envelnbin(fit.model) # a maioria dos pontos estao dentro da banda, diferente do que acontecia com o modelo anterior, ou seja, esse modelo se adequa melhor aos dados,
                     # sendo uma mudança muito significativa em relação à qualidade do ajuste

AICM2 <- AIC(fit.model)
BICM2 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(resultBN))$coef
rebetaM2 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
rphiM2 <- cbind(resultBN$theta,resultBN$SE.theta,resultBN$theta-ez*resultBN$SE.theta,resultBN$theta+ez*resultBN$SE.theta,0,0)
xtable(rebetaM2)
mrebeta <- rbind(rebetaM1,rbind(rebetaM2,rphiM2))
xtable(mrebeta)
### OBS > interessante notar q em relação aos valores da estimativas, ao erro-padrao, aos IC's e aos testes, nao mudou muita, coisa, nao teve nenhuma diferença discrepante


# Estimativas

xtable(mrebeta)
covbetaM2 <- vcov(fit.model)


# Estatísticas de ajuste

AICBIC <- rbind(cbind(AICM1,AICM2),cbind(BICM1,BICM2),cbind(desvioM1,desvioM2),cbind(pvdesvM1,pvdesvM2))
colnames(AICBIC) <- c("Poisson","Binomial negativo")
rownames(AICBIC) <- c("AIC","BIC","Desvio","p-valor desvio")
xtable(t(AICBIC))


# Comparação entre gráficos de diagnóstico

par(mfrow=c(1,2))
diagPoisson(resultM1)
abline(0,1)
diagnbin(resultBN)


# Comparação entre os envelopes

par(mfrow=c(1,2))
envelPoisson(resultM1,"log")
title("Modelo de Poisson",cex=1.2)
envelnbin(resultBN)
title("Modelo binomial negativo",cex=1.2)
### OBS > aqui fica claro que o modelo BN tem um ajuste melhor, e confirma que o verdadeiro problema era realmente a superdispersao,
###       pois a gente estava sendo muito restritivo ao modelar a media e a variancia pelo mesmo parametro


pred <- predict(resultBN,type=c("response"),se.fit = TRUE)
mupredM2 <- unique(pred$fit)
semupredM2 <- unique(pred$se.fit)
liICM1 = mupredM2-ez*semupredM2
lsICM1 = mupredM2+ez*semupredM2

resultpred <- cbind(c(mupredM1),c(semupredM1),c(liICM1),c(lsICM1))
xtable(resultpred)


# Análise preditiva

fi <- resultBN$theta
ypredBN <- rnegbin(nt,mupredM2,fi) # gerando os valores preditos da binomial negativa

# Comparação entre os modelos 

par(mfrow=c(1,2))
boxplot(cbind(v.x,ypred[anofac=="1961"],ypredBN[anofac=="1961"]),cex=1.3,cex.lab=1.3,names=c("Observado","Poisson","Binomial Negativo"),ylim=c(0,60))
title("Ano de 1961",cex=1.2)
boxplot(cbind(v.y,ypred[anofac=="1962"],ypredBN[anofac=="1962"]),cex=1.3,cex.lab=1.3,names=c("Observado","Poisson","Binomial Negativo"),ylim=c(0,60))
title("Ano de 1962",cex=1.2)