######################################
######################################
#          Regressão de Poisson
######################################
######################################
######################################


# Pacotes:
library(xtable)
library(plotrix)
library(plyr)
op <- options()


# Sumário com alguns comandos relativos à função glm "result<- glm(y~x,family=nomedadistribuição(link="nomedolink"))",
# em que y é a resposta e x a matriz com as covariáveis: > coef(result): extrai as estimativas dos coeficientes do modelo (beta)
#                                                        > vcov(result): extrai a matriz de variância e covariância de beta
#                                                        > fitted(result): valores ajustados (F(eta))
#                                                        > predict(result): valores preditos (eta)


# Diagnóstico e envelope para o modelo Poisson
source("https://www.ime.unicamp.br/~cnaber/diag_pois.r")
source("https://www.ime.unicamp.br/~cnaber/envel_pois.r")



########################################
########################################
# Dados de sobrevivência das bactérias
########################################
########################################


tempo = c(1,2,3,4,5,6,7,8,9,10,11,12) # variavel tempo
auxt <- tempo-mean(tempo) # variavel tempo corrgida pela media
auxt2 <- auxt^2 # variavel tempo corrigida pela media ao quadrado
nbac = c(175, 108, 95, 82 ,71 ,50 ,49, 31, 28, 17 ,16 ,11) # numero de bacterias
n = length(tempo) # tamanho da amostra


# Dados
tabela <- rbind(nbac,tempo)
xtable(tabela)


# gráfico de dispersão
plot(tempo,nbac,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="tempo de exposição",ylab="número de bactérias sobreviventes",pch=19)


# Modelo 1
resultM1 <- fit.model <- result <- glm(nbac~auxt,family=poisson(link="log"))
summary(result)
xtable(summary(result))

desvioM1 <- deviance(result)
p <- ncol(model.matrix(result))
pvdesvM1 <- 1-pchisq(desvioM1,df=n-p)

diagPoisson(fit.model)
par(mfrow=c(1,1))
envelPoisson(fit.model,"log")

AICM1 <- AIC(fit.model)
BICM1 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaM1 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM1)
mrebeta <- rbind(rebetaM1)
covbetaM1 <- vcov(fit.model)

pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredM1 <- pred$fit
semupredM1 <- pred$se.fit
liICM1 = mupredM1-ez*semupredM1
lsICM1 = mupredM1+ez*semupredM1


# Modelo 2
resultM2 <- fit.model <- result <- glm(nbac~auxt+auxt2,family=poisson(link="log"))
summary(result)
xtable(summary(result))

desvioM2 <- deviance(result)
p <- ncol(model.matrix(result))
pvdesvM2 <- 1-pchisq(desvioM2,df=n-p) 

diagPoisson(fit.model)
par(mfrow=c(1,1))
envelPoisson(fit.model,"log")

AICM2 <- AIC(fit.model)
BICM2 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaM2 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM2)
mrebeta <- rbind(mrebeta,rebetaM2)
covbetaM2 <- vcov(fit.model)

pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredM2 <- pred$fit
semupredM2 <- pred$se.fit
liICM2 = mupredM2-ez*semupredM2
lsICM2 = mupredM2+ez*semupredM2



# Estatísticas de comparação de modelos

AICBICD <- rbind(cbind(AICM1,AICM2),cbind(BICM1,BICM2),cbind(desvioM1,desvioM2),cbind(pvdesvM1,pvdesvM2))
colnames(AICBICD) <- c("Modelo 1","Modelo 2")
rownames(AICBICD) <- c("AIC","BIC","Desvio","p-valor desvio")
xtable(t(AICBICD))
### OBS: Modelo 1 preferível ao Modelo 2 pois, apesar deles terem equivalencia na qualidade do ajuste, o Modelo 1 é mais parcimonioso e obteve menores valores de AIC, BIC



# Estimativas dos parâmetros do modelo

xtable(mrebeta)



# Gráficos de diagnóstico

diagPoisson(resultM1)
diagPoisson(resultM2)

# Gráficos de envelope

par(mfrow=c(1,2))
envelPoisson(resultM1,"log")
title("Modelo 1",cex=1.2)
envelPoisson(resultM2,"log")
title("Modelo 2",cex=1.2)

# Valores preditos

par(mfrow=c(1,2))
plot(tempo,nbac,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="tempo de exposição",ylab="número de bactérias sobreviventes",pch=17)
legend(6,140,c("observado","predito"),col=c(1,2),pch=c(17,19),bty="n",cex=1.3)
title("Modelo 1",cex=1.2)
plotCI(tempo,mupredM1,li=liICM1,ui=lsICM1,add=TRUE,cex=1.3,cex.axis=1.3,cex.lab=1.3,col=2,pch=19)

plot(tempo,nbac,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="tempo de exposição",ylab="número de bactérias sobreviventes",pch=17)
legend(6,140,c("observado","predito"),col=c(1,2),pch=c(17,19),bty="n",cex=1.3)
title("Modelo 2",cex=1.2)
plotCI(tempo,mupredM2,li=liICM2,ui=lsICM2,add=TRUE,cex=1.3,cex.axis=1.3,cex.lab=1.3,col=2,pch=19)



# Estimativas da taxa de diminução exp(beta(11))

taup <- exp(rebetaM1[2,1])
auxta <-  rbind(as.numeric(c(0,taup)))
eptaup <-  sqrt(auxta%*%covbetaM1%*%t(auxta))
ictaup <- c(taup-ez*eptaup,taup+ez*eptaup)



# Modelo de regressão segmentada: 
# (isso esta sendo feito em razao da dificuldade em se escolher 
# qual o melhor modelo, devido a equivalencia dos modelos propostos)

mXseg <- rbind(cbind(1,auxt[1:3],0,0),cbind(0,0,1,auxt[4:12])) # nota-se que dividiu-se a matriz de deliamento em blocos, a primeira parte da obs 1 ate a 2, e a segunda da obs 4 ate a 12
xtable(mXseg)

resultM3 <- fit.model <- result <- glm(nbac~-1+mXseg,family=poisson(link="log"))
summary(result) # > os dois primeiros coeficientes são referentes aos betas das tres primeiras observações, ou seja, da primeira parte da matriz de delinemento
                # > os dois ultimos coeficientes são referentes aos betas das tos restante das observações, ou seja, da segunda parte da matriz de delinemento
xtable(summary(result))

desvioM3 <- deviance(result)
p <- ncol(model.matrix(result))
pvdesvM3 <- 1-pchisq(desvioM3,df=n-p)

diagPoisson(fit.model)
par(mfrow=c(1,1))
envelPoisson(fit.model,"log")

AICM3 <- AIC(fit.model)
BICM3 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaM3 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM3)
mrebeta <- rbind(mrebeta,rebetaM3)
covbetaM3 <- vcov(fit.model)


# Gr[aficos das estimativas

par(mfrow=c(1,2))
plotCI(rebeta[c(1,3),1],li=rebetaM3[c(1,3),3],ui=rebetaM3[c(1,3),4],cex=1.3,cex.axis=1.3,cex.lab=1.3,col=1,pch=19,axes=F,xlab="parâmetro",ylab="estimativa")
#axis(2,at=c(-0.40,4))
axis(2)
#axis(1,at=1:4,labels=c("beta01","beta11","beta02","beta12"))
#axis(1,at=1:4,labels=c(expression(paste(beta,"01")),expression(paste(beta,"11")),expression(paste(beta,"02")),expression(paste(beta,"12"))))
axis(1,at=1:2,labels=c(expression(paste(beta,"01")),expression(paste(beta,"02"))))
#
plotCI(rebeta[c(2,4),1],li=rebetaM3[c(2,4),3],ui=rebetaM3[c(2,4),4],cex=1.3,cex.axis=1.3,cex.lab=1.3,col=1,pch=19,axes=F,xlab="parâmetro",ylab="estimativa")
#axis(2,at=c(-0.40,4))
axis(2)
#axis(1,at=1:4,labels=c("beta01","beta11","beta02","beta12"))
#axis(1,at=1:4,labels=c(expression(paste(beta,"01")),expression(paste(beta,"11")),expression(paste(beta,"02")),expression(paste(beta,"12"))))
axis(1,at=1:2,labels=c(expression(paste(beta,"11")),expression(paste(beta,"12"))))

### OBS > nota-se que tanto no intercepto quanto no coeficiente angular da primeira parte possuem um nivel de incerteza maior em relação a segunda parte, em razao de ter menos dados



pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredM3 <- pred$fit
semupredM3 <- pred$se.fit
liICM3 = mupredM3-ez*semupredM3
lsICM3= mupredM3+ez*semupredM3



# Estatísticas de comparação de modelos

AICBICD <- rbind(cbind(AICM3),cbind(BICM3),cbind(desvioM3),cbind(pvdesvM3))
colnames(AICBICD) <- c("Modelo Segmentado")
rownames(AICBICD) <- c("AIC","BIC","Desvio","p-valor desvio")
xtable(t(AICBICD))
### OBS > o modelo 1 continua sendo preferivel, por ser o mais parcimonioso, ou seja, tem menores valores de AIC E BIC
###     > nota-se uma queda no desvio, ou seja, nota-se um ajuste melhor no modelo segmentado, mas essa melhora no ajuste nao compensa o aumento no numero de parametros



# Gráficos de diagnóstico

diagPoisson(resultM3)

# Gráficos de envelope

par(mfrow=c(1,1))
envelPoisson(resultM3,"log")


plot(tempo,nbac,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="tempo de exposição",ylab="número de bactérias sobreviventes",pch=17,ylim=c(0,200))
title("Modelo 2",cex=1.2)
plotCI(tempo,mupredM3,li=liICM3,ui=lsICM3,add=TRUE,cex=1.3,cex.axis=1.3,cex.lab=1.3,col=2,pch=19)
plotCI(tempo,mupredM1,li=liICM1,ui=lsICM1,add=TRUE,cex=1.3,cex.axis=1.3,cex.lab=1.3,col=3,pch=19)
plotCI(tempo,mupredM2,li=liICM2,ui=lsICM2,add=TRUE,cex=1.3,cex.axis=1.3,cex.lab=1.3,col=4,pch=19)
legend(6,140,c("observado","predito - modelo 3","predito - modelo 1","predito - modelo 2"),col=c(1,2,3,4),pch=c(17,19,19,19),bty="n",cex=1.3)