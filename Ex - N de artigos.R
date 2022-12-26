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
library(pscl)

op <- options()

# Diagnóstico e envelope para o modelo Poisson
source("https://www.ime.unicamp.br/~cnaber/diag_pois.r")
source("https://www.ime.unicamp.br/~cnaber/envel_pois.r")

# Diagnóstico e envelope para o modelo binomial negativo
source("https://www.ime.unicamp.br/~cnaber/diag_nbin.r")
source("https://www.ime.unicamp.br/~cnaber/envel_nbin.r")

# Teste CB=M Binomial Negativa
source("https://www.ime.unicamp.br/~cnaber/TestCBMBinNeg.r")

# Diagnóstico e envelope para o modelos de contagem inflacionada
source("https://www.ime.unicamp.br/~cnaber/diag_cont_inf.r")
source("https://www.ime.unicamp.br/~cnaber/envel_cont_inf.r")


# Definindo variaveis do banco:
m.dados <- bioChemists
n <- nrow(m.dados) # tamanho da amostra
vamos <- sample(seq(1:915),n,replace=FALSE) # pegando uma amostra aleatoria de tamanho n
m.dados <- m.dados[vamos,] # atualizando o banco ; transformando o banco original na amostra aleatoria

nart <- m.dados$art # número de artigos ; variavel resposta
gen <- m.dados$fem  # genero ; variavel explicativa
gen <- revalue(gen,c("Men"="masculino","Women"="feminino")) # alterando a variavel genero
ecivil <- m.dados$mar # estado civil ; variavel explicativa
ecivil <- revalue(ecivil,c("Single"="solteiro","Married"="casado")) # alterando a variavel estado civil
nfil <- m.dados$kid5 # numero de filhos ; variavel explicativa
nartment <- m.dados$ment # quantidade de artigos produzidos no ultimos 3 anos pelo orientador ; variavel explicativa
ephd <- m.dados$phd # prestigio do departamento aonde o aluno desenvolveu seus estudos (entre 0 e 5) ; variavel explicativa

n <- length(nart) # nota-se que, na maioria dos dados, tem-se uma grande concentração de zeros e um's
                  # essa grande concentração de zeros indica que, provavelmente, um modelo poisson ou binomial negativo não consiga captar essa freq. de zeros
hist(nart) # aqui no no histograma tbm podemos notar essa concentração de zeros
 
## centralizando nas médias
nfilc  <-  nfil - mean(nfil)
nartmentc  <-  nartment - mean(nartment)
nephc  <-  ephd - mean(ephd)



# dispersões
par(mfrow=c(2,2))
vgrupo <- interaction(gen,ecivil)
plot(nfil,nart,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de filhos",ylab="número de artigos",pch=19)
plot(ephd,nart,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="prestígio de departamento",ylab="número de artigos",pch=19)
plot(nartment,nart,cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de artigos do mentor",ylab="número de artigos",pch=19)
### OBS > em relação ao numero de filhos, nota-se uma tendencia de redução, a medida que aumenta o numero de filhos, o numero de artigos diminui
###     > em relação ao prestigio de departamento, no inicio percebe-se uma tendencia de aumento, mas depois parece que ha uma estabilização ; alem disso, nota-se valores mais altos de numero de artigos em valores de prestigios menores (proximo do prestigio 2)
###     > em relação ao numero de artigos do mentor, nao vemos muita associação ; importante notar a concentração na parte inferior

# dispersões por grupo
par(mfrow=c(1,2))
#plot(nfil[vgrupo=="Men.Single"],nart[vgrupo=="Men.Single"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de filhos",ylab="número de artigos",pch=19)
plot(nfil[vgrupo=="masculino.casado"],nart[vgrupo=="masculino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de filhos",ylab="número de artigos",pch=19)
title("homens casados",cex=1.2)
#plot(nfil[vgrupo=="feminino.solteiro"],nart[vgrupo=="feminino.solteiro"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de filhos",ylab="número de artigos",pch=19)
plot(nfil[vgrupo=="feminino.casado"],nart[vgrupo=="feminino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de filhos",ylab="número de artigos",pch=19)
title("mulheres casadas",cex=1.2)
### OBS > nota-se que o nuemro de filhos abaixo de cinco anos vao pesar mais pras mulheres, vao produzir menos artigos
###     > porem, em ambos notamos um tendencia de de redução com o aumento do numero de filhos abaixo de cinco anos

par(mfrow=c(2,2))
plot(nartment[vgrupo=="masculino.casado"],nart[vgrupo=="masculino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de artigos do mentor",ylab="número de artigos",pch=19)
title("homens solteiros",cex=1.2)
plot(nartment[vgrupo=="masculino.casado"],nart[vgrupo=="masculino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de artigos do mentor",ylab="número de artigos",pch=19)
title("homens casados",cex=1.2)
plot(nartment[vgrupo=="feminino.solteiro"],nart[vgrupo=="feminino.solteiro"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de artigos do mentor",ylab="número de artigos",pch=19)
title("mulheres solteiras",cex=1.2)
plot(nartment[vgrupo=="feminino.casado"],nart[vgrupo=="feminino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="número de artigos do mentor",ylab="número de artigos",pch=19)
title("mulheres casadas",cex=1.2)
### OBS > nao mudou muita coisa entre os homens solteiros e casados, e as mulheres solteiras e casadas em relação ao numero de artigos do mentor, nada muito relevante

par(mfrow=c(2,2))
plot(ephd[vgrupo=="masculino.casado"],nart[vgrupo=="masculino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="prestígio de departamento",ylab="número de artigos",pch=19)
title("homens solteiros",cex=1.2)
plot(ephd[vgrupo=="masculino.casado"],nart[vgrupo=="masculino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="prestígio de departamento",ylab="número de artigos",pch=19)
title("homens casados",cex=1.2)
plot(ephd[vgrupo=="feminino.solteiro"],nart[vgrupo=="feminino.solteiro"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="prestígio de departamento",ylab="número de artigos",pch=19)
title("mulheres solteiras",cex=1.2)
plot(ephd[vgrupo=="feminino.casado"],nart[vgrupo=="feminino.casado"],cex=1.3,cex.axis=1.3,cex.lab=1.3,xlab="prestígio de departamento",ylab="número de artigos",pch=19)
title("mulheres casadas",cex=1.2)
### OBS > nao mudou muita coisa entre os homens solteiros e casados, e as mulheres solteiras e casadas em relação ao prestigio de departamento, nada muito relevante


# Box plot
par(mfrow=c(2,2))
boxplot(nart~gen,names=c("homem","mulher"),cex=1.2,cex.lab=1.2,cex.main=1.2)
boxplot(nart~ecivil,names=c("casado","solteiro"),cex=1.2,cex.lab=1.2,cex.main=1.2)
boxplot(nart~vgrupo,names=c("homem solteiro","mulher solteira","homem casado","mulher casada"),cex=1.2,cex.lab=1.2,cex.main=1.2)
### OBS > em relação genero, tem-se uma diferença entre os sexos, mas nao chega a ser discrepante ; a concentração de ambos é baixa ; nota-se valores discrepante mais altos pros homens, os valores de maximo e 3quartil aumentam ; mas em questao de mediana, media e 1 quartil nao ha tanta diferença
###     > em relação ao estado civil, nota-se que nao ha uma diferença muito relevante, apenas quando observamos os valores extremos, que sao maiores para os solteiros
###     > agora analisando todas as combinações, nota-se que na parte inferior nao ha tanta difenrença ; a diferença se concentra mais na parte superior, em relação ao 3quartil, ao valor maximo e aos valores extremos


# Medidas resumo

# por gênero
bdgen <- data.frame(gen,nart)
rvp1 <- ddply(bdgen,.(gen),summarise,media=mean(nart),dp=sqrt(var(nart)),vari=var(nart),cv=100*((sqrt(var(nart))/mean(nart))),min=min(nart),max=max(nart),n=length(nart))
xtable(rvp1)

# por estado civil
bdeciv <- data.frame(ecivil,nart)
rvp2 <- ddply(bdeciv,.(ecivil),summarise,media=mean(nart),dp=sqrt(var(nart)),vari=var(nart),cv=100*((sqrt(var(nart))/mean(nart))),min=min(nart),max=max(nart),n=length(nart))
xtable(rvp2)

# por gênero x estado civil
begeneciv <- data.frame(gen,ecivil,nart)
rvp3 <- ddply(begeneciv,.(gen,ecivil),summarise,media=mean(nart),dp=sqrt(var(nart)),vari=var(nart),cv=100*((sqrt(var(nart))/mean(nart))),min=min(nart),max=max(nart),n=length(nart))
xtable(rvp3)

# Gráficos de perfis médios
par(mfrow=c(1,1))
ez <- qnorm(0.975)
plotCI(rvp3$media[rvp3$gen=="masculino"],uiw=ez*rvp3$dp[rvp3$gen=="masculino"]/sqrt(rvp3$n[rvp3$gen=="masculino"]),liw=ez*rvp3$dp[rvp3$gen=="masculino"]/sqrt(rvp3$n[rvp3$gen=="masculino"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="estado civil",ylab="número de artigos",pch=19,col=1,ylim=c(1,2.5))
lines(rvp3$media[rvp3$gen=="masculino"],lwd=2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("solteiro","casado"),cex.axis=1.2)
plotCI(rvp3$media[rvp3$gen=="feminino"],uiw=ez*rvp3$dp[rvp3$gen=="feminino"]/sqrt(rvp3$n[rvp3$gen=="feminino"]),liw=ez*rvp3$dp[rvp3$gen=="feminino"]/sqrt(rvp3$n[rvp3$gen=="feminino"]),axes=FALSE,cex.lab=1.2,cex.axis=1.2,cex=1.2,xlab="estado civil",ylab="número de artigos",pch=17,col=2,ylim=c(0,350),add=TRUE)
lines(rvp3$media[rvp3$gen=="feminino"],lwd=2,col=2)
legend(1.5,2.5,c("masculino","feminino"),col=c(1,2),lwd=c(2,2),pch=c(19,17),bty="n",cex=1.2)
### OBS > nota-se que, quando mudamos de estado civil (de solteiro pra casado), nota-se um aumento do sexo feminino e uma redução do sexo masculino, porem nao muito discrepantes



# Ajuste do modelo completo de Poisson (sem interação)

#result<-fit.model<-glm(ncli~ndomc+rendac+idadec+discc+dislc,family=poisson("log"))
result <- fit.model <- glm(nart~nfilc+nartmentc+nephc+gen+ecivil,family=poisson("log"))
summary(result)
xtable(summary(result))

desvioM1 <- deviance(result)
p <- ncol(model.matrix(result))
pvdesvM1 <- 1-pchisq(desvioM1,df=n-p) # p-valor muito baixo, ou seja, rejeita a hipotese nula, ou seja, rejeita o ajuste

diagPoisson(fit.model) # no primeiro, nota-se uma quantidade consideravel de pontos fora das bandas
                       # no segundo, nota-se uma concentração dos pontos na parte inferior
                       # nota-se um histograma pouco assimetrico, um concentração maior na parte inferior, abaixo de zero
                       # nota-se, no inicio do grafico, pontos muito acima da reta
#abline(0,1)
par(mfrow=c(1,1))
envelPoisson(fit.model,"log") # nota-se os pontos saindo bastante das bandas

AICM1 <- AIC(fit.model)
BICM1 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaM1 <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaM1)
mrebeta <- rbind(rebetaM1)
covbetaM1 <- vcov(fit.model)

estcomp <- cbind(AICM1,BICM1,desvioM1,pvdesvM1)


# Ajuste do modelo completo binomial negativo (sem interação)

#result<-fit.model<-glm(ncli~ndomc+rendac+idadec+discc+dislc,family=poisson("log"))
result <- fit.model <- glm.nb(nart~nfilc+nartmentc+nephc+gen+ecivil,link="log")
summary(result)
xtable(summary(result))
 
desvioM2 <- deviance(result)
p <- ncol(model.matrix(result))
pvdesvM2 <- 1-pchisq(desvioM2,df=n-p) # p-valor continua baixo, ou seja, continua rejeitando a hipotese nula, ou seja, continua rejeitando o ajuste

diagnbin(fit.model) # no primeiro grafico, nota-se uma quantidade menor de pontos fora das bandas ; a escala tbm diminui
                    # no segundo grafico ainda nota-se uma concentração na parte inferior ; a escala tbm diminui
                    # no histograma, nota-se um comportamento melhor dos residuos, mais proximo de uma simetria em torno do zero, apesar de ainda haver uma concentração forte abaixo de zero
                    # nota-se, no inicio do grafico, pontos muito acima da reta
                    # de maneira geral, pela analise grafica esse modelo se ajustou melhor aos dados, pois ele consegue corrigir u pouco dessa questão dos valores dicrepantes,
#                     pois esse modelo tem mais flexibilidade para modelar a dispersao dos dados, pois o modelo poisson tem uma restrição forte, ja que é um parametro a media e variancia
#abline(0,1)
par(mfrow=c(1,1))
ligacaonbin <- "log"
envelnbin(fit.model) # nota-se um ajuste bem melhor, praticamente todos os pontos dentro das bandas

AICM2 <- AIC(fit.model)
BICM2 <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
retheta <- result$theta
reeptheta <- result$SE.theta
rebetaM2 <- rbind(cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4]),
                  cbind(retheta,reeptheta,retheta-ez*reeptheta,retheta+ez*reeptheta,0,0))
xtable(rebetaM2)
mrebeta <- rbind(rebetaM2)
covbetaM2 <- vcov(fit.model)

estcomp <- rbind(estcomp,cbind(AICM2,BICM2,desvioM2,pvdesvM2))


# Ajuste do modelo Poisson inflacionado de zeros (modelando apenas a média
# da contagem original)

result <- fit.model <- zeroinfl(nart~nfilc+nartmentc+nephc+gen+ecivil|1,dist=c("poisson"))
summary(result)
#xtable(summary(result)) # nao da certo para a função zeroinfl

#desvioM3 <- deviance(result) # nao da certo
#p <- ncol(model.matrix(result)) # nao da certo
#pvdesvM3 <- 1-pchisq(desvioM3,df=n-p) # nao da certo

#diagnbin(fit.model) # nao da certo
diagcontinf(fit.model,0,nart) # no primeiro grafico nota-se uma correção da parte inferior, mas ha muitos pontos fora das bandas na parte de cima
                              # no segundo grafico ainda nota-se uma concentração na parte inferior
                              # no histograma há assimetria ainda
                              # no quarto grafico, os pontos ainda estao bem longe do que a gente espera, estao longe da reta
abline(0,1)                   # de maneira geral, ha uma melhora, mas nao é uma grande melhora, ainda notamos problemas que nao indicam um bom ajuste
par(mfrow=c(1,1))
ligacaonbin <- "log"
envelcontinf(fit.model,1,0) # alguns pontos estao dentro das bandas no inicio, mas notamos muitos mais pontos fora das bandas depois
                            # porem, comparando com o primeiro modelo há uma melhora nitida do ajuste

AICM3 <- AIC(fit.model)
BICM3 <- AIC(fit.model,k=log(n))

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef$count
repifr <- (summary(result))$coef$zero
repif <- exp(repifr[1])/(1+exp(repifr[1]))
reepif <- repifr[2]*(1/(1+exp(repif[1]))) # erro-padrão muito pequeno
rebetaM3 <- rbind(cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4]),
                  cbind(repif,reepif,repif-ez*reepif,repif+ez*reepif,0,0))
#cbind(retheta,reeptheta,retheta-ez*reeptheta,retheta+ez*reeptheta,0,0))
xtable(rebetaM3)
mrebeta <- rbind(rebetaM3)
covbetaM3 <- vcov(fit.model)

estcomp <- rbind(estcomp,cbind(AICM3,BICM3,0,0))


# Ajuste binomial negativo inflacionado de zeros (modelando apenas a média
# da contagem original)

result <- fit.model <- zeroinfl(nart~nfilc+nartmentc+nephc+gen+ecivil|1,dist=c("negbin"))
summary(result)
#xtable(summary(result)) # não da certo

#desvioM4 <- deviance(result) # não da certo
#p <- ncol(model.matrix(result)) # não da certo
#pvdesvM4 <- 1-pchisq(desvioM4,df=n-p) # não da certo

#diagnbin(fit.model) # não da certo
diagcontinf(fit.model,0,nart) # no primeiro grafico nota-se uma correção da parte inferior, mas ha muitos pontos fora das bandas na parte de cima
                              # no segundo grafico ainda nota-se uma concentração na parte inferior
                              # no histograma há uma concentração melhor em torno do zero, pois o modelo 2 nao conseguia captar a frequencia de zeros, entao havia uma concentração abaixo de zero consequentemente ; mas ainda nota-se uma concentração na parte positiva
                              # no quarto grafico, os pontos ainda estao bem longe do que a gente espera, estao longe da reta
                              # de maneira geral, comparando com o modelo 2, notamos ganho apenas em relaçao ao histograma, ou seja, teve muito mais perda em relação ao ajuste
abline(0,1)
par(mfrow=c(1,1))
ligacaonbin <- "log"
envelcontinf(fit.model,2,0) # nota-se q o modelo se ajusta bem aos dados, a maioria dos pontos estao dentro das bandas

AICM4 <- AIC(fit.model)
BICM4 <- AIC(fit.model,k=log(n))

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef$count
repifr <- (summary(result))$coef$zero
repif <- exp(repifr[1])/(1+exp(repifr[1]))
reepif <- repifr[2]*(1/(1+exp(repif[1]))) # erro-padrão muito pequeno
retheta <- result$theta
reeptheta <- sqrt(exp(retheta))*result$SE.logtheta
rebetaM4 <- rbind(cbind(rebeta[-7,1],rebeta[-7,2],rebeta[-7,1]-ez*rebeta[-7,2],rebeta[-7,1]+ez*rebeta[-7,2],rebeta[-7,3],rebeta[-7,4]),
                  cbind(retheta,reeptheta,retheta-ez*reeptheta,retheta+ez*reeptheta,0,0),
                  cbind(repif,reepif,0,0,0,0))
#cbind(retheta,reeptheta,retheta-ez*reeptheta,retheta+ez*reeptheta,0,0))
xtable(rebetaM4)
mrebeta <- rbind(rebetaM4)
covbetaM4 <- vcov(fit.model)

estcomp <- rbind(estcomp,cbind(AICM4,BICM4,0,0))
xtable(estcomp)

### OBS > nota-se que os modelos da binomial negativa se ajustaram melhor aos dados, obtendo menores valores de AIC E BIC
###     > nota-se, pelos valores do desvio, que os modelos possion rejeitaram o ajuste