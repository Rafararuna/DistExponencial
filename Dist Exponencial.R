### Exemplo - Distribui??o de Poisson

# Defininfo as vari?veis:

n <- 100 # tamanho da amostra
x <- runif(n,1,10)  # variavel explicativa
eta <- 2-.6*x # preditor linear ; onde 2 ? o intercepto (b0) e -0,6 ? o b1 ; estrutura linear
mu <- exp(eta) # fun??o de liga??o do log, ent?o a m?dia vai ser o exponencial do preditor linear
y <- rpois(n,mu) # vari?vel resposta
plot(x, y, pch=20)

X <- matrix(c(rep(1,n),x),nrow=n) # Matriz do modelo;



# Agora, vamos implementar o algoritmo Escore de Fisher:

N <- 10 # numero de intera??es
beta <- matrix(0, nrow=N, ncol=2)
epsilon <- numeric()
### beta vai guardar as estimativas de b0 e b1 a cada itera??o
### Epsilon vai guardar a medida de diferen?a quadr?tica entre estimativas obtidas para duas itera??es sucessivas.

## Vamos realizar 10 itera??es

for(i in 1:N){
  
  if(i==1) {mu <- y+.5 
  # Se for a primeira itera??o, ent?o usamos nossos chutes iniciais.
  eta <- log(mu)}
  
  if (i!=1) {
    eta <- X %*% beta[i-1,] # Preditor linear calculado com as estimativas do passo anterior
    mu <- exp(eta) # a partir da segunda itera??o, usamos os resultados do passo anterior para obter estimativas no passo atual.
  }
  
  vmu <- mu
  glinhamu <- 1/mu
  # Fun??o de variancia e derivada da fun??o de liga??o avaliadas em mu.
  
  z <- eta+(y-mu)*glinhamu
  W <- diag(as.numeric((vmu*(glinhamu**2))**(-1)))
  # Vetor z e matriz diagonal avaliadas em mu.
  
  beta[i,] <- solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%z) # Solu??o de minimos quadrados ponderados para o vetor beta no passo i.
  if(i>1) epsilon[i-1]=sum(((beta[i,]-beta[i-1,])/beta[i-1,])**2)
  
}

beta # nota-se que o valores ficaram proximos dos betas reais


# Plotando o ajuste:

plot(x,y,pch=20,cex=1.5,las=1)
ajuste = function(x) exp(beta[10,1]+beta[10,2]*x)
curve(ajuste,0,40,add=T)


## Vamos tentar o ajuste de uma regress?o linear simples

ajuste2 <- lm(y ~ x)
abline(coefficients(ajuste2), col='red')
### OBS: fica muito ruim, pois n?o pega a curva inicial e nem a estabiliza??o final


## Agora, vamos ajustar o modelo declarando a log-verossimilhan?a a um otimizador do R

require(bbmle)
logvero=function(b0,b1)
  -sum(dpois(y,exp(b0+b1*x), log=T))
### OBS: logvero armazena a fun??o de log verossimilhan?a (-)

ylinha <- y+.5
ajuste <- lm(log(ylinha)~x) # Regressando g(y+.5) em fun??o de x para obter valores iniciais para b0 e b1.
est2 <- mle2(logvero,start=list(b0=ajuste$coefficients[1],b1=ajuste$coefficients[2]))
est2 # praticamente os mesmos valores que encontramos no escore de fisher



## Finalmente, usemos a fun??o glm

ajuste <- glm(y~x,family=poisson(link = "log"))
names(ajuste)
coef(ajuste) # Estimativas dos coeficientes do modelo.
fitted(ajuste) # M?dias ajustadas pelo modelo para cada valor de x na amostra.
predict(ajuste,newdata=data.frame(x=c(3.5,5.5,7.5))) # Estimativas (na escala do preditor) para x=3.5; x=5.5 e x=7.5.
predict(ajuste,newdata=data.frame(x=c(3.5,5.5,7.5)),type='response') # Estimativas (na escala da resposta - m?dias estimadas) para x=3.5; x=5.5 e x=7.5.
summary(ajuste) # Resumo do modelo ajustado contendo, dentre outras coisas, as estimativas dos betas e os correspondnetes erros padr?es.
vcov(ajuste) # Matriz de vari?ncias e covari?ncias estimada.