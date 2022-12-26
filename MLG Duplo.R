# Pacotes

require(dglm)

source("https://www.ime.unicamp.br/~cnaber/envel_norm.r")
source("C:\\Users\\Leandro\\Desktop\\Arquivos HD antigo\\Aulas\\MLG\\Aula 22\\envel_gama_dglm.r")
# essa funlção acima é adaptação do código do prof. Gilberto Paula
# "https://www.ime.usp.br/~giapaula/envel_norm_dglm"



# Analises

n <- 50 # tamanho da amostra
x <- runif(n,10,20) # para compor o preditor linear da media
z <- runif(n,1,5) # parametro de dispersao

ls <- -4+ 1.3*z # estrutura do preditor linear para o parametro de dispersao
sig2 <- exp(ls) # para obter as variancias

y <- 2-.6*x+rnorm(n,0,sqrt(sig2)) # estrutura linear, supondo um modelo normal

for(i in 1:n){y[i]<-2-.6*x[i]+rnorm(1,0,sqrt(sig2[i]))}

X11()
par(mfrow=c(1,1))
plot(x,y,pch=16) # tendencia linear de redução
plot(z,y,pch=16) # o que a gente espera é q, a medida q z aumenta, a variabilidade aumenta

fit <- lm(y~x) # modelo homocedastico
summary(fit)
envelnorm(fit) # a heterocedasticidade atrapalha um pouco o ajuste, por isso temos pontos pra fora

fit2 <- dglm(y~x,~z,family = gaussian) # modelo heterocedastico, modelagem para o parametro de dispersao
summary(fit2)
envel_norm_dglm(fit2) # nota-se que teve uma melhora do ajuste