############################################################################
#                   Autor - Larangeiras 29/09/2015                         #
############################################################################

library(scatterplot3d)


############################################################################
#                             FUNCTION 1                                   #
############################################################################
# Essas imagens estao no formato enxuto, ou seja, pussui 3 bandas organiza #
# das em um array de dimensao 9. As 3 primeiras dimensoes do array contem  #
# as intensidades, ja as outras 3 a parte real dos numeros complexos e as 3#
# ultimas a parte imaginaria.                                              #
# Essa funcao tem como entradas a imagem e as cordenadas do pixel, dados a #
# imagem e as coordenadas do pixel retorna a matriz de covariancia complexa#
# que o representa                                                         #
############################################################################

Matriz_Cov_Complex <- function( G , x , y )
{
  
  if( (x>=1 && x<=dim(G)[1])  && (y>=1 && y<=dim(G)[2]) )
  {
    
    intensidades <- G[x,y,  1:3]
    Parte_real <- G[x,y,  4:6]
    Parte_imaginaria <- G[x,y,  7:9]
    Triangular_superior <- complex( real = Parte_real , imaginary = Parte_imaginaria )
    Triangular_inferior <- Conj(Triangular_superior)
    M <- matrix(0 , nrow = 3, ncol = 3)
    M[1,1] <- intensidades[1]
    M[2,2] <- intensidades[2]
    M[3,3] <- intensidades[3]
    M[1,2] <- Triangular_superior[1]
    M[1,3] <- Triangular_superior[2]
    M[2,3] <- Triangular_superior[3]
    M[2,1] <- Triangular_inferior[1]
    M[3,1] <- Triangular_inferior[2]
    M[3,2] <- Triangular_inferior[3]
    
    
    return( M )
    
  }
  else{ print( "Voce digitou coordenadas fora da dimensao da imagem" ) }
  
}

############################################################################
#                             FUNCTION 2                                   #
############################################################################
# A ideia desta funcao e selecionar uma parte da imagem e retornar as matri#
# zes de covariancia complexas de cada pixel. Para selecionar a parte dese #
# jada da imagem indique a linha inicial e final e tambem a coluna inicial #
# e final, dai a funcao retornara as matrizes de covariancia complexas de  #
# cada pixel.
############################################################################

Matrizes_Cov_Complex = function( G  , linha_inicial , linha_final ,  coluna_inicial , coluna_final )
{
  
  if(    (linha_inicial>=1 && linha_inicial<=dim(G)[1] )  &&  ( linha_final>=1 && linha_final<=dim(G)[1] ) && 
           ( coluna_inicial>=1 && coluna_inicial<=dim(G)[2] ) && ( coluna_final>=1 && coluna_final<=dim(G)[2] )  )
  {
    
    resultados = list(  ) 
    contador = 1
    for( i in 1:dim(G)[1] )
    {
      for( j in 1:dim(G)[2] )
      {
        if( (i>=linha_inicial && i<=linha_final ) && ( j>=coluna_inicial && j<=coluna_final ) )
        {
          resultados[[contador]] = Matriz_Cov_Complex(G,i,j)
          contador = contador+1
        }
      }
    }
    
    return(resultados)
    
  }
  else{ print( "Voce digitou coordenadas fora da dimensao da imagem" ) }
  
}

############################################################################
#                             FUNCTION 3                                   #
############################################################################
# A ideia desta funcao e verificar se uma matriz complexa e simetrica, ou se
# seja, se A = A^{*}, onde A^{*}= (\overline{A})^{t}. 
# Nota: A funcao identical( ) verifica se dois objetos sao iguais retornando 
# TRUE ou FALSE
############################################################################

# Matriz_Complex_simet <- function( M )
# {
#   M_estrela <- t( Conj(M) )
#   if( identical( M , M_estrela ) )
#   {
#     print("A matriz e simetrica")
#   }
#   else{ print("A matriz nao e simetrica") }
#   
# }

############################################################################
#                             FUNCTION 4                                   #
############################################################################
# A ideia desta funcao e calcular o determinante da matriz complexa que sao
# pixel das imagens PolSAR.
# Verificar com mais calma se as contas estao certas ou se existe algum paco
# te do R que calcule determinante de matrizes complexas 
############################################################################

det_Cov_Complex <- function(G)
{
  I1 <- G[1,1]
  I2 <- G[2,2]
  I3 <- G[3,3]
  
  det_Cov_Complex <- (I1*I2*I3)-( I1* ( G[2,3]*Conj(G[2,3]) ) )-( I3* ( G[1,2]*Conj(G[1,2]) ) )+
    ( G[1,2]*G[3,1]*G[2,3] )+(G[2,1]*G[1,3]*G[3,2])-(I2*( G[1,3]*Conj(G[1,3]) ) )
  
  
  # Peguei apenas a parte real dos dados porque o determinante dessas matrizes
  # sao numeros reais
  return( Re(det_Cov_Complex) )
}

############################################################################
#                             FUNCTION 5                                   #
############################################################################
# A ideia desta funcao e generalizar a funcao anterior "det_Cov_Complex( )"
# Observe que essse M e uma lista de matrizes, por isso quando acesso estou
# usando M[[i]]
############################################################################

det_Cov_Complex_de_varias_matrizes <- function( M )
{
  det_Cov_Complex_de_varias_matrizes = c( )# estou preenchendo num vetor porque sao numeros e nao matrizes
  contador=1
  for( i in 1:length(M) )
  {
    det_Cov_Complex_de_varias_matrizes[contador] = det_Cov_Complex(M[[i]])
    contador=contador+1
  }
  
  
  return( det_Cov_Complex_de_varias_matrizes )
}

############################################################################
#                             FUNCTION 6                                   #
############################################################################
# A variancia efetiva e uma estatistica escalar calculada da seguinte maneira
# determinante da matriz de covariancia elevado ao inverso da dimensao ou a
# raiz enesima do determinante, onde "n" e a dimensao da matriz.
# A funcao de forma geral teria que receber 2 parametros a matriz de covari-
# ancia e a dimensao da matriz, em particular, nossa matriz sempre e de 3X3.
############################################################################

Variancia_efetiva <- function(G)
{
  Variancia_efetiva <- (det_Cov_Complex(G))^(1/3)
  
  return(Variancia_efetiva)
}

############################################################################
#                             FUNCTION 7                                   #
############################################################################
# A ideia desta funcao e generalizar a funcao anterior "Variancia_efetiva()"
# Observe que essse M e uma lista de matrizes, por isso quando acesso estou
# usando M[[i]]
############################################################################

 
Variancia_efetiva_de_varias_matrizes <- function( M )
{
  Variancia_efetiva_de_varias_matrizes = c( )# estou preenchendo num vetor porque sao numeros e nao matrizes
  contador=1
  for( i in 1:length(M) )
  {
    Variancia_efetiva_de_varias_matrizes[contador] = Variancia_efetiva(M[[i]])
    contador=contador+1
  }
  
  
  return( Variancia_efetiva_de_varias_matrizes )
}

############################################################################
#                             FUNCTION 8                                   #
############################################################################
# A ideia desta funcao e plotar um histograma das "Variancia_efetiva()"
############################################################################

Histograma.VarianciaEfetiva <- function(dados)
{
  dados_para_saida <- Variancia_efetiva_de_varias_matrizes(dados)
  Histograma.VarianciaEfetiva <- hist(dados_para_saida, breaks = "FD", probability = T,
                                       xlab ="Effective Variance", main = "")
}
  
############################################################################
#                             FUNCTION 9                                   #
############################################################################
# A dependencia efetiva e uma estatistica escalar calculada da seguinte maneira
# 1 - determinante da matriz de correlacao elevado ao inverso da dimensao ou
# a raiz enesima do determinante, onde "n" e a dimensao da matriz.
# A funcao de forma geral teria que receber 2 parametros a matriz de covari-
# ancia e a dimensao da matriz, em particular, nossa matriz sempre e de 3X3.
############################################################################

Dependencia_efetiva <- function(G)
{
  det_Corr_Complex <- (det_Cov_Complex(G))^2
  Dependencia_efetiva <- 1-( (det_Corr_Complex)^(1/3) )
  
  return(Dependencia_efetiva)
}

############################################################################
#                             FUNCTION 10                                  #
############################################################################
# A ideia desta funcao e generalizar a funcao anterior "Dependencia_efetiva()"
# Observe que essse M e uma lista de matrizes, por isso quando acesso estou
# usando M[[i]]
############################################################################


Dependencia_efetiva_de_varias_matrizes <- function( M )
{
  Dependencia_efetiva_de_varias_matrizes = c( )# estou preenchendo num vetor porque sao numeros e nao matrizes
  contador=1
  for( i in 1:length(M) )
  {
    Dependencia_efetiva_de_varias_matrizes[contador] = Dependencia_efetiva(M[[i]])
    contador=contador+1
  }
  
  
  return( Dependencia_efetiva_de_varias_matrizes )
}

############################################################################
#                             FUNCTION 11                                  #
############################################################################
# A ideia desta funcao e plotar um histograma das "Dependencia_efetiva()"
############################################################################

Histograma.DependenciaEfetiva <- function(dados)
{
  dados_para_saida <- Dependencia_efetiva_de_varias_matrizes(dados)
  Histograma.DependenciaEfetiva <- hist(dados_para_saida, breaks = "FD", probability = T,
                                         xlab ="Effective Dependence", main = "")
}

############################################################################
#                             FUNCTION 12                                  #
############################################################################
# A idea dessa funcao e calcular o logaritmo do desvio padrao da matriz de
# covarinacia
############################################################################

log_DP_1 <- function(M)
{
  # verificar qual log usar: natural, base 10, base 2 etc.
  log_DP_1 <- log10(sqrt(Re(M[1,1]))) 
  
  return(log_DP_1)
}

############################################################################
#                             FUNCTION 13                                  #
############################################################################
# A idea dessa funcao e generalizar "log_DP( )", ou seja, calcular para uma
# lista de matrizes
############################################################################

log_DP_1_varias_matrizes <- function(M)
{
  log_DP_1_varias_matrizes <- c( )
  contador=1
  for( i in 1:length(M) )
  {
    log_DP_1_varias_matrizes[contador] = log_DP_1(M[[i]])
    contador=contador+1
  }
  
 return( log_DP_1_varias_matrizes )
}

############################################################################
#                             FUNCTION 14                                  #
############################################################################
# A ideia dessa funcao e plotar o log do desvio padrao de todas as matrizes
# escolhidas 
############################################################################

Histograma.Variancia <- function(dados)
{
  dados_para_saida <- log_DP_1_varias_matrizes(dados)
  Histograma.Variancia <- hist(dados_para_saida, breaks = "FD", probability = T, 
                              xlab = expression(log~sigma[1]), main = "")
}

############################################################################
#                             FUNCTION 15                                  #
############################################################################
# A ideia dessa funcao e calcular a correlacao marginal $\rho_{12}$ 
############################################################################

Corr_marginal_12 <- function(S)
{
  Corr_marginal_12 <- (Re(S[1,2])/sqrt(Re(S[1,1])*Re(S[2,2])))
  
  return(Corr_marginal_12)
}

############################################################################
#                             FUNCTION 16                                  #
############################################################################
# A idea dessa funcao e generalizar "Corr_marginal_12( )", ou seja, calcular
# para uma lista de matrizes
############################################################################

Corr_marginal_12_varias_matrizes <- function(S)
{
  Corr_marginal_12_varias_matrizes <- c( )
  contador=1
  for( i in 1:length(S) )
  {
    Corr_marginal_12_varias_matrizes[contador] = Corr_marginal_12(S[[i]])
    contador=contador+1
  }
  
  
  return( Corr_marginal_12_varias_matrizes )
}

############################################################################
#                             FUNCTION 17                                  #
############################################################################
# A ideia dessa funcao e plotar o log do desvio padrao de todas as matrizes
# escolhidas 
############################################################################

Histograma.Corr <- function(dados)
{
  dados_para_saida <- Corr_marginal_12_varias_matrizes(dados)
  Histograma.Corr <- hist(dados_para_saida, breaks ="FD", probability = T, xlab = expression(rho[12]) , main = "", xlim=c(-1,1))
  lines(density(dados_para_saida), col="green")
}

############################################################################
#                             FUNCTION 18                                  #
############################################################################
# A ideia dessa funcao e calcular a matriz de correlacao R dada a matriz de 
# Covariancia S 
# Note que essa matriz e simetrica
############################################################################

Matriz_Corr_Complex <- function(S)
{
  R <- matrix( , 3, 3)
  R[1,1] <- R[2,2] <- R[3,3] <- 1
  R[1,2] <- Re(  S[1,2]/sqrt(S[1,1]*S[2,2]) )
  R[1,3] <- Re(  S[1,3]/sqrt(S[1,1]*S[3,3]) )
  R[2,1] <- Re(  S[2,1]/sqrt(S[2,2]*S[1,1]) )
  R[2,3] <- Re(  S[2,3]/sqrt(S[2,2]*S[3,3]) )
  R[3,1] <- Re(  S[3,1]/sqrt(S[3,3]*S[1,1]) )
  R[3,2] <- Re(  S[3,2]/sqrt(S[3,3]*S[2,2]) )
  
  
  return(R)
}

############################################################################
#                             FUNCTION 19                                  #
############################################################################
# A idea dessa funcao e generalizar "Matriz_Corr_Complex( )", ou seja, calcular
# para uma lista de matrizes S
############################################################################

Matrizes_Corr_Complex <- function(S)
{
  Matrizes_Corr_Complex <- list( )
  contador = 1
  for(i in 1:length(S) )
  {
    Matrizes_Corr_Complex[[contador]] <- Matriz_Corr_Complex(S[[i]]) 
    contador = contador+1
  }
  
  
  
  return( Matrizes_Corr_Complex )
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Corr_marginal_13_varias_matrizes <- function(R)
{
  Corr_marginal_13_varias_matrizes <- c( )
  contador = 1
  for(i in 1: length(R))
  {
    Corr_marginal_13_varias_matrizes[contador] <- R[[i]][2,3]
    contador = contador+1
  }
  return(Re(Corr_marginal_13_varias_matrizes))  
}    
    

############################################################################
#                             FUNCTION 20                                  #
############################################################################
# A idea dessa funcao e calcular o logaritmo do desvio padrao da matriz de
# covarinacia
############################################################################

log_DP_2 <- function(S)
{
  # verificar qual log usar: natural, base 10, base 2 etc.
  log_DP_2 <- log10(sqrt(Re(S[2,2]))) 
  
  return(log_DP_2)
}

############################################################################
#                             FUNCTION 21                                  #
############################################################################
# A idea dessa funcao e generalizar "log_DP_2( )", ou seja, calcular para uma
# lista de matrizes
############################################################################

log_DP_2_varias_matrizes <- function(S)
{
  log_DP_2_varias_matrizes <- c( )
  contador=1
  for( i in 1:length(S) )
  {
    log_DP_2_varias_matrizes[contador] = log_DP_2(S[[i]])
    contador=contador+1
  }
  
  
  return( log_DP_2_varias_matrizes )
}

############################################################################
#                             FUNCTION 22                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Grafico2D.s1_s2 <- function(dados)
{
  y <- log_DP_1_varias_matrizes(dados)
  x <- log_DP_2_varias_matrizes(dados)
  Grafico2D.s1_s2 <- plot( x, y, xlab = expression(log~sigma[2]), ylab = expression(log~sigma[1]), pch = 20,
                                      lab = c(2,2,4) )
}

############################################################################
#                             FUNCTION 23                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Grafico.s1_r12 <- function(dados)
{
  y <- log_DP_1_varias_matrizes(dados)
  x <- Corr_marginal_12_varias_matrizes(dados)
  Plot_Dispersao_s1_r12 <- plot( x, y, xlab = expression(rho[12]), ylab = expression(log~sigma[1]), pch = 20,
                                       lab = c(2,2,4), xlim = c(-1,1) )
}

############################################################################
#                             FUNCTION 24                                  #
############################################################################
# A ideia dessa funcao e calcular a correlacao marginal $\rho_{23}$ 
############################################################################

Corr_marginal_23 <- function(S)
{
  Corr_marginal_23 <- (Re(S[2,3])/sqrt(Re(S[2,2])*Re(S[3,3])))
  
  return(Corr_marginal_23)
}

############################################################################
#                             FUNCTION 25                                  #
############################################################################
# A idea dessa funcao e generalizar "Corr_marginal_23( )", ou seja, calcular
# para uma lista de matrizes
############################################################################

Corr_marginal_23_varias_matrizes <- function(S)
{
  Corr_marginal_23_varias_matrizes <- c( )
  contador=1
  for( i in 1:length(S) )
  {
    Corr_marginal_23_varias_matrizes[contador] = Corr_marginal_23(S[[i]])
    contador=contador+1
  }

  return( Corr_marginal_23_varias_matrizes )
}

############################################################################
#                             FUNCTION 26                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Grafico.s1_r23 <- function(dados)
{
  y <- log_DP_1_varias_matrizes(dados)
  x <- Corr_marginal_23_varias_matrizes(dados)
  Plot_Dispersao_s1_r23 <- plot( x, y, xlab = expression(rho[23]), ylab = expression(log~sigma[1]), pch = 20, lab = c(2,2,4), xlim = c(-1,1) )
}

############################################################################
#                             FUNCTION 27                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Grafico2D.r12_r23 <- function(dados)
{
  y <- Corr_marginal_12_varias_matrizes(dados)
  x <- Corr_marginal_23_varias_matrizes(dados)
  Plot_Dispersao_r12_r23 <- plot( x, y, xlab = expression(rho[23]), ylab = expression(rho[12]), pch = 20, lab = c(2,2,4),
                                        xlim = c(-1,1), ylim = c(-1,1) )
}


############################################################################
#                             FUNCTION 28                                  #
############################################################################
# A ideia dessa funcao e exibir um plot tridimensional
############################################################################


Corr3D = function(S)
{
  H = Matrizes_Corr_Complex(S)
  J = list()
  contador = 1
  for(i in 1:length(S))
  {
    J[[contador]] = H[[i]]
    contador = contador +1
    x = c()
    y = c()
    z = c()
    contador2=1
    for(j in 1: length(J))
    {
      x[contador2] = J[[j]][1,2] 
      y[contador2] = J[[j]][1,3]
      z[contador2] = J[[j]][2,3]
      contador2=contador2+1
    }
    
  }
  #return(x)
  op = par(no.readonly = TRUE)
  scatterplot3d(x,y,z, xlab = expression(rho[12]), ylab = expression(rho[23]), zlab = expression(rho[13]),
                highlight.3d = FALSE, color = "grey35", angle = 24, pch=20, cex.symbols = .5,
                mar = op$mar, cex.axis = .5, x.ticklabs = c("-1","0","1"), y.ticklabs = c(-1,0,1),
                z.ticklabs = c(-1,0,1), xlim = c(-1,1), ylim = c(-1,1), zlim = c(-1,1),
                lab = c(2,2,3) )
  
}

############################################################################
#                             FUNCTION 29                                  #
############################################################################

# A ideia da funcao e: plotar todos os graficos em um unico plot. Mas ainda
# esta em construcao. O que quero e colocar um plot de 4 linhas e 4 colunas
# ainda nao sei como fazer isto.
############################################################################
VisCovMatComplex <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f )
{
  I <- Matrizes_Cov_Complex(Imagem,linha_i,linha_f,coluna_i,coluna_f)
  size.label = 0.9 
  ppcex = 0.6   ## Pointsize
  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), 3, 4, byrow = TRUE))
  par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), mai = c(0.75,0.60,0,0), pty = "s")
  
  graficos = list(Histograma.Variancia(I), Grafico2D.s1_s2(I), Histograma.VarianciaEfetiva(I),
                  Histograma.DependenciaEfetiva(I), Histograma.Corr(I),
                  Grafico2D.r12_r23(I), Grafico.s1_r12(I), Grafico.s1_r23(I),
                  Corr3D(I))
  # return(graficos)
}
############################################################################


############################################################################
#                             FUNCTION 30                                  #
############################################################################

#                   Selecionando exibicao de graficos                      
############################################################################

VisCovMatComplex.selecionada <- function(entrada, lista.imagens )
{
  n.imagens = length(lista.imagens)
  
  graficos = c("Histograma.Variancia", "Grafico2D.s1_s2", "Histograma.VarianciaEfetiva",
               "Histograma.DependenciaEfetiva", "Histograma.Corr",
               "Grafico2D.r12_r23", "Grafico.s1_r12", "Grafico.s1_r23",
               "Corr3D")


  if(n.imagens <= 0)
  {
    print("E necessario conter imagens no formato enxuto para executar a funcao.")
  }
  if(n.imagens == 1)
  {
    print("O objetivo dessa funcao e compara 2 a 4 imagens.")
  }
  if(n.imagens >= 5)
  {
    print("Nao e adequado comparar mais de 4 graficos.")
  }



  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(n.imagens == 2)
  { # Para 2 imagens

    size.label = 0.9
    ppcex = 0.6   ## Pointsize
    layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), mai = c(0.75,0.90,0,0), pty = "s")

    I1 = Matrizes_Cov_Complex(lista.imagens[[1]][[1]],lista.imagens[[1]][[2]],
                              lista.imagens[[1]][[3]],lista.imagens[[1]][[4]],
                              lista.imagens[[1]][[5]])
    I2 = Matrizes_Cov_Complex(lista.imagens[[2]][[1]],lista.imagens[[2]][[2]],
                              lista.imagens[[2]][[3]],lista.imagens[[2]][[4]],
                              lista.imagens[[2]][[5]])
    
    
    dados1 = list(Variancia_efetiva_de_varias_matrizes(I1), Dependencia_efetiva_de_varias_matrizes(I1),
                  log_DP_1_varias_matrizes(I1), Corr_marginal_12_varias_matrizes(I1), Matrizes_Corr_Complex(I1),
                  Corr_marginal_13_varias_matrizes(I1), log_DP_2_varias_matrizes(I1),
                  Corr_marginal_23_varias_matrizes(I1) 
                  )
    dados2 = list(Variancia_efetiva_de_varias_matrizes(I2), Dependencia_efetiva_de_varias_matrizes(I2),
                  log_DP_1_varias_matrizes(I2), Corr_marginal_12_varias_matrizes(I2), Matrizes_Corr_Complex(I2),
                  Corr_marginal_13_varias_matrizes(I2), log_DP_2_varias_matrizes(I2),
                  Corr_marginal_23_varias_matrizes(I2)
                  )
    

    if(entrada == graficos[1])
    {
      graficos = list(Histograma.Variancia(I1),
                      Histograma.Variancia(I2) )
    }
    if(entrada == graficos[2])
    {
      graficos = list(Grafico2D.s1_s2(I1),
                      Grafico2D.s1_s2(I2) )
    }
    if(entrada == graficos[3])
    {
      graficos = list(Histograma.VarianciaEfetiva(I1),
                      Histograma.VarianciaEfetiva(I2) )
    }
    if(entrada == graficos[4])
    {
      graficos = list(Histograma.DependenciaEfetiva(I1),
                      Histograma.DependenciaEfetiva(I2) )
    }
    if(entrada == graficos[5])
    {
      graficos = list(Histograma.Corr(I1),
                      Histograma.Corr(I2)
      )
    }
    if(entrada == graficos[6])
    {
      graficos = list(Grafico2D.r12_r23(I1), Grafico2D.r12_r23(I2))
    }
    if(entrada == graficos[7])
    {
      graficos = list(Grafico.s1_r12(I1),
                      Grafico.s1_r12(I2) )
    }
    if(entrada == graficos[8])
    {
      graficos = list(Grafico.s1_r23(I1),
                      Grafico.s1_r23(I2) )
    }
    if(entrada == graficos[9])
    {
      graficos = list(Corr3D(I1), Corr3D(I2))
    }
  }

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(n.imagens == 3)
  { # Para 3 imagens

    size.label = 0.9
    ppcex = 0.6   ## Pointsize
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
    par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), mai = c(0.75,0.60,0,0), pty = "s")

    
    I1 = Matrizes_Cov_Complex(lista.imagens[[1]][[1]],lista.imagens[[1]][[2]],
                              lista.imagens[[1]][[3]],lista.imagens[[1]][[4]],
                              lista.imagens[[1]][[5]])
    I2 = Matrizes_Cov_Complex(lista.imagens[[2]][[1]],lista.imagens[[2]][[2]],
                              lista.imagens[[2]][[3]],lista.imagens[[2]][[4]],
                              lista.imagens[[2]][[5]])
    I3 = Matrizes_Cov_Complex(lista.imagens[[3]][[1]],lista.imagens[[3]][[2]],
                              lista.imagens[[3]][[3]],lista.imagens[[3]][[4]],
                              lista.imagens[[3]][[5]])
    
    
    dados1 = list(Variancia_efetiva_de_varias_matrizes(I1), Dependencia_efetiva_de_varias_matrizes(I1),
                  log_DP_1_varias_matrizes(I1), Corr_marginal_12_varias_matrizes(I1), Matrizes_Corr_Complex(I1),
                  Corr_marginal_13_varias_matrizes(I1), log_DP_2_varias_matrizes(I1),
                  Corr_marginal_23_varias_matrizes(I1)
                  )
    dados2 = list(Variancia_efetiva_de_varias_matrizes(I2), Dependencia_efetiva_de_varias_matrizes(I2),
                  log_DP_1_varias_matrizes(I2), Corr_marginal_12_varias_matrizes(I2), Matrizes_Corr_Complex(I2),
                  Corr_marginal_13_varias_matrizes(I2), log_DP_2_varias_matrizes(I2),
                  Corr_marginal_23_varias_matrizes(I2)
                  )
    dados3 = list(Variancia_efetiva_de_varias_matrizes(I3), Dependencia_efetiva_de_varias_matrizes(I3),
                  log_DP_1_varias_matrizes(I3), Corr_marginal_12_varias_matrizes(I3), Matrizes_Corr_Complex(I3),
                  Corr_marginal_13_varias_matrizes(I3), log_DP_2_varias_matrizes(I3),
                  Corr_marginal_23_varias_matrizes(I3)
                  )
    
    

    if(entrada == graficos[1])
    {
      graficos = list(Histograma.Variancia(I1),
                      Histograma.Variancia(I2),
                      Histograma.Variancia(I3) )
    }
    if(entrada == graficos[2])
    {
      graficos = list(Grafico2D.s1_s2(I1),
                      Grafico2D.s1_s2(I2),
                      Grafico2D.s1_s2(I3) )
    }
    if(entrada == graficos[3])
    {
      graficos = list(Histograma.VarianciaEfetiva(I1),
                      Histograma.VarianciaEfetiva(I2),
                      Histograma.VarianciaEfetiva(I3) )
    }
    if(entrada == graficos[4])
    {
      graficos = list(Histograma.DependenciaEfetiva(I1),
                      Histograma.DependenciaEfetiva(I2),
                      Histograma.DependenciaEfetiva(I3) )
    }
    if(entrada == graficos[5])
    {
      graficos = list(Histograma.Corr(I1),
                      Histograma.Corr(I2),
                      Histograma.Corr(I3) )
    }
    if(entrada == graficos[6])
    {
      graficos = list(Grafico2D.r12_r23(I1), Grafico2D.r12_r23(I2), Grafico2D.r12_r23(I3))
    }
    if(entrada == graficos[7])
    {
      graficos = list(Grafico.s1_r12(I1),
                      Grafico.s1_r12(I2),
                      Grafico.s1_r12(I3) )
    }
    if(entrada == graficos[8])
    {
      graficos = list(Grafico.s1_r23(I1),
                      Grafico.s1_r23(I2),
                      Grafico.s1_r23(I3) )
    }
    if(entrada == graficos[9])
    {
      graficos = list(Corr3D(I1), Corr3D(I2), Corr3D(I3))
    }
  }

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(n.imagens == 4)
  { # Para 4 imagens

    size.label = 0.9
    ppcex = 0.6   ## Pointsize
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE))
    par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), mai = c(0.75,0.60,0,0), pty = "s")

    I1 = Matrizes_Cov_Complex(lista.imagens[[1]][[1]],lista.imagens[[1]][[2]],
                              lista.imagens[[1]][[3]],lista.imagens[[1]][[4]],
                              lista.imagens[[1]][[5]])
    I2 = Matrizes_Cov_Complex(lista.imagens[[2]][[1]],lista.imagens[[2]][[2]],
                              lista.imagens[[2]][[3]],lista.imagens[[2]][[4]],
                              lista.imagens[[2]][[5]])
    I3 = Matrizes_Cov_Complex(lista.imagens[[3]][[1]],lista.imagens[[3]][[2]],
                              lista.imagens[[3]][[3]],lista.imagens[[3]][[4]],
                              lista.imagens[[3]][[5]])
    I4 = Matrizes_Cov_Complex(lista.imagens[[4]][[1]],lista.imagens[[4]][[2]],
                              lista.imagens[[4]][[3]],lista.imagens[[4]][[4]],
                              lista.imagens[[4]][[5]])
    
    dados1 = list(Variancia_efetiva_de_varias_matrizes(I1), Dependencia_efetiva_de_varias_matrizes(I1),
                  log_DP_1_varias_matrizes(I1), Corr_marginal_12_varias_matrizes(I1), Matrizes_Corr_Complex(I1),
                  Corr_marginal_13_varias_matrizes(I1), log_DP_2_varias_matrizes(I1),
                  Corr_marginal_23_varias_matrizes(I1)
                  )
    dados2 = list(Variancia_efetiva_de_varias_matrizes(I2), Dependencia_efetiva_de_varias_matrizes(I2),
                  log_DP_1_varias_matrizes(I2), Corr_marginal_12_varias_matrizes(I2), Matrizes_Corr_Complex(I2),
                  Corr_marginal_13_varias_matrizes(I2), log_DP_2_varias_matrizes(I2),
                  Corr_marginal_23_varias_matrizes(I2)
                  )
    dados3 = list(Variancia_efetiva_de_varias_matrizes(I3), Dependencia_efetiva_de_varias_matrizes(I3),
                  log_DP_1_varias_matrizes(I3), Corr_marginal_12_varias_matrizes(I3), Matrizes_Corr_Complex(I3),
                  Corr_marginal_13_varias_matrizes(I3), log_DP_2_varias_matrizes(I3),
                  Corr_marginal_23_varias_matrizes(I3)
                  )
    dados4 = list(Variancia_efetiva_de_varias_matrizes(I4), Dependencia_efetiva_de_varias_matrizes(I4),
                  log_DP_1_varias_matrizes(I4), Corr_marginal_12_varias_matrizes(I4), Matrizes_Corr_Complex(I4),
                  Corr_marginal_13_varias_matrizes(I4), log_DP_2_varias_matrizes(I4),
                  Corr_marginal_23_varias_matrizes(I4)
                  )
    
    
    
    
    if(entrada == graficos[1])
    {
      graficos = list(Histograma.Variancia(I1),
                      Histograma.Variancia(I2),
                      Histograma.Variancia(I3),
                      Histograma.Variancia(I4) )
    }
    if(entrada == graficos[2])
    {
      graficos = list(Grafico2D.s1_s2(I1),
                      Grafico2D.s1_s2(I2),
                      Grafico2D.s1_s2(I3),
                      Grafico2D.s1_s2(I4) )
    }
    if(entrada == graficos[3])
    {
      graficos = list(Histograma.VarianciaEfetiva(I1),
                      Histograma.VarianciaEfetiva(I2),
                      Histograma.VarianciaEfetiva(I3),
                      Histograma.VarianciaEfetiva(I4) )
    }
    if(entrada == graficos[4])
    {
      graficos = list(Histograma.DependenciaEfetiva(I1),
                      Histograma.DependenciaEfetiva(I2),
                      Histograma.DependenciaEfetiva(I3),
                      Histograma.DependenciaEfetiva(I4) )
    }
    if(entrada == graficos[5])
    {
      graficos = list(Histograma.Corr(I1),
                      Histograma.Corr(I2),
                      Histograma.Corr(I3),
                      Histograma.Corr(I4) )
    }
    if(entrada == graficos[6])
    {
      graficos = list(Grafico2D.r12_r23(I1), Grafico2D.r12_r23(I2),
                      Grafico2D.r12_r23(I3),Grafico2D.r12_r23(I4))
    }
    if(entrada == graficos[7])
    {
      graficos = list(Grafico.s1_r12(I1),
                      Grafico.s1_r12(I2),
                      Grafico.s1_r12(I3),
                      Grafico.s1_r12(I4) )
    }
    if(entrada == graficos[8])
    {
      graficos = list(Grafico.s1_r23(I1),
                      Grafico.s1_r23(I2),
                      Grafico.s1_r23(I3),
                      Grafico.s1_r23(I4) )
    }
    if(entrada == graficos[9])
    {
      graficos = list(Corr3D(I1), Corr3D(I2),
                      Corr3D(I3),Corr3D(I4))
    }
  }

}

#########################################################################################





######################################################################################################################################
#                                                     Medidas resumo básicas                                                         #
######################################################################################################################################

fivenum_log_DP_1 <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_log_DP_1 <- fivenum( log_DP_1_varias_matrizes(  Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f) ) )
  return(fivenum_log_DP_1)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fivenum_log_DP_2 <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_log_DP_2 <- fivenum( log_DP_2_varias_matrizes(  Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f) ) )
  return(fivenum_log_DP_2)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fivenum_Variancia_Efetiva <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_Variancia_Efetiva <- fivenum( Variancia_efetiva_de_varias_matrizes( Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)  ) )
  return(fivenum_Variancia_Efetiva)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fivenum_Dependencia_efetiva <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_Dependencia_efetiva <- fivenum( Dependencia_efetiva_de_varias_matrizes( Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)  ) )
  return(fivenum_Dependencia_efetiva)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fivenum_Corr_12 <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_Corr_12 <- fivenum( Corr_marginal_12_varias_matrizes(  Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f) ) )
  return(fivenum_Corr_12)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fivenum_Corr_23 <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_Corr_23 <- fivenum( Corr_marginal_23_varias_matrizes(  Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f) ) )
  return(fivenum_Corr_23)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fivenum_Corr_13 <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  fivenum_Corr_13 <- fivenum( Corr_marginal_13_varias_matrizes(  Matrizes_Cov_Complex(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f) ) )
  return(fivenum_Corr_13)
}




