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
# A ideia da funcao e montar uma matriz diagonal, onde os elementos da diago
# nal sao os autovalores.
# Para isso estou calculando os autovalores com a funcao do R, "eigen()" e 
# preencho a diagonal de uma matriz nula com esses autovalores calculados 
############################################################################

# Diag_mat_Complex <- function(M)
# {
#   autovalores <- eigen(M)$values
#   Diag_mat_Complex <- diag( autovalores ) # diag( ) vai colocar os autovalores na diagonal da matriz
# #   Diag_mat_Complex <- matrix(0, nrow = 3, ncol = 3)
# #   for(i in 1:3)
# #   {
# #     Diag_mat_Complex[i,i] <- autovalores[i]
# #   }
#   
#   return(Diag_mat_Complex) 
# }

############################################################################
#                             FUNCTION 7                                   #
############################################################################
# A ideia desta funcao e calcular o determinante da matriz complexa, mas dia-
# gonalizando antes.
############################################################################

# det_Cov_Complex_Diag <- function(H)
# {
#   det_Cov_Complex_Diag <- det( Diag_mat_Complex(H) )
#   
#   return(det_Cov_Complex_Diag)  
# }

############################################################################
#                             FUNCTION 8                                   #
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
#                             FUNCTION 9                                   #
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
#                             FUNCTION 10                                  #
############################################################################
# A ideia desta funcao e plotar um histograma das "Variancia_efetiva()"
############################################################################

Histograma_Variancia_Efetiva <- function(dados)
{
  dados_para_saida <- Variancia_efetiva_de_varias_matrizes(dados)
  Histograma_Variancia_Efetiva <- hist(dados_para_saida, breaks = "FD", probability = T,
                                       xlab ="Effective Variance", main = "")
}
  
############################################################################
#                             FUNCTION 11                                  #
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
#                             FUNCTION 12                                  #
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
#                             FUNCTION 13                                  #
############################################################################
# A ideia desta funcao e plotar um histograma das "Dependencia_efetiva()"
############################################################################

Histograma_Dependencia_Efetiva <- function(dados)
{
  dados_para_saida <- Dependencia_efetiva_de_varias_matrizes(dados)
  Histograma_Dependencia_Efetiva <- hist(dados_para_saida, breaks = "FD", probability = T,
                                         xlab ="Effective Dependence", main = "")
}

############################################################################
#                             FUNCTION 14                                  #
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
#                             FUNCTION 15                                  #
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
#                             FUNCTION 16                                  #
############################################################################
# A ideia dessa funcao e plotar o log do desvio padrao de todas as matrizes
# escolhidas 
############################################################################

Histograma_log_DP_1 <- function(dados)
{
  dados_para_saida <- log_DP_1_varias_matrizes(dados)
  Histograma_log_DP_1 <- hist(dados_para_saida, breaks = "FD", probability = T,
                              xlab = expression(log~sigma[1]), main = "")
}

############################################################################
#                             FUNCTION 17                                  #
############################################################################
# A ideia dessa funcao e calcular a correlacao marginal $\rho_{12}$ 
############################################################################

Corr_marginal_12 <- function(S)
{
  Corr_marginal_12 <- (Re(S[1,2])/sqrt(Re(S[1,1])*Re(S[2,2])))
  
  return(Corr_marginal_12)
}

############################################################################
#                             FUNCTION 18                                  #
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
#                             FUNCTION 19                                  #
############################################################################
# A ideia dessa funcao e plotar o log do desvio padrao de todas as matrizes
# escolhidas 
############################################################################

Histograma_Corr_marginal_12 <- function(dados)
{
  dados_para_saida <- Corr_marginal_12_varias_matrizes(dados)
  Histograma_Corr_marginal_12 <- hist(dados_para_saida, breaks ="FD", probability = T, xlab = expression(rho[12]) , main = "", xlim=c(-1,1))
  lines(density(dados_para_saida), col="green")
}

############################################################################
#                             FUNCTION 20                                  #
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
#                             FUNCTION 21                                  #
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
#                             FUNCTION 22                                  #
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
#                             FUNCTION 23                                  #
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
#                             FUNCTION 24                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Plot_Dispersao_s1_s2 <- function(dados)
{
  y <- log_DP_1_varias_matrizes(dados)
  x <- log_DP_2_varias_matrizes(dados)
  Plot_Dispersao_s1_s2 <- plot( x, y, xlab = expression(log~sigma[2]), ylab = expression(log~sigma[1]), pch = 20,
                                      lab = c(2,2,4) )
}

############################################################################
#                             FUNCTION 25                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Plot_Dispersao_s1_r_12 <- function(dados)
{
  y <- log_DP_1_varias_matrizes(dados)
  x <- Corr_marginal_12_varias_matrizes(dados)
  Plot_Dispersao_s1_r12 <- plot( x, y, xlab = expression(rho[12]), ylab = expression(log~sigma[1]), pch = 20,
                                       lab = c(2,2,4), xlim = c(-1,1) )
}

############################################################################
#                             FUNCTION 26                                  #
############################################################################
# A ideia dessa funcao e calcular a correlacao marginal $\rho_{23}$ 
############################################################################

Corr_marginal_23 <- function(S)
{
  Corr_marginal_23 <- (Re(S[2,3])/sqrt(Re(S[2,2])*Re(S[3,3])))
  
  return(Corr_marginal_23)
}

############################################################################
#                             FUNCTION 27                                  #
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
#                             FUNCTION 28                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Plot_Dispersao_s1_r_23 <- function(dados)
{
  y <- log_DP_1_varias_matrizes(dados)
  x <- Corr_marginal_23_varias_matrizes(dados)
  Plot_Dispersao_s1_r23 <- plot( x, y, xlab = expression(rho[23]), ylab = expression(log~sigma[1]), pch = 20, lab = c(2,2,4), xlim = c(-1,1) )
}

############################################################################
#                             FUNCTION 29                                  #
############################################################################
# A ideia dessa funcao e exibir um plot de dispersao entre o s1 e o s2
############################################################################

Plot_Dispersao_r12_r_23 <- function(dados)
{
  y <- Corr_marginal_12_varias_matrizes(dados)
  x <- Corr_marginal_23_varias_matrizes(dados)
  Plot_Dispersao_r12_r23 <- plot( x, y, xlab = expression(rho[23]), ylab = expression(rho[12]), pch = 20, lab = c(2,2,4),
                                        xlim = c(-1,1), ylim = c(-1,1) )
}


############################################################################
#                             FUNCTION 30                                  #
############################################################################
# A ideia dessa funcao e exibir um plot tridimensional
############################################################################


CORRELACAO_ESPACIAL = function(S)
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
#                             FUNCTION 31                                  #
############################################################################

# A ideia da funcao e: plotar todos os graficos em um unico plot. Mas ainda
# esta em construcao. O que quero e colocar um plot de 4 linhas e 4 colunas
# ainda nao sei como fazer isto.
############################################################################
VisCovMatComplex <- function(Imagem  , linha_i , linha_f ,  coluna_i , coluna_f)
{
  I <- Matrizes_Cov_Complex(Imagem,linha_i,linha_f,coluna_i,coluna_f)
  size.label = 0.9 
  ppcex = 0.6   ## Pointsize
  layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), 3, 4, byrow = TRUE))
  par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), mai = c(0.75,0.60,0,0), pty = "s")
  
  graficos = list(Histograma_log_DP_1(I), Plot_Dispersao_s1_s2(I), Histograma_Variancia_Efetiva(I),
               Histograma_Dependencia_Efetiva(I), Histograma_Corr_marginal_12(I),
               Plot_Dispersao_r12_r_23(I), Plot_Dispersao_s1_r_12(I), Plot_Dispersao_s1_r_23(I),
               CORRELACAO_ESPACIAL(I))
  
  VisCovMatComplex = list( )
  contador =1
  for(i in 1: length(graficos))
  {
    VisCovMatComplex[[contador]] = graficos[[i]]
    contador = contador+1
  }
}



##############################################################################################
# obter.selecao <- function (nome)
# {
#     graficos = c("Histograma_log_DP_1", "Plot_Dispersao_s1_s2", "Histograma_Variancia_Efetiva",
#                  "Histograma_Dependencia_Efetiva", "Histograma_Corr_marginal_12",
#                  "Plot_Dispersao_r12_r_23", "Plot_Dispersao_s1_r_12", "Plot_Dispersao_s1_r_23",
#                  "CORRELACAO_ESPACIAL")
#   numero.selecionado = which(graficos == nome)
#   if (length(numero.selecionado) == 0){
#     stop(cat("Wrong specification of conditions!\n"))
#   }
#   return(numero.selecionado)
# }
##############################################################################################


##############################################################################################
painel.selecionado = function(grafico.selecionado, lista.viscovmatcomplex)
{
  lista.viscovmatcomplex = list(VisCovMatComplex1, VisCovMatComplex2, VisCovMatComplex3)
  n.viscovmatcomplex = length(lista.viscovmatcomplex)
  layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
  # par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), mai = c(0.75,0.60,0,0), pty = "s")
  
  graficos = c("Histograma_log_DP_1", "Plot_Dispersao_s1_s2", "Histograma_Variancia_Efetiva",
               "Histograma_Dependencia_Efetiva", "Histograma_Corr_marginal_12",
               "Plot_Dispersao_r12_r_23", "Plot_Dispersao_s1_r_12", "Plot_Dispersao_s1_r_23",
               "CORRELACAO_ESPACIAL")
  
  for(i in 1:length(graficos))
  {
    for(j in 1:n.viscovmatcomplex)
    {
      if(grafico.selecionado[i] == graficos[1])
      {
        Histograma_log_DP_1(I[j])
      }
      if(grafico.selecionado[i] == graficos[2])
      {
        Plot_Dispersao_s1_s2(I[j])
      }
      if(grafico.selecionado[i] == graficos[3])
      {
        Histograma_Variancia_Efetiva(I[j])
      }
      if(grafico.selecionado[i] == graficos[4])
      {
        Histograma_Dependencia_Efetiva(I[j])
      }
      if(grafico.selecionado[i] == graficos[5])
      {
        Histograma_Corr_marginal_12(I[j])
      }
      if(grafico.selecionado[i] == graficos[6])
      {
        Plot_Dispersao_r12_r_23(I[j])
      }
      if(grafico.selecionado[i] == graficos[7])
      {
        Plot_Dispersao_s1_r_12(I[j])
      }
      if(grafico.selecionado[i] == graficos[8])
      {
        Plot_Dispersao_s1_r_23(I[j])
      }
      if(grafico.selecionado[i] == graficos[9])
      {
        CORRELACAO_ESPACIAL(I[j])
      }
    }
  }
}
##############################################################################################












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



