#= Algoritmos desenvolvidos em Análise Numérica I =#

function elim_gauss(A, b)
    A, n = [A b], size(A,1)
    A = float(A)
    if A[n,n] == 0 
        print(" ∄ unique solution ")
    else
        for k = 1:n-1
            for i = k+1:n
                m = A[i,k]/A[k,k]
                A[i,k] = 0
                for j = k+1:n+1
                    A[i,j] = A[i,j] - m*A[k,j]
                end
            end
        end
    end
    return A
end

function subst_rever(A)
    n = size(A,1)
    x = zeros(n)
    x[n] = A[n,n+1]/A[n,n]
    for i = n-1:-1:1
        s = A[i,n+1]
        for j = i+1:n
            s = s - A[i,j]*x[j]
        end
        x[i] = s/A[i,i]
    end
    return x
end

function subst_prog(A)
    n = size(A,1)
    x = zeros(n)
    x[1] = A[1,n+1]/A[1,n]
    for i = 1:n
        s = 0
        for j = 1:i-1
            s += A[i,j]*x[j]
        end
        x[i] = (A[i,n+1] - s)/A[i,i]
    end
    return x
end

function elim_gauss_pivot_parcial(A,b)
    A, n = [A b], size(A,1)
    A = float(A)
    for k = 1:n-1
        aux = abs(A[k,k])
        for j = k:n
            if abs(A[j,k]) > aux
                aux = abs(A[j,k])
                s = j
                A[k,:], A[s,:] = A[s,:], A[k,:]
            end
        end
        for i = k+1:n
            m = A[i,k]/A[k,k]
            A[i,k] = 0
            for j = k+1:n+1
                A[i,j] = A[i,j] - m*A[k,j]
            end
        end
    end
    return A
end
            
function fator_lu(A)
    n = size(A,1)
    A = float(A)
    for i = 1:n-1
        for j = i+1:n
            m = A[j,i]/A[i,i]
            A[j,i] = m
            for k = i+1:n
                A[j,k] = A[j,k] - m*A[i,k]
            end
        end
    end
    return A
end

function decomp_lu(A)
    n = size(A,1)
    L, U = zeros(n,n), zeros(n,n)
    for j = 1:n
        for i = 1:n
            if i > j
                L[i,j] = A[i,j]
            elseif i <= j
                U[i,j] = A[i,j]
            end
            L[i,i] = 1
        end
    end
    return L, U
end

function fator_lu_pivot(A)
    A = float(A)
    n = size(A,1)
    p = zeros(n)
    for i = 1:n
        p[i] = i
    end
    for k = 1:n
        piv = abs(A[k,k])
        r = k
        for i = k+1:n
            if abs(A[i,k]) > piv
                piv = abs(A[i,k])
                r = i
            end
        end
        if piv == 0 
            print("A é singular")
            break
        end
        if r != k
            aux = p[k]
            p[k] = p[r]
            p[r] = aux 
        end
        for j = 1:n
            aux = A[k,j]
            A[k,j] = A[r,j]
            A[r,j] = aux
        end
        for i = k+1:n
            m = A[i,k]/A[k,k]
            A[i,k] = m
            for j = k+1:n
                A[i,j] = A[i,j] - m*A[k,j]
            end
        end
    end
    return A
end

function fator_cholesky(A)
    A = float(A)
    n = size(A,1)
    L = zeros(n,n)
    if isposdef(A) == true
        for i = 1:n
            L[i,i] = sqrt(A[i,i])
            for j = i+1:n
                L[j,i] = A[j,i]/L[i,i]
                for k = i+1:j
                    A[j,k] = A[j,k] - L[j,i]*L[k,i]
                end
            end
        end
    else
        @error("Matrix is not pos. def.")
    end
    return L       
end

function ldu(A)
    n = size(A,1)
    L = zeros(n,n)
    D = zeros(n,n)
    U = zeros(n,n)
    for i = 1:n
        D[i,i] = A[i,i]
        for j = 1:n
            if i < j
                U[i,j] = A[i,j]
            elseif i > j
                L[i,j] = A[i,j]
            end
        end
    end
    return L, D, U
end

function jacobi(L, D, U, b, sol)
    L, D, U = float(L), float(D), float(U)
    T, c = -inv(D)*(L+U), inv(D)*b
    iter, x = 1, [sol T*sol+c]
    while iter < 100
        if norm(x[:,iter+1] - x[:,iter], Inf)/norm(x[:,iter+1], Inf) <= 1e-4
            return x[:,iter], iter
            break
        end
        x = [x T*x[:,iter+1]+c]
        iter += 1
    end             
    return x[:,iter], iter
end

function gauss_seidel(L, D, U, b, sol)
    L, D, U = float(L), float(D), float(U)
    T, c = -inv(D+L)*U, inv(D+L)*b
    iter, x = 1, [sol T*sol+c]
    while iter < 100
        if norm(x[:,iter+1] - x[:,iter], Inf)/norm(x[:,iter+1], Inf) <= 1e-4
            return x[:,iter], iter
            break
        end
        x = [x T*x[:,iter+1]+c]
        iter += 1
    end
    return x[:,iter], iter
end

function sor(L, D, U, b, sol, ω = 1.22)
    L, D, U = float(L), float(D), float(U)
    T, c = -inv(D+L)*U, inv(D+L)*b
    iter, x = 1, [sol (1 - ω)*sol + ω*(T*sol + c) ]
    while iter < 100
        if norm(x[:,iter+1] - x[:,iter], Inf)/norm(x[:,iter+1], Inf) <= 1e-4
            return x[:,iter], iter
            break
        end
        tv = (1 - ω)*x[:,iter+1] + ω*(T*x[:,iter+1] + c) 
        x = [x tv]
        iter += 1
    end
    return x[:,iter], iter
end

function h(n)
    H = zeros(n,n)
    for i = 1:n
        for j = 1:n
            H[i,j] = 1/(i + j - 1)
        end
    end
    return H
end


#== Projeto 1 ==#
using LinearAlgebra

#=
a) 

Para n = 4:
  b = h(4)*ones(4)
  palp = [0.882; 0.216; 0.731; 0.817]

- Fatoração LU:

  L, U = decomp_lu(fator_lu(h(4)))
  y = subst_prog([L b])
  x = subst_rever([U y])

  Ou seja, 
  x = [1.0000000000000058
       0.9999999999999312
       1.0000000000001652
       0.9999999999998932].

  Erro relativo:
  norm(ones(4) - x)/norm(x) = 1.04e-13

  Valor do Resíduo:
  norm(h(4)*x - b) = 0

- Fatoração de Cholesky:

  L = fator_cholesky(h(4))
  y = subst_prog([L b])
  x = subst_rever([L' y])

  Ou seja,
  x = [1.0000000000000127
       0.9999999999998604
       1.0000000000003235
       0.9999999999997958].

  Erro relativo:
  norm(ones(4) - x)/norm(x) = 2.03e-13
   
  Valor do Resíduo:
  norm(h(4)*x - b) = 0

- Método iterativo de Jacobi:

  L, D, U = ldu(h(4))
  jacobi(L, D, U, b, palp)
  
  Que estoura com palpite inicial 
  diferente da solução.
  x = [1.0730154686196049e40
       2.1474601258584184e40
       2.7426054310094827e40
       3.130767656863611e40],

  iter = maxiter = 100

  Erro relativo:
  norm(ones(4) - x)/norm(x) = 1

  Valor do Resíduo:
  norm(h(4)*x - b) = 5.27e40

- Metodo iterativo de Gauss-Seidel:

  L, D, U = ldu(h(4))
  gauss_seidel(L, D, U, b, palp)

  Que converge para
  x = [1.0136268960842676 
       0.9418404443889442
       1.0549674239623932
       0.993447648418551],

  iter = 12.

  Erro relativo:
  norm(ones(4) - x)/norm(x) = 0.0406

  Valor do Resíduo:
  norm(h(4)*x - b) = 0.00124

- Método SOR (ω = 1.22):
  
  L, D, U = ldu(h(4))
  sor(L, D, U, b, palp)

  Que converge para
  x = [1.0019675140001632 
       0.9817234852512118 
       1.039685766823361 
       0.9758439098540963],

  iter = 71.

  Erro relativo:
  norm(ones(4) - x)/norm(x) = 0.025

  Valor do Resíduo:
  norm(h(4)*x - b) = 2.626e-5


Para n = 10:
  b = h(10)*ones(10)
  palp = [0.6426320683292519; 0.8980540353126407; 0.6765257582390736; 0.10276221103749394; 0.7879925034971433; 0.551805814105925; 0.6476579742426902; 0.21412847698232884; 0.1797099189242406; 0.4442178522435265]

- Fatoração LU:
  
  L, U = decomp_lu(fator_lu(h(10)))
  y = subst_prog([L b])
  x = subst_rever([U y])
  
  Ou seja, 
  x = [0.9999999986378105
       1.0000001160832035
       0.9999975518008497
       1.000022090797208
       0.9998952497304442
       1.0002865848469427
       0.9995316657163531
       1.0004510716338584
       0.9997638743427951
       1.0000517972798855].

  Erro relativo:
  norm(ones(10) - x)/norm(x) = 0.00024

  Valor do Resíduo:
  norm(h(10)*x - b) = 7.2e-16


- Fatoração de Cholesky:

  L = fator_cholesky(h(10))
  y = subst_prog([L b])
  x = subst_rever([L' y])
  
  Ou seja, 
  x = [1.000000000480556
       0.9999999617615873
       1.0000007598529697
       0.9999934922526903
       1.0000294668952292
       0.9999226288348547
       1.0001218624860921
       0.9998864693536743
       1.0000576642113488
       0.999987693759108]

  Erro relativo:
  norm(ones(10) - x)/norm(x) = 6.17e-5

  Valor do Resíduo:
  norm(h(10)*x - b) = 5.98e-16
     

- Método iterativo de Jacobi:

  L, D, U = ldu(h(10))
  jacobi(L, D, U, b, palp)

  x = [1.8806243776632377e87 
        4.400184459781355e87
        6.159434840334995e87
        7.494713256032087e87
        8.554287053558874e87
        9.42035738957688e87
        1.0143838974049734e88
        1.0758528720123265e88
        1.128797925903644e88
        1.1749216633807762e88],

  iter = maxiter = 100.

  Erro relativo:
  norm(ones(10) - x)/norm(x) = 1

  Valor do Resíduo:
  norm(h(10)*x - b) = 3.033e88

- Metodo iterativo de Gauss-Seidel:

  L, D, U = ldu(h(10))
  gauss_seidel(L, D, U, b, palp)

  x = [1.0208104336462795
       0.9410075787965178
       1.0284008008666843
       0.6391421501849809
       1.4027398135145672
       1.1842128920800488
       1.2656900871966894
       0.801914189020182
       0.7302830489704283
       0.9554683486903366],

  iter = 74.

  Erro relativo:
  norm(ones(10) - x)/norm(x) = 0.222

  Valor do Resíduo:
  norm(h(10)*x - b) = 0.001

- Método SOR (ω = 1.22):

  L, D, U = ldu(h(10))
  sor(L, D, U, b, palp)

  Que converge para
  x = [ 1.0017392370010025
        1.0366942047279193
        0.9811255762324658
        0.5839808699741541
        1.3679312603796379
        1.1716355923730668
        1.2710946742367693
        0.8204324868465193
        0.7578187387153918
        0.9888906016538977],

  iter = maxiter = 100.

  Erro relativo:
  norm(ones(10) - x)/norm(x) = 0.219

  Valor do Resíduo:
  norm(h(10)*x - b) = 0.0004


Para  n = 20:
  b = h(20)*ones(20)
  palp = [0.9949754299792115; 0.15579008325848775; 0.17339790755815487; 0.5625497178604368; 0.549852344548744; 0.5222005499577691; 0.5952166662508702; 0.8164207312851908; 0.9918719969153333; 0.7681494805378093; 0.20636291545042207; 0.27124237698305054; 0.46767567791063325; 0.6423547702717181; 0.2339783599219334; 0.9596824374462323; 0.6125221087422925; 0.32637723891238446; 0.9908583844281353; 0.36700013183978375]

- Fatoração LU:

  L, U = decomp_lu(fator_lu(h(20)))
  y = subst_prog([L b])
  x = subst_rever([U y])

  Ou seja, 
  x = [1.0000005589224557
       0.999917857960547
       1.002931164162727
       0.9563134455258951
       1.3282125376864222
      -0.29529866358027745
       3.2585921488966028
       2.0696421854628992
      -9.399109111284226
       9.212338183503137
      27.984046204292646
     -65.19362202748488
      57.459785845944836
     -22.07845185963078
      18.62675954068388
       1.5450070428946734
     -58.107340989513474
      85.20186941082892
     -46.51983575910282
      10.94824265431891]

  Erro relativo:
  norm(ones(20) - x)/norm(x) = 0.9995

  Valor do Resíduo:
  norm(h(20)*x - b) = 2.45e-15


- Fatoração de Cholesky:

  Note que para n = 20, h não é positiva
  definida, logo, o método não funciona.

- Método iterativo de Jacobi:

  L, D, U = ldu(h(20))
  jacobi(L, D, U, b, palp)

  x = [2.086221569933349e119
       5.309416883039617e119
       7.847942155171201e119
       9.950147971625833e119
       1.1738888329531537e120
       1.3288825879961668e120
       1.4650021478368624e120
       1.585812534252713e120
       1.6939601993163838e120
       1.791471014673531e120
       1.879933242438657e120
       1.9606161867008885e120
       2.0345504524900745e120
       2.1025840964713976e120
       2.1654230112739103e120
       2.223660649099877e120
       2.2778003327335996e120
       2.328272288529837e120
       2.3754468435289605e120
       2.4196447845810397e120],

  iter = maxiter = 100.

  Erro relativo:
  norm(ones(20) - x)/norm(x) = 1

  Valor do Resíduo:
  norm(h(20)*x - b) = 8.74e120

- Metodo iterativo de Gauss-Seidel:

  L, D, U = ldu(h(20))
  gauss_seidel(L, D, U, b, palp)

  x = [0.9743860770995272
       1.2548596711787772
       0.6583235597365864
       0.9283245897887146
       0.9246411703451116
       0.932475123407934
       1.0380481110397102
       1.2816598887441226
       1.469122601421961
       1.248743772780494
       0.6836614811887469
       0.7404072077020234
       0.9253215887519203
       1.086216521081417
       0.6626325635130255
       1.3723281006142405
       1.0088122946676565
       0.7062919234862035
       1.3546092906549823
       0.7149599147611585],
    
  iter = 63.

  Erro relativo:
  norm(ones(20) - x)/norm(x) = 0.2462

  Valor do Resíduo:
  norm(h(20)*x - b) = 0.001

- Método SOR (ω = 1.22):

  L, D, U = ldu(h(20))
  sor(L, D, U, b, palp)

  x =  [0.9690829022895874
        1.2203360377877304
        0.7594417847410774
        0.9584679259625388
        0.9070341397425795
        0.895200542147438
        0.9973609499343907
        1.2454601362132798
        1.4406480588172907
        1.228828997691549
        0.6719782410409058
        0.7361198378996933
        0.9274241926623361
        1.0937002450379607
        0.6745629768999494
        1.3878756173997084
        1.0272586603573701
        0.7270250942925948
        1.377113257701962
        0.7388028959877653],
    
  iter = maxiter = 100.

  Erro relativo:
  norm(ones(20) - x)/norm(x) = 0.235

  Valor do Resíduo:
  norm(h(20)*x - b) = 0.0004

b)

cond(h(4)) = norm(h(4))*norm(inv(h(4))) = 15513.73

cond(h(10)) = norm(h(10))*norm(inv(h(10))) = 1.602e13

cond(h(20)) = norm(h(20))*norm(inv(h(20))) = 1.314e18

A matriz de Hilbert para n ≥ 3 é mal condicionada.

c)

Para matrizes de Hilbert, os métodos de fatoração LU 
e decomposição de Cholesky se apresentam mais confiáveis.
Considerando as matrizes em estudo, o número de condicionamento 
delas explode conforme n aumenta, ou seja, em geral são matrizes 
mal condicionadas.


#== Exercício 2 ==#

M1 = [3 0 4; 7 4 2; -1 1 2] ambos divergem (ρ(t) > 1)
M2 = [-3 3 -6; -4 7 -8; 5 7 -9] J converge G diverge
M3 = [4 1 1; 2 -9 0; 0 -8 -6] ambos convergem
M4 = [7 6 9; 4 5 -4; -7 -3 8] ambos convergem

a)
O número máximo de iterações escolhido foi 100.

-Método de Jacobi:

 M1: x = [10064.44
          33805.36
         -23249.73],

     com iter = maxiter = 100.

     Resíduo: 172622.94.

 M2: x = [-0.0409
          -0.0444
          -0.3684],

     com iter = 32

     Resíduo: 0.0004628.

 M3: x = [0.254706
          0.056581
         -0.075509],

     com iter = 11.

     Resíduo: 0.0004512.

 M4: x = [0.073903
         -0.010552
          0.060640],

     com iter = 18

     Resíduo: 0.0006527.

-Método de Gauss-Seidel:

 M1: x = [8.55e18
         -1.17e19
          1.01e19],

     com iter = maxiter = 100.

     Resíduo: 7.4e19.

 M2: x = [-6359.55
          -0.44444
          -3533.43],

     com iter = maxiter = 100.

     Resíduo: 67128.2.

 M3: x = [0.254629
          0.056584
         -0.075445],

     com iter = 3.

     Resíduo: 0.0003429.

 M4: x = [0.073833
         -0.010490
          0.060670],

     com iter = 31

     Resíduo: 0.0002134.

b)

- Método de Jacobi:

  M1: ρ(Tj) = 1.125

  M2: ρ(Tj) = 0.813
 
  M3: ρ(Tj) = 0.443

  M4: ρ(Tj) = 0.641

- Método de Gauss_Seidel:

  M1: ρ(Tg) = 1.583

  M2: ρ(Tg) = 1.111

  M3: ρ(Tg) = 0.018

  M4: ρ(tg) = 0.774

c)

Para as matrizes M3 e M4 ambos os métodos con-
vergem, dado que o raio espectral de suas ma-
trizes de iteração (Tj e Tg) são menores que 1.

d)

Ocorre apenas para M2.

e)

Para nenhuma das matrizes Mi, i =1:4, o método 
de Gauss Seidel converge enquanto o de Jacobi
diverge.

f)

Comparando somente o número de iterações entre M3 
e M4 entre os dois métodos, para M3 têm-se que o 
método de Gauss Seidel converge mais rápido, mas 
para M4 o contrário ocorre.

Utilizando somente as quatro amostras de raio es-
pectral, concluí-se que, quanto menor for o raio
espectral da matriz de iteração do método iterati-
vo, menos iterações serão necessárias para atingir
a convergência, neste caso considerada como
||xᵏ⁺¹ - xᵏ||₂ ≤ 10⁻⁴.

Veja também:

    
     |  iter(Tj)  |  ρ(Tj)  |  iter(Tg)  |  ρ(Tg)  |
====================================================
  M3 |    11      |  0.443  |     03     |  0.018  |
====================================================
  M4 |    18      |  0.641  |     31     |  0.774  |
====================================================

#== Exercício 3 ==#

E = 20V
R1 = 10Ω
R2 = R3 = R4 = R5 = 100Ω

a)
Sendo as equações obtidas através da Lei de Kirchoff:
I1R1 + I4R4 - E = 0
I1R1 + I5R5 - I2R2 = 0
I5R5 + I3R3 - I4R4 = 0
I6 = I1 + I2 = I5 + I4 + I2
I1 = I5 + I4
I3 = I2 + I5

As resistências são dadas, portanto
as incógnitas são as correntes.
Organizando as equações, temos:

I1*R1 +  0*I2 +  0*I3 + I4*R4 +  0*I5 + 0*I6  = E
I1*R1 - I2*R2 +  0*I3 +  0*I4 + I5*R5 + 0*I6 = 0
 0*I1 +  0*I2 + I3*R3 - I4*R4 + I5*R5 + 0*I6 = 0
  -I1 -    I2 +  0*I3 +  0*I4 +  0*I5 +   I6 = 0
   I1 +  0*I2 +  0*I3 -    I4 -    I5 + 0*I6 = 0
 0*I1 +    I2 -    I3 +  0*I4 +    I5 + 0*I6 = 0

[10    0   0  100   0  0] * [I1] = [20]
[10 -100   0    0 100  0]   [I2]   [0]
[ 0    0 100 -100 100  0]   [I3]   [0]
[-1   -1   0    0   0  1]   [I4]   [0]
[ 1    0   0   -1  -1  0]   [I5]   [0]
[ 0    1  -1    0   1  0]   [I6]   [0]

Ou seja, queremos resolver Ax = b, em que
A = [10    0   0  100   0  0 
     10 -100   0    0 100  0  
      0    0 100 -100 100  0  
     -1   -1   0    0   0  1  
      1    0   0   -1  -1  0  
      0    1  -1    0   1  0]

e

b = [20 0 0 0 0 0]'


b) Não. Note que na primeira linha, temos que 
|a₁₁| < |(0 + 0 + 100 + 0 + 0)|.
Portanto a matriz dos coeficientes não é estrita- 
mente diagonalmente dominante.

c)

Fazendo permutações para que tenhamos 
Ã = [1    0   0   -1  -1  0
     0    1  -1    0   1  0
     0    0 100 -100 100  0 
    10    0   0  100   0  0 
    10 -100   0    0 100  0  
    -1   -1   0    0   0  1]
 
e

b̃ = [0 0 0 20 0 0]', e ainda,

usando como palpite inicial para todos os métodos
x⁰ = [0 0 0 0 0 0]'

- Método de Jacobi:
    x = [-1.5928036184242602e7
          1.0965266021180712e7 
          1.3507762466757648e7 
         -313328.587091269 
         -1.4490636409214836e7 
         -1.1044019751210878e7]

    iter = maxiter = 100.

- Método de Gauss-Seidel:
    x = [0.22890581747549732
         0.07412593324575897
         0.12540231600617097
         0.17710941825245027
         0.05123535149820923 
         0.3030317507212563]

    iter = maxiter = 100

- Método SOR (ω = 1.22):
    x = [-3.192409243730852e9
         -1.5991852832825754e9
          2.2671581091926174e9
          3.1924092457308525e8
         -1.27994435890949e9
         -4.791594527013428e9]

    iter = maxiter = 100.

d)
Não. O que chegou mais perto de convergir foi o método
de Gauss Seidel, como pode se observar os resíduos:
- Método de Jacobi: 2.712e9
- Método de Gauss Seidel: 0.0471
- Método SOR: 6.703e10

  =#
