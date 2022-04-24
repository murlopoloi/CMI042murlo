function tspawn(t=1e-10)
    for i = 2:11
        t = [t 10.0^(-11+i)]
    end
    return t
end

function aspawn(t)
    A = zeros(11,6)
    for i = 1:11
        for j = 1:6
            A[i,j] = t[i]^(j-1)
        end
    end
    return A
end 

b = A*ones(6)

btil = b + ((0.5 .- rand(11,11)).*1.0e-10)

function delbspawn(b)
    tmp = zeros(11,11)
    delb = (0.5 .- rand(11,11))*b
    for j = 1:11
        tmp[:,j] = 10.0^(-j).*delb
    end
    return tmp
end

function letrae(A, b, delb)
    sols = zeros(6,11)
    for j = 1:11
        tmp = pinv(A)*(b + delb[:,j])
        sols[:,j] = tmp
    end
    return sols
end

function errspawn(sol,delb)
    erra, errb = zeros(11), zeros(11)
    for j = 1:11
        erra[j] = norm(sol[:,j] - ones(6), 2)
        errb[j] = norm(delb[:,j], 2)
    end
    return erra, errb
end

function pwrmtd(A, palp, maxiter=10000, tol1=1e-3, tol2=1e-4)
    iter = 1
    y = A * palp
    x = [palp y / norm(y, 2)]
    mu = palp' * A * palp
    while iter < maxiter
        y = [y A * x[:, iter]]
        x = [x y[:, iter] / norm(y[:, iter], 2)]
        mu = [mu x[:, iter]' * A * x[:, iter]]
        iter += 1
        if norm(A * x[:, iter] - mu[iter] * x[:, iter], 2) < tol2
            return mu[iter], iter
            break
        end
    end
    return mu[iter], iter
end


#=== Projeto 2 ===#
#=

#=== Exercício 2 ===#
a) A matriz do sistema é:

A =  [1.0  1.0e-10  1.0e-20  1.0e-30  1.0e-40  1.0e-50;
      1.0  1.0e-9   1.0e-18  1.0e-27  1.0e-36  1.0e-45;
      1.0  1.0e-8   1.0e-16  1.0e-24  1.0e-32  1.0e-40;
      1.0  1.0e-7   1.0e-14  1.0e-21  1.0e-28  1.0e-35;
      1.0  1.0e-6   1.0e-12  1.0e-18  1.0e-24  1.0e-30;
      1.0  1.0e-5   1.0e-10  1.0e-15  1.0e-20  1.0e-25;
      1.0  0.0001   1.0e-8   1.0e-12  1.0e-16  1.0e-20;
      1.0  0.001    1.0e-6   1.0e-9   1.0e-12  1.0e-15;
      1.0  0.01     0.0001   1.0e-6   1.0e-8   1.0e-10;
      1.0  0.1      0.01     0.001    0.0001   1.0e-5;
      1.0  1.0      1.0      1.0      1.0      1.0]

Note que A é uma matriz 11x6 e portanto
não admite inversa, portanto usamos sua 
pseudo-inversa, A⁺.
cond(A) = ||A||*||A⁺|| = 7.3e10;
Novamente usamos a pseudo-inversa de AᵀA pois a matriz
original tem det(AᵀA) = 0.
cond(AᵀA) = ||AᵀA||*||(AᵀA)⁺|| = 2.83e13

b) Para resolver o problema no sentido dos QM
com equações normais, basta resolver o sistema
(AᵀA)a̅ = Aᵀb̃  
para o b̃ = [1.000000000319151
            1.0000000010940948
            1.0000000100438697
            1.0000001000518761
            1.000001000272732
            1.0000100000700827
            1.0001000096769226
            1.0010010006549706
            1.0101010100019265
            1.1111100001486796
            5.999999999962692]
    obtido.
    Não conseguimos resolver o sistema
    via equações normais devido a sin- 
    gularidade de AᵀA.

c) Fazendo svd() para obter a decomposição, para
resolver o sistema basta fazer x̂ = V*S⁻¹*Uᵀ*b̃
    ̂x = [1.0000000001192366
         0.9999949719524042
         1.0050638281037132
         0.494040471784181
         6.009054742040767
        -3.5081554863265048]

    ou ainda, resolvendo via ̂x = A⁺b̃, temos:
    ̂x = [1.0000000001192364
         0.9999949719523978
         1.0050638299677228
         0.4940403525748913
         6.009056651252046
        -3.508155486328322]

        
 d)  Para a letra c), obtemos para o 
 primeiro método de resolução:
   Erroₐ = 6.758
   e para o segundo:
   Erroₐ = 6.758

   Usando como critério o erro calculado,
considerando todas as casas decimais possíveis,
têm-se que o primeiro método SVD utilizado fornece
um erro menor para este caso.

e) Para este item será usada a resolução via 
pseudoinversa da matriz A.
            
            δb⁽¹⁾          δb⁽²⁾        δb⁽³⁾          δb⁽⁴⁾
   δb = [-0.0061288   -0.00061288   -6.1288e-5    -6.1288e-6
          0.166411     0.0166411     0.00166411    0.000166411
          0.0564344    0.00564344    0.000564344   5.64344e-5
          0.151056     0.0151056     0.00151056    0.000151056
          0.012015     0.0012015     0.00012015    1.2015e-5
          0.0887275    0.00887275    0.000887275   8.87275e-5
         -0.0580855   -0.00580855   -0.000580855  -5.80855e-5
          0.0417439    0.00417439    0.000417439   4.17439e-5
         -0.00243264  -0.000243264  -2.43264e-5   -2.43264e-6
         -0.121766    -0.0121766    -0.00121766   -0.000121766
          0.148113     0.0148113     0.00148113    0.000148113
          
             δb⁽⁵⁾        δb⁽⁶⁾         δb⁽⁷⁾        δb⁽⁸⁾
          -6.1288e-7   -6.1288e-8   -6.1288e-9   -6.1288e-10
          1.66411e-5   1.66411e-6   1.66411e-7   1.66411e-8
          5.64344e-6   5.64344e-7   5.64344e-8   5.64344e-9
          1.51056e-5   1.51056e-6   1.51056e-7   1.51056e-8
          1.2015e-6    1.2015e-7    1.2015e-8    1.2015e-9
          8.87275e-6   8.87275e-7   8.87275e-8   8.87275e-9
         -5.80855e-6  -5.80855e-7  -5.80855e-8  -5.80855e-9
          4.17439e-6   4.17439e-7   4.17439e-8   4.17439e-9
         -2.43264e-7  -2.43264e-8  -2.43264e-9  -2.43264e-10
         -1.21766e-5  -1.21766e-6  -1.21766e-7  -1.21766e-8
          1.48113e-5   1.48113e-6   1.48113e-7   1.48113e-8
          
              δb⁽⁹⁾         δb⁽¹⁰⁾        δb⁽¹¹⁾
          -6.1288e-11   -6.1288e-12   -6.1288e-13;
           1.66411e-9    1.66411e-10   1.66411e-11;
           5.64344e-10   5.64344e-11   5.64344e-12;
           1.51056e-9    1.51056e-10   1.51056e-11;
           1.2015e-10    1.2015e-11    1.2015e-12;
           8.87275e-10   8.87275e-11   8.87275e-12;
          -5.80855e-10  -5.80855e-11  -5.80855e-12;
           4.17439e-10   4.17439e-11   4.17439e-12;
          -2.43264e-11  -2.43264e-12  -2.43264e-13;
          -1.21766e-9   -1.21766e-10  -1.21766e-11;
           1.48113e-9    1.48113e-10   1.48113e-11]

    que nos forneceu as seguines soluções para cada respectivo
    δb:
    ̂x = [1.06096        1.0061        1.00061     1.00006
        -348.244       -33.9244      -2.49244     0.650756
         3.16851e5      31686.0       3169.5      317.85
        -3.11491e7     -3.11491e6    -3.1149e5   -31148.1
         2.77074e8      2.77074e7     2.77074e6   2.77075e5
         3.0818e7       3.0818e6      3.0818e5    30818.1
         ...
         1.00001       1.0         1.0        1.0
         0.965076      0.996508    0.999651   0.999965
         32.685        4.1685      1.31685    1.03169
        -3113.91      -310.492    -30.1501   -2.11589
         27708.5       2771.84     278.172    28.8051
         3081.92       308.301     30.9388    3.20263
         ...
         1.0       1.0       1.0;
         0.999997  1.0       1.0;
         1.00317   1.00032   1.00003;
         0.687523  0.967866  0.9959;
         3.86842   1.37475   1.12538;
         0.429013  0.151651  0.123915]

f) Erroₐ = [2.8051781524286914e8
            2.8051781524288066e7
            2.8051781524267476e6
            280517.8152430002
            28051.781538019328
            2805.178291137437
            280.51920778846653
            28.06572632535562
            2.9413427562552807
            0.9279903713583947
            0.8850209635039678]
        
   Erro_b = [0.32192622721071773
              0.03219262272107177
              0.003219262272107177
              0.00032192622721071774
              3.219262272107177e-5
              3.219262272107177e-6
              3.2192622721071765e-7
              3.219262272107177e-8
              3.2192622721071776e-9
              3.2192622721071766e-10
              3.219262272107177e-11]

g) Sendo os elementos da diagonal da matriz S, da decomposição SVD:
    3.4315598794145354
    2.056162079008163
    0.08243666808507952
    0.0007244300719504885
    6.640573520304083e-7
    5.4804937345545984e-11,
   nota-se que há três valores singulares próximos de zero (para uma tolerância
   de 1e-3).
   Levando-se em consideração também o número de condição da matriz, espera-se que
   as soluções não sejam boas soluções.


#=== Exercício 3 ===#
a)
=#
function pwrmtd(A, palp, maxiter = 10000, tol1 = 1e-3, tol2 = 1e-4)
    iter = 1
    y = A*palp
    x = [palp y/norm(y, 2)]
    mu = palp'*A*palp
    while iter < maxiter
        y = [y A*x[:,iter]]
        x = [x y[:,iter]/norm(y[:,iter], 2)]
        mu = [mu x[:,iter]'*A*x[:,iter]]
        iter += 1
        #if norm(A*x[:,iter] - mu[iter]*x[:,iter], 2) < tol2
        if norm(x[:,iter] - x[:,iter-1], 2) < tol2
            return mu[iter], iter
            break
        end
    end
    return mu[iter], iter
end
#=
b) A = [1 1 1; -1 9 2; 0 -1 2]
   B = [1 2 1; 2 4 -1; 0 0 6]
   Para a matriz A:
   O autovalor dominante da matriz A é 8.5844.
   O método das potências implementado encontrou o valor aproximado do
   autovalor dominante como sendo 8.4 em três iterações, para dado cri-
   tério de convergência.

   Para a matriz B:
   O autovalor dominante da matriz B é 6.
   O método das potências implementado encontrou o valor aproximado do
   autovalor dominante como sendo 5.2727 em três iterações, para dado cri-
   tério de convergência.

   Para o critério de convergência estabelecido no enunciado, o método conver-
   giu no mesmo número de iterações para ambas matrizes A e B. Porém, a matriz
   que se aproximou melhor da solução foi a matriz A.