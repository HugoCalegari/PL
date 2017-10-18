function [xot, fot, h] = funcsimplex(m,n,A,b,c)
    
    %Parte de incialização
    x = zeros(length(c),1);
    
    %Como se trabalha com o problema na forma Ax >= b, o C é definifo como os índices da particição básica inicial. 
    %É a identidade. (Observe que C é o vetor com os índices das colunas nas quais a matriz A é a identidade.)
    C = ((n-m+1):n);

    %Parte da iteração
    max_it = 20;        %Realiza-se no máximo no iterações; caso seja necessário muda-se para outro valor plausível
    for h = 1:max_it    %Loop principal

        % Cria-se a partição básica inicial e calcula a solução básica inicial pela sua inversa.
        B=A(:,C); %Partição básica inicial
        invB = inv(B); %Utilização do comando "inv()" para se calcular a matriz inversa de B
        x(C) = invB*b; %Cálculo da solução básica inicial

        %Cálculo dos custos reduzidos

        c_B = c(C);         %Custos básicos associados à partição básica
        p = (c_B'*invB)';   %Cálculo do vetor multiplicador simplex
        c_r = c' - p'*A;    %Vetor de custos relativos

        %Verificação se a solução é ótima, por meio da identificação dos custos relativos menores que zero; caso seja ótima
        %simplesmente usar a função "return"
        
        H = find(c_r < 0); % Identificação dos índices que possuem o custo relativo menor que zero

        if (isempty(H)) %Verifica se o vetor J definido como aquele que possui os custos relativos menores que zero
                        %é vazio ou não. Caso seja vazio, estamos na solução ótima
            fot = c'*x; %Cálculo da função objetivo na solução ótima
            xot = x; %Valor da solução ótima
            return;
        end

        %Variável que entrará na base por possuir o custo menor que zero      
        j_suporte = (find(min(c_r(H)) == c_r)); %Aqui determina-se o custo mínimo; nesta igualdade identifica-se o índice
                                         %cujo custo é o menor possível
        
        j_entra = j_suporte(1); %Caso houver custos relativos mínimos iguais, pegar a primeira posição com custo relativo
                             %mínimo

        %Cálculo da direção simplex

        y = invB*A(:,j_entra);
        I = find(y > 0); %Determina os índices que possuem as direções positivas para posteriormente calcular o tamanho
                         %do passo        

        if (isempty(I)) %Verifica se o vetor com o conjunto de índices é vazio ou não. Caso for vazio a solução é ilimitada
          
           fot = -inf;  %A função objetivo tende ao menos infinito
           xot = [];        %Vetor vazio; como o problema é ilimitado, não existe um valor fixo de "x" que minimize
                              %a função objetivo. 
           fprintf('Problema ilimitado!\n'); %Mensagem de PL ilimitado
           return         
        end

        %Determinação do tamanho do passo
        epsilon = min(x(C(I))./y(I));

        L = find(x(C)./y == epsilon); %Encontra o índice da solução básica que é igual a epsilon, ou seja, a variável que 
                                      %vai sair da base
        
        %Determinação da variável que sai da base

        l_sai = L(1);

        %Estratégia simplex; alteração de como é a variável básica (SlideAula7 - número 26)
        x(C) = x(C) - epsilon*y;
        
        %O valor da variável que entra é epsilon
        x(j_entra) = epsilon;

        %Atualização dos índices da matriz da partição básica
        C(l_sai) = j_entra;

        %Verifica se o n´umero de iterações máximo foi atingido. Caso tenha alterar o valor de n_max
        if(h == max_it)
            fprintf('Número máximo de iterações executado\n\n');
            return
        end
    end  %Fim do loop principal
end  %Fim da função
