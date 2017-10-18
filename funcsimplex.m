function [xot, fot, h] = funcsimplex(m,n,A,b,c)
    
    %Parte de inicialização
    x = zeros(length(c),1); %A solução x é iniciada como um vetor coluna de zeros
    
    %Como se trabalha com o problema na forma Ax <= b, o C é definido como os índices da particição básica inicial. 
    %É a identidade. (Observe que C é o vetor com os índices das colunas nas quais a matriz A é a identidade.)
    C = ((n-m+1):n);

    %Parte da iteração
    max_it = 20;        %Realiza-se no máximo 20 iterações; caso seja necessário muda-se para outro valor plausível
    for h = 1:max_it    %Loop principal

        % Cria-se a partição básica inicial e calcula a solução básica inicial pela sua inversa.
        B=A(:,C); %Partição básica inicial
        invB = inv(B); %Utilização do comando "inv()" para se calcular a matriz inversa de B
        x(C) = invB*b; %Cálculo da solução básica inicial

        %Cálculo dos custos relativos

        c_B = c(C);         %Custos básicos associados à partição básica
        lambda = (c_B'*invB)';   %Cálculo do vetor multiplicador simplex
        c_r = c' - lambda'*A;    %Vetor de custos relativos
        
        H = find(c_r < 0); % Identificação dos índices que possuem o custo relativo menor que zero
		
		%Verificação se a solução é ótima, por meio da identificação dos custos relativos maiores que zero; caso sejam
        %a função condicional a seguir retornará verdadeiro e, assim, simplesmente usar a função "return". Já estamos
        %na solução ótima

        if (isempty(H)) %Verifica se o vetor H definido como aquele que possui os custos relativos menores que zero
                        %é vazio ou não. Caso seja vazio, estamos na solução ótima
            fot = c'*x; %Cálculo da função objetivo na solução ótima
            xot = x; %Valor da solução ótima
            return;
        end

        %Variável que entrará na base por possuir o custo menor que zero      
        j_suporte = (find(min(c_r(H)) == c_r)); %Aqui determina-se o custo relativo mínimo; nesta igualdade identifica-se o índice
	                                            %cujo custo é o menor possível; regra de Dantzig
        
        j_entra = j_suporte(1); %Caso houver custos relativos mínimos iguais, pegar a primeira posição com custo relativo
                                %mínimo

        y = invB*A(:,j_entra); %Cálculo da direção simplex
        I = find(y > 0); %Determina os índices que possuem as direções positivas para posteriormente calcular o tamanho
                         %do passo        

        if (isempty(I)) %Verifica se o vetor com o conjunto de índices é vazio ou não. Caso for vazio a solução é ilimitada
          
           fot = -inf;  %A função objetivo tende ao menos infinito
           xot = [];    %Vetor vazio; como o problema é ilimitado, não existe um valor fixo de "x" que minimize
                        %a função objetivo. 
           fprintf('Problema ilimitado!\n'); %Mensagem de PL ilimitado
           return         
        end

        epsilon = min(x(C(I))./y(I)); %Determinação do tamanho do passo

        L = find(x(C)./y == epsilon); %Encontra o índice da solução básica que é igual a epsilon, ou seja, a variável que 
                                      %vai sair da base

        l_sai = L(1);         %Determinação da variável que sai da base

        x(C) = x(C) - epsilon*y; %Estratégia simplex; alteração de como é a variável básica (SlideAula7 - número 26)
        
        x(j_entra) = epsilon; %O valor da variável que entra é epsilon (perturbação)

        C(l_sai) = j_entra; %Atualização dos índices da matriz da partição básica

        %Verifica se o número de iterações máximo foi atingido. Caso alcançado alterar o valor de max_it
        if(h == max_it)
            fprintf('Número máximo de iterações executado\n\n');
            return
        end
    end  %Fim do loop principal
end  %Fim da função
