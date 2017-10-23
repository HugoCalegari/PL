function [xot, fot, h] = funcsimplex(m,n,A,b,c)
    %A saída da função funcsimplex é (xot, fot) que é a representação da solução ótima; e existe o argumento h que
    %representa o número de iterações necessários para a obtenção da solução ótima ou identificação da solução ilimitada 
    
    %Parte de inicialização
    x = zeros(length(c),1); %A solução x é iniciada como um vetor coluna de zeros
    
    %Como se trabalha com o problema na forma Ax <= b, o C é definido como os índices da particição básica inicial. 
    %É a identidade. (Observe que C é o vetor com os índices das colunas nas quais a matriz A é a identidade.)
    index = ((n-m+1):n);

    %Parte da iteração
    maxit = 50;        %Realiza-se no máximo 50 iterações; caso seja necessário muda-se para outro valor plausível
    for h = 1:maxit    %Loop principal

        % Cria-se a partição básica inicial e calcula-se a solução básica inicial com o uso da inversa
        B=A(:,index); %Partição básica inicial
        invB = inv(B); %Utilização do comando "inv()" para se calcular a matriz inversa de B
        x(index) = invB*b; %Cálculo da solução básica inicial

        %Cálculo dos custos relativos

        cB = c(index);         %Custos básicos associados à partição básica
        lambda = (cB'*invB)';   %Cálculo do vetor multiplicador simplex
        cr = c' - lambda'*A;    %Vetor de custos relativos
        
        H = find(cr < 0); % Identificação dos índices que possuem o custo relativo menor que zero
		
		    %Verificação se a solução é ótima, por meio da identificação dos índices cujos custos relativos são maiores que zero; caso sejam
        %a função condicional a seguir retornará verdadeiro e, assim, simplesmente usar a função "return". Já estamos
        %na solução ótima

        if (isempty(H)) %Verifica se o vetor H definido como aquele que possui os índices dos custos relativos menores que zero
                        %é vazio ou não. Caso seja vazio, estamos na solução ótima
            fot = c'*x; %Cálculo da função objetivo na solução ótima
            xot = x; %Valor da solução ótima
            return;
        end

        %Variável que entrará na base por possuir o custo menor que zero      
        jsuporte = (find(min(cr(H)) == cr)); %Aqui determina-se o índice no qual o custo relativo é mínimo; regra de Dantzig
        
        jentra = jsuporte(1); %Caso houver custos relativos mínimos iguais, simplesmente utilizar a primeira posição
                              %da variável jsuporte, definida na variável jentra

        y = invB*A(:,jentra); %Cálculo da direção simplex
        I = find(y > 0); %Determina os índices que possuem as direções positivas para posteriormente calcular o tamanho
                         %do passo        

        if (isempty(I)) %Verifica se o vetor com o conjunto de índices da direção simplex é vazio ou não. 
                        %Caso for vazio a solução é ilimitada
          
           fot = -inf;  %A função objetivo tende ao menos infinito
           xot = [];    %Vetor vazio; como o problema é ilimitado, não existe um valor fixo de "x" que minimize
                        %a função objetivo. 
           fprintf('Problema ilimitado!\n'); %Mensagem de PL ilimitado
           return         
        end

        epsilon = min(x(index(I))./y(I)); %Determinação do tamanho do passo

        L = find(x(index)./y == epsilon); %Encontra o índice da solução básica que é igual a epsilon, ou seja, a variável que 
                                      %vai sair da base

        lsai = L(1);         %Determinação da variável que sai da base

        x(index) = x(index) - epsilon*y; %Estratégia simplex; alteração de como é a variável básica (SlideAula7 - número 26)
        
        x(jentra) = epsilon; %O valor da variável que entra é epsilon

        index(lsai) = jentra; %Atualização dos índices da matriz da partição básica

        %Verifica se o número de iterações máximo foi atingido. Caso alcançado alterar o valor de maxit
        if(h == maxit)
            fprintf('Número máximo de iterações executado\n\n');
            return
        end
    end
end
