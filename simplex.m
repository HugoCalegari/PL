%% Implementation of the revised Simplex. Solves a linear
% programming problem of the form
%
%   min c'*x
%   s.t. Ax  = b
%         x >= 0
%
% The function input parameters are the following:
%     A: The constraint matrix 
%     b: The rhs vector 
%     c: The vector of cost coefficients 
%     C: The indices of the basic variables corresponding to an
%        initial basic feasible solution
% 
% The function returns:
%     x_opt: Decision variable values at the optimal solution  
%     f_opt: Objective function value at the optimal solution
%
% Usage: [x_opt, f_opt] = S12345X(A,b,c,C)
%  NOTE: Replace 12345X with your own student number 
%        and rename the file accordingly

function [x_opt, f_opt] = teste(A,b,c,C)
    %% Initialization phase 
    % Initialize the vector of decision variables
    x = zeros(length(c),1);

    % Create the initial Basis matrix, compute its inverse and    
    % compute the inital basic feasible solution
    B=A(:,C);
    invB = inv(B);
    x(C) = invB*b;


    %% Iteration phase
    n_max = 10;        % At most n_max iterations
    for n = 1:n_max    % Main loop

        % Compute the vector of reduced costs c_r 

        c_B = c(C);         % Basic variable costs
        p = (c_B'*invB)';   % Dual variables
        c_r = c' - p'*A;    % Vector of reduced costs

        % Check if the solution is optimal. If optimal, use 
        % 'return' to break from the function, e.g.

        J = find(c_r < 0); % Find indices with negative reduced costs

        if (isempty(J))
            f_opt = c'*x;
            x_opt = x;
            return;
        end

        % Choose the entering variable
        j_in = J(1);

        % Compute the vector u (i.e., the reverse of the basic directions) 

        u = invB*A(:,j_in);
        I = find(u > 0);

        if (isempty(I))
           f_opt = -inf;  % Optimal objective function cost = -inf
           x_opt = [];        % Produce empty vector []
           return         % Break from the function
        end

        % Compute the optimal step length theta

        theta = min(x(C(I))./u(I));

        L = find(x(C)./u == theta); % Find all indices with ratio theta

        % Select the exiting variable

        l_out = L(1);

        % Move to the adjacent solution 

        x(C) = x(C) - theta*u;
        % Value of the entering variable is theta
        x(j_in) = theta;


        % Update the set of basic indices C 

        C(l_out) = j_in;

        % Compute the new inverse basis B^-1 by performing elementary row 
        % operations on [B^-1 u] (pivot row index is l_out). The vector u is trans-
        % formed into a unit vector with u(l_out) = 1 and u(i) = 0 for
        % other i.

        M=horzcat(u, invB);
        [f g]=size(M);
        if (theta~=0)
        M(l_out,:)=M(l_out,:)/M(l_out,1); % Copy row l_out, normalizing M(l_out,1) to 1
        end
        for k = 1:f % For all matrix rows
           if (k ~= l_out) % Other then l_out
           M(k,:)=M(k,:)-M(k,1)*M(l_out,:); % Set them equal to the original matrix Minus a multiple of normalized row l_out, making R(k,j_in)=0
        end
        end
        invB=M(1:3,2:end);


        % Check if too many iterations are performed (increase n_max to 
        % allow more iterations)
        if(n == n_max)
            fprintf('Max number of iterations performed!\n\n');
            return
        end
    end  % End for (the main iteration loop)
end  % End function

%% Example 3.5 from the book (A small test problem)
% Data in standard form:
 A = [1 2 2 1 0 0;
     2 1 2 0 1 0;
     2 2 1 0 0 1];
 b = [20 20 20]';
 c = [-10 -12 -12 0 0 0]';
% C = [4 5 6];           % Indices of the basic variables of 
%                        % the initial basic feasible solution
% 
% The optimal solution
% x_opt = [4 4 4 0 0 0]' % Optimal decision variable values
% f_opt = -136           % Optimal objective function cost


Exemplos tirados das listas de exercício
A=[-1 1 1 0; 2 -1 0 1];
b = [2 6]';
c = [-1 -1 0 0]';
m = 2; n = 4;

x_ot = [8 10 0 0]'
f_ot = -18



A = [1 1 1 0 0; 1 0 0 1 0; 0 1 0 0 1];
b = [4 3 7/2]';
c = [-2 -1 0 0 0]';
m = 3; n = 5;

x_opt = [3.0 1.0 0.0 0.0 2.5] f_opt = -7


A = [1 1 -1 0 1 0; -1 1 0 -1 0 1];
b = [2 1]';
c = [0 0 0 0 1 1]';
m = 2; n = 6;

x_ot = [0.5 1.5 0.0 0.0 0.0 0.0]' f_ot = 0


Exemplo de problema ilimitado
A = [-1 -1 1 0; -3 -5 0 1];
b = [8 30]';
c = [-4 -5 0 0]';
m = 2; n = 4;

Problema ilimitado!
x_ot = [](0x0) com função objetivo: f_ot = -Inf