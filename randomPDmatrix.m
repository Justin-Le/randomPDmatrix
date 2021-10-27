%% Example

clc; clearvars;

seed = 0;
n = 10;
condition_num = 1e3;
rand_eigvals_flag = 1;
P = generate_randomPDmatrix(seed, n, condition_num, rand_eigvals_flag);

% b = ones(n, 1);
% [sol, fval] = quadprog(P, b); % solve a matrix equation involving P

function P = generate_randomPDmatrix(seed, n, condition_num, rand_eigvals_flag)
    % Generates a random positive-definite matrix with a user-specified
    % condition number and with the smallest eigenvalue being 1.
    
    % Inputs: the RNG seed; the size of the matrix; the condition number; 
    % a flag indicating whether the eigenvalues are 
    % to be random or uniformly spaced.
    
    % Outputs: the positive-definite matrix.

    rng(seed);
    P = rand(n);
    [U, S, V] = svd(P);

    if rand_eigvals_flag == 1
        S(S~=0) = [1; (sqrt(condition_num) - 1)*rand(n-2, 1) + 1; sqrt(condition_num)];
    else
        S(S~=0) = linspace(sqrt(condition_num), 1, n);
    end

    P = U*S*V';
    P = P*P';

    fprintf('Condition number: %f.\n\n', cond(P));
    disp('Eigenvalues: '); fprintf('%.3f ', eig(P)); fprintf('\n');
end
