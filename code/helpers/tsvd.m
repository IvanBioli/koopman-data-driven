function [U,Sigma,V] = tsvd(A,tol,k)
%   [U,SIGMA,V] = TSVD(A,TOL,K) calcola la Truncated Singular Value
%   Decomposition della matrice A, eliminando i vettori/valori singolari 
%   relativi ai i valori singolari sigma_r tali che r > k o 
%   sigma_r/sigma_1 < TOL.
%   INPUT:
%       - A: matrice di cui calcolare la TSVD
%       - TOL: tolleranza relativa per il troncamento dei valori singolari
%       - k: numero massimo di valori singolari mantenuti
%   OUTPUT:
%       - U,Sigma,V: matrici tali che Atilde = U*Sigma*V' fornisce
%         l'approssimazione di A cercata

%Gestione nel caso in cui k<=0 o k>min(m,n)
[m,n] = size(A);
if (k <= 0 || k > min(m,n))
    disp('Warning: k <= 0 o k > min(m,n), impostato di default k = min(m,n)');
    k = min(m,n);
end

%Determinazione dell'indice per il troncamento 
[U,Sigma,V] = svd(A,'econ');
sigma1 = Sigma(1,1);
r = find((diag(Sigma)/sigma1 < tol),1);
if ~isempty(r)
    k = min(k,r);
end

%Esecuzione del troncamento
U = U(:,1:k);
Sigma = Sigma(1:k,1:k);
V = V(:,1:k);
end

