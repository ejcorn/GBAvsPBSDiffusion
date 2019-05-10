addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/BCT/'))
%% linear systems practice
% goal is to speed up damped oscillation of one node in a linear system
% independently

N = 500;         % N: number of nodes
t = 1:100;      % t: time
pct_skew = -1;  % pct_skew: amount of skewness to add

%%
i = randperm(N,1);      % node whose time course you want to alter
%Cij = makeevenCIJ(N,64,3);
Cij = ones(N);
A = rand(N) .* Cij;
A = A - pct_skew*A';    % skew A
A(logical(eye(length(A)))) = 0;
A = A / (max(real(eig(A))));   % normalize by max eigenvalue
disp(['Max eigenvalue of A is ',num2str(max(eig(A)))]); % max eigenvalue should be 1
c = 1;
A = A;% - c*eye(length(A));
disp(['Max eigenvalue of A is ',num2str(max(eig(A)))]); % max eigenvalue should be v. close to 1 - c;
%%
Xo = zeros(N,1);
Xo(i) = 1;
Xt = zeros(N,length(t));
for j = 1:length(t)
    Xt(:,j) = expm(A*t(j))*Xo;
end

figure; subplot(1,2,1); imagesc(A);
subplot(1,2,2); plot(Xt(:,t)');

%% steady state behavior vs. largest eigenmode
[V,D] = eig(A);
[~,which_max] = max(diag(D));

corr(Xt(:,end),V(:,which_max))
