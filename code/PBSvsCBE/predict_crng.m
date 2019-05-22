basedir = '~/Dropbox/Neurodegeneration/GBAvsPBSDiffusion/';
outputdir = 'GBAvsPBSasyndiffusionCMax60';
savedir = fullfile(basedir,outputdir,'PBSvsCBE');

% load in structural connectivity
load(fullfile(basedir,outputdir,'processed','W.mat'));

t = 1; % time point to testp
n_c = 2000; % number of time constants to test
c_rng = linspace(0,40,n_c); % set range of time constants
n_regions = size(W,1);
iCP = 17; % index of injection site
Xo = zeros(n_regions,1); Xo(iCP) = 1; % make initial conditions

% Where i is row element and j is column element
% Wij is a connection from region i to region j
% convention is the opposite, so without transposing W
% I am capturing "retrograde" connectivity

W = W .* ~eye(n_regions); % get rid of diagonal
in_deg = sum(W,1);
out_deg = sum(W,2);
L_out = diag(out_deg) - W;% outdegree laplacian

Xt_c = zeros(n_regions,n_c);
for i = 1:n_c
	%A = L_out*c_rng(i)*t;
	%[V,D] = eig(A);
	%A_expm = V*diag(exp(diag(D)))/V;
	Xt_c(:,i) = expm(-L_out*c_rng(i)*t)*Xo;
end

save(fullfile(savedir,'RegionalTrajectoriesCRange.mat'),'c_rng','Xt_c');

thrsh = [0 20 50 60 70 80];
for cutoff = thrsh
	f=figure; plot(c_rng(c_rng > cutoff),Xt_c(:,c_rng>cutoff)');
	xlabel('c'); ylabel('X(t)');
	f.Units = 'centimeters';
	f.PaperSize = [8 8];
	f.PaperPosition = [0 0 8 8];
	saveas(f,fullfile(savedir,['RegionalTrajectories',num2str(cutoff),'-110.pdf']))
end

%%
% similarity in predicted values... regional distribution doesn't change that much over time, despite changes in absolute value of pathology
nticks = 5; tick_space = round(linspace(1,n_c,nticks),-2); tick_space(1) = 1;
f=figure; 
subplot(1,2,1);
imagesc(corr(Xt_c)); caxis([-1 1]); colorbar; axis square
title('Correlation between points in X(t)')
xlabel('c'); ylabel('c'); 
xticks(tick_space); xticklabels(round(c_rng(tick_space)));
yticks(tick_space); yticklabels(round(c_rng(tick_space)));
subplot(1,2,2);
imagesc(squareform(pdist(Xt_c'))); colorbar; axis square
xlabel('c'); ylabel('c');
xticks(tick_space); xticklabels(round(c_rng(tick_space)));
yticks(tick_space); yticklabels(round(c_rng(tick_space)));
title('Euclidean distance between points in X(t)');
f.Units = 'centimeters';
f.PaperSize = [18 9];
f.PaperPosition = [0 0 18 9];
saveas(f,fullfile(savedir,['XtCorrelationDistance.pdf']))
