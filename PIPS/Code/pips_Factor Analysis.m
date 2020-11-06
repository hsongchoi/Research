addpath('C:\Users\haesong\Desktop\PIPS')
%data=load('C:\Users\haesong\Desktop\PIPS\Total.csv');
data=load('C:\Users\haesong\Desktop\factor_student.csv');
%% none
%  delete = [6 7];
%  data(:,delete) = []

[Loadings1,specVar2,T,stats] = factoran(data,2,'rotate','none');
figure(1)
biplot(Loadings1, 'varlabels',num2str((1:24)'));
title('Unrotated Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');

%% Varimax rotation
[Loadings2,specVar2,T,stats] = factoran(data,2);
figure(2)
biplot(Loadings2, 'varlabels',num2str((1:24)'));
title('Varimax Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');

%% Promax rotation
[Loadings3,specVar3,rotation3] = factoran(data,2,'rotate','promax');

figure(3)
biplot(Loadings3, 'varlabels',num2str((1:24)'));
title('Promax Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');

%% 'equamax'

[Loadings4,specVar4,rotation4] = factoran(data,2,'rotate','equamax');

figure(4)
biplot(Loadings4, 'varlabels',num2str((1:24)'));
title('Equamax Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');


%% 'orthomax'
[Loadings5,specVar5,rotation5] = factoran(data,2,'rotate','orthomax');

figure(5)
biplot(Loadings5, 'varlabels',num2str((1:24)'));
title('Orthomax Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');

%% 'quartimax'
[Loadings6,specVar6,rotation6] = factoran(data,2,'rotate','quartimax');

figure(6)
biplot(Loadings6, 'varlabels',num2str((1:24)'));
title('quartimax Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');