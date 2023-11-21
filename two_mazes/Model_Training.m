addpath(genpath('/home/kemerelab/Downloads/PGPLVM/LMT-master'));warning off

%% Train model

% Define training data from neural activity during running periods in two mazes

load('maze1_run_v1_3_v2_5_500ms.mat');
figure;plot(pos_linear)
s1 = 401;
e1 = 800;
portion = s1:e1;
xx1 = pos_linear(portion)';
yy1 = double(spikes(:,portion)');
tgrid1 = double(idx(portion)');

load('maze2_run_v1_3_v2_5_500ms.mat');
figure;plot(pos_linear)
s2 = 401;
e2 = 800;
portion = s2:e2;
xx2 = pos_linear(portion)';
plot(xx2)
yy2 = double(spikes(:,portion)');
tgrid2 = double(idx(portion)');

xx=[xx1;xx2];
yy=[yy1;yy2];
tgrid=[tgrid1;tgrid2+tgrid1(end)];

name = ['m1_',num2str(s1),'_',num2str(e1),'_m2_',num2str(s2),'_',num2str(e2),'_result_la.mat']
clearvars -except xx yy tgrid name


% Train the model

yy = yy(:,stable_id); % select stable cell ids if needed
nf = 3; % number of latent space dimensions
niter = 30; % number of iterations

figure; % show latent variables in every iteration
[result_la,setopt]=run_pgplvm(xx,yy,tgrid,nf,niter,[]);
figure
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'line_only',0)


% Discover number of separate manifolds

[tbl_id, comp_id,k_log,n_debris_log] = find_components(result_la.xxsamp, 0, 1);
x1=result_la.xxsamp(comp_id==1,:);x2=result_la.xxsamp(comp_id==2,:);
std_all=sqrt(mean([vecnorm(x1-mean(x1),2,2).^2;vecnorm(x2-mean(x2),2,2).^2]));

save(name,'result_la','xx','setopt','std_all')