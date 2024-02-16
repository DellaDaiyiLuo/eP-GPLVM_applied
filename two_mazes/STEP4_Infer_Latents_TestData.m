addpath(genpath('/home/kemerelab/Downloads/PGPLVM/LMT-master'));warning off

%% Infer latents of test data - run

% Define training data from neural activity during running

load('/home/kemerelab/Downloads/new_mazes_eeg/maze1_run_v1_3_v2_5_500ms.mat');
figure;plot(pos_linear)
s1 = 508; e1 = 667; %maze1 test data
% s2=719; e2=982; %maze2 test data
portion = s1:e1;
xx0 = pos_linear(portion)';
yy0 = double(spikes(:,portion)');
tgrid = double(idx(portion)');
d = [0,find(diff(tgrid')>1),numel(tgrid)];


% Load trained P-GPLVM model

load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat')
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)}; % for visualization


% Infer latent variables of test data

len = std_all/5;
k=6;
niter=20;
figure
[result_ori, setopt, fftc, ~, logpy_ori] = run_projection(xx0,yy0,tgrid,result_la,niter,d);


% Estimate metrics

[curv_ori,stepdis_ori] = estimate_smoothness(result_ori.xxsamp,d,1);
[knndis_ori, ~,~] = projected_xx_measures_segs(result_ori.xxsamp,{result_la.xxsamp},d,0,k);
consist_seg_ori = estimate_consistency(result_ori.xxsamp,xbackground,d,k,len);
consist_seg1_ori = consist_seg(1,:);
consist_seg2_ori = consist_seg(2,:);


% Visualize - gray xbackground & colored test data

figure
show_latent_variable(result_ori.xxsamp,xx0,xbackground,setopt.tgrid,'scatter_only',0)
save('ori_.mat','xx0','yy0','result_ori', 'setopt', 'fftc','curv_ori','logpy_ori','consist_seg1_ori','consist_seg2_ori','consist_seg_ori','d')

% Visualize - colored xbackground & black test data

figure
scatter3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),50,xx,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',0.2)
hold on
show_latent_variable(result_ori.xxsamp,xx0,[],setopt.tgrid,'scatter_only',1)
show_latent_variable(result_ori.xxsamp,xx0,[],setopt.tgrid,'line_only',1,'line_color',[.7,.7,.7])


% shuffle

consist_seg_all = [];
curv_all=[];
stepdis_all=[];
results=[];
logpy_all = {[],[],[]};
knndis_all = [];

type = {'cell','segtime','cir','ind_time'};
filename = {'cellshuf_pool.mat','segtimeshuf_pool.mat','localcirshuf_pool.mat','indtimeshuf_pool.mat'};

shuftype=3;
load(filename{shuftype})
load([type{shuftype} 'results.mat'])
tic
for j=1:24
    for i=1:19
        display(['batch' num2str(j) ', shuf iter ' num2str(i)])
        [result, setopt, fftc, order, logpy] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
            'shuffle',type{shuftype},'draw',0,'TCgrid',0);
        results=[results,result.xxsamp];
        [curv,stepdis] = estimate_smoothness(result.xxsamp,d,1);
        curv_all=[curv_all;curv];
        stepdis_all=[stepdis_all;stepdis];
        [knndis_seg, ~,~] = projected_xx_measures_segs(result.xxsamp,xbackground,d,0,k);
        knndis_all=[knndis_all;knndis_seg];

        consist_seg_all = [consist_seg_all;estimate_consistency(result.xxsamp,xbackground,d,k,len)]; % return same row number as cell number in xbackground
        for t=1:3
            logpy_all{t} = [logpy_all{t};logpy(t,:)];
        end
    end
    save(filename{shuftype},'logpy_all','curv_all','stepdis_all','consist_seg_all','knndis_all');
    save([type{shuftype} 'results.mat'],'results')
end
toc
figure
show_latent_variable(result.xxsamp,xx,xbackground,setopt.tgrid,'scatter_only',0)


% ---- additional metrics ---- %

clear
load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat')
len = std_all/5;
k=6;
load('ori_.mat')
d = [0,find(round(diff(setopt.tgrid'))>1),numel(setopt.tgrid)];
load('cellresults.mat')
consist_seg_all = [];
for i=1:500
    consist_seg = estimate_consistency(results(:,i*3-2:i*3),xbackground,d,k,len);
    consist_seg_all = [consist_seg_all;consist_seg];
end
consist_seg_ori = estimate_consistency(result_ori.xxsamp,xbackground,d,k,len);
clearvars -except consist_seg_all


%% Infer latents of test data - PBE

% Select 100 PBEs in maze1 exploration as test data

tb = 17;
load(['/home/kemerelab/Downloads/new_mazes_eeg/pbe_maze1_' num2str(tb) 'ms.mat'])
scaler=tb/500;

% --- select randomly --- %
id_event = randsample(size(event_edge,1),100);
id_event = sort(id_event);
portion = [];
for i=1:100
    portion = [portion,event_edge(id_event(i),1)+1:event_edge(id_event(i),2)+1];
end
% % --- select certain range --- %
% n_event_range = [400,450]; 
% portion = event_edge(n_event_range(1),1)+1:event_edge(n_event_range(2),2)+1;

yy0 = double(spikes(:,portion))';
xx0 = zeros(size(yy0,1),1);
tgrid = (time(portion)-time(portion(1)))/tb*1000;
tgrid = tgrid';
event_length=[diff(event_edge(:,1));diff(event_edge(end,:))+1];
d = double([0;cumsum(event_length(id_event))]);

clearvars -except xx0 yy0 tgrid scaler id_event d
save('select_pbe_maze2_17ms.mat')



% Infer latent variables of PBEs

load('/home/kemerelab/Downloads/new_mazes_eeg/review/select_100_pbe_maze1n2_17ms.mat')
load('/home/kemerelab/Downloads/new_mazes_eeg/run_run_1maze/m1_101_500_result_la.mat')
tgrid_tc = setopt.tgrid;
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};
len = std_all/5;
k=6;


% project all PBEs at once

niter=20;
figure
[result_ori, setopt, fftc, ~,logpy_ori] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
    'newdatatype','pbe','scaler',scaler);


% % project one PBE at a time

% xxsamp = [];
% ffmat = [];
% logpy_ori = [];
% niter=10;
% for i=1:100
%     por = d(i)+1:d(i+1);
%     [result_ori, setopt, fftc, ~,logpy] = run_projection(xx0(por),yy0(por,:),tgrid(por),result_la,niter,d(i:i+1)-d(i),'xbackground',xbackground,...
%         'tgrid_tc',tgrid_tc,'newdatatype','pbe','scaler',scaler,'draw',0);
%     xxsamp = [xxsamp;result_ori.xxsamp];
%     ffmat = [ffmat;result_ori.ffmat];
%     logpy_ori = [logpy_ori,logpy];
% end
% result_ori.xxsamp = xxsamp;
% result_ori.ffmat = ffmat;

[curv_ori,stepdis_ori] = estimate_smoothness(result_ori.xxsamp,d,1);
[knndis_ori, ~,~] = projected_xx_measures_segs(result_ori.xxsamp,{result_la.xxsamp},d,0,k);
consist_seg_ori = estimate_consistency(result_ori.xxsamp,xbackground,double(d),k,len);


save('/home/kemerelab/Downloads/new_mazes_eeg/review/ori_pbe.mat','result_ori','curv_ori','stepdis_ori','logpy_ori','consist_seg_ori','knndis_ori')

figure
show_latent_variable(result_ori.xxsamp,xx0,xbackground,setopt.tgrid,'line_only',1)
save('ori_17ms.mat','result_ori', 'setopt', 'fftc','s_ori','logpy_ori','consist_seg1_ori','consist_seg2_ori','consist_seg_ori')


%% shuffle PBE

load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat','result_la','std_all','setopt')
% load('select_pbe_maze2.mat')
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};
% tgrid_tc = setopt.tgrid;
tb = 17;
scaler=tb/500;

type = {'cell','segtime','cir','ind_time'};
filename = {'cellshuf_pool_runs_to_maze1.mat','segtimeshuf_pool_runs_to_maze2.mat.mat','localcirshuf_pool','indtimeshuf_pool.mat'};
shuftype=2;
len = std_all/5;
consist_seg1_all = [];
consist_seg2_all = [];
curv_all=[];
results=[];
logpy_all = {[],[],[]};
% consist_seg_all = [];
stepdis_all = [];
knndis_all = [];

figure
for j=1:5
    for i=1:10
        [result, setopt, fftc, order, logpy] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
            'newdatatype','pbe','scaler',scaler,'shuffle',type{shuftype},'draw',0);
%         results=[results,result.xxsamp];
        [curv,stepdis] = estimate_smoothness(result.xxsamp,d,1);
        curv_all=[curv_all;curv];
        stepdis_all=[stepdis_all;stepdis];
        [knndis, ~,~] = projected_xx_measures_segs(result.xxsamp,{result_la.xxsamp},d,0,k);
        knndis_all = [knndis_all;knndis];
        consist_seg = estimate_consistency(result.xxsamp,xbackground,d,k,len);
        consist_seg1_all = [consist_seg1_all;consist_seg(1,:)];
%         consist_seg2_all = [consist_seg2_all;consist_seg(2,:)];
%         consist_seg_all = [consist_seg_all;estimate_consistency(result.xxsamp,{result_la.xxsamp},d,k,len)];
        for t=1:3
            logpy_all{t} = [logpy_all{t};logpy(t,:)];
        end
    end
    save(filename{shuftype},'consist_seg1_all','logpy_all','knndis_all','stepdis_all','curv_all');
%     save([type{shuftype} 'results.mat'],'results')
end


