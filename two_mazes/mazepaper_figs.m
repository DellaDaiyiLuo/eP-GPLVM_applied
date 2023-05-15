addpath(genpath('/home/kemerelab/Downloads/PGPLVM/LMT-master'));warning off

load('maze1_run_v1_3_v2_5_500ms.mat');
figure;plot(pos_linear)
s1 = 101;
e1 = 500;
portion = s1:e1;
xx1 = pos_linear(portion)';
yy1 = double(spikes(:,portion)');
tgrid1 = double(idx(portion)');

load('maze2_run_v1_3_v2_5_500ms.mat');
figure;plot(pos_linear)
s2 = 301;
e2 = 700;
portion = s2:e2;
xx2 = pos_linear(portion)';
plot(xx2)
yy2 = double(spikes(:,portion)');
tgrid2 = double(idx(portion)');

xx=[xx1;xx2];
yy=[yy1;yy2];
tgrid=[tgrid1;tgrid2+tgrid1(end)];

name = ['m1_',num2str(s1),'_',num2str(e1),'_m2_',num2str(s2),'_',num2str(e2),'_result_la.mat'];
clearvars -except xx yy tgrid name

nf = 3;
niter = 30;
figure;
%plot(1:n1,xx(1:n1));hold on;plot(n1+1:n1+n2,xx(n1+1:end));scatter(select,xx(select),'.')
[result_la,setopt]=run_pgplvm(xx,yy,tgrid,nf,niter,[]);
figure
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'line_only',0)


[tbl_id, comp_id,k_log,n_debris_log] = find_components(result_la.xxsamp, 0, 1);
x1=result_la.xxsamp(comp_id==1,:);x2=result_la.xxsamp(comp_id==2,:);
std_all=sqrt(mean([vecnorm(x1-mean(x1),2,2).^2;vecnorm(x2-mean(x2),2,2).^2]));

save(name,'result_la','xx','setopt','std_all')

%% projection
load('/home/kemerelab/Downloads/new_mazes_eeg/maze1_run_v1_3_v2_5_500ms.mat');
figure;plot(pos_linear)
s1 = 508;
e1 = 667;
portion = s1:e1;
xx0 = pos_linear(portion)';
yy0 = double(spikes(:,portion)');
tgrid = double(idx(portion)');
d = [0,find(diff(tgrid')>1),numel(tgrid)];

load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat','result_la')
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};

len = std_all/5;
k=6;
niter=30;
figure
[result_ori, setopt, fftc, ~, logpy_ori] = run_projection(xx0,yy0,tgrid,...
    result_la,niter,d,'givemeasure',0);

[curv_ori,stepdis_ori] = estimate_smoothness(result_ori.xxsamp,d,1);
[knndis_ori, ~,~] = projected_xx_measures_segs(result_ori.xxsamp,{result_la.xxsamp},d,0,k);
% logpy_ori = estimate_py(result_ori, fftc, result_la.xxsamp, yy0,tgrid,d);
consist_seg_ori = estimate_consistency(result_ori.xxsamp,{result_la.xxsamp},d,k,len);
consist_seg = estimate_consistency(result_ori.xxsamp,xbackground,d,k,len);
consist_seg1_ori = consist_seg(1,:);
consist_seg2_ori = consist_seg(2,:);

figure
show_latent_variable(result_ori.xxsamp,xx0,xbackground,setopt.tgrid,'scatter_only',0)
save('ori_.mat','xx0','yy0','result_ori', 'setopt', 'fftc','curv_ori','logpy_ori','consist_seg1_ori','consist_seg2_ori','consist_seg_ori','d')

figure
scatter3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),50,xx,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',0.2)
hold on
show_latent_variable(result_ori.xxsamp,xx0,[],setopt.tgrid,'scatter_only',1)
show_latent_variable(result_ori.xxsamp,xx0,[],setopt.tgrid,'line_only',1,'line_color',[.7,.7,.7])


% shuffle
consist_seg1_all = [];
consist_seg2_all = [];
consist_seg_all = [];
curv_all=[];
stepdis_all=[];
results=[];
logpy_all = {[],[],[]};
knndis_all = [];

type = {'cell','segtime','cir','ind_time'};
filename = {'cellshuf_pool.mat','segtimeshuf_pool.mat','localcirshuf_pool.mat','indtimeshuf_pool.mat'};
% for shuftype=1:20
shuftype=3;
load(filename{shuftype})
load([type{shuftype} 'results.mat'])
tic
for j=1:24
    for i=1:19
        display(['batch' num2str(j) ', shuf iter ' num2str(i)])
        [result, setopt, ~, ~,logpy] = run_projection(xx0,yy0,tgrid,...
            result_la,niter,d,'givemeasure',0,'TCgrid',0,'shuffle',type{shuftype},'draw',0);
        results=[results,result.xxsamp];
        [curv,stepdis] = estimate_smoothness(result.xxsamp,d,1);
        curv_all=[curv_all;curv];
        stepdis_all=[stepdis_all;stepdis];
        [knndis_seg, ~,~] = projected_xx_measures_segs(result.xxsamp,{result_la.xxsamp},d,0,k);
        knndis_all=[knndis_all;knndis_seg];
        
%         consist_seg = estimate_consistency(result.xxsamp,xbackground,d,k,len);
%         consist_seg1_all = [consist_seg1_all;consist_seg(1,:)];
%         consist_seg2_all = [consist_seg2_all;consist_seg(2,:)];
        consist_seg_all = [consist_seg_all;estimate_consistency(result.xxsamp,{result_la.xxsamp},d,k,len)];
%         switch shuftype
%             case 1
%                 logpy = estimate_py(result, fftc, result_la.xxsamp, yy0(:,order),tgrid,d);
%             case 2
%                 logpy = estimate_py(result, fftc, result_la.xxsamp, yy0(order,:),tgrid,d);
%         end
        for t=1:3
            logpy_all{t} = [logpy_all{t};logpy(t,:)];
        end
    end
    save(filename{shuftype},'logpy_all','curv_all','stepdis_all','consist_seg_all','knndis_all');
%     save(filename{shuftype},'consist_seg1_all','consist_seg2_all','logpy_all','curv_all','stepdis_all','consist_seg_all','knndis_all');
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
% consist_seg_all = [];
c_all = [];
for i=1:500
    c_all=[c_all;estimate_smoothness(results(:,i*3-2:i*3),d)];
%     consist_seg = estimate_consistency(results(:,i*3-2:i*3),{result_la.xxsamp},d,k,len);
%     consist_seg_all = [consist_seg_all;consist_seg];
end
consist_seg_ori = estimate_consistency(result_ori.xxsamp,{result_la.xxsamp},d,k,len);
c_ori = estimate_smoothness(result_ori.xxsamp,d);
% clearvars -except consist_seg_all
clearvars -except c_all


%------------------------------%

clear
clc
vars = {'consist_seg1','consist_seg2','stepdis','logpy'};
data_time = load('segtimeshuf_pool.mat');
data_cell = load('cellshuf_pool.mat');
data_cir = load('localcirshuf_pool.mat');
data_ori = load('ori_.mat');
d = [0,find(diff(data_ori.setopt.tgrid')>1),numel(data_ori.setopt.tgrid)];
seglen = diff(d);

%% projection PBE
% firing rate scaling param
tb = 17;
load(['/home/kemerelab/Downloads/new_mazes_eeg/pbe_post1_' num2str(tb) 'ms.mat'],'spikes')
totalspikes_pbe = sum(spikes,1);
clearvars -except xx0 yy0 tgrid scaler id_event d
save('select_pbe_maze1_17ms.mat')

% load('/home/kemerelab/Downloads/new_mazes_eeg/maze2_run_v1_3_v2_5_500ms.mat','spikes');
% totalspikes_run = sum(spikes,1);
% figure
% h1=histogram(totalspikes_pbe);
% hold on
% h2=histogram(totalspikes_run);
% [~,I]=max(h1.Values);
% p_pbe=h1.BinEdges(I)+h1.BinWidth/2;
% [~,I]=max(h2.Values);
% p_run = h2.BinEdges(I)+h2.BinWidth/2;
% scaler=p_pbe/p_run;

% project PBE
load(['/home/kemerelab/Downloads/new_mazes_eeg/pbe_maze1_' num2str(tb) 'ms.mat'])
scaler=tb/500;
id_event = randsample(size(event_edge,1),100);
id_event = sort(id_event);
% n_event_range = [400,450]; 
% portion = event_edge(n_event_range(1),1)+1:event_edge(n_event_range(2),2)+1;
portion = [];
for i=1:100
    portion = [portion,event_edge(id_event(i),1)+1:event_edge(id_event(i),2)+1];
end
yy0 = double(spikes(:,portion))';
xx0 = zeros(size(yy0,1),1);
tgrid = (time(portion)-time(portion(1)))/tb*1000;
tgrid = tgrid';
event_length=[diff(event_edge(:,1));diff(event_edge(end,:))+1];
d = double([0;cumsum(event_length(id_event))]);

clearvars -except xx0 yy0 tgrid scaler id_event d
save('select_pbe_maze2_17ms.mat')

load('/home/kemerelab/Downloads/new_mazes_eeg/run_pbe_2mazes_m1_17ms/select_pbe_maze1_17ms.mat')
load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat')
tgrid_tc = setopt.tgrid;
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};
len = std_all/5;
k=6;

% project all
niter=20;
figure
[result_ori, setopt, fftc, ~,logpy_ori] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
    'newdatatype','pbe','scaler',scaler);

% project single PBEs
xxsamp = [];
ffmat = [];
logpy_ori = [];
niter=10;
for i=1:100
    por = d(i)+1:d(i+1);
    [result_ori, setopt, fftc, ~,logpy] = run_projection(xx0(por),yy0(por,:),tgrid(por),result_la,niter,d(i:i+1),'xbackground',xbackground,...
        'tgrid_tc',tgrid_tc,'newdatatype','pbe','scaler',scaler,'draw',0);
    xxsamp = [xxsamp;result_ori.xxsamp];
    ffmat = [ffmat;result_ori.ffmat];
    logpy_ori = [logpy_ori,logpy];
end
result_ori.xxsamp = xxsamp;
result_ori.ffmat = ffmat;
s_ori = estimate_smoothness(result_ori.xxsamp,d);
consist_seg = estimate_consistency(result_ori.xxsamp,xbackground,double(d),k,len);
consist_seg1_ori = consist_seg(1,:);
consist_seg2_ori = consist_seg(2,:);
consist_seg_ori = estimate_consistency(result_ori.xxsamp,{result_la.xxsamp},double(d),k,len);

save('ori_1.mat','result_ori', 'setopt', 'fftc','s_ori','logpy_ori','consist_seg1_ori','consist_seg2_ori','consist_seg_ori')

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
filename = {'cellshuf_pool.mat','segtimeshuf_pool.mat','localcirshuf_pool','indtimeshuf_pool.mat'};
shuftype=3;
len = std_all/5;
consist_seg1_all = [];
consist_seg2_all = [];
curv_all=[];
results=[];
logpy_all = {[],[],[]};
consist_seg_all = [];
stepdis_all = [];
knndis_all = [];


for j=1:9
    for i=1:20
        [result, setopt, fftc, order, logpy] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
            'newdatatype','pbe','scaler',scaler,'shuffle',type{shuftype},'draw',0);
        results=[results,result.xxsamp];
        [curv,stepdis] = estimate_smoothness(result.xxsamp,d,1);
        curv_all=[curv_all;curv];
        stepdis_all=[stepdis_all;stepdis];
        [knndis, ~,~] = projected_xx_measures_segs(result.xxsamp,{result_la.xxsamp},d,0,k);
        knndis_all = [knndis_all;knndis];
        consist_seg = estimate_consistency(result.xxsamp,xbackground,d,k,len);
        consist_seg1_all = [consist_seg1_all;consist_seg(1,:)];
        consist_seg2_all = [consist_seg2_all;consist_seg(2,:)];
        consist_seg_all = [consist_seg_all;estimate_consistency(result.xxsamp,{result_la.xxsamp},d,k,len)];
        for t=1:3
            logpy_all{t} = [logpy_all{t};logpy(t,:)];
        end
        plot(sum(logpy),'color',[.6,.6,.6])
    end
    save(filename{shuftype},'consist_seg_all','consist_seg1_all','consist_seg2_all','logpy_all','knndis_all','stepdis_all','curv_all');
    save([type{shuftype} 'results.mat'],'results')
end

for i=1:200
    consist_seg = estimate_consistency(results(:,i*3-2:i*3),xbackground,d,k,len);
    consist_seg1_all = [consist_seg1_all;consist_seg(1,:)];
    consist_seg2_all = [consist_seg2_all;consist_seg(2,:)];
end
figure
show_latent_variable(result.xxsamp,xx0,xbackground,setopt.tgrid,'line_only',1)
