load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat')
tgrid_tc = setopt.tgrid;
xx_tc = xx;

%% Position vs time overview

% maze1
load('/home/kemerelab/Downloads/new_mazes_eeg/maze1_run_v1_3_v2_5_500ms.mat', 'pos_origin', 'pos_linear', 'idx')
cmap1 = parula(256);
cdata1 = interp1(linspace(min(xx_tc(1:400)),max(xx_tc(1:400)),256)',cmap1,pos_linear);
% maze2
load('/home/kemerelab/Downloads/new_mazes_eeg/maze2_run_v1_3_v2_5_500ms.mat', 'pos_origin', 'pos_linear', 'idx')
xx_tc(401:800) = xx_tc(401:800)+145;
pos_linear = pos_linear+145;
cmap2 = copper(256);
cdata2 = interp1(linspace(min(xx_tc(401:800)),max(xx_tc(401:800)),256)',cmap2,pos_linear);

%% TC in latent space
figure
show_latent_variable(result_la.xxsamp,xx,[],tgrid_tc,'line_only',1,'line_color',[.7,.7,.7])
scatter3(result_la.xxsamp(1:400,1),result_la.xxsamp(1:400,2),result_la.xxsamp(1:400,3),120,cdata1(portion_tc1,:),'Marker','.')
scatter3(result_la.xxsamp(401:800,1),result_la.xxsamp(401:800,2),result_la.xxsamp(401:800,3),120,cdata2(portion_tc2,:),'Marker','.')


%% metrics
data_time = load('segtimeshuf_pool.mat');
data_cell = load('cellshuf_pool.mat');
data_cir = load('localcirshuf_pool.mat');
data_ori = load('ori_.mat');
load('/home/kemerelab/Downloads/new_mazes_eeg/run_pbe_2mazes_m1_17ms/select_pbe_maze1_17ms.mat', 'd')

data_cell.logpy=data_cell.logpy_all;
data_cir.logpy=data_cir.logpy_all{1}+data_cir.logpy_all{2}+data_cir.logpy_all{3};
data_time.logpy=data_time.logpy_all;
data_ori.logpy = sum(data_ori.logpy_ori);

% to pbe_plot ---------------------


% % z-score
% vars = {'consist_seg1','consist_seg2','consist_seg','s','logpy'};
% 
% for j=4 %1:numel(vars)-1
% %     m_cell = mean(data_cell.([vars{j},'_all']),1);
% %     std_cell = std(data_cell.([vars{j},'_all']),1);
% %     z_score_cell.([vars{j}])=(data_ori.([vars{j},'_ori'])-m_cell)./std_cell;
%     
%     m_time = mean(data_time.([vars{j},'_all']),1);
%     std_time = std(data_time.([vars{j},'_all']),1);
%     z_score_time.([vars{j}])=(data_ori.([vars{j},'_ori'])-m_time)./std_time;
% end
% figure
% scatter3(z_score_cell.s,z_score_cell.consist_seg,data_ori.consist_seg_ori,'*')
% xlabel('z-scored s')
% ylabel('z_scored consist')
% zlabel('consist')

% for j=1:3
%     m_cell = mean(data_cell.('logpy_all'){j},1);
%     std_cell = std(data_cell.('logpy_all'){j},1);
% %     z_score_cell.(['logpy',num2str(j)])=(data_ori.('logpy_ori')(j,:)-max(data_cell.('logpy_all'){j},[],1))./std_cell;
%     z_score_cell.(['logpy',num2str(j)])=(data_ori.('logpy_ori')(j,:)-m_cell)./std_cell;
% end
% figure
% scatter3(z_score_cell.logpy1,z_score_cell.logpy2,z_score_cell.logpy3)
% xlabel('log p(Y|F)')
% ylabel('log p(F|X)')
% zlabel('log p(X)')
% title({'z-scored among test-TC scale cell-shuf','m1 pbe, TC mazes m2'})

cc_time=sum(data_time.logpy_all<data_ori.logpy_ori,1)/size(data_time.logpy_all,1); % log_py value percentile
cc_cell=sum(data_cell.logpy_all<data_ori.logpy_ori,1)/size(data_cell.logpy_all,1);
find(cc_cell<1)
numel(find(cc_time>=1))
numel(find(cc_time>.95))

% threshold consistency value of PBEs
consis_seg_thr = prctile(data_cell.consist_seg_all(:),99);
figure
histogram(data_cell.consist_seg_all(:),'EdgeColor','None','Normalization','probability');
hold on
plot(ones(2,1)*consis_seg_thr,ylim,'--')
ylabel('relative frequency')
xlabel('consistency value')

% metric and selected PBEs
PBEid = [24,30,96];%[39,54,23]; % maze1: [46,33,96]; %maze2:
markershape = {'o','square','^'};
% 1-------------
length=diff(d);
color = zeros(100,3);
color(cc_time>0.95,1)=0.7;
figure
hold on
scatter(z_score_time.s,data_ori.consist_seg_ori,length*1.5,color,'filled','MarkerEdgeColor',[1,1,1])
% scatter(z_score_time.s,data_ori.consist_seg_ori,diff(d),cc_time,'o','filled')
% % scatter(z_score_time.s,data_ori.consist_seg_ori,100,diff(d),'.')
% hold on
plot(xlim,ones(1,2)*consis_seg_thr,'--','Color',[.65,.65,.65])
labels = [];
for k=1:3
    scatter(z_score_time.s(PBEid(k)),data_ori.consist_seg_ori(PBEid(k)),80,markershape{k},'MarkerEdgeColor',[0.93,0.69,0.13])
    labels = [labels,{['No.',num2str(PBEid(k))]}];
end
ylabel('consistency value')
xlabel('z-scored smoothness')
legend([{'n.s. PBEs','sig. PBEs','threshold'},labels])
% colormap copper
% cb = colorbar;
% cb.Label.String='LLH percentile (time-shuf)';


% 2-------------
validid=data_ori.consist_seg_ori>consis_seg_thr;

figure
hold on
scatter(data_ori.consist_seg2_ori(not(validid)),data_ori.consist_seg1_ori(not(validid)),length(not(validid))*1.5,'filled','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[.6,.6,.6])
scatter(data_ori.consist_seg2_ori(validid),data_ori.consist_seg1_ori(validid),length(validid)*1.5,color(validid,:),'filled','MarkerEdgeColor',[1,1,1])
axis([0,2.5,0,2.5])
xticks(0:0.5:2.5)
labels = [];
for k=1:3
    scatter(data_ori.consist_seg2_ori(PBEid(k)),data_ori.consist_seg1_ori(PBEid(k)),80,markershape{k},'MarkerEdgeColor',[0.93,0.69,0.13])
end
xlabel('consistency to maze2')
ylabel('consistency to maze1')
plot([0,2.5],[0,2.5],'--','Color',[.65,.65,.65])
grid on

% 3-----------
figure
scatter(length*0.017,cc_time,length*1.5,color,'o','filled','MarkerEdgeColor',[1,1,1])
xlabel('PBE duration (s)')
ylabel('LLH percentile (time-shuf)')
% colormap copper
axis([0,1,-0.01,1.01])
hold on
for k=1:3
    scatter(0.016*length(PBEid(k)),cc_time(PBEid(k)),80,markershape{k},'MarkerEdgeColor',[0.93,0.69,0.13])
end
grid on

% 4---------
data_cell.logpy=data_cell.logpy_all{1}+data_cell.logpy_all{2}+data_cell.logpy_all{3};
time_llh=sum(data_time.logpy_all<sum(data_ori.logpy_ori))/size(data_time.logpy_all,1);
cell_llh=sum(data_cell.logpy<sum(data_ori.logpy_ori))/size(data_cell.s_all,1);
figure;scatter(time_llh,cell_llh)
hold on
plot(xlim,[.95,.95])

%% Projection in latent space

d = [0,find(round(diff(data_ori.setopt.tgrid'))>1),numel(data_ori.setopt.tgrid)];
line_color = [0,0,0];
i=PBEid(3);
show_latent_variable(data_ori.setopt.xplds,data_ori.setopt.tgrid,[],data_ori.setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',[1,.58,0.66])
show_latent_variable(data_ori.result_ori.xxsamp,data_ori.setopt.tgrid,[],data_ori.setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',line_color)
title(i)
scatter3(data_ori.result_ori.xxsamp(d(i)+1,1),data_ori.result_ori.xxsamp(d(i)+1,2),data_ori.result_ori.xxsamp(d(i)+1,3),50,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])

for k=1:50
    show_latent_variable(results(:,k*3-2:k*3),data_ori.setopt.tgrid,[],data_ori.setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',[1,0,0])
end


