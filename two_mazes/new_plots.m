%% add metrics
load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat', 'result_la')
% load('m1_101_500_result_la.mat', 'result_la')
load('segtimeresults.mat')
load('segtimeshuf_pool.mat')

k=6;
% stepdis_all = [];
knndis_all = [];
for i=1:size(curv_all,1)
%     [~,stepdis] = estimate_smoothness(results(:,i*3-2:i*3),d,1);
%     stepdis_all=[stepdis_all;stepdis];
    [knndis_seg, ~,~] = projected_xx_measures_segs(results(:,i*3-2:i*3),{result_la.xxsamp},d,0,k);
    knndis_all=[knndis_all;knndis_seg];
end
save('cellshuf_pool.mat','consist_seg1_all','consist_seg2_all','logpy_all','s_all','d','curv_all','stepdis_all','consist_seg_all','knndis_all');




%% plots
data_time = load('segtimeshuf_pool.mat');
data_cir = load('localcirshuf_pool.mat');
data_cell = load('cellshuf_pool.mat');
data_ori = load('ori_.mat');

i=3;
if size(data_ori.logpy_ori,1)>1
    data_ori.logpy_ori=sum(data_ori.logpy_ori);
end
if numel(data_time.logpy_all)==3
    data_time.logpy_all=data_time.logpy_all{1}+data_time.logpy_all{2}+data_time.logpy_all{3};
end
if numel(data_cir.logpy_all)==3
    data_cir.logpy_all=data_cir.logpy_all{1}+data_cir.logpy_all{2}+data_cir.logpy_all{3};
end
if numel(data_cell.logpy_all)==3
    data_cell.logpy_all=data_cell.logpy_all{1}+data_cell.logpy_all{2}+data_cell.logpy_all{3};
end

[size(data_ori.logpy_ori);size(data_time.logpy_all);size(data_cir.logpy_all);size(data_cell.logpy_all)]
%% logpy histogram

figure
hold on
vall = [data_cell.logpy_all(:,i);data_cir.logpy_all(:,i);data_time.logpy_all(:,i);ones(2,1)*data_ori.logpy_ori(i)];
edge = linspace(min(vall),max(vall),100);
h1=histogram(data_cell.logpy_all(:,i),edge,'EdgeColor','None','FaceColor',[0.00,0.45,0.74],'Normalization','probability');
h2=histogram(data_cir.logpy_all(:,i),edge,'EdgeColor','None','FaceColor',[0.47,0.67,0.19],'Normalization','probability');
h3=histogram(data_time.logpy_all(:,i),edge,'EdgeColor','None','FaceColor',[0.85,0.67,0.19],'Normalization','probability');
h4=plot(ones(2,1)*data_ori.logpy_ori(i),ylim,'--','Color',[.7,0,0],'LineWidth',2);
legend([h1,h2,h3,h4],{'cell-ID-shuf','ind-time-shuf','time-shuf','ori'})
xlabel('log likelihood of test data')
ylabel('relative frequency')


figure 
ex=[1,5];
scatter3(data_cell.stepdis_all(:,i),data_cell.consist_seg1_all(:,i),data_cell.knndis_all(:,i),'filled','MarkerFaceColor',[0.00,0.45,0.74],'MarkerFaceAlpha',.3)
hold on
scatter3(data_cir.stepdis_all(:,i),data_cir.consist_seg_all(:,i),data_cir.knndis_all(:,i),'filled','MarkerFaceColor',[0.47,0.67,0.19],'MarkerFaceAlpha',.3)
scatter3(data_time.stepdis_all(:,i),data_time.consist_seg1_all(:,i),data_time.knndis_all(:,i),'filled','MarkerFaceColor',[0.85,0.67,0.19],'MarkerFaceAlpha',.3)
scatter3(data_ori.stepdis_ori(i),data_ori.consist_seg_ori(i),data_ori.knndis_ori(i),'x','MarkerEdgeColor',[.7,0,0])
title(num2str(i))
% examples
scatter(data_cell.s_all(ex(1),i),data_cell.consist_seg_all(ex(1),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.00,0.45,0.74])
scatter(data_time.s_all(ex(2),i),data_time.consist_seg_all(ex(2),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.47,0.67,0.19])
xlabel('stepdis')
ylabel('consistency')
zlabel('knndis')
legend({'cell-shuf','ind-time-shuf','time-shuf','original','example1','example2'})


figure
scatter(data_cell.consist_seg2_all(:,i),data_cell.consist_seg1_all(:,i),'filled','MarkerFaceColor',[0.00,0.45,0.74],'MarkerFaceAlpha',.3)
hold on
scatter(data_cir.consist_seg2_all(:,i),data_cir.consist_seg1_all(:,i),'filled','MarkerFaceColor',[0.47,0.67,0.19],'MarkerFaceAlpha',.3)
scatter(data_time.consist_seg2_all(:,i),data_time.consist_seg1_all(:,i),'filled','MarkerFaceColor',[0.85,0.67,0.19],'MarkerFaceAlpha',.3)
scatter(data_ori.consist_seg2_ori(i),data_ori.consist_seg1_ori(i),'x','MarkerEdgeColor',[.7,0,0])
title(num2str(i))
xlabel('consist to maze2')
ylabel('consist to maze1')
axis([0,3,0,3])
grid 
plot([0,3],[0,3])