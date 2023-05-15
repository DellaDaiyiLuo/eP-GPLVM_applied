load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat')
tgrid_tc = setopt.tgrid;
xx_tc = xx;

%% Position vs time overview

figure
hold on
% maze1
load('/home/kemerelab/Downloads/new_mazes_eeg/maze1_run_v1_3_v2_5_500ms.mat', 'pos_origin', 'pos_linear', 'idx')
cmap1 = parula(256);
cdata1 = interp1(linspace(min(xx_tc(1:400)),max(xx_tc(1:400)),256)',cmap1,pos_linear);
portion_tc1 = 101:500;
s1=scatter(1:numel(pos_linear),pos_linear,30,cdata1,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3);
d_tc = [0,find(diff(idx(portion_tc1))>1),numel(portion_tc1)];
for i=1:numel(d_tc)-1
    seg = portion_tc1(d_tc(i)+1:d_tc(i+1));
    s2=plot(seg,pos_linear(seg),'LineWidth',2,'Color',[.85,.3,.1]);
    s6=scatter(seg(1),pos_linear(seg(1)),50,'MarkerFaceColor',[.85,.3,.1],'MarkerEdgeColor','None');
end
% maze2
load('/home/kemerelab/Downloads/new_mazes_eeg/maze2_run_v1_3_v2_5_500ms.mat', 'pos_origin', 'pos_linear', 'idx')
xx_tc(401:800) = xx_tc(401:800)+145;
pos_linear = pos_linear+145;
cmap2 = copper(256);
cdata2 = interp1(linspace(min(xx_tc(401:800)),max(xx_tc(401:800)),256)',cmap2,pos_linear);
portion_tc2 = 301:700;
s3=scatter([1:numel(pos_linear)]+size(cdata1,1),pos_linear,30,cdata2,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3);
d_tc = [0,find(diff(idx(portion_tc2))>1),numel(portion_tc2)];
for i=1:numel(d_tc)-1
    seg = portion_tc2(d_tc(i)+1:d_tc(i+1));
    s4=plot(seg+size(cdata1,1),pos_linear(seg),'LineWidth',2,'Color',[0.00,0.45,0.74]);
    s7=scatter(seg(1)+size(cdata1,1),pos_linear(seg(1)),50,'MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','None');
end

s5=plot(part,xx0(d(i)+1:d(i+1)),'LineWidth',2,'Color',line_color);
scatter(part(1),xx0(d(i)+1),30,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
plot([1300,1340],[400,400],'k')
legend([s1,s2,s3,s4,s5,s6,s7],{'maze1 all','maze1 train','maze2 all','maze2 train', 'test data','start','start2'})
ylim([-9,457])
ylabel('animal position')

%% TC in latent space
figure
show_latent_variable(result_la.xxsamp,xx,[],tgrid_tc,'line_only',1,'line_color',[.7,.7,.7])
scatter3(result_la.xxsamp(1:400,1),result_la.xxsamp(1:400,2),result_la.xxsamp(1:400,3),120,cdata1(portion_tc1,:),'Marker','.')
scatter3(result_la.xxsamp(401:800,1),result_la.xxsamp(401:800,2),result_la.xxsamp(401:800,3),120,cdata2(portion_tc2,:),'Marker','.')
%% Projection in latent space
load('ori_.mat')
portion = 508:667;
d = [0,find(diff(setopt.tgrid')>1),numel(setopt.tgrid)];
i=3;
line_color = [0,0,0];
part = portion(d(i)+1:d(i+1));

figure
scatter3(result_la.xxsamp(1:400,1),result_la.xxsamp(1:400,2),result_la.xxsamp(1:400,3),25,cdata1(portion_tc1,:),...
    'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3)
hold on
scatter3(result_la.xxsamp(401:800,1),result_la.xxsamp(401:800,2),result_la.xxsamp(401:800,3),25,cdata2(portion_tc2,:),...
    'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3)
show_latent_variable(setopt.xplds,xx0,[],setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',[1,.58,0.66])
show_latent_variable(result_ori.xxsamp,xx0,[],setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',line_color)
scatter3(result_ori.xxsamp(d(i)+1,1),result_ori.xxsamp(d(i)+1,2),result_ori.xxsamp(d(i)+1,3),50','MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
set(gca,'CameraPosition',[-1100,700,66]);
set(gca,'CameraPosition',[490,850,235]);
exportgraphics(gcf,'fig_seg3_timeshuf_latent1.pdf','BackgroundColor','none','ContentType','vector')
%% 2D pos

figure
scatter(pos_origin(1,portion_tc1),pos_origin(2,portion_tc1),36,cdata1(portion_tc1,:),'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3)
hold on
xx_origin = pos_origin(:,part);
plot(xx_origin(1,:),xx_origin(2,:),'LineWidth',2,'Color',line_color)
scatter(xx_origin(1,1),xx_origin(2,1),50,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
load('/home/kemerelab/Downloads/new_mazes_eeg/maze2_run_v1_3_v2_5_500ms.mat', 'pos_origin')
scatter(pos_origin(1,portion_tc2)+200,pos_origin(2,portion_tc2),36,cdata2(portion_tc2,:),'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3)
axis([-100,120,-90,100])

%% metrics
data_time = load('segtimeshuf_pool.mat');
data_cell = load('cellshuf_pool.mat');
data_cir=load('localcirshuf_pool.mat');
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

zmedllh_time = (data_ori.logpy_ori-median(data_time.logpy_all))./mad(data_time.logpy_all,1);
zmedllh_cir = (data_ori.logpy_ori-median(data_cir.logpy_all))./mad(data_cir.logpy_all,1);

length_seg=data_ori.d(i+1)-data_ori.d(i);
colors=[1,0,0;0.49,0.18,0.56;0.00,0.45,0.74;0.47,0.67,0.19];

% logpy histogram
figure
hold on
vall = length_seg*[min(data_cell.logpy_all(:,i)),max(data_cell.logpy_all(:,i)),min(data_cir.logpy_all(:,i)),max(data_cir.logpy_all(:,i)),min(data_ori.logpy_ori(:,i)),max(data_ori.logpy_ori(:,i))];
edge = linspace(-4800,max(vall),70);
% edge = linspace(min(vall),max(vall),70);
h1=histogram(data_cell.logpy_all(:,i)*length_seg,edge,'EdgeColor','None','FaceColor',colors(2,:),'Normalization','probability');
h2=histogram(data_cir.logpy_all(:,i)*length_seg,edge,'EdgeColor','None','FaceColor',colors(3,:),'Normalization','probability');
h3=histogram(data_time.logpy_all(:,i)*length_seg,edge,'EdgeColor','None','FaceColor',colors(4,:),'Normalization','probability');
h4=plot(ones(2,1)*data_ori.logpy_ori(i)*length_seg,ylim,'--','Color',colors(1,:),'LineWidth',2);
legend([h1,h2,h3,h4],{'cell ID-shuf','cir-shuf','time-shuf','ori'})
xlabel('log likelihood of test data')
ylabel('relative frequency')

% metric stepdis vs. consist
figure 
ex=[1,80,430]; % single maze
scatter(data_cell.stepdis_all(:,i),data_cell.consist_seg_all(:,i),'filled','MarkerFaceColor',colors(2,:),'MarkerFaceAlpha',.3)
hold on
scatter(data_cir.stepdis_all(:,i),data_cir.consist_seg_all(:,i),'filled','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',.3)
scatter(data_time.stepdis_all(:,i),data_time.consist_seg_all(:,i),'filled','MarkerFaceColor',colors(4,:),'MarkerFaceAlpha',.3)
scatter(data_ori.stepdis_ori(i),data_ori.consist_seg_ori(i),70,'+','MarkerEdgeColor',colors(1,:))
title(num2str(i))
% examples
scatter(data_cell.stepdis_all(ex(1),i),data_cell.consist_seg_all(ex(1),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(2,:))
scatter(data_cir.stepdis_all(ex(2),i),data_cir.consist_seg_all(ex(2),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(3,:))
scatter(data_time.stepdis_all(ex(3),i),data_time.consist_seg_all(ex(3),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(4,:))
xlabel('step distance')
ylabel('consistency')
legend({'cell ID-shuf','cir-shuf','time-shuf','original','example1','example2','example3'})

figure 
scatter(data_cell.consist_seg2_all(:,i),data_cell.consist_seg1_all(:,i),'filled','MarkerFaceColor',colors(2,:),'MarkerFaceAlpha',.3)
hold on
scatter(data_cir.consist_seg2_all(:,i),data_cir.consist_seg1_all(:,i),'filled','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',.3)
scatter(data_time.consist_seg2_all(:,i),data_time.consist_seg1_all(:,i),'filled','MarkerFaceColor',colors(4,:),'MarkerFaceAlpha',.3)
scatter(data_ori.consist_seg2_ori(i),data_ori.consist_seg1_ori(i),70,'+','MarkerEdgeColor',colors(1,:))
% examples
scatter(data_cell.consist_seg2_all(ex(1),i),data_cell.consist_seg1_all(ex(1),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(2,:))
scatter(data_cir.consist_seg2_all(ex(2),i),data_cir.consist_seg1_all(ex(2),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(3,:))
scatter(data_time.consist_seg2_all(ex(3),i),data_time.consist_seg1_all(ex(3),i),20,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(4,:))
xlabel('consistency to maze2')
ylabel('consistency to maze1')
axis([0,3,0,3])
plot([0,3],[0,3],'--','color',[.6,.6,.6])
grid on
legend({'cell ID-shuf','cir-shuf','time-shuf','original','example1','example2','example3'})

% shuf examples
figure 
load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat','result_la','xx')
load('ori_.mat', 'setopt','xx0')
d = [0,find(diff(setopt.tgrid')>1),numel(setopt.tgrid)];
shuf_type = {'cell','cir','segtime'};
type = 2;
line_color=[0,0,0];
load([shuf_type{type},'results.mat'])
% scatter3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),36,xx,...
%     'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.3)
hold on
show_latent_variable(results(:,ex(type)*3-2:ex(type)*3),xx0,[],setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',line_color)
scatter3(results(d(i)+1,ex(type)*3-2),results(d(i)+1,ex(type)*3-1),results(d(i)+1,ex(type)*3),50','MarkerFaceColor',[0.6,0,0],'MarkerEdgeColor',[0.6,0,0])
xlabel('PC1');ylabel('PC2');zlabel('PC3');
title(['example' num2str(type) '. individual circular shuffled'])

set(gca,'CameraPosition',[-360,800,300]);
exportgraphics(gcf,['fig_seg3_' shuf_type{type} 'shuf_latent1.pdf'],'BackgroundColor','none','ContentType','vector')
set(gca,'CameraPosition',[490,850,235]);
exportgraphics(gcf,['fig_seg3_' shuf_type{type} 'shuf_latent2.pdf'],'BackgroundColor','none','ContentType','vector')