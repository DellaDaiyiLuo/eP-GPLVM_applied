% data load using two_mazes_run.m
load('select_pbe_maze1_17ms.mat', 'd')

zmedllh_time = (data_ori.logpy-median(data_time.logpy))./mad(data_time.logpy,1);
zmedllh_cir = (data_ori.logpy-median(data_cir.logpy))./mad(data_cir.logpy,1);
sum(zmedllh_cir>3)
%-------------------------------------%
% cc_time=sum(data_time.logpy<data_ori.logpy,1)/size(data_time.logpy,1); % log_py value percentile
% cc_cir=sum(data_cir.logpy<data_ori.logpy)/size(data_cir.logpy,1);
% sum(cc_time>.95)
% sum(cc_cir>.95)

length=diff(d);

% figure
% scatter3(cc_time,cc_cir,1:100)
% xlabel('time-shuf')
% ylabel('ind-time-shuf')
% zlabel('id')

% figure
% hold on
% scatter(cc_time(not(validid)),cc_cir(not(validid)),length(not(validid))*1.5,'MarkerFaceColor',[.7,.7,.7],'MarkerEdgeColor',[1,1,1])
% scatter(cc_time(validid),cc_cir(validid),length(validid)*1.5,[zeros(sum(validid),2),cc_time(validid)'],'filled','MarkerEdgeColor',[1,1,1])
% xlabel('percentile among time-shuf')
% ylabel('percentile among ind-time-shuf')


% fig1. stepsize-consist
validid=(zmedllh_cir>3)&(data_ori.logpy_ori>max(data_cir.logpy_all));
figure
hold on
scatter(data_ori.stepdis_ori(not(validid)),data_ori.consist_seg_ori(not(validid)),length(not(validid))*1.5,'MarkerFaceColor',[.7,.7,.7],'MarkerEdgeColor',[1,1,1])
scatter(data_ori.stepdis_ori(validid),data_ori.consist_seg_ori(validid),length(validid)*1.5,zmedllh_cir(validid),'filled','MarkerEdgeColor',[1,1,1]);
caxis([3,10])
ylabel('consistency')
xlabel('step distance')

figure
hold on
scatter(data_ori.consist_seg2_ori(not(validid)),data_ori.consist_seg1_ori(not(validid)),length(not(validid))*1.5,'MarkerFaceColor',[.7,.7,.7],'MarkerEdgeColor',[1,1,1])
scatter(data_ori.consist_seg2_ori(validid),data_ori.consist_seg1_ori(validid),length(validid)*1.5,[cc_time(validid)',zeros(sum(validid),2)],'filled','MarkerEdgeColor',[1,1,1])
ylabel('consist to m1')
xlabel('consist to m2')


% ------------ single PBE
d = [0,find(round(diff(setopt.tgrid'))>1),numel(setopt.tgrid)];
line_color = [0,0,0];
i=33;
show_latent_variable(data_ori.setopt.xplds,data_ori.setopt.tgrid,[],data_ori.setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',[1,.58,0.66])
show_latent_variable(data_ori.result_ori.xxsamp,data_ori.setopt.tgrid,[],data_ori.setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',line_color)
title(i)
scatter3(data_ori.result_ori.xxsamp(d(i)+1,1),data_ori.result_ori.xxsamp(d(i)+1,2),data_ori.result_ori.xxsamp(d(i)+1,3),50,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
for k=1:50
    show_latent_variable(results(:,k*3-2:k*3),setopt.tgrid,[],setopt.tgrid,'line_only',1,'part',d(i)+1:d(i+1),'line_color',[1,0,0])
end