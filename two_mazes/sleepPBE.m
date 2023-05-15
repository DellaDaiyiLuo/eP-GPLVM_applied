%% Load PBEs
NAMES = {'pre','maze2','post1','post2'};

for nameID = 1:4
    session = NAMES{nameID};
    tb = 17;
    load(['/home/kemerelab/Downloads/new_mazes_eeg/pbe_' session '_' num2str(tb) 'ms.mat'])
    n_event_range = [1,size(event_edge,1)];
    scaler=tb/500; 
    portion = event_edge(n_event_range(1),1)+1:event_edge(n_event_range(2),2)+1;
    d = double([event_edge(n_event_range(1):n_event_range(2),1);portion(end)]-event_edge(n_event_range(1),1));
    yy0 = double(spikes(:,portion))';
    xx0 = zeros(size(yy0,1),1);
    tgrid = (time(portion)-time(portion(1)))/tb*1000;
    tgrid = tgrid';

%     load('/home/kemerelab/Downloads/new_mazes_eeg/pyr_id.mat')
%     yy0 = yy0(:,pyrid);
%     clearvars -except xx0 yy0 tgrid scaler n_event_range d session

    %% project PBE individually
    load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat')
    xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};
    len = std_all/5;
    k=6;
    niter=10;

    xxsamp = [];
    ffmat = [];
    logpy_ori = [];
    tic
    for i=1:numel(d)-1
        por = d(i)+1:d(i+1);
        [result_ori, setopt, fftc, ~,logpy] = run_projection(xx0(por),yy0(por,:),tgrid(por),result_la,niter,[0,numel(por)],'xbackground',xbackground,...
            'newdatatype','pbe','scaler',scaler,'draw',0);
        xxsamp = [xxsamp;result_ori.xxsamp];
        ffmat = [ffmat;result_ori.ffmat];
        logpy_ori = [logpy_ori,logpy];
        disp(num2str(i))
    end
    toc
    result_ori.xxsamp = xxsamp;
    result_ori.ffmat = ffmat;
    setopt.tgrid = tgrid;
    s_ori = estimate_smoothness(result_ori.xxsamp,d);
    consist_seg = estimate_consistency(result_ori.xxsamp,xbackground,double(d),k,len);
    consist_seg1_ori = consist_seg(1,:);
    consist_seg2_ori = consist_seg(2,:);
    consist_seg_ori = estimate_consistency(result_ori.xxsamp,{result_la.xxsamp},double(d),k,len);

    save(['ori_' session '.mat'],'result_ori', 'setopt','s_ori','logpy_ori','consist_seg1_ori','consist_seg2_ori','consist_seg_ori')
end
show_latent_variable(result_ori.xxsamp,xx0,xbackground,tgrid,'line_only',1)

%% shuffle PBE
NAMES = {'pre','maze1','maze2','post1','post2'};
session = NAMES{3};

load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat','result_la','std_all')
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};
tb = 17;
scaler=tb/500;

load(['/home/kemerelab/Downloads/new_mazes_eeg/pbe_' session '_' num2str(tb) 'ms.mat'])
type = {'cell','segtime','cir','ind_time'};
filename = {'cellshuf_pool.mat','segtimeshuf_pool.mat','localcirshuf_pool.mat','indtimeshuf_pool.mat'};
shuftype=3;
len = std_all/5;
k=6;
niter=20;
batchsize=50;
result_ori.xxsamp = [];
result_ori.ffmat = [];
setopt.tgrid = [];
setopt.xplds = [];
logpy_ori = [];
consist_seg1_all = zeros(100,size(event_edge,1));
consist_seg2_all = zeros(100,size(event_edge,1));
logpy_all = {zeros(100,size(event_edge,1)),zeros(100,size(event_edge,1)),zeros(100,size(event_edge,1))};
consist_seg_all = zeros(100,size(event_edge,1));
stepdis_all = zeros(100,size(event_edge,1));
curv_all = zeros(100,size(event_edge,1));
knndis_all = zeros(100,size(event_edge,1));

n_event_range=0:batchsize:size(event_edge,1);
n_event_range(end)=size(event_edge,1);

for i_batch = 1:16 % FOR
    display(['batch' num2str(i_batch)])
    load(['/home/kemerelab/Downloads/new_mazes_eeg/pbe_' session '_' num2str(tb) 'ms.mat'])
    portion = event_edge(n_event_range(i_batch)+1,1)+1:event_edge(n_event_range(i_batch+1),2)+1;
    d = double([event_edge(n_event_range(i_batch)+1:n_event_range(i_batch+1),1);portion(end)]-event_edge(n_event_range(i_batch)+1,1));
    yy0 = double(spikes(:,portion))';
    xx0 = zeros(size(yy0,1),1);
    tgrid = (time(portion)-time(portion(1)))/tb*1000;
    tgrid = tgrid';
    clear time spikes
    r=n_event_range(i_batch)+1:n_event_range(i_batch+1);
    
    [r_ori, s, ~, ~,l_ori] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
                'newdatatype','pbe','scaler',scaler,'draw',0);
    result_ori.xxsamp = [result_ori.xxsamp;r_ori.xxsamp];
    result_ori.ffmat = [result_ori.ffmat;r_ori.ffmat];
    setopt.tgrid = [setopt.tgrid;s.tgrid];
    setopt.xplds = [setopt.xplds;s.xplds];
    logpy_ori = [logpy_ori,l_ori];
    save(['ori_' session '.mat'],'result_ori','setopt','logpy_ori');

    for i=1:100
        display(['batch' num2str(i_batch) ', shuf iter ' num2str(i)])
        [result,~, ~, ~,logpy] = run_projection(xx0,yy0,tgrid,result_la,niter,d,'xbackground',xbackground,...
            'newdatatype','pbe','scaler',scaler,'draw',0,'shuffle',type{shuftype});
        [curv,stepdis] = estimate_smoothness(result.xxsamp,d,1);
        curv_all(i,r)=curv;
        stepdis_all(i,r)=stepdis;
        [knndis, ~,~] = projected_xx_measures_segs(result.xxsamp,{result_la.xxsamp},d,0,k);
        knndis_all(i,r)=knndis;
        consist_seg = estimate_consistency(result.xxsamp,xbackground,d,k,len);
        consist_seg1_all(i,r)=consist_seg(1,:);
        consist_seg2_all(i,r)=consist_seg(2,:);
        consist_seg_all(i,r)= estimate_consistency(result.xxsamp,{result_la.xxsamp},d,k,len);
        for t=1:3
            logpy_all{t}(i,r)= logpy(t,:);
        end
        save(['/home/kemerelab/Downloads/new_mazes_eeg/sleep/' session '_' filename{shuftype}],'consist_seg_all','consist_seg1_all','consist_seg2_all','logpy_all','knndis_all','stepdis_all','curv_all');
    end

end

n_event_range = [1,size(event_edge,1)];
portion = event_edge(n_event_range(1),1)+1:event_edge(n_event_range(2),2)+1;
d = double([event_edge(n_event_range(1):n_event_range(2),1);portion(end)]-event_edge(n_event_range(1),1));
consist_seg = estimate_consistency(result_ori.xxsamp,xbackground,double(d),k,len);
consist_seg1_ori = consist_seg(1,:);
consist_seg2_ori = consist_seg(2,:);
consist_seg_ori = estimate_consistency(result_ori.xxsamp,{result_la.xxsamp},double(d),k,len);
[~,stepdis_ori] = estimate_smoothness(result_ori.xxsamp,d,1);
save(['ori_' session '.mat'],'result_ori','d','setopt','logpy_ori','consist_seg1_ori','consist_seg2_ori','consist_seg_ori','stepdis_ori');

%% plottings
% single session
logpyjoint=logpy_all{1}+logpy_all{2}+logpy_all{3};
logpyjoint=logpyjoint(:,1:2700);
logpyjoint_ori=sum(logpy_ori);
logpyjoint_ori=logpyjoint_ori(1:2700);
length=diff(d);
zmedllh_cir = (logpyjoint_ori-median(logpyjoint))./mad(logpyjoint,1);
validid=(zmedllh_cir>3)&(logpyjoint_ori>max(logpyjoint));
sum(validid)

% consist m1 vs. m2 with PBE ID 
figure
scatter3(consist_seg2_ori,consist_seg1_ori,1:2700,length(1:2700)*1.5)
hold on
scatter3(consist_seg2_ori(validid),consist_seg1_ori(validid),find(validid),length(validid)*1.5)

% consist m1 vs. m2
figure
scatter(consist_seg2_ori,consist_seg1_ori,length*1.5,'filled','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[1,1,1])
hold on
scatter(consist_seg2_ori(validid),consist_seg1_ori(validid),length(validid)*1.5,zmedllh_cir(validid),'filled')
title('PBE maze1, z-scored by cir-shuf')
c=colorbar;
c.Label.String='modified z-scored LLH';
axis([0,3,0,3])
xlabel('consist to m2')
ylabel('consist to m1')
grid

% consist ratio
cratio=max(consist_seg1_ori,consist_seg2_ori)./consist_seg_ori;
figure
scatter(consist_seg2_ori(validid),consist_seg1_ori(validid),30,cratio(validid))
title({'maze1, all sig PBEs','color: consist max(to m1, to m2) / to all'})

% step vs. consist
figure
hold on
scatter(stepdis_ori(not(validid)),consist_seg_ori(not(validid)),length(not(validid))*1.5,'MarkerFaceColor',[.7,.7,.7],'MarkerEdgeColor',[1,1,1])
scatter(stepdis_ori(validid),consist_seg_ori(validid),length(validid)*1.5,zmedllh_cir(validid),'filled','MarkerEdgeColor',[1,1,1]);
caxis([3,9])
c=colorbar;
c.Label.String='z-scored LLH';
ylabel('consistency')
xlabel('step distance')
title('PBE maze1, z-scored by cir-shuf')


% PBE ID vs. consist (m1/m2/all)
figure
plot(consist_seg1_ori(1:2700),'x')
hold on
plot(consist_seg2_ori(1:2700),'+')
plot(consist_seg_ori(1:2700))
plot(find(validid),consist_seg_ori(validid),'o')

% show all sig PBEs
figure
load('/home/kemerelab/Downloads/new_mazes_eeg/m1_101_500_m2_301_700_result_la.mat','result_la')
xbackground = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:end,:)};
i=find(validid,1);
show_latent_variable(result_ori.xxsamp,setopt.tgrid,xbackground,setopt.tgrid,'line_only',1,'part',d(i)+1:d(1+i))
for i=find(validid)
    show_latent_variable(result_ori.xxsamp,setopt.tgrid,[],setopt.tgrid,'line_only',1,'part',d(i)+1:d(1+i),'line_color',[min(cratio(i),1),0,0])
end
title({'maze1 PBE sig','color: consist max(to m1, to m2) / to all'})

figure
i=437;
show_latent_variable(result_ori.xxsamp,setopt.tgrid,xbackground,setopt.tgrid,'line_only',1,'part',d(i)+1:d(1+i))
show_latent_variable(setopt.xplds,setopt.tgrid,[],setopt.tgrid,'line_only',1,'part',d(i)+1:d(1+i),'line_color',[1,.5,.5])
title({['maze2 PBE No.',num2str(i)],['z-scored LLH among cir: ' num2str(zmedllh_cir(i))]})

% categorize
labels = {'n.s.','maze1','maze2','mixed'};
mixid = cratio<.9;
context1id = consist_seg1_ori>consist_seg2_ori;
p = [sum(not(validid)),sum(validid&context1id&not(mixid)),sum(validid&not(context1id)&not(mixid)),sum(validid&mixid)];
figure;pie(p,{num2str(p(1)),num2str(p(2)),num2str(p(3)),num2str(p(4))});legend(labels);colormap pink
title(['maze1, ' num2str(numel(validid)) ' PBEs'])

% all sessions
sessions = {'pre','maze1','post1','maze2','post2'};
axes = [];
for j=1:numel(sessions)
    load(['ori_' sessions{j} '.mat'])
%     figure(1)
%     subplot(2,3,j)
%     hold on
%     scatter(consist_seg1_ori,consist_seg2_ori,30,[1,find(round(diff(setopt.tgrid'))~=1)+1],'.')
% %     scatter(consist_seg1_ori,consist_seg2_ori,30,logpy_ori(1,:)+logpy_ori(2,:),'.')
%     xlabel('consist m1');ylabel('consist m2');
%     grid on
%     title(sessions{j})
%     axis([0,3,0,3])

%     figure(2)
%     for i=1:3
%     subplot(1,3,i)
%     hold on
%     histogram(logpy_ori(i,:),'DisplayStyle','stairs','Normalization','probability','BinWidth',1)
% %     histogram(logpy_ori(i,:),'DisplayStyle','stairs','Normalization','cdf','BinWidth',0.1)
%     end
    
    figure(3)
    ax = subplot(5,1,j);
    axes = [axes,ax];
    hold on
    scatter(1:size(consist_seg1_ori,2),consist_seg1_ori,5,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.4)
    scatter(1:size(consist_seg1_ori,2),consist_seg2_ori,5,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.4)
    plot(smoothdata(consist_seg1_ori,'movmedian',50),'LineWidth',2,'color',[0.00,0.45,0.74]);
    plot(smoothdata(consist_seg2_ori,'movmedian',50),'LineWidth',2,'color',[0.85,0.33,0.10]);
    title(sessions{j})
    ylim([0,3])

%     hold on
%     scatter3(logpy_ori(1,:),logpy_ori(2,:),logpy_ori(3,:),'.')
end
linkaxes(axes,'x')
xlabel('PBE ID')
legend(sessions)

i=i+1
plot(yy(:,i))