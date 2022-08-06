%% learn tc from maze data
load('maze1_run_v1_3_v2_5_500ms.mat')
idx=idx+1;
plot(pos_linear(idx));
position = pos_linear(idx);
spikes = spikes(:,idx);
portion = [1:891];% 671:820];% 901:1050];
xx = position(portion);
yy = double(spikes(:,portion))';
tgrid = double(idx(portion))';
    

%% spiking profiles (run vs pbes)
load('pbe_maze2_14ms.mat')
spikes_pbe = spikes;
time_pbe = time;
load('maze2_run_v1_3_v2_5_500ms_new.mat') % run
load('maze2_run_v1_3_v2_5_500ms.mat') %behav

% exclude pbe time bins in behavior data %
event_edge = event_edge - event_edge(1,1) +1;
time_behav_bin_st = time-0.25; % 0.5s time bin / 2
exclude_behav_bin = [];
for i=1:size(event_edge,1)
    st = find(time_behav_bin_st<time_pbe(event_edge(i,1)),1,'last');
    endidx = find(time_behav_bin_st>time_pbe(event_edge(i,2)),1);
    exclude_behav_bin = [exclude_behav_bin, st:endidx-1];
end
select = ones(size(pos_linear),'logical');
select(exclude_behav_bin) = 0;
spikes = spikes(:,select); % behavior excluding pbe bins
% spikes = spikes(:,idx); % run
%----------------------------------------%

figure(1)
h = histogram(sum(spikes,1),'Normalization','pdf');
hold on
figure;plot(h.Values)
hold on
sm_hist = smoothdata(h.Values,'gaussian',5);
plot(sm_hist)
[~,I]=max(sm_hist);
peak_spkcnt = mean(h.BinEdges(I:I+1))

figure(1)
h_pbe = histogram(sum(spikes_pbe,1),'Normalization','pdf')
figure(2);plot(h_pbe.Values)
hold on
sm_hist_pbe = smoothdata(h_pbe.Values,'gaussian',4);
plot(sm_hist_pbe)
[~,I]=max(sm_hist_pbe);
peak_spkcnt_pbe = mean(h_pbe.BinEdges(I:I+1));

scaler = peak_spkcnt_pbe/peak_spkcnt


%% scaler of mazes
load('maze1_run_v1_5_v2_10.mat')
spikes_run = spikes;
load('maze2_run_v1_5_v2_10.mat')
spikes_run = [spikes_run,spikes];
mu_run = mean(spikes_run,2);
load('pbe_maze1_6ms.mat')
mu_pbe = mean(spikes,2);

scaler = mu_pbe./mu_run;
scaler(mu_run==0) = 1;
plot(scaler)
hold on
plot([0,110],[0.05, 0.05])

save('scaler_maze12_run_maze1_pbe_6ms.mat','scaler')

%% maze1 four ends
end1=find(xx>423);
end2=find((xx>287)&(xx<320));
end3=find((xx>137)&(xx<180));
end4=find(xx<41);

figure
scatter3(pos_origin(1,idx),pos_origin(2,idx),xx,3,[0.7,0.7,0.7])
hold on
scatter3(pos_origin(1,idx(end1)),pos_origin(2,idx(end1)),xx(end1),3,[0.85,0.325,0.098])
scatter3(pos_origin(1,idx(end2)),pos_origin(2,idx(end2)),xx(end2),3,[0.929,0.694,0.125])
scatter3(pos_origin(1,idx(end3)),pos_origin(2,idx(end3)),xx(end3),3,[0.494,0.184,0.556])
scatter3(pos_origin(1,idx(end4)),pos_origin(2,idx(end4)),xx(end4),3,[0.466,0.674,0.188])


figure
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'line_only',1,'line_color',[0.7,0.7,0.7])
% plot3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),'Color',[0.7,0.7,0.7])
hold on 
scatter3(result_la.xxsamp(end1,1),result_la.xxsamp(end1,2),result_la.xxsamp(end1,3),3,[0.85,0.325,0.098])
scatter3(result_la.xxsamp(end2,1),result_la.xxsamp(end2,2),result_la.xxsamp(end2,3),3,[0.929,0.694,0.125])
scatter3(result_la.xxsamp(end3,1),result_la.xxsamp(end3,2),result_la.xxsamp(end3,3),3,[0.494,0.184,0.556])
scatter3(result_la.xxsamp(end4,1),result_la.xxsamp(end4,2),result_la.xxsamp(end4,3),3,[0.466,0.674,0.188])

%% bouts
figure;plot(xx)
% bp = [0,311,472,680,886,1055,1241,1414,1597,1795,2042,2198,2371,2484];
% bp = [0,114,165,249,322,363,405,444,622,660,707,787,860,964,1035];
bp = [0,66,144,181,304,372,414,450,482,533,591,682,751,853,891];
% figure(1);
% hold on
figure(2);
hold on
for i=1:numel(bp)-1
%     figure(1);
%     plot(bp(i)+1:bp(i+1),xx(bp(i)+1:bp(i+1)),'Color',[1,0,0]/numel(bp)*i)
    figure(2);
%     show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'line_only',1,'line_color',[1,0,0]/numel(bp)*i,'part',bp(i)+1:bp(i+1))
    plot(result_la.xxsamp(bp(i)+1:bp(i+1),1),result_la.xxsamp(bp(i)+1:bp(i+1),2),'Color',[1,0,0]/numel(bp)*i);
%     plot3(result_la.xxsamp(bp(i)+1:bp(i+1),1),result_la.xxsamp(bp(i)+1:bp(i+1),2),result_la.xxsamp(bp(i)+1:bp(i+1),3),'Color',[1,0,0]/numel(bp)*i);
%     i=i+1;
end

%% direction
fwd = [];
bwd = [];
for i=3:numel(bp)-1
    if rem(i,2)>0
        fwd = [fwd bp(i)+1:bp(i+1)];
    else
        bwd = [bwd bp(i)+1:bp(i+1)];
    end
end
figure;
plot(xx);hold on
scatter(fwd,xx(fwd),'.');scatter(bwd,xx(bwd),'.')
figure
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'line_only',1,'line_color',[.7,.7,.7])
hold on
scatter3(result_la.xxsamp(fwd,1),result_la.xxsamp(fwd,2),result_la.xxsamp(fwd,3),3,[0.85,0.325,0.098])
scatter3(result_la.xxsamp(bwd,1),result_la.xxsamp(bwd,2),result_la.xxsamp(bwd,3),3,[0.929,0.694,0.125])

%% bout routes
figure;plot(xx)
bp = [0,68,119,144,181,304,372,414,450,482,533,591,682,751,853,891];
hold on; plot(bp(2:end),xx(bp(2:end)),'*')

interp_traj = cell(numel(bp)-1,1);
nf = size(result_la.xxsamp,2);
for i=1:numel(bp)-1
    xseg = xx(bp(i)+1:bp(i+1));
    xgrid = round(min(xseg)):round(max(xseg));
    xxinterp = xgrid';
    for j=1:nf
        curve = fit(xseg,result_la.xxsamp(bp(i)+1:bp(i+1),j),'smoothingspline','SmoothingParam',0.5);
        y = feval(curve,xgrid);
        xxinterp = [xxinterp,y];
    end
    interp_traj{i}=xxinterp;
end
figure
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'scatter_only',1)
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'line_only',1,'part',bp(i)+1:bp(i+1))
plot3(interp_traj{i}(:,2),interp_traj{i}(:,3),interp_traj{i}(:,4))
i=i+1;

figure;hold on
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'scatter_only',1)
for i=1:numel(bp)-1
    plot(interp_traj{i}(:,2),interp_traj{i}(:,3),'Color',[1,0,0]/numel(bp)*i)
%     plot3(interp_traj{i}(:,2),interp_traj{i}(:,3),interp_traj{i}(:,4),'Color',[1,0,0]/numel(bp)*i)
    pause
end