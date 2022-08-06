load('maze2_run_v1_3_v2_5_500ms_new.mat')
plot(pos_linear);
position = pos_linear(310:end); % maze1: 300:end, maze2: 310:end
spikes = double(spikes(:,310:end));

plot(position)

% manual_seg = [0,55,95,137,185,238,290,347,401,483,516]; % maze1 500ms run
manual_seg = [0,30,46,58,103,120,147,162,197,259,299,340,375,415,427,450,459,503,539,555,581,627,660,690,703,720,737,783,811,847]; %maze2
fwd = [];
bwd = [];
for i=1:numel(manual_seg)-1
    if rem(i,2)==1
        fwd = [fwd,manual_seg(i)+1:manual_seg(i+1)];
    else
        bwd = [bwd,manual_seg(i)+1:manual_seg(i+1)];
    end
end
figure;hold on;plot(position)
scatter(fwd,position(fwd),5)
scatter(bwd,position(bwd),5)

%-----------repeat for fwd and bwd -----------------%
pos_dir = position(bwd);
spk_dir = spikes(:,bwd);

pp = round(pos_dir/10);
tbl = tabulate(pp);
tc = zeros(size(spk_dir,1),size(tbl,1));
for i=1:size(tbl,1)
    if tbl(i,2)>0
        tc(:,i)=mean(spk_dir(:,pp==tbl(i,1)),2);
    end
end

tc = tc(:,tbl(:,2)>0);
tbl_ = tbl(tbl(:,2)>0,:);

% plot
[M, I] = max(tc,[],2);
[~, I_cell] = sort(I);
tc_n = tc./repmat(M,1,size(tc,2));
imshow(tc_n(I_cell,:))


tc_sm = smoothdata(tc,2,'gaussian',5);
figure;plot(tc(6,:));
hold on;plot(tc_sm(6,:))
title({'tc smoothed with Gaussian of window 10'})
R=corrcoef(tc_sm);
figure;image(R,'CDataMapping','scaled')

%------------------------%
tc_sm_fwd = tc_sm;
tbl_fwd = tbl_;

tc_sm_bwd = tc_sm;
tbl_bwd = tbl_;

save('bayes_tc_maze2_500ms.mat','fwd','bwd','position','tc_sm_fwd','tc_sm_bwd','tbl_fwd','tbl_bwd')
clearvars -except fwd bwd position tc_sm_fwd tc_sm_bwd tbl_fwd tbl_bwd

%% Decode PBE and line fitting

load('bayes_tc_maze2_500ms.mat')
tbl = [-1*tbl_bwd(end:-1:1,:);tbl_fwd];
tc_sm = [tc_sm_bwd(:,end:-1:1),tc_sm_fwd];

% two direction interleave2
[tbl,od] = sort(abs(tbl(:,1)));
tc_sm = tc_sm(:,od);

pbe_time_bin = 14;
load(['pbe_maze2_' num2str(pbe_time_bin) 'ms.mat'])

nt = size(spikes,2);
spikes = double(spikes);
scaler = pbe_time_bin/500; %0.06; %
tc_sc = tc_sm.*scaler; %tc = tc_sm.*scaler_cell+0.0001;
tc_sc = tc_sc+min(nonzeros(tc_sc))/10;

% tc_scc = tc.*scaler_ratio+0.00001;
% tc_sc = smoothdata(tc_scc,2,'gaussian',10);

loglikelihood = -repmat(sum(tc_sc',2)',nt,1) + spikes'*log(tc_sc);
[Ls, xinitidx] = max(loglikelihood,[],2);
xinit_ratio = tbl(xinitidx,1);

save('bayes_decoded_maze2maze2pbe.mat','xinit_ratio')

% figure
% plot(abs(xinit_ratio),'.')
% hold on;
% plot(event_edge(:,2),abs(xinit_ratio(event_edge(:,2))),'*')
% xlim([4739,4793]);
% xlim([4296,4350]);
%-------

matrix = exp(loglikelihood');
% matrix = exp(loglikelihood(4739:4793,:)'); 
matrix_n = zeros(size(matrix));
for i=1:size(matrix,2)
    matrix_n(:,i) = matrix(:,i)/sum(matrix(:,i));
%     matrix_n(:,i) = (matrix(:,i)-min(matrix(:,i)))/(max(matrix(:,i))-min(matrix(:,i)));
end
figure;image(1:size(matrix,2),tbl(:,1)*2,matrix_n,'CDataMapping','scaled')
c = gray;
c = flipud(c);
colormap(c);
set(gca,'YDir','normal')
xlabel('time bin')
ylabel('Inbound      Outbound')
xlim(event_edge(829,:))
cb = colorbar;
cb.Label.String='Probability';
title('PBE example2')


i=141;%2452; %
pberange = event_edge(i,1):event_edge(i,2);
x = zeros(numel(pberange),1000);
y = zeros(numel(pberange),1000);
k = 1;
for l=event_edge(i,1):event_edge(i,2)
    y(k,:) = randsample(tbl,1000,true,matrix_n(:,l));
    x(k,:) = k;
    k = k+1;
end
mdl = fitlm(x(:),y(:));
plot(mdl)
mdl.Rsquared.Ordinary
% [p,S] = polyfit(x(:),y(:),1);
% [yy,delta] = polyval(p,x(:,1),S);
% figure;plot(x(:),y(:),'.')
% hold on
% plot(x(:,1),yy)
xlabel(['time bin (' num2str(pbe_time_bin) 'ms)'])
ylabel('spatial bin (4cm)')
title(['error: ', num2str(error(i),'%.2f'), ' bins'])
legend('sample draw from posterior','fitting line')
figure;image(1:k-1,tbl(:,1),matrix_n(:,pberange),'CDataMapping','scaled')
c = gray;
c = flipud(c);
colormap(c);
set(gca,'YDir','normal')
title({['error: ', num2str(error(i),'%.2f'), ' bins'],'posterior, two dirs interleaved'})
xlabel(['time bin (' num2str(pbe_time_bin) 'ms)'])
ylabel('spatial bin (4cm)')

save(['linefit_error_' num2str(pbe_time_bin) 'ms.mat'],'error')



%% shuffle
% time shuffle
event_edge = event_edge+1-event_edge(1,1);
error = zeros(size(event_edge,1),1);
error_timeshuf = zeros(size(event_edge,1),1000);
percentile_timeshuf = zeros(size(event_edge,1),1);
llh_lst = [];
slope_lst = [];
p_lst = [];
for i=1:size(event_edge,1)
    pberange = event_edge(i,1):event_edge(i,2);
    pbelen = numel(pberange);
    x = zeros(pbelen,1000);
    y = zeros(pbelen,1000);
    k = 1;
    for l=pberange
        y(k,:) = randsample(tbl,1000,true,matrix_n(:,l));
        x(k,:) = k;
        k = k+1;
    end
    mdl = fitlm(x(:),y(:));
    error(i) = mdl.Rsquared.Ordinary;
    slope_lst = [slope_lst, mdl.Coefficients.Estimate(2)];
    p_lst = [p_lst, mdl.Coefficients.pValue(2)];
    
    for j=1:1000
        yshu = y(randperm(pbelen),:);
        mdl = fitlm(x(:),yshu(:));
        error_timeshuf(i,j) = mdl.Rsquared.Ordinary;
    end
    percentile_timeshuf(i) = numel(find(error_timeshuf(i,:)<error(i)))/1000;
    llh_lst = [llh_lst,sum(max(loglikelihood(pberange,:),[],2))];
    
    [i,percentile_timeshuf(i),p_lst(i)]
end

numel(find(percentile_timeshuf>0.95))/size(event_edge,1)
histogram(percentile_timeshuf,100)
title({['tc: maze1; PBEs: ' num2str(pbe_time_bin) 'ms time bin, maze1 (' num2str(size(event_edge,1)) ' events)'],'20.28% PBEs significant than time-shuffled'})
hold on;plot([0.95 0.95],[0,90])
xlabel('real R^2 percentile among time-shuffled versions')
save(['maze' num2str(maze) 'maze' num2str(maze) 'pbe_linefit_R2_' num2str(pbe_time_bin) 'ms_scaler28.mat'],'error','error_timeshuf','percentile_timeshuf','scaler','llh_lst')
clear error_timeshuf percentile_timeshuf

significant_28 = percentile_timeshuf>0.95;
[sum(significant_28) sum(significant_60) sum(significant_28 & significant_60)]
[sum(llh_lst_60(significant_60&significant_28)),sum(llh_lst_28(significant_60&significant_28))]
% cell_ID shuffle - not working, all heronzontal lines

segid = 1:size(event_edge,1);
sig_pbe_maze2_id = segid(significant_60&significant_28);
%%
pos_list = unique(tbl);
matrix_comb = zeros(numel(pos_list),size(matrix,2)); % sum of posterior prob of two dirs
for i=1:numel(pos_list)
    matrix_comb(i,:)=sum(matrix(tbl==pos_list(i),:),1);
end


