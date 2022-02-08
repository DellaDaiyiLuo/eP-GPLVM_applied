%% bayesian model during run
% estimaate tc
load('track1_100ms.mat')
spikes = spikes(:,idx);
position = position(idx);

spikes = spikes(:,1000:2000);
position = position(1000:2000);

pp = round(position/2);
tbl = tabulate(pp);
tc = zeros(size(spikes,1),size(tbl,1));
for i=1:size(tbl,1)
    if tbl(i,2)>0
        tc(:,i)=mean(spikes(:,pp==tbl(i,1)),2);
    end
end

tc = tc(:,tbl(:,2)>0);
tbl_ = tbl(tbl(:,2)>0,:);

% plot
[M, I] = max(tc,[],2);
[~, I_cell] = sort(I);
tc_n = tc./repmat(M,1,size(tc,2));
imshow(tc_n(I_cell,:))

save('t1_bayesian_tc_run_1000_2000.mat', 'tc', 'pp','tbl');

%% predict location of PBEs
load('t1_bayesian_tc_run_1000_2000.mat')
tc_sm = smoothdata(tc,2,'gaussian',10);

load('track1_PBEs_4ms.mat')

nt = size(spikes,2);
spikes = double(spikes);

% % same scaler across cells
% L = [];
% PosDist = [];
% norms = [];
% scalerlist= [];
% for scaler = 1:50
%     scalerlist=[scalerlist scaler];
%     tc = tc_sm/scaler+0.0001;
%     loglikelihood = -repmat(sum(tc',2)',nt,1) + spikes'*log(tc);
%     [Ls, xinitidx] = max(loglikelihood,[],2);
%     L = [L sum(Ls)];
%     xinit = tbl(xinitidx,1);
%     [N, edges] = histcounts(xinit,20,'BinLimits',[tbl(1,1),tbl(end,1)]);
%     PosDist = [PosDist;N];
%     norms = [norms norm(N)];
% end




% different scalers for different cells
% iterations: best decoded position -> best match scaler for each cell
% init position by scaler=18 (maximize log likelihood when cells have same
% scaler)

tc = tc_sm+0.0001;
loglikelihood = -repmat(sum(tc',2)',nt,1) + spikes'*log(tc);
[Ls, xinitidx] = max(loglikelihood,[],2);
    
scalerlist=1./[1:50]';
niters=20;
scalercells=zeros(117,niters);
L=[];
xdiff=[];
PosDist=[];
norms=[];

for i=1:niters
    xinitidx0=xinitidx;
    
    % update scaler for each cell
    for cell=1:117
        tc=tc_sm(cell,:)+0.0001;
        lambdas=scalerlist*tc(xinitidx);
        llh=repmat(spikes(cell,:),numel(scalerlist),1).*log(lambdas)-lambdas;
        [~,m]=max(sum(llh,2));
        scalercells(cell,i)=scalerlist(m);
    end
    
    % update decoded position
    tc = tc_sm.*scalercells(:,i)+0.0001;
    loglikelihood = -repmat(sum(tc',2)',nt,1) + spikes'*log(tc);
    [Ls, xinitidx] = max(loglikelihood,[],2);
    L = [L sum(Ls)];
    xdiff=[xdiff norm(xinitidx-xinitidx0)];
    
    xinit = tbl(xinitidx,1);
    [N, edges] = histcounts(xinit,20,'BinLimits',[tbl(1,1),tbl(end,1)]);
    PosDist = [PosDist;N];
    norms = [norms norm(N)];
end

figure;histogram(scalercells(:,end))
title('histogram of cell scalers')
xlabel('scaler (rate=rate*scaler)')
ylabel('cell counts')

figure;
plot(xdiff)
title('Decoded position change in each iteration')
xlabel('iteration')

figure; % plot log likelihood
plot(L)
title('Log likelihood of PBE data')
xlabel('iteration')

figure; % plot norm of position histogram
plot(norms)
title('2-norm of position histogram')
xlabel('iteration')

figure; % plot histogram at some scaler
hold on;
plot(0.5*(edges(1:end-1)+edges(2:end)),PosDist(i,:))
i=i+10;
title('histogram of PBE decoded position')
xlabel('position')
ylabel('counts')

figure;plot(idx,xinit,'.');