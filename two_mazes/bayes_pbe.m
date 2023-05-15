maze_id = 1;
load(['/home/kemerelab/Downloads/new_mazes_eeg/maze',num2str(maze_id),'_run_v1_3_v2_5_500ms.mat'])
switch maze_id
    case 1
        tc_portion = 101:500;
        manual_seg = [0,56,111,150,193,227,289,345,400];
    case 2
        tc_portion = 301:700;
        manual_seg = [0,5,31,58,72,100,131,169,210,247,271,284,298,310,342,379,400];
end
plot(pos_linear(tc_portion));
position = pos_linear(tc_portion);
spikes = double(spikes(:,tc_portion));

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
vars = {'fwd','bwd'};
for var_id = 1:2
    pos_dir = position(eval(vars{var_id}));
    spk_dir = spikes(:,eval(vars{var_id}));

    pos_bin = floor((max(pos_dir)-min(pos_dir))/29);
    pp = round((pos_dir-min(pos_dir))/pos_bin);
    tbl = tabulate(pp);
    tc = zeros(size(spk_dir,1),size(tbl,1));
    for i=1:size(tbl,1)
        if tbl(i,2)>0
            tc(:,i)=mean(spk_dir(:,pp==tbl(i,1)),2);
        end
    end

    tc = tc(:,tbl(:,2)>0);
    tbl_ = tbl(tbl(:,2)>0,:);

    tc_sm = smoothdata(tc,2,'gaussian',5);
    tc_sm_all.(vars{var_id}) = tc_sm;
    tbl_all.(vars{var_id}) = tbl_*pos_bin;
end

% plot
[M, I] = max(tc,[],2);
[~, I_cell] = sort(I);
tc_n = tc./repmat(M,1,size(tc,2));
imshow(tc_n(I_cell,:))
% 
% 
% tc_sm = smoothdata(tc,2,'gaussian',5);
% figure;plot(tc(6,:));
% hold on;plot(tc_sm(6,:))
% title({'tc smoothed with Gaussian of window 10'})
% R=corrcoef(tc_sm);
% figure;image(R,'CDataMapping','scaled')


save(['/home/kemerelab/Downloads/new_mazes_eeg/bayes_tc_maze',num2str(maze_id),'_500ms.mat'],'tc_sm_all','tbl_all')

%% Decode PBE and line fitting
clear
tbl = [];
tc_sm = [];

for maze_id = 1:2
    load(['/home/kemerelab/Downloads/new_mazes_eeg/bayes_tc_maze',num2str(maze_id),'_500ms.mat'])
    tbl = [tbl;[-1*tbl_all.bwd(end:-1:1,:);tbl_all.fwd]];
    tc_sm = [tc_sm,[tc_sm_all.bwd(:,end:-1:1),tc_sm_all.fwd]];
end

% % two direction interleave2
% [tbl,od] = sort(abs(tbl(:,1)));
% tc_sm = tc_sm(:,od);

pbe_time_bin = 17;
load('select_pbe_maze1_17ms.mat')

yy0 = yy0';
nt = size(yy0,2);
scaler = pbe_time_bin/500; %
tc_sc = tc_sm.*scaler; %tc = tc_sm.*scaler_cell+0.0001;
tc_sc = tc_sc+min(nonzeros(tc_sc))/10;

% tc_scc = tc.*scaler_ratio+0.00001;
% tc_sc = smoothdata(tc_scc,2,'gaussian',10);

loglikelihood = -repmat(sum(tc_sc',2)',nt,1) + yy0'*log(tc_sc);
% [Ls, xinitidx] = max(loglikelihood,[],2);
% xinit_ratio = tbl(xinitidx,1);

save(['bayes_decoded_maze',num2str(maze_id),'pbe_17.mat'],'loglikelihood','tbl')



tbl_plot=abs(tbl([31:60,91:120],1));
seg_id = 42;
matrix = exp(loglikelihood(d(seg_id)+1:d(seg_id+1),:)');
% matrix = exp(loglikelihood(4739:4793,:)'); 
matrix_n = zeros(size(matrix));
for i=1:size(matrix,2)
    matrix_n(:,i) = matrix(:,i)/sum(matrix(:,i));
%     matrix_n(:,i) = (matrix(:,i)-min(matrix(:,i)))/(max(matrix(:,i))-min(matrix(:,i)));
end
matrix_plot = matrix_n([30:-1:1,90:-1:61],:)+matrix_n([31:60,91:120],:);

figure;
image(1:size(matrix,2),1:size(matrix_plot,1),matrix_plot,'CDataMapping','scaled')
hold on
for j=1:size(matrix_plot,1)/30-1
    plot(xlim,ones(1,2)*j*30+0.5,'--','Color',[0.65,0.65,0.65])
end
c = gray;
c = flipud(c);
colormap(c);
set(gca,'YDir','normal')
xlabel('time bin')
yticks([15,45]);
yticklabels({'maze1','maze2'})
cb = colorbar;
cb.Label.String='Probability';
title(['PBE No.',num2str(seg_id)])
cmap1 = parula(30);
scatter(ones(30,1)*.5,1:30,30,cmap1,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.7)
cmap2 = copper(30);
scatter(ones(30,1)*.5,31:60,30,cmap2,'filled','MarkerEdgeColor','None','MarkerFaceAlpha',.7)
l=xlim;
xlim([0.1,l(2)])
caxis([0,1])

