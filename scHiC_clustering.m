% %% Basic analysis (Scree, rank 1 matrices)
% [evec,eval] = eig(H);
% [~,I] = sort(diag(eval),'descend');
% evec = evec(:,I);
% eval = diag(eval);
% eval = eval(I);
% figure
% plot(eval,'k.')
% 
% [U,S,V] = svd(H);
% [~,I] = sort(diag(S),'descend');
% [U, S, V] = svd(H);
% [~,Isin] = sort(diag(S),'descend');
% U = U(:,Isin);
% S = diag(S);
% S = S(Isin);
% V = V(:,Isin);
% 
% figure
% subplot(1,4,1)
% imagesc(mylog2_neg_inf(H)), indika_figure_style, erez_imagesc
% subplot(1,4,2)
% imagesc(U(:,1)*S(1)*V(:,1)'), indika_figure_style, erez_imagesc
% subplot(1,4,3)
% imagesc(U(:,2)*S(2)*V(:,2)'), indika_figure_style, erez_imagesc
% subplot(1,4,4)
% imagesc(U(:,3)*S(3)*V(:,3)'), indika_figure_style, erez_imagesc
% 
% figure
% for i = 1:3
%     [evec_indiv{i},eval_indiv{i}] = eig(H_indiv(:,:,i));
%     [~,Ieig_indiv{i}] = sort(diag(eval_indiv{i}),'descend');
%     evec_indiv{i} = evec_indiv{i}(:,Ieig_indiv{i});
%     eval_indiv{i} = diag(eval_indiv{i});
%     eval_indiv{i} = eval_indiv{i}(Ieig_indiv{i});
%     
%     % Scree plots for each cell
%     subplot(1,3,i)
%     plot(eval_indiv{i},'k.'), axis square
%     
%     [U_indiv{i}, S_indiv{i}, V_indiv{i}] = svd(H_indiv(:,:,i));
%     [~,Isin_indiv{i}] = sort(diag(S_indiv{i}),'descend');
%     U_indiv{i} = U_indiv{i}(:,Isin_indiv{i});
%     S_indiv{i} = diag(S_indiv{i});
%     S_indiv{i} = S_indiv{i}(Isin_indiv{i});
%     V_indiv{i} = V_indiv{i}(:,Isin_indiv{i});
% end
% 
% linkaxesInFigure
% ylim_temp = get(gca,'YLim');
% set(gca,'YLim',[0, ylim_temp(2)])

%% Vectorizing individual cell Hi-C matrices and trying to cluster

% all_H_vec = zeros(4850055,16);
% for i = 1:16
%     At = H_indiv(:,:,i).';
%     m  = (1:size(At,1)).' >= (1:size(At,2));
%     all_H_vec(:,i)  = At(m);
% end

figure
imagesc(mylog2_neg_inf(H_indiv(:,:,1)))

temp2 = [];
for i = 1:50
    temp2 = [temp2; diag(H_indiv(:,:,1),i)];
end
    
temp = zeros(size(temp2,1),16);
for j = 2:16
    temp2 = [];
    for i = 1:50
        temp2 = [temp2; diag(H_indiv(:,:,j),i)];
    end
    temp(:,j) = temp2;
end

all_H_vec = temp(find(sum(temp,2) > 0),:);

[u, s, v] = svd(all_H_vec);


%% UMAP of vectorized scHi-C
tic
[reduction, umap, clusterIdentifiers] = run_umap(all_H_vec');
toc

figure
scatter(reduction(:,1), reduction(:,2))%,10,cmap(labels+1,:),'filled')
box on
xlabel('UMAP X')
ylabel('UMAP Y')
title('UMAP of scHi-C')
set(gca,'LineWidth',1.5)

%% tSNE of vectorized scHi-C
restoredefaultpath
tic
ydata = tsne(all_H_vec')%,[],2,size(all_H_vec,1));
toc

figure
scatter(ydata(:,1), ydata(:,2))%,10,cmap(labels+1,:),'filled')
box on
xlabel('Component 1')
ylabel('Component 2')
title('t-SNE of Hi-C')
set(gca,'LineWidth',1.5)

%% PCA of vectorized scHi-C
tic
[mappedX, mapping] = pca(all_H_vec);
toc

figure
scatter(mappedX(:,1), mappedX(:,2),10)%,cmap(labels+1,:),'filled')
box on
xlabel('PC1')
ylabel('PC2')
title('PCA of Hi-C')
set(gca,'LineWidth',1.5)
