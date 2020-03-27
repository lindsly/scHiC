addpath(genpath('\\172.17.109.24\internal_4DN\tools\matlab\scott_stephen_SF_functions'))

if (exist('H_indiv_cells_1_114','var')) == 0
    load('H_indiv_cells_1_114.mat')
end

h = H_indiv_cells_1_114(:,:,1);

% Get loci that are missing for all cells (mostly centromeres)
if 1 == 0
    h_comb_all = sum(H_indiv,3);

    figure
    imagesc(mylog2_neg_inf(h_comb_all))
    erez_imagesc

    [h_comb_all_trim, bad_locs_all] = hic_trim(h_comb_all);

    figure
    imagesc(mylog2_neg_inf(h_comb_all_trim))
    erez_imagesc
end

H_trim = H_indiv_cells_1_114(~bad_locs_all,~bad_locs_all,:);

H_trim_reg = .95*H_trim + .05*ones(size(H_trim,1),size(H_trim,2),size(H_trim,3));

%%
% h = H_trim_reg(:,:,1);
% h = H_comb_all_trim;
h = H_indiv_cells_1_114(1:chrStart(2),1:chrStart(2),2);
[h, bad_loc_chr1] = hic_trim(h,1);
% h = .95*h + .05*ones(size(h,1),size(h,2),size(h,3));

rna_chr1 = rna_compatible(1:chrStart(2));
rna_chr1 = rna_chr1(~bad_loc_chr1);
D = diag(sum(h,1));
L = D - h;

L_hat = D^(-.5)*L*D^(-.5);

[evec, eval] = eig(L);%,5,'smallestreal');

evec_test = evec(:,2);

% evec_test(evec_test>0) = 1;
% evec_test(evec_test<0) = -1;

figure('Position', [825 42 659 1074])
subplot(5,1,1)
    bar(rna_chr1)
subplot(5,1,2)
    bar(evec_test)
subplot(5,1,3:5)
    imagesc(mylog2_neg_inf(h))
    axis square


h2 = H_comb_all(1:chrStart(2),1:chrStart(2));
h2 = ToepNorm2(h2(~bad_loc_chr1,~bad_loc_chr1,:));

D = diag(sum(h2,1));
L = D - h2;

L_hat = D^(-.5)*L*D^(-.5);

[evec, eval] = eig(L);%,5,'smallestreal');

evec_test = evec(:,2);

evec_test(evec_test>0) = 1;
evec_test(evec_test<0) = -1;

figure('Position', [667 42 817 1074])
subplot(4,1,1:3)
    imagesc(mylog2_neg_inf(h2))
    axis square
subplot(4,1,4)
    bar(evec_test)


% figure
% hist(h(h>0),50)