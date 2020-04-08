%% Load data and prep
addpath(genpath('\\172.17.109.24\internal_4DN\tools\matlab\scott_stephen_SF_functions'))

% Only run this chunk if necessary (loading in .mat files)
if 1 == 0
    load('H_indiv_cells_1_114.mat')
    load('H_comb_all.mat')
    load('chrStart.mat')
    load('H_bulk_R1.mat')
    load('R_indiv_cells.mat')
end

for i = 1:114
    tot(i) = sum(sum(H_indiv_cells_1_114(:,:,i)));
end

% Cell 4 has ~87,000 contacts
h = H_indiv_cells_1_114(1:chrStart(24)-1,1:chrStart(24)-1,4);

% Set up matrices to be trimmed
h_trim = h;
H_comb_all_trim = H_comb_all(1:chrStart(24)-1,1:chrStart(24)-1);
H_bulk_R1_trim = H_bulk_R1;

% Gather all loci with a zero on the diagonal
trimLocs = ismember(diag(h), 0);

% Gather all loci that would be trimmed in bulk due to few contacts
[~, bad_locs_bulk] = hic_trim(H_bulk_R1,1,.1);

% Find union of both trimming location variables 
% (using OR operator on logical vectors)
trimLocsComb = trimLocs | bad_locs_bulk';

h_trim(:,trimLocsComb) =  [];
h_trim(trimLocsComb,:) =  [];
H_comb_all_trim(:,trimLocsComb) =  [];
H_comb_all_trim(trimLocsComb,:) =  [];
H_bulk_R1_trim(:,trimLocsComb) =  [];
H_bulk_R1_trim(trimLocsComb,:) =  [];

figure
imagesc(h_trim)

clear tot

%% Python imputation call

% Save into txt file
fid = fopen('Mymatrix.txt','wt');
A = h_trim;
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end

commandStr = 'py "C:\Users\lindsly\Google Drive\UMICH\Research\Single Cell Hi-C\Matlab\HiCluster\HiCluster\schicluster\single_cell.py" Mymatrix.txt';
[status, commandOut] = system(commandStr);
if status==0
    disp('Imputation complete');
else
    disp('Error: check commandOut');
end
% Return 'Q' imputed matrix
load('imputed_matrix.mat')

% Q not symmetric
issymmetric(Q);
% Use (Q*Q')/2 to make this matrix symmetric as discussed with Can
Q_sym = (Q+Q')/2;
issymmetric(Q_sym);

fclose('all');
clear status commandOut commandStr A fid ii

%% Test plots
% figure
% imagesc(mylog2_neg_inf(1000*(Q)))
% 
% figure
% imagesc(log(Q))
% 
% figure
% imagesc((Q*1000))
% erez_imagesc
% 
% figure
% imagesc(Q*1000-diag(diag(Q*1000)))
% erez_imagesc

%% Isolate one chromosome
chrStartTrim = chrStart;

trimLocIdx = sort(find(trimLocsComb),'descend');
for iTrimLocs = 1:length(trimLocIdx)
    chrStartTrim(chrStartTrim >= trimLocIdx(iTrimLocs)) =...
        chrStartTrim(chrStartTrim >= trimLocIdx(iTrimLocs))-1;
end

h_trim_chr14 = h_trim(chrStartTrim(14)+1:chrStartTrim(15),chrStartTrim(14)+1:chrStartTrim(15));

h_chr14 = Q_sym(chrStartTrim(14)+1:chrStartTrim(15),chrStartTrim(14)+1:chrStartTrim(15));

H_comb_all_chr14 = H_comb_all_trim(chrStartTrim(14)+1:chrStartTrim(15),chrStartTrim(14)+1:chrStartTrim(15));

H_chr14 = H_bulk_R1_trim(chrStartTrim(14)+1:chrStartTrim(15),chrStartTrim(14)+1:chrStartTrim(15));

rna_neg_trim = rna_neg_compatible(~trimLocsComb,:);

% Cell 2 in ALDH- scRNA-seq has more counts than 1
r_chr14 = rna_neg_trim(chrStartTrim(14)+1:chrStartTrim(15),2);

% figure
% imagesc(mylog2_neg_inf(100*h_chr14))

clear iTrimLocs chrEnd startPos

%% Fiedler number/vector

% Single cell (with imputation)
D = diag(sum(h_chr14,1));
L = D - h_chr14;
L_sym = D^(-.5)*L*D^(-.5);
[f_vec, f_val] = eigs(L_sym,2,'sm');

f_vec = f_vec(:,2);
f_val = f_val(2,2);

figure('Position', [661 135 560 825])
subplot(6,1,1)
bar(r_chr14)
title('RNA-seq')
subplot(6,1,2)
bar(f_vec)
title('Fiedler Vector')
subplot(6,1,3:6)
imagesc(mylog2_neg_inf(100*h_chr14))
axis square
title('Imputed scHi-C')

% Sum of single cells (no imputation)
D = diag(sum(H_comb_all_chr14,1));
L = D - H_comb_all_chr14;
L_sym = D^(-.5)*L*D^(-.5);
[f_vec, f_val] = eigs(L_sym,2,'sm');
f_vec = f_vec(:,2);
f_val = f_val(2,2);

figure('Position', [661 135 560 825])
% subplot(6,1,1)
% bar(r_chr14)
% title('RNA-seq')
subplot(6,1,2)
bar(-f_vec)
title('Fiedler Vector')
subplot(6,1,3:6)
imagesc(mylog2_neg_inf(H_comb_all_chr14))
axis square
title('Sum of scHi-C')

% Bulk 
D = diag(sum(H_chr14,1));
L = D - H_chr14;
L_sym = D^(-.5)*L*D^(-.5);
[f_vec, f_val] = eigs(L_sym,2,'sm');
f_vec = f_vec(:,2);
f_val = f_val(2,2);

figure('Position', [661 135 560 825])
% subplot(6,1,1)
% bar(r_chr14)
% title('RNA-seq')
subplot(6,1,2)
bar(-f_vec)
title('Fiedler Vector')
subplot(6,1,3:6)
imagesc(mylog2_neg_inf(H_chr14))
axis square
title('Bulk Hi-C')

% create temp variables for sharing
RNA_SC_chr14 = r_chr14;
HiC_SC_raw_chr14 = h_trim_chr14;
HiC_SumSC_raw_chr14 = H_comb_all_chr14;
HiC_SC_impute_chr14 = h_chr14;
HiC_Bulk_raw_chr14 = H_chr14;


%% OLD CODE
% Get loci that are missing for all cells (mostly centromeres)
% load('bad_locs_all.mat') % Calculated using all 314 cells
% if 1 == 0
% %     h_comb_all = sum(H_indiv,3);
% 
%     figure
%     imagesc(mylog2_neg_inf(H_comb_all))
%     erez_imagesc
% 
%     [h_comb_all_trim, bad_locs_all] = hic_trim(H_comb_all);
% %     [h_comb_test, bad_locs_test] = hic_trim(H_indiv_cells_1_114);
% 
%     figure
%     imagesc(mylog2_neg_inf(h_comb_all_trim))
%     erez_imagesc
% end
% 
% H_trim = H_indiv_cells_1_114(~bad_locs_all,~bad_locs_all,:);

