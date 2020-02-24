% Convenient access to lines I use a lot

opengl software
opengl hardware


w = warning ('on','all');

w = warning ('off','all');

clearvars -except H Bru_AS_NonAS_FPKMtable_Avg RNA_AS_NonAS_FPKMtable_Avg R alt_parent%variables

clearvars -except H Bru_AS_NonAS_FPKMtable_Avg RNA_AS_NonAS_FPKMtable_Avg de_genes_in_hic AS_and_Cycle_DESpec


hicpath = '\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Wenlong_phased\Phased_Hi-C_Files\';

restoredefaultpath
addpath(genpath('\\172.17.109.24\internal_4DN\tools\matlab\scott_stephen_SF_functions'))

if (exist('T','var')) == 0
    load('Master_Table.mat')
end

if (exist('H','var')) == 0
    load('AS_HiC_RNASeq_BruSeq_intra_inter.mat', 'H')
end

if (exist('AS_and_Cycle_DESpec','var')) == 0
    load('AS_and_Cycle_DESpec.mat')
end

if (exist('R','var')) == 0
    load('AS_HiC_RNASeq_BruSeq_intra_inter.mat', 'R')
end

if (exist('B','var')) == 0
    load('AS_HiC_RNASeq_BruSeq_intra_inter.mat', 'B')
end
% 
% if (exist('RNA_AS_NonAS_FPKMtable_Avg','var')) == 0
% 	('RNA_BRU_All_Gene_Phased_mutRatioAS_5050NONAS_FIXED.mat')
% end

load('HumanTfB.mat')


if (exist('RNA_AS_NonAS_FPKMtable_Avg','var')) == 0
	load('AS_NonAS_RNASeq_BruSeq_FPKM.mat')
end

sample_names = {'mat_g1' 'mat_s' 'mat_g2' 'pat_g1' 'pat_s' 'pat_g2'};
nice_sample_names = {'Maternal G1' 'Maternal S' 'Maternal G2' 'Paternal G1' 'Paternal S' 'Paternal G2'};
inform_names =  {'MG1vsPG1', 'MSvsPS', 'MG2vsPG2', 'MG1vsMS', 'MSvsMG2', 'MG1vsMG2', 'PG1vsPS', 'PSvsPG2', 'PG1vsPG2'};
nine_set = [1 4; 2 5; 3 6; 1 2; 2 3; 1 3; 4 5; 5 6; 4 6];


chr_nums = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X'}';

nice_chr_nums = {'Chromosome 1' 'Chromosome 2' 'Chromosome 3' 'Chromosome 4' 'Chromosome 5'...
    'Chromosome 6' 'Chromosome 7' 'Chromosome 8' 'Chromosome 9' 'Chromosome 10' 'Chromosome 11'...
    'Chromosome 12' 'Chromosome 13' 'Chromosome 14' 'Chromosome 15' 'Chromosome 16' 'Chromosome 17' ...
    'Chromosome 18' 'Chromosome 19' 'Chromosome 20' 'Chromosome 21' 'Chromosome 22' 'Chromosome X'}';

color_set_chrs = [rgb('Black');rgb('DarkSlateGray');rgb('Red');rgb('Salmon');rgb('Crimson');...
    rgb('DarkRed');rgb('Pink');rgb('HotPink');rgb('DeepPink');rgb('Orange');...
    rgb('Tomato');rgb('OrangeRed');rgb('Gold');rgb('DarkKhaki');rgb('Tan');rgb('RosyBrown');...
    rgb('Chocolate');rgb('Green');rgb('LightGreen');rgb('YellowGreen');rgb('Olive');rgb('Teal');...
    rgb('Blue');rgb('Cyan');rgb('DeepSkyBlue');rgb('Purple');rgb('BlueViolet');rgb('Indigo')];


writetable(HiC_RNA_Change,['HiC_RNA_Change_90pctsparse_',phase,'.csv'])

save('filename.mat','var1','var2','var3')

temp = hic(:,:,1);
imagesc(temp, [prctile(temp(:),1), prctile(temp(:),99)])

% LATEX FONT
set(0,'defaulttextinterpreter','latex')
set(0,'defaulttextinterpreter','none')
set(0,'defaultAxesFontName','Arial')
set(0,'defaultTextFontName','Arial')


%    ['\Large{\textbf{','text','}}']

figure('units','normalized','outerposition',[0 0 1 1])

Q(isnan(Q))=0;

linkaxes([ax1,ax2,ax3],'xy')

% WHEN USING KMEANS, GET RID OF SIJIA'S VERSOIN
rmpath('\\172.17.109.24\internal_4dn\tools\matlab\functions\scott\sijia_centrality_codes\tensorlab_2016-03-28')

