tic
warning('off','all')
% close all
% clear
restoredefaultpath
addpath(genpath('\\172.17.109.24\internal_4DN\tools\matlab\'))
% addpath(genpath('\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Wenlong_phased\Phased_Hi-C_final\'))
% % addpath(genpath('\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Phased_Hi-C_intra_inter_2\'))
% addpath(genpath('\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Preliminary_analysis\RNAseq'))
addpath(genpath(pwd))

binSize = 1E6;
chrSizes = readtable(sprintf('%s.chrom.sizes','hg19'),'filetype','text');
chrStart = [1;cumsum(ceil(chrSizes{:,2}./binSize))+1];
chr1mbLength = diff(chrStart,1);
chrStart = chrStart(1:24);

%% Hi-C Loading

calc_chr = 0;
calc_whole_gen = 1;

hic_base = '\\172.17.109.24\internal_4dn\projects\scHiC_VARI-068_487-HC\processed\all_valid_libraries_cells\';
% sample = {'Sample_487-HC-1_AACCGT\'};
% mkviif = {'Sample_487-HC-1_AACCGT_MKVII_F_3\','Sample_487-HC-1_AACCGT_MKVII_F_7\','Sample_487-HC-1_AACCGT_MKVII_F_9\'};

libraries = struct2table(dir(hic_base));
libraries = natsortfiles(libraries.name(3:end));

cells =[];
for i = 1:size(libraries,1)
    cells_temp = struct2table(dir(['\\172.17.109.24\internal_4dn\projects\scHiC_VARI-068_487-HC\processed\all_valid_libraries_cells\', libraries{i}]));
    cells = [cells; natsortfiles(cells_temp.name(3:end))];
end

H_indiv = [];

if calc_chr == 1;
    for i = 1:size(cells,1)
        H_indiv = padconcatenation_sr(H_indiv, hic2mat('observed','KR',...
                  [hic_base,cells{i},'\aligned\inter_30.hic'],...
                  '2','2','bp',binSize,1,0),3);
    end
end

%% Load all 314 single cells into matlab from .hic files DO ONCE
% H_indiv_1_114 = [];
% if calc_whole_gen == 1
%     for i = 1:114 %size(cells,1)
%         library_temp =  [extractBefore(cells{i},'_MK'),'\'];
% 
%         H_indiv_1_114 = padconcatenation_sr(H_indiv_1_114, ...
%             extractHicGenomewide([hic_base,library_temp,cells{i},'\aligned\inter_30.hic'],...
%             'hg19',1E6,'BP','NONE','observed'),3);
%     end
% end
% save('H_indiv_cells_1_114.mat', 'H_indiv_1_114', '-v7.3')
% 
% H_indiv_115_214 = [];
% if calc_whole_gen == 1
%     for i = 115:214 %size(cells,1)
%         i
%         library_temp =  [extractBefore(cells{i},'_MK'),'\'];
% 
%         H_indiv_115_214 = padconcatenation_sr(H_indiv_115_214, ...
%             extractHicGenomewide([hic_base,library_temp,cells{i},'\aligned\inter_30.hic'],...
%             'hg19',1E6,'BP','NONE','observed'),3);
%     end
% end
% save('H_indiv_cells_115_214.mat', 'H_indiv_115_214', '-v7.3')
% clear H_indiv_115_214
% 
% 
% H_indiv_215_314 = [];
% if calc_whole_gen == 1
%     for i = 215:size(cells,1)
%         i
%         library_temp =  [extractBefore(cells{i},'_MK'),'\'];
% 
%         H_indiv_215_314 = padconcatenation_sr(H_indiv_215_314, ...
%             extractHicGenomewide([hic_base,library_temp,cells{i},'\aligned\inter_30.hic'],...
%             'hg19',1E6,'BP','NONE','observed'),3);
%     end
% end
% save('H_indiv_cells_215_314.mat', 'H_indiv_215_314', '-v7.3')
% clear H_indiv_215_314
% H_indiv(isnan(H_indiv)) = 0;
% 
% H_indiv = cat(3,H_indiv,H_indiv_215_314);



H = sum(H_indiv,3); 
% H_indiv(:,:,1)+H_indiv(:,:,2)+H_indiv(:,:,3);
% figure
%     imagesc(H)
%     erez_imagesc

%% Plot 3 example single cell Hi-C matrices
to_plot = 1;

if to_plot == 1
    figure
    subplot(2,4,1)
        imagesc(H_indiv(:,:,1))
        erez_imagesc
        title('Cell 3')
    subplot(2,4,2)
        imagesc(H_indiv(:,:,2))
        erez_imagesc
        title('Cell 5')
    subplot(2,4,3)
        imagesc(H_indiv(:,:,3))
        erez_imagesc
        title('Cell 6')
    subplot(2,4,4)
        imagesc(H)
        erez_imagesc
        title('Sum')
    subplot(2,4,5)
        imagesc(mylog2_neg_inf(H_indiv(:,:,1)))
        erez_imagesc
        title('log(Cell 3)')
    subplot(2,4,6)
        imagesc(mylog2_neg_inf(H_indiv(:,:,2)))
        erez_imagesc
        title('log(Cell 5)')
    subplot(2,4,7)
        imagesc(mylog2_neg_inf(H_indiv(:,:,3)))
        erez_imagesc
        title('log(Cell 6)')
    subplot(2,4,8)
        imagesc(mylog2_neg_inf(H))
        erez_imagesc
        title('log(Sum)')
end
    
%% RNA-Seq loading
rna_seq_pos = readtable('\\172.17.109.24\internal_4dn\projects\SciHi-C_VARI068\processed\rnaseq\114097_VARI068_ALDHpos_DGE.txt');
rna_seq_neg = readtable('\\172.17.109.24\internal_4dn\projects\SciHi-C_VARI068\processed\rnaseq\114098_VARI068_ALDHneg_DGE.txt');

rna_temp = rna_seq_pos(:,1:2);

% Find gene locations with Biomart file
biomart = readtable('mart_export_ensembl_hg37_info.txt');
all_genes_temp = biomart(find(ismember(biomart.HGNCSymbol,rna_temp.GENE)),:);
[C,IA,IC] = unique(all_genes_temp.HGNCSymbol);
all_genes = all_genes_temp(IA,:);
all_genes = movevars(all_genes,'HGNCSymbol','Before','GeneStableID');
all_genes = movevars(all_genes,'Chromosome_scaffoldName','After','HGNCSymbol');
all_genes = movevars(all_genes,'GeneStart_bp_','After','Chromosome_scaffoldName');
all_genes = movevars(all_genes,'GeneEnd_bp_','After','GeneStart_bp_');
all_genes = all_genes(:,1:4);
all_genes.Properties.VariableNames{'HGNCSymbol'} = 'genename';
all_genes.Properties.VariableNames{'Chromosome_scaffoldName'} = 'chr';
all_genes.Properties.VariableNames{'GeneStart_bp_'} = 'start';
all_genes.Properties.VariableNames{'GeneEnd_bp_'} = 'stop';
all_genes.chr = str2double(all_genes.chr);
all_genes.chr(find(isnan(all_genes.chr))) = 23;

idx = find(ismember(rna_temp.GENE,all_genes.genename));

% sum(strcmp(rna_temp.GENE(idx), all_genes.genename))

all_genes.expr = rna_temp{idx,2};

regions_hic = [1:binSize:(size(H,1))*binSize];
regions_hic = [regions_hic; binSize:binSize:(size(H,1))*binSize];

chr_genes = all_genes(find(ismember(all_genes.chr, 3)),:);
clear rna_compatible
for i = 1:length(regions_hic)
    rna_compatible(i) = length(find(chr_genes.start > regions_hic(1,i) & chr_genes.start < regions_hic(2,i)));
end

%% Plotting of chr 3
figure('Position', [1094 42 826 1074])
subplot(4,1,1)
    bar(rna_compatible,'k')
    ylabel('Reads')
subplot(4,1,2:4)
    imagesc(mylog2_neg_inf(H))
    erez_imagesc
    
%% Trimming Hi-C and RNA-Seq
trim = find(diag(H,0) == 0);

H(trim,:) = [];
H(:,trim) = [];
rna_compatible(trim) = [];


D = diag(sum(H,1));
L = D - H;
L_sym = D^(-1/2)*L*D^(-1/2);

[u, s, v] = svd(L_sym);

%% Basic analysis (Scree, rank 1 matrices)
[evec,eval] = eig(H);
[~,I] = sort(diag(eval),'descend');
evec = evec(:,I);
eval = diag(eval);
eval = eval(I);
figure
plot(eval,'k.')

[U,S,V] = svd(H);
[~,I] = sort(diag(S),'descend');
[U, S, V] = svd(H);
[~,Isin] = sort(diag(S),'descend');
U = U(:,Isin);
S = diag(S);
S = S(Isin);
V = V(:,Isin);

figure
subplot(1,4,1)
imagesc(mylog2_neg_inf(H)), indika_figure_style, erez_imagesc
subplot(1,4,2)
imagesc(U(:,1)*S(1)*V(:,1)'), indika_figure_style, erez_imagesc
subplot(1,4,3)
imagesc(U(:,2)*S(2)*V(:,2)'), indika_figure_style, erez_imagesc
subplot(1,4,4)
imagesc(U(:,3)*S(3)*V(:,3)'), indika_figure_style, erez_imagesc

figure
for i = 1:3
    [evec_indiv{i},eval_indiv{i}] = eig(H_indiv(:,:,i));
    [~,Ieig_indiv{i}] = sort(diag(eval_indiv{i}),'descend');
    evec_indiv{i} = evec_indiv{i}(:,Ieig_indiv{i});
    eval_indiv{i} = diag(eval_indiv{i});
    eval_indiv{i} = eval_indiv{i}(Ieig_indiv{i});
    
    % Scree plots for each cell
    subplot(1,3,i)
    plot(eval_indiv{i},'k.'), axis square
    
    [U_indiv{i}, S_indiv{i}, V_indiv{i}] = svd(H_indiv(:,:,i));
    [~,Isin_indiv{i}] = sort(diag(S_indiv{i}),'descend');
    U_indiv{i} = U_indiv{i}(:,Isin_indiv{i});
    S_indiv{i} = diag(S_indiv{i});
    S_indiv{i} = S_indiv{i}(Isin_indiv{i});
    V_indiv{i} = V_indiv{i}(:,Isin_indiv{i});
end

linkaxesInFigure
ylim_temp = get(gca,'YLim');
set(gca,'YLim',[0, ylim_temp(2)])

%% OLD CODE
% 
% type = {'intra', 'inter'};
% sample = {'Sample_487-HC-1_AACCGT\'};
% mkviif = {'Sample_487-HC-1_AACCGT_MKVII_F_3\'};
% res_name = {'s1mb', 's100kb'};
% norm_type = {'obs','oe','kr','oekr'};
% res = [1E6, 1E5];
%   
%     for n = 1 %1:4
%         n
%         for h = 1 % 1mb 100kb
%             h
%             for chr_select = 1:23 % 1:23
%                 if chr_select == 23
%                     chr_str = 'X'
%                 else
%                     chr_str = string(chr_select)
%                 end
% 
%                 for i = 1 % intra inter
%                    par = [];
%                    for j = 1:2 % mat pat
%                       temp = [];
%                       for k = 1:3 % G1 S G2
%                          ext = [type{i}, '\', sample{j}, '\', mkviif{k}, '\'];
%         %                  ext =  strcat(strcat(itype, iparent), iphase);
%                          hic_path = [hic_base, ext];
%                          if strcmp(norm_type{n},'kr')
%                              if j == 2 && k == 3
%                                 obs_none = juicer2mat(juicer_tools_dump_mat('observed','NONE',...
%                                      [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
%                                  
%                                 badLocs_temp = sum(obs_none)' <= 0 | diag(obs_none) <=0;
%                                 obs_none_no_zero = obs_none;
%                                 obs_none_no_zero(badLocs_temp,:) = [];
%                                 obs_none_no_zero(:,badLocs_temp) = [];
%                                 
%                                 [x,~] = KR_NORM_bnewt(obs_none_no_zero);
%                                 
%                                 obs_none_no_zero_kr = diag(x)*obs_none_no_zero*diag(x);
%                                                                
%                                 obs_none_kr = zeros(size(obs_none));
%                                 obs_none_kr(~badLocs_temp,~badLocs_temp) = obs_none_no_zero_kr;
% 
%                                 obs_none_kr = obs_none_kr*(nanmean(obs_none(:))/nanmean(obs_none_kr(:)));
%                                  
%                                 temp = padconcatenation_sr(temp, obs_none_kr,3);
%                              else
%                                 temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('observed','KR',...
%                                      [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
%                              end
%                          elseif strcmp(norm_type{n},'oekr')
%                              if j == 2 && k == 3
%                                 obs_none = juicer2mat(juicer_tools_dump_mat('observed','NONE',...
%                                      [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
%                                  
% %                                 figure, imagesc(obs_none), title('obs none')
% 
%                                 badLocs_temp = sum(obs_none)' <= 0 | diag(obs_none) <=0;
%                                 obs_none_no_zero = obs_none;
%                                 obs_none_no_zero(badLocs_temp,:) = [];
%                                 obs_none_no_zero(:,badLocs_temp) = [];
% 
% %                                 figure, imagesc(obs_none_no_zero), title('obs none no zero')
%                                 
%                                 [x,~] = KR_NORM_bnewt(obs_none_no_zero);
%                                 
%                                 obs_none_no_zero_kr = diag(x)*obs_none_no_zero*diag(x);
%                                 
% %                                 figure, imagesc(obs_none_no_zero_kr), title('obs none no zero kr')
%                                 
%                                 obs_none_kr = zeros(size(obs_none));
%                                 obs_none_kr(~badLocs_temp,~badLocs_temp) = obs_none_no_zero_kr;
% 
%                                 obs_none_kr = obs_none_kr*(nanmean(obs_none(:))/nanmean(obs_none_kr(:)));
% 
%                                 
% %                                 figure, imagesc(log(obs_none_kr)), title('obs none kr scaled')
%                                 
%                                 mat_g2_kr = juicer2mat(juicer_tools_dump_mat('observed','KR',...
%                                      [hic_base,'intra\maternal\g2\inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
%                                 mat_g2_oe_kr = juicer2mat(juicer_tools_dump_mat('OE','KR',...
%                                      [hic_base,'intra\maternal\g2\inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
%                                  
%                                 expected = mat_g2_kr./mat_g2_oe_kr;
%                                 expec_vec = [];
%                                 for d = 1:size(expected,1)
%                                    expec_vec(d) = nanmean(diag(expected,d-1)); 
%                                    if isnan(expec_vec(d))
%                                        expec_vec(d) = expec_vec(d-1);
%                                    end
%                                 end
%                                 
%                                 while length(expec_vec) < size(obs_none_kr,1)
%                                     expec_vec = [expec_vec expec_vec(end)];
%                                 end
%                                 
%                                 expected_no_nan = toeplitz(expec_vec);
%                                 obs_oe_kr = obs_none_kr./expected_no_nan;
%                                 
% %                                 figure, imagesc(obs_oe_kr), title('obs oe kr')
% 
%                                 temp = padconcatenation_sr(temp, obs_oe_kr,3);
%                              else
%                                 temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('OE','KR',...
%                                      [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
%                              end
%                          elseif strcmp(norm_type{n},'oe')
%                              temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('OE','NONE',...
%                                      [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
%                          else
%                              temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('observed','NONE',...
%                                      [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
%                          end
%                          
%                       end
% 
%                       % Dimension mismatch
% %                       if j == 2
%                           if size(temp,1) ~= size(par,1)
%                               par = padconcatenation_sr(temp, par,3);
%                           else
%                               par(:,:,(1:3)+3*(j-1)) = temp;
%                           end
% %                       else
% %                           par(:,:,(1:3)+3*(j-1)) = temp;
% %                       end
% 
%                       clear temp
%                    end
%                    if strcmp(norm_type{n},'kr')
%                        H.(res_name{h}).(type{i}).kr{chr_select,1} = par;
%                        [H.(res_name{h}).(type{i}).kr_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
%                    elseif strcmp(norm_type{n},'oekr')
%                        H.(res_name{h}).(type{i}).oekr{chr_select,1} = par;
%                        [H.(res_name{h}).(type{i}).oekr_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
%                    elseif strcmp(norm_type{n},'oe')
%                        H.(res_name{h}).(type{i}).oe{chr_select,1} = par;
%                        [H.(res_name{h}).(type{i}).oe_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
%                    else
%                        H.(res_name{h}).(type{i}).obs{chr_select,1} = par;
%                        [H.(res_name{h}).(type{i}).obs_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
%                    end
%                    
% %                    figure
% %                    subplot(1,2,1)
% %                    imagesc(H.s1mb.intra.full{1}(:,:,3))
% %                    axis square
% %                    subplot(1,2,2)
% %                    imagesc(H.s1mb.intra.full{1}(:,:,6))
% %                    axis square
%                    clear par
%                 end
%             end
%         end
%     end
%     for i = 1:23
%         chr_bin_lengths_via_hic(i) = size(H.s100kb.intra.oekr{i},1);
%         chr_bin_lengths_via_hic_1mb(i) = size(H.s1mb.intra.oekr{i},1);
%     end
% end
% 
% % MATLAB numerical error fix for Pat G2
% for chr = 1:23
%     H.s1mb.intra.oekr{chr}(:,:,6) = (H.s1mb.intra.oekr{chr}(:,:,6)+H.s1mb.intra.oekr{chr}(:,:,6)')/2;
%     H.s1mb.intra.oekr_trim_10pct{chr}(:,:,6) = (H.s1mb.intra.oekr_trim_10pct{chr}(:,:,6)+H.s1mb.intra.oekr_trim_10pct{chr}(:,:,6)')/2;
% end
% 
% disp(['Main Hi-C complete: ', num2str(toc), ' seconds'])
