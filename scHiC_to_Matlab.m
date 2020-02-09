tic
warning('off','all')
% close all
% clear
restoredefaultpath
addpath(genpath('\\172.17.109.24\internal_4DN\tools\matlab\'))
addpath(genpath('\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Wenlong_phased\Phased_Hi-C_final\'))
% addpath(genpath('\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Phased_Hi-C_intra_inter_2\'))
addpath(genpath('\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Preliminary_analysis\RNAseq'))

binSize = 1E6;
chrSizes = readtable(sprintf('%s.chrom.sizes','hg19'),'filetype','text');
chrStart = [1;cumsum(ceil(chrSizes{:,2}./binSize))+1];
chr1mbLength = diff(chrStart,1);
chrStart = chrStart(1:24);

calc_chr = 1;
calc_whole_gen = 1;


%% Hi-C Loading
hic_base = '\\172.17.109.24\internal_4dn\projects\W50K_cell_cycle\Wenlong_phased\Phased_Hi-C_final\';
% hic_base = '\\172.17.109.24\internal_4dn\projects\W50K_cell_cycle\Wenlong_phased\Phased_Hi-C_intra_inter_2\';

type = {'intra', 'inter'};
parent = {'maternal', 'paternal'};
phase = {'g1', 's', 'g2'};
res_name = {'s1mb', 's100kb'};
norm_type = {'obs','oe','kr','oekr'};
res = [1E6, 1E5];
  
if calc_chr == 1
%     hicpath = '\\172.17.109.24\internal_4DN\projects\W50K_cell_cycle\Wenlong_phased\Phased_Hi-C_Files\';
    % Hi-C Matrix is MatG1, MatS, MatG2, PatG1, PatS, PatG2
    for n = 1 %1:4
        n
        for h = 1 % 1mb 100kb
            h
            for chr_select = 1:23 % 1:23
                if chr_select == 23
                    chr_str = 'X'
                else
                    chr_str = string(chr_select)
                end

                for i = 1 % intra inter
                   par = [];
                   for j = 1:2 % mat pat
                      temp = [];
                      for k = 1:3 % G1 S G2
                         ext = [type{i}, '\', parent{j}, '\', phase{k}, '\'];
        %                  ext =  strcat(strcat(itype, iparent), iphase);
                         hic_path = [hic_base, ext];
                         if strcmp(norm_type{n},'kr')
                             if j == 2 && k == 3
                                obs_none = juicer2mat(juicer_tools_dump_mat('observed','NONE',...
                                     [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
                                 
                                badLocs_temp = sum(obs_none)' <= 0 | diag(obs_none) <=0;
                                obs_none_no_zero = obs_none;
                                obs_none_no_zero(badLocs_temp,:) = [];
                                obs_none_no_zero(:,badLocs_temp) = [];
                                
                                [x,~] = KR_NORM_bnewt(obs_none_no_zero);
                                
                                obs_none_no_zero_kr = diag(x)*obs_none_no_zero*diag(x);
                                                               
                                obs_none_kr = zeros(size(obs_none));
                                obs_none_kr(~badLocs_temp,~badLocs_temp) = obs_none_no_zero_kr;

                                obs_none_kr = obs_none_kr*(nanmean(obs_none(:))/nanmean(obs_none_kr(:)));
                                 
                                temp = padconcatenation_sr(temp, obs_none_kr,3);
                             else
                                temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('observed','KR',...
                                     [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
                             end
                         elseif strcmp(norm_type{n},'oekr')
                             if j == 2 && k == 3
                                obs_none = juicer2mat(juicer_tools_dump_mat('observed','NONE',...
                                     [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
                                 
%                                 figure, imagesc(obs_none), title('obs none')

                                badLocs_temp = sum(obs_none)' <= 0 | diag(obs_none) <=0;
                                obs_none_no_zero = obs_none;
                                obs_none_no_zero(badLocs_temp,:) = [];
                                obs_none_no_zero(:,badLocs_temp) = [];

%                                 figure, imagesc(obs_none_no_zero), title('obs none no zero')
                                
                                [x,~] = KR_NORM_bnewt(obs_none_no_zero);
                                
                                obs_none_no_zero_kr = diag(x)*obs_none_no_zero*diag(x);
                                
%                                 figure, imagesc(obs_none_no_zero_kr), title('obs none no zero kr')
                                
                                obs_none_kr = zeros(size(obs_none));
                                obs_none_kr(~badLocs_temp,~badLocs_temp) = obs_none_no_zero_kr;

                                obs_none_kr = obs_none_kr*(nanmean(obs_none(:))/nanmean(obs_none_kr(:)));

                                
%                                 figure, imagesc(log(obs_none_kr)), title('obs none kr scaled')
                                
                                mat_g2_kr = juicer2mat(juicer_tools_dump_mat('observed','KR',...
                                     [hic_base,'intra\maternal\g2\inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
                                mat_g2_oe_kr = juicer2mat(juicer_tools_dump_mat('OE','KR',...
                                     [hic_base,'intra\maternal\g2\inter_30.hic'],chr_str,chr_str,'BP',res(h)),1);
                                 
                                expected = mat_g2_kr./mat_g2_oe_kr;
                                expec_vec = [];
                                for d = 1:size(expected,1)
                                   expec_vec(d) = nanmean(diag(expected,d-1)); 
                                   if isnan(expec_vec(d))
                                       expec_vec(d) = expec_vec(d-1);
                                   end
                                end
                                
                                while length(expec_vec) < size(obs_none_kr,1)
                                    expec_vec = [expec_vec expec_vec(end)];
                                end
                                
                                expected_no_nan = toeplitz(expec_vec);
                                obs_oe_kr = obs_none_kr./expected_no_nan;
                                
%                                 figure, imagesc(obs_oe_kr), title('obs oe kr')

                                temp = padconcatenation_sr(temp, obs_oe_kr,3);
                             else
                                temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('OE','KR',...
                                     [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
                             end
                         elseif strcmp(norm_type{n},'oe')
                             temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('OE','NONE',...
                                     [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
                         else
                             temp = padconcatenation_sr(temp, juicer2mat(juicer_tools_dump_mat('observed','NONE',...
                                     [hic_path,'inter_30.hic'],chr_str,chr_str,'BP',res(h)),1),3);
                         end
                         
                      end

                      % Dimension mismatch
%                       if j == 2
                          if size(temp,1) ~= size(par,1)
                              par = padconcatenation_sr(temp, par,3);
                          else
                              par(:,:,(1:3)+3*(j-1)) = temp;
                          end
%                       else
%                           par(:,:,(1:3)+3*(j-1)) = temp;
%                       end

                      clear temp
                   end
                   if strcmp(norm_type{n},'kr')
                       H.(res_name{h}).(type{i}).kr{chr_select,1} = par;
                       [H.(res_name{h}).(type{i}).kr_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
                   elseif strcmp(norm_type{n},'oekr')
                       H.(res_name{h}).(type{i}).oekr{chr_select,1} = par;
                       [H.(res_name{h}).(type{i}).oekr_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
                   elseif strcmp(norm_type{n},'oe')
                       H.(res_name{h}).(type{i}).oe{chr_select,1} = par;
                       [H.(res_name{h}).(type{i}).oe_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
                   else
                       H.(res_name{h}).(type{i}).obs{chr_select,1} = par;
                       [H.(res_name{h}).(type{i}).obs_trim_10pct{chr_select,1}, bad_locs.(res_name{h}){chr_select,1}] = hic_trim(par,1,.1);
                   end
                   
%                    figure
%                    subplot(1,2,1)
%                    imagesc(H.s1mb.intra.full{1}(:,:,3))
%                    axis square
%                    subplot(1,2,2)
%                    imagesc(H.s1mb.intra.full{1}(:,:,6))
%                    axis square
                   clear par
                end
            end
        end
    end
    for i = 1:23
        chr_bin_lengths_via_hic(i) = size(H.s100kb.intra.oekr{i},1);
        chr_bin_lengths_via_hic_1mb(i) = size(H.s1mb.intra.oekr{i},1);
    end
end

% MATLAB numerical error fix for Pat G2
for chr = 1:23
    H.s1mb.intra.oekr{chr}(:,:,6) = (H.s1mb.intra.oekr{chr}(:,:,6)+H.s1mb.intra.oekr{chr}(:,:,6)')/2;
    H.s1mb.intra.oekr_trim_10pct{chr}(:,:,6) = (H.s1mb.intra.oekr_trim_10pct{chr}(:,:,6)+H.s1mb.intra.oekr_trim_10pct{chr}(:,:,6)')/2;
end

disp(['Main Hi-C complete: ', num2str(toc), ' seconds'])
