%% A/B compartments

% Filtering
B = abs(medfilt2(mylog2_neg_inf(H)));
figure
subplot(1,2,1)
    imagesc(mylog2_neg_inf( H))
    erez_imagesc
subplot(1,2,2)
    imagesc(B)
    erez_imagesc

abcomp_mat = hic_abcomp(H, 'fiedler', rna_compatible, 'yes', 'no');
arm_split = find(diff(sign(abcomp_mat)))+1;

abcomp_mat1 = hic_abcomp(H(abcomp_mat<0,abcomp_mat<0), 'fiedler', rna_compatible(abcomp_mat<0), 'yes', 'yes');
abcomp_mat2 = hic_abcomp(H(abcomp_mat>0,abcomp_mat>0), 'fiedler', rna_compatible(abcomp_mat>0), 'yes', 'yes');

abcomp_mat = [abcomp_mat1; abcomp_mat2];

% ab_intervals = sort([1; find(diff(sign(abcomp_mat)))+1; arm_split; size(abcomp_mat,1)+1],'ascend');
ab_intervals = sort([1; find(diff(sign(abcomp_mat)))+1; size(abcomp_mat,1)+1],'ascend');


figure('Position', [1272 42 648 1074])
subplot(5,1,1)
    bar(rna_compatible,'k')
    ylabel('Reads')
subplot(5,1,2)
    bar(abcomp_mat,'k')
    ylabel('A/B')
    hold on
    plot([arm_split arm_split],get(gca,'YLim'),'r--')
subplot(5,1,3:5)
    imagesc(mylog2_neg_inf( H))
    erez_imagesc
    hold on
    plot_TADs(ab_intervals,0) % not TADs, actually A/B
      
%% TADs
sig0 = .95;
ms0 = 3;

tad_intervals = TAD_Laplace_Sijia(H,sig0, ms0);

figure('Position', [1272 42 648 1074])
subplot(5,1,1)
    bar(rna_compatible,'k')
    ylabel('Reads')
subplot(5,1,2)
    bar(abcomp_mat,'k')
    ylabel('A/B')
subplot(5,1,3:5)
    imagesc(mylog2_neg_inf( H))
    erez_imagesc
    hold on
    plot_TADs(tad_intervals,0)
