figure
imagesc(H)
erez_imagesc

figure
imagesc(mylog2_neg_inf( H))
erez_imagesc


rna_seq_pos = readtable('\\172.17.109.24\internal_4dn\projects\SciHi-C_VARI068\processed\rnaseq\114097_VARI068_ALDHpos_DGE.txt');
rna_seq_neg = readtable('\\172.17.109.24\internal_4dn\projects\SciHi-C_VARI068\processed\rnaseq\114098_VARI068_ALDHneg_DGE.txt');

abcomp_mat = hic_abcomp(hic(:,:,1), 'fiedler', R.s100kb.FPKM_avg_trim_10pct{chr}(:,1), 'yes', 'no')

sum(isnan(H),'all')
sum(isinf(H),'all')
abcomp_mat = hic_abcomp(H,'fiedler');

figure
bar(rna_seq.TAAATATCCCGT)


for i = 1:length(regions_HiC_chr15)
    ctcf_compatible(i) = length(find(CTCF_chr15 > regions_HiC_chr15(1,i) & CTCF_chr15 < regions_HiC_chr15(2,i)));
end
