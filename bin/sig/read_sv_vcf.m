% note that sometimes the last field is shifted with respect to the header, needed to be corrected manually
FH=readtable('/xchip/beroukhimlab/ofer/DIPG/Unique_HGG_non_hyper_06_29_2017.csv','Delimiter',',');

clear SVTable
vcf_list=FH.snowman_somatic_vcf;
SVTable=table();
for c1=1:length(vcf_list)
    if ~isempty(vcf_list{c1}),
	 disp(['parsing ' vcf_list{c1}]);
         SVTable=read_vcf(SVTable,vcf_list(c1));
    end
end

writetable(SVTable,'uniuq_HGG_nonhyper_SVTable_05_29_2017_dup.csv');

[u_svid,ia_svid,ic_svid] = unique(SVTable.sv_id,'stable');
SVTable_dd=SVTable(ia_svid,:);
writetable(SVTable_dd,'uniuq_HGG_nonhyper_SVTable_05_29_2017.csv');

