nume=length(events00);

[sij1dy_l,opt1dy0_l,fval0_l,lij,intra_chr_l] = Lijoptim_to_length(events00,chsize,bins,CHR,R,mfull{1}+mfull{2}+mfull{3}+mfull{4},sij1dx,len_factor);

Rn=R;
Rn=2*Rn/sum(Rn);
p_mult = kron(Rn,Rn').*lij;
p_mult = 2*p_mult ./ sum(sum(p_mult));


%[qFDR_MM, pa_MM, pval_tophits_MM, mfull_pval_MM] = PValMCV(mfull{1}+mfull{2}+mfull{3}+mfull{4}, p_mult);

%[hitstable_MM,hitstableMM_lookup] = HitsTableCV(mfull_pval_MM,pa_MM, pval_tophits_MM, bins_event_tble, qFDR_MM, events, refgene_tble);
%TbyGene_MM=TophitsByGenes(hitstable_MM,hitstableMM_lookup,1e4,bins,refgene,refgene_tble,UTumor,CosmicCencus,uFusionTable,bins_annot);

