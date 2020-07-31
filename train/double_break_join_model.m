%global variables%
global firstbin
global lastbin
global intra_chr_a

nume=length(events00);

[sij1dy_l,opt1dy0_l,fval0_l,lij,intra_chr_l] = Lijoptim_to_length(events00,chsize,bins,CHR,R,mfull{1}+mfull{2}+mfull{3}+mfull{4},sij1dx,len_factor);

Rn=R;
Rn=2*Rn/sum(Rn);
p_mult = kron(Rn,Rn').*lij;

%change ic ratio here% 
ratio = (sum(sum((p_mult(firstbin:lastbin,lastbin+1:end)))) + sum(sum(p_mult(firstbin:lastbin,1:firstbin-1))))/sum(sum(p_mult));


%want ratio to be 1-  intra_chr_a
p_mult(firstbin:lastbin,lastbin+1:end) = ((1 - intra_chr_a)/ratio)*p_mult(firstbin:lastbin,lastbin+1:end);
p_mult(firstbin:lastbin,1:firstbin-1) = ((1 - intra_chr_a)/ratio)*p_mult(firstbin:lastbin, 1:firstbin-1);


p_mult = 2*p_mult ./ sum(sum(p_mult));

%recalculate ratio
(sum(sum((p_mult(firstbin:lastbin,lastbin+1:end)))) + sum(sum(p_mult(firstbin:lastbin,1:firstbin-1))))/sum(sum(p_mult));

%[qFDR_MM, pa_MM, pval_tophits_MM, mfull_pval_MM] = PValMCV(mfull{1}+mfull{2}+mfull{3}+mfull{4}, p_mult);

%[hitstable_MM,hitstableMM_lookup] = HitsTableCV(mfull_pval_MM,pa_MM, pval_tophits_MM, bins_event_tble, qFDR_MM, events, refgene_tble);
%TbyGene_MM=TophitsByGenes(hitstable_MM,hitstableMM_lookup,1e4,bins,refgene,refgene_tble,UTumor,CosmicCencus,uFusionTable,bins_annot);

