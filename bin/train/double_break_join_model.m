%global variables%
global firstbin
global lastbin
global intra_chr_a



nume=length(events00);

%[sij1dy_l,opt1dy0_l,fval0_l,lij,intra_chr_l] = Lijoptim_to_length_copy(events00,chsize,bins,CHR,R,mfull{1}+mfull{2}+mfull{3}+mfull{4},sij1dx,len_factor);
[sij1dy_l,opt1dy0_l,fval0_l,lij,intra_chr_l] = Lijoptim_to_length_test(events00,chsize,bins,CHR,R,mfull{1}+mfull{2}+mfull{3}+mfull{4},sij1dx);
%[sij1dy_l,opt1dy0_l,fval0_l,lij,intra_chr_l] = Lijoptim_to_length_mediandiag_interp(events00,chsize,bins,CHR,R,mfull{1}+mfull{2}+mfull{3}+mfull{4},sij1dx, bins_event_tble);


%normalize lij by event ratios 
[lij]=renormalize_tiles(mfull00, lij, events00, bins, CHR);

Rn=R;
Rn=2*Rn/sum(Rn);        
p_mult = kron(Rn,Rn').*lij;

%normalize pmult by event ratios 
[p_mult]=renormalize_tiles(mfull00, p_mult, events00, bins, CHR);
%make sure pmult sums to 2 
p_mult = 2*p_mult ./ sum(sum(p_mult));




