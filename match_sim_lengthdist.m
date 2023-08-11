
%%%%these variables taken at line 96 (before conditional prob stuff) 
%%%%of break_invasion_copy.mat
mixmodel_pcawg=load('/Users/shu/2d_results/20220315_sijdx_pcawg_events_greater1e5.mat', 'sij1dx', 'sij1dy', 'annot_frac');
pcawg_sij1dx=mixmodel_pcawg.sij1dx;
pcawg_sij1dy=mixmodel_pcawg.sij1dy./sum(mixmodel_pcawg.sij1dy, 'all');
pcawg_annot_frac=mixmodel_pcawg.annot_frac;

% for ca=1:4
%     pcawg_sij1dy(ca,:) = pcawg_sij1dy(ca,:) * pcawg_annot_frac(ca);
% end
pcawg_sij1dy = sum(pcawg_sij1dy,2);     % summing across types of events (4)
pcawg_sij1dy = pcawg_sij1dy/sum(pcawg_sij1dy);  


n_pcawg=length(pcawg_sij1dy);
pcawg_sij1dy_interp=interp1(pcawg_sij1dx,pcawg_sij1dy,pcawg_sij1dx,'pchip');
plot(log10(pcawg_sij1dx),log10(pcawg_sij1dy),'*c',log10(pcawg_sij1dx(1:(n_pcawg-1))),log10(pcawg_sij1dy_interp(1:(n_pcawg-1))),'-c');

hold on

%plot simulation pcawg events
sij1dy_current=sum(sij1dy(:,1,:), 3);
sij1dy_current=sij1dy_current./sum(sij1dy_current, 'all');

n=length(sij1dy_current);
ratio=n/n_pcawg;
sij1dy_interp=interp1(sij1dx,sij1dy_current,sij1dx,'pchip');
plot(log10(sij1dx),log10(sij1dy_current),'*k',log10(sij1dx(1:(n-1))),log10(sij1dy_interp(1:(n-1))),'-k');



%%%%%%%%%%%%%%%%
%if need to calculate annot frac 
events=events00;
nume=length(events)
annot_frac(1) = sum(events(:,3)==1&events(:,6)==1)/nume;
annot_frac(2) = sum(events(:,3)==1&events(:,6)==2)/nume;
annot_frac(3) = sum(events(:,3)==2&events(:,6)==1)/nume;
annot_frac(4) = sum(events(:,3)==2&events(:,6)==2)/nume;

