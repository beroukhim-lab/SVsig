

mixmodel_pcawg=load('/Users/shu/2d_results/pcawg_noshort_sijdx.mat', 'sij1dx', 'sij1dy');
pcawg_sij1dx=mixmodel_pcawg.sij1dx;
pcawg_sij1dy=mixmodel_pcawg.sij1dy./sum(mixmodel_pcawg.sij1dy, 'all');

sijdy_norm=sij1dy/sum(sij1dy, 'all');

pcawg_sij1dy_interp=interp1(pcawg_sij1dx,pcawg_sij1dy,pcawg_sij1dx,'pchip');
sij1dy_interp=interp1(sij1dx,sijdy_norm,sij1dx,'pchip');


plot(log10(pcawg_sij1dx),log10(pcawg_sij1dy),'*c',log10(pcawg_sij1dx(1:99)),log10(pcawg_sij1dy_interp(1:99)),'-c');
hold on
plot(log10(sij1dx),log10(sijdy_norm),'*k',log10(sij1dx(1:99)),log10(sij1dy_interp(1:99)),'-k');





mixmodel_pcawg=load('/Users/shu/2d_results/20211205_pcawg_500kb.mat', 'sij1dx', 'sij1dy1');
pcawg_sij1dx=mixmodel_pcawg.sij1dx;
pcawg_sij1dy=mixmodel_pcawg.sij1dy1./sum(mixmodel_pcawg.sij1dy1, 'all');

pcawg_sij1dy_interp=interp1(pcawg_sij1dx,pcawg_sij1dy,pcawg_sij1dx,'pchip');
plot(log10(pcawg_sij1dx),log10(pcawg_sij1dy),'*c',log10(pcawg_sij1dx(1:9)),log10(pcawg_sij1dy_interp(1:9)),'-c');


test=diff(pcawg_sij1dx);

plot(log10(pcawg_sij1dx(1:9)),log10((pcawg_sij1dy_interp(1:9)).* test),'-k');


hold on

sij1dy_current=sum(sij1dy(:,1,:), 3);
sij1dy_current=sij1dy_current./sum(sij1dy_current, 'all');

sij1dy_interp=interp1(sij1dx,sij1dy_current,sij1dx,'pchip');
plot(log10(sij1dx),log10(sij1dy_current),'*k',log10(sij1dx(1:9)),log10(sij1dy_interp(1:9)),'-k');


hold on 

plot(log10(sij1dx), log10(sij1dy_current), '*m')

%old_a1=load('20210917_mixmodel_a1.mat', 'sij1dx', 'sij1dy');
old_a1=load('20210623_mixmodel_a0.mat', 'sij1dx', 'sij1dy');

olda1_sij1dx=old_a1.sij1dx;
olda1_sij1dy=sum(old_a1.sij1dy(:,1,:), 3);
olda1_sij1dy=olda1_sij1dy./sum(olda1_sij1dy, 'all');

olda1_sij1dy_interp=interp1(olda1_sij1dx,olda1_sij1dy,olda1_sij1dx,'pchip');
plot(log10(olda1_sij1dx),log10(olda1_sij1dy),'*r',log10(olda1_sij1dx(1:9)),log10(olda1_sij1dy_interp(1:9)),'-r');


mixmodel_matrices=load('/Users/shu/2d_results/20211205_pcawg_500kb.mat', 'sij1', 'lij', 'p', 'p_mult');
p_diff=(mixmodel_matrices.p)./p;
pmult_diff=(mixmodel_matrices.p_mult)./p_mult;
sij_diff=(mixmodel_matrices.sij1)./sij1;
lij_diff=(mixmodel_matrices.lij)./lij;



idx=find(p_diff>15);
[rows, columns] = ind2sub(size(p_diff), idx);


dist1=[]
 [rowi, coli]=ind2sub(size(p_diff), idx);
 
 
 for i=1:length(rowi)
     
     if bins(rowi(i),1)==bins(coli(i),1)
        first_index=bins(rowi(i),2);
        second_index=bins(coli(i),2);
        if second_index>=first_index
            dist1=[dist1 (second_index-first_index)];
        end 
     end
 end 
    
    
    
    
    disp(rowi, coli)
    dist = [dist newval]
    
    
    
    
    




