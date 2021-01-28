function [bins, numbins] = SetBins_bylength(events,avg_bin_length,chsize,CHR)

disp('setting bins boundaries...');

max_no_events=1e5;

chr = [events(:,1);events(:,4)];
pos = [events(:,2);events(:,5)];

bc=1;
%for all the chrmosomes
for c1 = CHR
    
    num_chr=floor(chsize(c1)/avg_bin_length);
    avg_length_chr=floor(avg_bin_length + ((chsize(c1)-num_chr*avg_bin_length)/num_chr));
    
    chri = find(chr==c1);
    sorted_pos = sort(pos(chri));
    
    start1=1;
    for c2 = 1:num_chr        
        bins(bc,1) = c1;
        bins(bc,2)=start1;
        bins(bc,3) =start1+avg_length_chr-1;        
        breaks_per_bin=find(sorted_pos(:,1)<=(start1+avg_length_chr) & sorted_pos(:,1)>=start1);
        bins(bc,4) = length(breaks_per_bin);     
        bc=bc+1;
        start1=start1+avg_length_chr;
        
    end
        
   

end

numbins = length(bins);   

