% generate a table of bins with chr number, start and end position and
% number of events


function [bins, numbins] = SetBins(events,num_breakpoints_per_bin,chsize,CHR, min_bin_dist)

disp('setting bins boundaries...');

max_no_events=1e5;
range = round(num_breakpoints_per_bin/2); % range to look for maximum distance between bin-separating events

chr = [events(:,1);events(:,4)];
pos = [events(:,2);events(:,5)];

bc=1;
%for all the chrmosomes
for c1 = CHR

    posc=0;
    nextpos=0;
    chri = find(chr==c1);
    sorted_pos = sort(pos(chri));
    bins(bc,2)=1;
    
    while posc<length(chri)-num_breakpoints_per_bin-range
 
        [max_diff,max_diff_loc]=max(diff(sorted_pos(posc+num_breakpoints_per_bin-range:posc+num_breakpoints_per_bin+range)));

        while max_diff<min_bin_dist && posc<length(chri)-num_breakpoints_per_bin-3*range
            posc=posc+2*range;
            [max_diff,max_diff_loc]=max(diff(sorted_pos(posc+num_breakpoints_per_bin-range:posc+num_breakpoints_per_bin+range)));
        end
        
        bins(bc,1) = c1;
        next_event=posc+num_breakpoints_per_bin-range+max_diff_loc;
        bins(bc,3) = sorted_pos(next_event)-round(max_diff/2);
        bins(bc+1,2) = bins(bc,3)+1;
        bins(bc,4) = next_event - nextpos;
        posc = next_event;
        nextpos = posc;
        bc=bc+1;

    end
    
    if length(sorted_pos)-nextpos < num_breakpoints_per_bin
        bins(bc-1,3) = chsize(c1);
        bins(bc-1,4) = bins(bc-1,4) + length(sorted_pos)-nextpos;
    else
        bins(bc,1) = c1;
        bins(bc,3) = chsize(c1);
        bins(bc,4) = length(sorted_pos)-nextpos;
        bc=bc+1;
    end

end

if bins(end,4)==0,bins(end,:)=[];end
numbins = length(bins);   

