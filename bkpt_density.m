   

events=events0;
bins=bins0;
% generate a sparse matrix with number of events per bin 
bins_event_tble = zeros(length(events),3);

%loop through all events (table of pairs of breakpoints)
%for each row
for c1 = 1:length(events)   
    %find the indices of the bins in which the ith and jth breaks fall in        
    bins_event_tble(c1,1) = find(bins(:,1)==events(c1,1) & bins(:,2)<=events(c1,2) & bins(:,3)>=events(c1,2));
    bins_event_tble(c1,2) = find(bins(:,1)==events(c1,4) & bins(:,2)<=events(c1,5) & bins(:,3)>=events(c1,5));
    %bin i is the the "smallest" or leftmost bin if we think of the genome
    %as chr 1 -> chr 24
    %bin j is the "largest" or rightmost bin
    bini=min(bins_event_tble(c1,1),bins_event_tble(c1,2));
    binj=max(bins_event_tble(c1,1),bins_event_tble(c1,2));
    
    %sub2ind turn two dimensional indices to a unique one-dimensional tile
    %index
    bins_event_tble(c1,3) = sub2ind([length(bins) length(bins)],bini,binj); % note: assign values only to upper tria of the matrix
    %annot is the type of rearrangement as determined by the strand
    %orientation 
end

 
bins_event_tble10k=bins_event_tble;
