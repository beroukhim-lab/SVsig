                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        % generate a sparse matrix with number of events per bin 
% and a table with the starting and ending bin of each event  
function [mfull , bins_event_tble] = BuildMatrix(events,bins,num_annot)

disp('building event matrix...');

%why make num_annot > 1 ?
%artifact from Ofer's code for binning 4 types of rearrangments
%but I don't see where events are being subsetted

for ca=1:num_annot
    mfull{ca}=sparse(length(bins),length(bins));
end

%bins: column 1 is chromosome, column 2 is start position, column 3 is end
%position, column 4 is number of events

%events: column 1: chr_i, column 2: start_i, column 3: strand_i, column 4:
%chr_j, column 5: start_j, column 6: stand_j

%initialize bins_event_tble
bins_event_tble = zeros(length(events),5);

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
    annot=(events(c1,3)-1)*2+events(c1,6);
    %mark an event appropriate mfull matrix with 
    mfull{annot}(bini,binj) = mfull{annot}(bini,binj)+1*events(c1, 11) ;%multiply each event by its weight (events column 11)
end

%this is adding the matrix and its transpose
for ca=1:num_annot
    mfull{ca}(:,:) = mfull{ca}(:,:) + mfull{ca}(:,:)';
end

