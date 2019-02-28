% generate a sparse matrix with number of events per bin 
% and a table with the starting and ending bin of each event  
function [mfull , bins_event_tble] = BuildMatrix(events,bins,num_annot)

disp('building event matrix...');

for ca=1:num_annot,
    mfull{ca}=sparse(length(bins),length(bins));
end

%initialize bins_event_tble
bins_event_tble = zeros(length(events),5);
for c1 = 1:length(events),
    bins_event_tble(c1,1) = find(bins(:,1)==events(c1,1) & bins(:,2)<=events(c1,2) & bins(:,3)>=events(c1,2));
    bins_event_tble(c1,2) = find(bins(:,1)==events(c1,4) & bins(:,2)<=events(c1,5) & bins(:,3)>=events(c1,5));
    bini=min(bins_event_tble(c1,1),bins_event_tble(c1,2));
    binj=max(bins_event_tble(c1,1),bins_event_tble(c1,2));
    bins_event_tble(c1,3) = sub2ind([length(bins) length(bins)],bini,binj); % note: assign values only to upper tria of the matrix
    annot=(events(c1,3)-1)*2+events(c1,6);
    mfull{annot}(bini,binj) = mfull{annot}(bini,binj)+1;
end

for ca=1:num_annot
    mfull{ca}(:,:) = mfull{ca}(:,:) + mfull{ca}(:,:)';
end

