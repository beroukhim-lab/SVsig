function [bins_event_tble, bins, mfull, events, removed_events] = RemoveZeroVarSampleEvents(bins_event_tble, bins, mfull, events)

% remove events from more than one sample per bin
disp('removing bins with zero variance...');
[sort_bins,sort_bins_I]=sort(bins_event_tble(:,3)); % sort the linear indices of the bins events table
unique_nnz_bins = unique(sort_bins); % list of all unique values of the linear indices
sizem=size(mfull{1});
removed_events=[];
for c1=1:length(unique_nnz_bins),
    idxs = sort_bins_I(sort_bins==unique_nnz_bins(c1)); % lines in events table with bin number unique_nnz_bins(c1)
    std_i = std(events(idxs,2));
    std_j = std(events(idxs,5));
    [binj,bini] = ind2sub(sizem,unique_nnz_bins(c1));
    if length(idxs)>=2 && (std_i<10 || std_j<10),
        bins(bini,4)=bins(bini,4)-length(idxs);
        bins(binj,4)=bins(binj,4)-length(idxs);
%        annot=(events(line2remove,3)-1)*2+events(line2remove,6);
        mfull{1}(bini,binj) = 0;
        mfull{1}(binj,bini) = 0;
        mfull{2}(bini,binj) = 0;
        mfull{2}(binj,bini) = 0;
        mfull{3}(bini,binj) = 0;
        mfull{3}(binj,bini) = 0;
        mfull{4}(bini,binj) = 0;
        mfull{4}(binj,bini) = 0;
        removed_events = [removed_events;idxs];
%        disp(strcat('removed ',num2str(length(idxs)),' events from bin ',num2str(bini),',',num2str(binj)));
    end
end
events(removed_events,:) = [];
bins_event_tble(removed_events,:) = [];

disp(['...removed ' num2str(length(removed_events)) ' events'])