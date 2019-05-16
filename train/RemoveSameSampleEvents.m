function [bins_event_tble, bins, mfull, events, removed_events] = RemoveSameSampleEvents(bins_event_tble, bins, mfull, events,patient_column,flag)

%patient_column=8; % column in events array containing the sample information

% remove events from more than one sample per bin
disp('removing events from more than one sample per bin...');
[sort_bins,sort_bins_I]=sort(bins_event_tble(:,3)); % sort the linear indices of the bins events table
unique_nnz_bins = unique(sort_bins); % list of all unique values of the linear indices
removed_events=[];
%find rows in events table that are in the same bin%
for c1=1:length(unique_nnz_bins)
    idxs = sort_bins_I(sort_bins==unique_nnz_bins(c1)); % lines in events table with bin number unique_nnz_bins(c1)
    [sort_sample,sort_sample_I] = sort(events(idxs,patient_column)); % sorted sample numbers in the above bin 
%    [binj,bini] = ind2sub([numbins numbins],unique_nnz_bins(c1));
%for events in the same bin sort by sampleID and remove if two rows next to each other have the same sample ID   
for c2=2:length(idxs)
        if sort_sample(c2)==sort_sample(c2-1) 
            line2remove = idxs(sort_sample_I(c2)); % line in events table to remove
            bini = bins_event_tble(line2remove,1);
            binj = bins_event_tble(line2remove,2);
            if flag
                bins(bini,4)=bins(bini,4)-1;
                bins(binj,4)=bins(binj,4)-1;
            end
            annot=(events(line2remove,3)-1)*2+events(line2remove,6);
            %multiply removed events by their weights
            mfull{annot}(bini,binj) = mfull{annot}(bini,binj) - 1*events(line2remove, 11);
            mfull{annot}(binj,bini) = mfull{annot}(binj,bini) - 1*events(line2remove, 11);
            removed_events = [removed_events line2remove];
        end
    end
end
events(removed_events,:) = [];
bins_event_tble(removed_events,:) = [];

disp(['...removed ' num2str(length(removed_events)) ' events'])