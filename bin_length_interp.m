for c1=CHR
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin:lastbin; 
    
    c2=firstbin;
    %find midpoint of bin
    diag_bin = sum(bins(c2,2:3),2)/2;
    diag_bin_size=(bins(c2,3)-bins(c2,2));
    %find distance of other bins in row from diag bin
    upper_diag_bins=sum(bins(c2+1:lastbin,2:3),2)/2-diag_bin;
    %bp size of each bin in the row divided by 2 
    half_size=(bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2))/2;
    %does cubic interpolation for start and end of each bin (based on bp
    %distance away from diag bin), averages these starts and ends
    upper_diag=(interp1(sij1dx',squeeze(sij1dy(:,c1,:)),upper_diag_bins-half_size,'pchip')+interp1(sij1dx',squeeze(sij1dy(:,c1,:)),upper_diag_bins+half_size,'pchip'))/2;
    %sij(c2,c2+1:lastbin,:) = bsxfun(@times,upper_diag,bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2));
    sij(c2,c2+1:lastbin,:) = upper_diag;
    last_diag = find(sij1dx<diag_bin_size,1,'last');
    % add similar correction to bin 2 and 3...



median_interp_tble = zeros(length(events),3);

%loop through all events (table of pairs of breakpoints)
%for each row

intra_events=events(events(:,1)==events(i,4));
for i = 1:length(intra_events)   
    %find the indices of the bins in which the ith and jth breaks fall in        
    bin1 = find(bins(:,1)==events(i,1) & bins(:,2)<=intra_events(i,2) & bins(:,3)>=intra_events(i,2));
    bin2 = find(bins(:,1)==events(i,4) & bins(:,2)<=intra_events(i,5) & bins(:,3)>=intra_events(i,5));
    %bin i is the the "smallest" or leftmost bin 
    %bin j is the "largest" or rightmost bin
    bini=min(bin1,bin2);
    binj=max(bin1,bin2);
    
    median_interp_tble(i,1) = (bins(binj,3)-bins(binj,2))/2 - (bins(bini,3)-bins(bini,2))/2;
    median_interp_tble(i,2) = abs(intra_events(i,4) - intra_events(i,1));

end


