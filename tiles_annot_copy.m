function annot_tiles=tiles_annot_copy(track,events,bins,CHR)
% columns organization:
% 1: all ones
% 2: short events 0-5e5
% 3: medium events 5e5-1e7
% 4: short events 1e7-
% 5: translocations
% 6: deletions (+-)
% 7: duplications (-+)
% 8: inversions (--/++)
% 9: TADs events
% 10 and on are annotations given by varargin
pad=1e4;
nume=length(events);
lenbins=length(bins);

short0=1e6;
long0=5e7;

if ischar(track),
    if strcmp(track,'length')
        annot_tiles=false(lenbins,lenbins,4);
        annot_tiles(:,:,4)=true(lenbins,lenbins);
        for c1=CHR,
            firstbin=find(bins(:,1)==c1,1);
            lastbin=find(bins(:,1)==c1,1,'last');
            chr_intra = firstbin:lastbin; 
            annot_tiles(chr_intra,chr_intra,4)=false;
            for c2=chr_intra
                diag_bin = sum(bins(c2,2:3),2)/2;
                diag_bin_size=(bins(c2,3)-bins(c2,2));
                upper_diag_bins=abs(sum(bins(chr_intra,2:3),2)/2-diag_bin);
                annot_tiles(c2,c2,1)=true;
                annot_tiles(c2,chr_intra,2)=upper_diag_bins>=0&upper_diag_bins<short0;
                annot_tiles(c2,c2,2)=false;
                annot_tiles(c2,chr_intra,3)=upper_diag_bins>=short0&upper_diag_bins<1e9;
%                annot_tiles(c2,chr_intra,3)=upper_diag_bins>=long0&upper_diag_bins<1e9;
            end
        end
    end    
end


%uncomment below to do based on percentage of events per tile

% [row1, col1] = find(annot_tiles(:,:,1) ==1);
% indices=cat(2,row1,col1);
% upper_i=(indices(:,1)>=indices(:,2));
% indices=indices(upper_i,:);
% 
% for x=1:length(indices)
%     
%     ibin=indices(x,1);
%     jbin=indices(x,2);
%     
%     ind1 = (events(:,1) == bins(ibin,1)) & (events(:,2) >bins(ibin,2)) & (events(:,2) <bins(ibin,3));
%     events_ij = events(ind1,:);
%     ind2 = (events_ij(:,4) == bins(jbin,1)) & (events_ij(:,5) >bins(jbin,2)) & (events_ij(:,5) <bins(jbin,3));
%     events_ij = events_ij(ind2,:);
%     short_events_ij=find(abs(events_ij(:,5) - events_ij(:,2)) <= 1e6);
%     
%     ind1 = (events(:,1) == bins(jbin,1)) & (events(:,2) >bins(jbin,2)) & (events(:,2) <bins(jbin,3));
%     events_ji = events(ind1,:);
%     ind2 = (events_ji(:,4) == bins(ibin,1)) & (events_ji(:,5) >bins(ibin,2)) & (events_ji(:,5) <bins(ibin,3));
%     events_ji = events_ji(ind2,:);
%     short_events_ji=find(abs(events_ji(:,5) - events_ji(:,2)) <= 1e6);
% 
%     short_events_x=cat(1, short_events_ij, short_events_ji);
%     
%   
%     
%     short_ratio=length(short_events_x)/(length(events_ij) + length(events_ji));
%     if (length(events_ij) + length(events_ji))>0
%         if short_ratio<=0.5 || isempty(short_events_x) 
%             annot_tiles(ibin,jbin,1)=false;
%             annot_tiles(jbin,ibin,1)=false;
%             annot_tiles(ibin,jbin,2)=true;
%             annot_tiles(jbin,ibin,2)=true;
% 
%         
%         end
%   
%     end 
% end 


    
     

return


