function annot_tiles=tiles_annot(track,events,bins,CHR)
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
        annot_tiles=false(lenbins,lenbins,3);
        annot_tiles(:,:,3)=true(lenbins,lenbins);
        for c1=CHR,
            firstbin=find(bins(:,1)==c1,1);
            lastbin=find(bins(:,1)==c1,1,'last');
            chr_intra = firstbin:lastbin; 
            annot_tiles(chr_intra,chr_intra,3)=false;
            for c2=chr_intra
                diag_bin = sum(bins(c2,2:3),2)/2;
                diag_bin_size=(bins(c2,3)-bins(c2,2));
                upper_diag_bins=abs(sum(bins(chr_intra,2:3),2)/2-diag_bin);
                annot_tiles(c2,chr_intra,1)=upper_diag_bins>=0&upper_diag_bins<short0;
                annot_tiles(c2,chr_intra,2)=upper_diag_bins>=short0&upper_diag_bins<1e9;
%                annot_tiles(c2,chr_intra,3)=upper_diag_bins>=long0&upper_diag_bins<1e9;
            end
        end
    end    
end

return


% for c1
%     annot_array(c1,2)=events(c1,1)==events(c1,4)&abs(events(c1,5)-events(c1,2))<short;
%     annot_array(c1,3)=events(c1,1)==events(c1,4)&(abs(events(c1,5)-events(c1,2))>=short&abs(events(c1,5)-events(c1,2))<long);
%     annot_array(c1,4)=events(c1,1)==events(c1,4)&abs(events(c1,5)-events(c1,2))>=long;
%     annot_array(c1,5)=events(c1,1)~=events(c1,4);
%     annot_array(c1,6)=events(c1,1)==events(c1,4)&events(c1,3)==1&events(c1,6)==2;
%     annot_array(c1,7)=events(c1,1)==events(c1,4)&events(c1,3)==2&events(c1,6)==1;
%     annot_array(c1,8)=events(c1,1)==events(c1,4)&events(c1,3)==events(c1,6);
%     c2=1;
%     annot_array(c1,9)=sum( (events(c1,1)==track(:,1)&events(c1,2)>=track(:,2)-pad&events(c1,2)<=track(:,3)+pad& ...
%             events(c1,4)==track(:,1)&events(c1,5)>=track(:,2)-pad&events(c1,5)<=track(:,3)+pad) )>0;
%     for c2=2:numa
%         annot_array(c1,8+c2)=sum( (events(c1,1)==track(:,1)&events(c1,2)>=track(:,2)-pad&events(c1,2)<=track(:,3)+pad) | ...
%             (events(c1,4)==track(:,1)&events(c1,5)>=track(:,2)-pad&events(c1,5)<=track(:,3)+pad)) >0;
%     end
% end

