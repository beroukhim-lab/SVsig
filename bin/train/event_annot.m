function annot_array=event_annot(events,varargin)
% columns organization:
% 1: all ones
% 2: short events 0-1e6
% 3: long events 1e6-
% 4: inter-chromosomal translocations
% 5: simple events
% 6: complex events
% 7: HR
% 8: MMEJ
% 9: NHEJ
% 10: MMBIR
% 11 and on are annotations given by varargin

% 6: deletions (+-)
% 7: duplications (-+)
% 8: inversions (--/++)
% 9: TADs events


%events matrix column organization
%12: topo
%13: topo_n
%14: mech
%15:homseq

short=1e6;
%long=1e7;

pad=1e4;
numa=nargin-1;
nume=length(events);

annot_array=zeros(nume,7+numa);
annot_array(:,1)=ones(nume,1);
annot_array(:,2)=events(:,1)==events(:,4)&abs(events(:,5)-events(:,2))<=short;
annot_array(:,3)=events(:,1)==events(:,4)&abs(events(:,5)-events(:,2))>short;
annot_array(:,4)=events(:,1)~=events(:,4);
annot_array(:,5)=events(:,13)==1;
annot_array(:,6)=events(:,13)>1;
annot_array(:,7)=events(:,14)==3|events(:,14)==5&events(:,15)<=2;
annot_array(:,8)=events(:,14)==3|events(:,14)==5&events(:,15)>2&events(:,15)<13;
annot_array(:,9)=events(:,14)>=3&events(:,14)<=6&events(:,15)>=13;
annot_array(:,10)=events(:,14)==2;
%annot_array(:,7)=events(:,14)==4|events(:,14)==6;
%annot_array(:,8)=events(:,14)==3;
%annot_array(:,9)=events(:,14)==5;
%annot_array(:,10)=events(:,14)==2;

for c1=1:nume 
    c2=1;
    annot_array(c1,11)=sum( (events(c1,1)==varargin{c2}(:,1)&events(c1,2)>=varargin{c2}(:,2)-pad&events(c1,2)<=varargin{c2}(:,3)+pad& ...
                events(c1,4)==varargin{c2}(:,1)&events(c1,5)>=varargin{c2}(:,2)-pad&events(c1,5)<=varargin{c2}(:,3)+pad) )>0;
    for c2=2:numa
        annot_array(c1,10+c2)=sum( (events(c1,1)==varargin{c2}(:,1)&events(c1,2)>=varargin{c2}(:,2)-pad&events(c1,2)<=varargin{c2}(:,3)+pad) | ...
            (events(c1,4)==varargin{c2}(:,1)&events(c1,5)>=varargin{c2}(:,2)-pad&events(c1,5)<=varargin{c2}(:,3)+pad)) >0;
    end
end
% for c1=1:nume,
%     
% 
% %    annot_array(c1,3)=events(c1,1)==events(c1,4)&abs(events(c1,5)-events(c1,2))>short;
% %    annot_array(c1,4)=events(c1,1)~=events(c1,4);
% %    annot_array(c1,5)=events(c1,1)==events(c1,4)&events(c1,3)==1&events(c1,6)==2;
% %    annot_array(c1,6)=events(c1,1)==events(c1,4)&events(c1,3)==2&events(c1,6)==1;
% %    annot_array(c1,7)=events(c1,1)==events(c1,4)&events(c1,3)==events(c1,6);
%     c2=1;
%     annot_array(c1,8)=sum( (events(c1,1)==varargin{c2}(:,1)&events(c1,2)>=varargin{c2}(:,2)-pad&events(c1,2)<=varargin{c2}(:,3)+pad& ...
%             events(c1,4)==varargin{c2}(:,1)&events(c1,5)>=varargin{c2}(:,2)-pad&events(c1,5)<=varargin{c2}(:,3)+pad) )>0;
%     for c2=2:numa
%         annot_array(c1,7+c2)=sum( (events(c1,1)==varargin{c2}(:,1)&events(c1,2)>=varargin{c2}(:,2)-pad&events(c1,2)<=varargin{c2}(:,3)+pad) | ...
%             (events(c1,4)==varargin{c2}(:,1)&events(c1,5)>=varargin{c2}(:,2)-pad&events(c1,5)<=varargin{c2}(:,3)+pad)) >0;
%     end
% end
