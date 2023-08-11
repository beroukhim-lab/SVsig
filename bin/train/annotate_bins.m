% generates a vector for each bin and each annoataion with 1 if:
% flag==0 if at least th=0.9 of events are of from the track
% flag==1 fortces 1 if any event from the track is in the bin
function bins_annot=annotate_bins(bins,bins_event_tble,events_annot,flag,varargin)

th=0.9;
pad=1e3;
if flag==1,
    numa=nargin-4;
else
    numa=length(events_annot(1,:));
end
nume=length(bins);
events_annot=logical(events_annot);

bins_annot=zeros(nume,numa);
if flag==1,
    for c1=1:nume,
        for c2=1:numa
            bins_annot(c1,c2)=sum(bins(c1,1)==varargin{c2}(:,1)&bins(c1,3)>=varargin{c2}(:,2)-pad&bins(c1,2)<=varargin{c2}(:,3)+pad)>0;
        end
    end
    bins_annot=[ones(nume,1) bins_annot];
else
    for c1=1:nume,
        for c2=1:numa
            bins_annot(c1,c2)=(sum(bins_event_tble(events_annot(:,c2),1)==c1)+sum(bins_event_tble(events_annot(:,c2),2)==c1))>bins(c1,4)*th;
        end
    end
end

