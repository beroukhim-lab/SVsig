function [bins0, events0, numbins] = remove_low_density_bins(bins0,events0)
% remove bins AND coresponding events with density lower than a threshold
% equal to: ???
% density represents the number of events per base pair in a particular bin
%density takes into account the variable bin size, otherwise if the bins
%were the same size you could just set a threshold for removing bins that have less than a certain number of events


%[sorted_desnity,idx]=sort(bins0(:,4)./(bins0(:,3)-bins0(:,2)));
%th=5e-5;
%20220628 original threshold is 5e-5
th=5e-5;
%th=5e-6;

%bins columns: 1 is for the chromosome, 2 is the start position, 3 is the
%end position, 4 is the number of events per bin
bins_to_remove=(bins0(:,4)./(bins0(:,3)-bins0(:,2))<th);

btr=find(bins_to_remove);
events_to_remove=zeros(length(events0),1);
for c1=1:length(btr)
    posi=bins0(btr(c1),1:3);
    events_to_remove = events_to_remove | events0(:,1)==posi(1)&events0(:,2)>=posi(2)&events0(:,2)<=posi(3)...
    | events0(:,4)==posi(1)&events0(:,5)>=posi(2)&events0(:,5)<=posi(3);
end

etr=find(events_to_remove);
for c1=1:sum(events_to_remove)
    bini=events0(etr(c1),1)==bins0(:,1)&events0(etr(c1),2)>=bins0(:,2)&events0(etr(c1),2)<=bins0(:,3);
    bins0(bini,4)=bins0(bini,4)-1;
    binj=events0(etr(c1),4)==bins0(:,1)&events0(etr(c1),5)>=bins0(:,2)&events0(etr(c1),5)<=bins0(:,3);
    bins0(binj,4)=bins0(binj,4)-1;
    if sum(bini)>1 || sum(binj)>1
        c1
    end
end

events0(events_to_remove,:)=[];
bins0(bins_to_remove,:)=[];

numbins=length(bins0);

disp(['removed ' num2str(sum(events_to_remove)) ' events in ' num2str(sum(bins_to_remove)) 'bins'])

end

