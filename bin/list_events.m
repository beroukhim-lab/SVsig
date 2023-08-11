function event_list=list_events(events,posi,posj,pad,UTumor, Uevent, Usample, Upatient,filename)

if ~isempty(pad),
    posi(2)=posi(2)-pad;
    posi(3)=posi(3)+pad;
    posj(2)=posj(2)-pad;
    posj(3)=posj(3)+pad;
end

if istable(events),
    if isempty(posj),
        event_list=events( events{:,1}==posi(1)&events{:,2}>=posi(2)&events{:,2}<posi(3)...
        | events{:,4}==posi(1)&events{:,5}>=posi(2)&events{:,5}<posi(3) ,:);
    else
        event_list=events( events{:,1}==posi(1)&events{:,2}>=posi(2)&events{:,2}<posi(3)&events{:,4}==posj(1)&events{:,5}>=posj(2)&events{:,5}<posj(3)...
        | events{:,1}==posj(1)&events{:,2}>=posj(2)&events{:,2}<posj(3)&events{:,4}==posi(1)&events{:,5}>=posi(2)&events{:,5}<posi(3) ,:);
    end
else
    if isempty(posj),
        event_list=events( events(:,1)==posi(1)&events(:,2)>=posi(2)&events(:,2)<posi(3)...
        | events(:,4)==posi(1)&events(:,5)>=posi(2)&events(:,5)<posi(3) ,:);
    else
        event_list=events( events(:,1)==posi(1)&events(:,2)>=posi(2)&events(:,2)<posi(3)&events(:,4)==posj(1)&events(:,5)>=posj(2)&events(:,5)<posj(3)...
        | events(:,1)==posj(1)&events(:,2)>=posj(2)&events(:,2)<posj(3)&events(:,4)==posi(1)&events(:,5)>=posi(2)&events(:,5)<posi(3) ,:);
    end

    if ~isempty(filename);
        chr_i=[num2str(event_list(:,1))];
        pos_i=[num2str(event_list(:,2))];
        chr_j=[num2str(event_list(:,4))];
        pos_j=[num2str(event_list(:,5))];
        tumor=UTumor{event_list(:,7),1};
        event=Uevent{event_list(:,8),1};
        sample=Usample{event_list(:,9),1};
        patient=Upatient{event_list(:,10),1};
        T=table(chr_i,pos_i,chr_j,pos_j,tumor,event,sample,patient);
        writetable(T,filename);
    end
end

return

    