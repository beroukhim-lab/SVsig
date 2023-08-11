function gen_bedfile( event_list )
%generate bed file for each sample and save it in a desktop library
% example of event_list: 
%event_list=annotated_table(annotated_table{:,end}==hit_num,:);

c_dir = strcat(pwd, '/bedfiles/')

[u_sample,ia,ic]=unique(event_list.sid);

for c1=1:length(u_sample),
    c_events=ic==c1;
    clear sample_event_list0 sample_event_list
    sample_event_list0=[event_list.seqnames(c_events) event_list.altchr(c_events) event_list.start(c_events) event_list.altpos(c_events)];
    ct=1;
    for c2=1:length(sample_event_list0(:,1)),
        if sample_event_list0(c2,1)==sample_event_list0(c2,2),
            sample_event_list(ct,:)=[sample_event_list0(c2,1) sample_event_list0(c2,3:4)];
        else
            sample_event_list(ct,:)=[sample_event_list0(c2,1) sample_event_list0(c2,3) sample_event_list0(c2,3)+100];
            ct=ct+1;
            sample_event_list(ct,:)=[sample_event_list0(c2,2) sample_event_list0(c2,4) sample_event_list0(c2,4)+100];
        end
        ct=ct+1;
    end
    writetable(array2table(sample_event_list),[c_dir u_sample{c1}],'WriteVariableName',0,'Delimiter','\t','QuoteStrings',false);
    movefile([c_dir u_sample{c1} '.txt'],[c_dir u_sample{c1} '.bed'])
end

