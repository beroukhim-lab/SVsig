function [tier1,tier2]=load_tumor_subtype(UTumor)
%%% I think this function takes in a table and returns a array of
%%% histologies that match the table of samples(?) 

%%preprocessed dictionary of histology classes%%
[num,txt,raw] =xlsread('tumor_subtype.xlsx');

%%% match samples with histology dictionary%
tier1=cell(height(UTumor),1);
tier2=cell(height(UTumor),1);
for c1=1:height(UTumor)
    for c2=2:length(txt),
        if cell2mat(regexp(txt{c2,5},UTumor{c1,1}))>0
            tier1(c1)=txt(c2,3);
            tier2(c1)=txt(c2,1);
        end
    end
end



