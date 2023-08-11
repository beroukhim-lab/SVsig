function d = plot_table(TbyGene,flag)

if flag==0,
    cnames={'hit','Gene i','pos i','Gene j','pos j','p-val','N','annot','Type-1','Type-2','Type-3'};
    column_width = {30 100 60 100 60 60 30 30 100 100 100};
    column_format = {[] [] [] [] [] 'numeric' 'numeric' [] [] [] []};

    f = figure('Position',[340 500 800 650]);
    t = uitable('ColumnName', cnames,'FontSize',12,'RowName',[],'ColumnWidth',column_width,'ColumnFormat',column_format);
    t.RowStriping='off';
    t.Position = [20 20 sum([column_width{:}]) 700];

    ct=1;
    for c1=1:length(TbyGene),
        
        % do not plot clusters with only one event
        if TbyGene(c1).num_events<=1,
            continue
        end
        
        if isempty(TbyGene(c1).tumor_1),
            tumor1=[]; tumor2=[]; tumor3=[]; 
        else
            tumor1=TbyGene(c1).tumor_1{1}; tumor2=TbyGene(c1).tumor_2{1}; tumor3=TbyGene(c1).tumor_3{1}; 
        end
        if ~isempty(TbyGene(c1).gene_i) & isempty(TbyGene(c1).gene_j),        
            d(ct,:) = {c1 TbyGene(c1).gene_i{1} strcat(TbyGene(c1).pos_i(1:6),'M') [] strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).annot tumor1 tumor2 tumor3};
        elseif isempty(TbyGene(c1).gene_i) & ~isempty(TbyGene(c1).gene_j),
            d(ct,:) = {c1 [] strcat(TbyGene(c1).pos_i(1:6),'M') TbyGene(c1).gene_j{1} strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).annot tumor1 tumor2 tumor3};
        elseif isempty(TbyGene(c1).gene_i) & isempty(TbyGene(c1).gene_j),
            d(ct,:) = {c1 [] strcat(TbyGene(c1).pos_i(1:6),'M') [] strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).annot tumor1 tumor2 tumor3};
        else
            d(ct,:) = {c1 TbyGene(c1).gene_i{1} strcat(TbyGene(c1).pos_i(1:6),'M') TbyGene(c1).gene_j{1} strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).annot tumor1 tumor2 tumor3};
        end
        ct=ct+1;
        len_i=length(TbyGene(c1).gene_i);
        len_j=length(TbyGene(c1).gene_j);
        if len_i>1 || len_j>1,
            for c2=2:max(length(TbyGene(c1).gene_i),length(TbyGene(c1).gene_j))
                if c2>len_i
                    d(ct,:) = {[] [] [] TbyGene(c1).gene_j{c2} [] [] [] [] [] [] []};
                elseif c2>len_j
                    d(ct,:) = {[] TbyGene(c1).gene_i{c2} [] [] [] [] [] [] [] [] []};
                else
                    d(ct,:) = {[] TbyGene(c1).gene_i{c2} [] TbyGene(c1).gene_j{c2} [] [] [] [] [] [] []};
                end
                ct=ct+1;
            end
        end
        
    end


    t.Data = d;
else
        cnames={'hit','Gene i','pos i','Gene j','pos j','p-val','N','Type-1','Type-2','Type-3'};
    column_width = {30 100 60 100 60 60 30 100 100 100};
    column_format = {[] [] [] [] [] 'numeric' 'numeric' [] [] []};

    f = figure('Position',[340 500 800 650]);
    t = uitable('ColumnName', cnames,'FontSize',12,'RowName',[],'ColumnWidth',column_width,'ColumnFormat',column_format);
    t.RowStriping='off';
    t.Position = [20 20 sum([column_width{:}]) 600];

    ct=1;
    for c1=1:length(TbyGene),
        if ~isempty(TbyGene(c1).gene_i) & isempty(TbyGene(c1).gene_j),        
            d(ct,:) = {c1 TbyGene(c1).gene_i{1} strcat(TbyGene(c1).pos_i(1:6),'M') [] strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).tumor_1{1} TbyGene(c1).tumor_2{1} TbyGene(c1).tumor_3{1}};
        elseif isempty(TbyGene(c1).gene_i) & ~isempty(TbyGene(c1).gene_j),
            d(ct,:) = {c1 [] strcat(TbyGene(c1).pos_i(1:6),'M') TbyGene(c1).gene_j{1} strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).tumor_1{1} TbyGene(c1).tumor_2{1} TbyGene(c1).tumor_3{1}};
        elseif isempty(TbyGene(c1).gene_i) & isempty(TbyGene(c1).gene_j),
            d(ct,:) = {c1 [] strcat(TbyGene(c1).pos_i(1:6),'M') [] strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).tumor_1{1} TbyGene(c1).tumor_2{1} TbyGene(c1).tumor_3{1}};
        else
            d(ct,:) = {c1 TbyGene(c1).gene_i{1} strcat(TbyGene(c1).pos_i(1:6),'M') TbyGene(c1).gene_j{1} strcat(TbyGene(c1).pos_j(1:6),'M') TbyGene(c1).p_val TbyGene(c1).num_events TbyGene(c1).tumor_1{1} TbyGene(c1).tumor_2{1} TbyGene(c1).tumor_3{1}};
        end
        ct=ct+1;
        len_i=length(TbyGene(c1).gene_i);
        len_j=length(TbyGene(c1).gene_j);
        if len_i>1 || len_j>1,
            for c2=2:max(length(TbyGene(c1).gene_i),length(TbyGene(c1).gene_j))
                if c2>len_i
                    d(ct,:) = {[] [] [] TbyGene(c1).gene_j{c2} [] [] [] [] [] []};
                elseif c2>len_j
                    d(ct,:) = {[] TbyGene(c1).gene_i{c2} [] [] [] [] [] [] [] []};
                else
                    d(ct,:) = {[] TbyGene(c1).gene_i{c2} [] TbyGene(c1).gene_j{c2} [] [] [] [] [] []};
                end
                ct=ct+1;
            end
        end
    end


    t.Data = d;
end

%t.BackgroundColor=[0 0 0.8];
%t.Position(3) = t.Extent(3);


            
