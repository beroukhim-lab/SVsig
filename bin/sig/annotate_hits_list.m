function annotated_table = annotate_hits_list( TbyGene_Table,SVTable,bins,hitstable_lookup,pa_mix,qFDR_mix)

if nargin < 6
  qFDR_mix = [];
end

if isempty(TbyGene_Table) || ...
   (istable(TbyGene_Table) && height(TbyGene_Table) == 0)
    error('annotate_hits_list:EmptyInput', ...
        'TbyGene_Table is empty: no significant hits were found to annotate.');
end

if isstruct(TbyGene_Table)

    annotated_table=table();
    chits=table(); ct=1;
    for c1=1:length(TbyGene_Table)
        for c2=1:length(TbyGene_Table(c1).bins)
            clear chits
            hitstable_bin=hitstable_lookup(TbyGene_Table(c1).bins(c2),1:2);

            chri=bins(hitstable_bin(1),1);
            posi=bins(hitstable_bin(1),2:3);
            chrj=bins(hitstable_bin(2),1);
            posj=bins(hitstable_bin(2),2:3);

            SVTable_lines = SVTable.seqnames==chri&SVTable.start>=posi(1)&SVTable.start<=posi(2) & ...
                            SVTable.altchr==chrj&SVTable.altpos>=posj(1)&SVTable.altpos<=posj(2) | ...
                            SVTable.altchr==chri&SVTable.altpos>=posi(1)&SVTable.altpos<=posi(2) & ...
                            SVTable.seqnames==chrj&SVTable.start>=posj(1)&SVTable.start<=posj(2);

            chits=SVTable(SVTable_lines,:);
if isempty(TbyGene_Table(c1).gene_i),
  TbyGene_Table(c1).gene_i = {'none'};
end
if isempty(TbyGene_Table(c1).gene_j),
  TbyGene_Table(c1).gene_j = {'none'};
end
            chits.gene_i = repmat(TbyGene_Table(c1).gene_i(1),height(chits),1);
            chits.gene_j = repmat(TbyGene_Table(c1).gene_j(1),height(chits),1);
            chits.nearby_genes_i = repmat({TbyGene_Table(c1).gene_i},height(chits),1);
            chits.nearby_genes_j = repmat({TbyGene_Table(c1).gene_j},height(chits),1);
            chits.hit_num=c1*ones(sum(SVTable_lines),1);
            chits.pval=str2double(TbyGene_Table(c1).p_val)*ones(sum(SVTable_lines),1);
            chits.tile_num = c2 * ones(sum(SVTable_lines),1);
            chits.u_tile_num = ct * ones(sum(SVTable_lines),1);
            ct=ct+1;
            tile_i = min(hitstable_bin);
            tile_j = max(hitstable_bin);
            chits.bin_i = tile_i * ones(sum(SVTable_lines),1);
            chits.bin_j = tile_j * ones(sum(SVTable_lines),1);
            qval_ij = get_tile_metric(qFDR_mix, tile_i, tile_j, hitstable_lookup, TbyGene_Table(c1).bins(c2));
            chits.tile_qval = qval_ij * ones(sum(SVTable_lines),1);

            pmix_ij = get_tile_metric(pa_mix, tile_i, tile_j, hitstable_lookup, TbyGene_Table(c1).bins(c2));
            chits.p_mix = pmix_ij * ones(sum(SVTable_lines),1);
            annotated_table=[annotated_table;chits];

        end
    end
    
else
    
    annotated_table=table();
    chits=table();
    for c1=1:length(TbyGene_Table)        
            clear chits
            hitstable_bin=TbyGene_Table(c1,:);

            chri=bins(hitstable_bin(1),1);
            posi=bins(hitstable_bin(1),2:3);
            chrj=bins(hitstable_bin(2),1);
            posj=bins(hitstable_bin(2),2:3);

            SVTable_lines = SVTable.seqnames==chri&SVTable.start>=posi(1)&SVTable.start<=posi(2) & ...
                            SVTable.altchr==chrj&SVTable.altpos>=posj(1)&SVTable.altpos<=posj(2) | ...
                            SVTable.altchr==chri&SVTable.altpos>=posi(1)&SVTable.altpos<=posi(2) & ...
                            SVTable.seqnames==chrj&SVTable.start>=posj(1)&SVTable.start<=posj(2);

            chits=SVTable(SVTable_lines,:);
            chits.hit_num = c1 * ones(sum(SVTable_lines),1);
            chits.pval = str2double(TbyGene_Table(c1).p_val) * ones(sum(SVTable_lines),1);
            chits.num_tiles = size(TbyGene_Table(c1).bins,2) * ones(sum(SVTable_lines),1);
            tile_i = min(hitstable_bin);
            tile_j = max(hitstable_bin);
            chits.bin_i = tile_i * ones(sum(SVTable_lines),1);
            chits.bin_j = tile_j * ones(sum(SVTable_lines),1);
            qval_ij = get_tile_metric(qFDR_mix, tile_i, tile_j, hitstable_lookup, []);
            chits.tile_qval = qval_ij * ones(sum(SVTable_lines),1);

            pmix_ij = get_tile_metric(pa_mix, tile_i, tile_j, hitstable_lookup, []);
            chits.p_mix = pmix_ij * ones(sum(SVTable_lines),1);
            annotated_table=[annotated_table;chits];
        
    end
    
end
end

function val = get_tile_metric(metric, tile_i, tile_j, hitstable_lookup, lookup_row)
val = nan;

if isempty(metric)
  return;
end

if isvector(metric)
  if ~isempty(lookup_row) && lookup_row >= 1 && lookup_row <= numel(metric)
    val = metric(lookup_row);
    return;
  end

  if ~isempty(hitstable_lookup) && size(hitstable_lookup,2) >= 2
    is_match = (hitstable_lookup(:,1) == tile_i & hitstable_lookup(:,2) == tile_j) | ...
           (hitstable_lookup(:,1) == tile_j & hitstable_lookup(:,2) == tile_i);
    idx = find(is_match, 1, 'first');
    if ~isempty(idx) && idx <= numel(metric)
      val = metric(idx);
    end
  end
  return;
end

if ndims(metric) >= 2
  if tile_i <= size(metric,1) && tile_j <= size(metric,2)
    val = metric(tile_i, tile_j);
  elseif tile_j <= size(metric,1) && tile_i <= size(metric,2)
    val = metric(tile_j, tile_i);
  end
end

if isempty(val) || numel(val) ~= 1
  val = nan;
end
end
