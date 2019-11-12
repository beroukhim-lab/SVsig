function [output, Uevent, Usample, Upatient, UTumor, Ustrand1, Ustrand2] = GenerateSVarray(input,length_th,CHR,Tumor_column,Event_column,Sample_column,Patient_column,Weights_column)
global samp_num

%Tumor_column = 7;
%Sample_column = 10;

%reading in SV file and filling in output events array with appropriate
%columns
%need to have an extra column with the weight of each arrangement 
%ultimately multiply binned events by the weights?


disp('generating events array...');

snowman_table=readtable(input);



[Ustrand1, ia_strand1, ic_strand1]=unique(snowman_table(:,3));
output(:,3)=ic_strand1;
[Ustrand2, ia_strand2, ic_strand2]=unique(snowman_table(:,6));
output(:,6)=ic_strand2;
output(:,2)=table2array(snowman_table(:,2));
output(:,5)=table2array(snowman_table(:,5));

if ~isempty(Tumor_column)
    [UTumor, ia_code, ic_code]=unique(snowman_table(:,Tumor_column));
    output(:,7)=ic_code;
end
if ~isempty(Event_column)
    [Uevent, ia_event, ic_event]=unique(snowman_table(:,Event_column));
    output(:,8)=ic_event;
end
if ~isempty(Sample_column)
    [Usample, ia_sample, ic_sample]=unique(snowman_table(:,Sample_column));
    output(:,9)=ic_sample;
end
if ~isempty(Patient_column)
    [Upatient, ia_patient, ic_patient]=unique(snowman_table(:,Patient_column));
    output(:,10)=ic_patient;
end
if ~isempty(Weights_column)
    %[UWeights, ia_weights, ic_weights] = unique(snowman_table(:,Weights_column));
    %output(:, 11) = ic_weights;
    %we want weights to be actual fractions rather than coded by an integer
    output(:,11) = snowman_table.weights;
end

if size(snowman_table,2)>11
    [Utopo, ia_topo, ic_topo]=unique(snowman_table.topo);
    output(:,12)=ic_topo;
    output(:,13)=snowman_table.topo_n;
    [Umech, ia_mech, ic_mech]=unique(snowman_table.mech);
    output(:,14)=ic_mech;
    output(:,15)=snowman_table.homseq;
end

%choose samples to include if subsetting and only keep those in output

if ~isempty(samp_num)
 %randomly select patients WITHOUT REPLACEMENT
idx = randsample(size(Upatient, 1), samp_num, false);
%subset the output and snowman_table    
idxl = ismember(output(:,10), idx); %find rows
output = output(idxl, :);
snowman_table = snowman_table(idxl, :);
end 

if 0 % changed 10/30/2016
    snowman=table2struct(snowman_table);
    output(:,1) = chr_str2num({snowman.seqnames});
    output(:,4) = chr_str2num({snowman.altchr});
else
%    output(:,1) = chr_str2num(snowman_table{:,1});
%    output(:,4) = chr_str2num(snowman_table{:,4});
    output(:,1) = snowman_table{:,1};
    output(:,4) = snowman_table{:,4};
end

%remove lines with events smaller than length_th bp
len_thr_idx=(output(:,1)==output(:,4))&(abs(output(:,2)-output(:,5))<length_th);
output(len_thr_idx,:)=[];

%remove lines with chromosomes out of range
for c1=1:24
    if sum(CHR==c1)==0
        chr_idx=(output(:,1)==c1 | output(:,4)==c1);
        output(chr_idx,:)=[];
    end
    
    
 
end

%save icgc_snowman output Ucode Ustrand1 Ustrand2
    