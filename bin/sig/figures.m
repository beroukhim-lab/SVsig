uniqu_HGG = readtable('/Volumes/xchip_beroukhimlab/ofer/DIPG/matlab/uniuq_HGG_nonhyper_SVTable_05_29_2017.csv');

% frequency of reaarangements per sample
[u_sid, ia_sid, ic_sid] = unique(uniqu_HGG.sid);
num_sv_per_sample = histc(ic_sid,[1:max(ic_sid)]);
[sv_num_val, sv_num_idx] = sort(num_sv_per_sample);
bar(log10(sv_num_val))
set(gca,'XTick',[1:max(ic_sid)],'XTickLabel',u_sid(sv_num_idx),'XTickLabelRotation',90)

% generate morpheus data
sample_annot = table();
sample_annot.sid = u_sid;
sample_annot.numSV = num_sv_per_sample;
