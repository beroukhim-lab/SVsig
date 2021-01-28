function [normalized_mat] = renormalize_tiles(mat, short_ratio, long_ratio, inter_ratio, events, bins, CHR)

%ic ratio is ratio of inter:intra events
%intra ratio is ratio of short:long events 

annot_tiles=tiles_annot('length',events,bins,CHR);

normalized_mat=zeros(size(mat,1),size(mat,2));

short_annot = (sum(sum(mat(annot_tiles(:,:,1)))))/(sum(sum(mat)));
long_annot = (sum(sum(mat(annot_tiles(:,:,2)))))/(sum(sum(mat)));
inter_annot= (sum(sum(mat(annot_tiles(:,:,3)))))/(sum(sum(mat)));
%sum_mat=sum(sum(mat));

%short_ratio=short_annot/(short_annot + long_annot);
%long_ratio=short_annot/(short_annot + long_annot);
%inter_ratio=short_annot/(short_annot + long_annot);

%
%normalized_mat(annot_tiles(:,:,1)) = (short_ratio/short_annot)*mat(annot_tiles(:,:,1));
%normalized_mat(annot_tiles(:,:,2)) = (long_ratio/long_annot)*mat(annot_tiles(:,:,2));
%normalized_mat(annot_tiles(:,:,3)) = (inter_ratio/inter_annot)*mat(annot_tiles(:,:,3));

%mat(firstbin:lastbin,1:firstbin-1) = ((1 - intra_ratio)/ratio)*p_mult(firstbin:lastbin, 1:firstbin-1);

%new_short = (sum(sum(normalized_mat(annot_tiles(:,:,1)))))/(sum(sum(normalized_mat)));
%new_long = (sum(sum(normalized_mat(annot_tiles(:,:,2)))))/(sum(sum(normalized_mat)));
%new_inter = (sum(sum(normalized_mat(annot_tiles(:,:,3)))))/(sum(sum(normalized_mat)));

%disp(['old short, long, inter ratios = ' num2str(short_annot) ' , ' num2str(long_annot) ' , ' num2str(inter_annot)]);

%disp(['new short, long, inter ratios = ' num2str(new_short) ' , ' num2str(new_long) ' , ' num2str(new_inter)]);

normalized_mat(annot_tiles(:,:,1)) = (short_ratio/short_annot)*mat(annot_tiles(:,:,1));
normalized_mat(annot_tiles(:,:,2)) = (long_ratio/long_annot)*mat(annot_tiles(:,:,2));
normalized_mat(annot_tiles(:,:,3)) = (inter_ratio/inter_annot)*mat(annot_tiles(:,:,3));


new_short = (sum(sum(normalized_mat(annot_tiles(:,:,1)))))/(sum(sum(normalized_mat)));
new_long = (sum(sum(normalized_mat(annot_tiles(:,:,2)))))/(sum(sum(normalized_mat)));
new_inter = (sum(sum(normalized_mat(annot_tiles(:,:,3)))))/(sum(sum(normalized_mat)));


disp(['old short, long, inter ratios = ' num2str(short_annot) ' , ' num2str(long_annot) ' , ' num2str(inter_annot)]);

disp(['new short, long, inter ratios = ' num2str(new_short) ' , ' num2str(new_long), ' , ' num2str(new_inter)]);






