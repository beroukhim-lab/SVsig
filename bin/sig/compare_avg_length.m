%function [h, p_ks, n] = compare_avg_length( tile, bins, events , plot_flag )
function [p_mw, n] = compare_avg_length( tile, bins, events , plot_flag )

% function returns the denisty of events in a SINGLE tile (support)


tile_pos_x=[]; tile_pos_y=[]; tile_rand_x=[]; tile_rand_y=[];
for c1=1:length(tile(:,1))
    tile1(c1)=min(tile(c1,:));
    tile2(c1)=max(tile(c1,:));
    pos1=bins(tile1(c1),1:3);
    pos2=bins(tile2(c1),1:3);
    tile_events=list_events(events(:,1:6),bins(tile1(c1),1:3),bins(tile2(c1),1:3),0,0,[],[],[],[]);
    
       rng('default')
   
       tile_pos_x=[tile_pos_x;tile_events(:,2)];
       tile_pos_y=[tile_pos_y;tile_events(:,5)];      
       %need to subtract from each the index of bin start
       tile_pos_x=tile_pos_x-pos1(2); tile_pos_y= tile_pos_y-pos2(2);
       %real_diff=abs(tile_pos_y - tile_pos_x);
       real_tiles= [tile_pos_x tile_pos_y];
      
       
       tile_rand_x=[tile_rand_x;randi([pos1(2) pos1(3)], 20, 1)];
       tile_rand_y=[tile_rand_y;randi([pos2(2) pos2(3)], 20, 1)];
       tile_rand_x=tile_rand_x-pos1(2); tile_rand_y= tile_rand_y-pos2(2);

       %rand_diff=abs(tile_rand_x - tile_rand_y);
       rand_tiles= [tile_rand_x tile_rand_y];
       
       v = 1:length(real_tiles);
       C = nchoosek(v,2);
       
       real_dist=zeros(length(C), 1);
       for i =1:length(C)-1
            a=C(i,1);
            b=C(i,2);
            real_set=[real_tiles(a,:); real_tiles(b,:)];
            real_dist(i,1)=pdist(real_set,'euclidean');
  
            
       end
       
       v1 = 1:length(rand_tiles);
       C1 = nchoosek(v1,2);
       rand_dist=zeros(length(C1), 1);
       for i =1:length(C1)-1
            a=C1(i,1);
            b=C1(i,2);
            rand_set=[rand_tiles(a,:); rand_tiles(b,:)];
            rand_dist(i,1)=pdist(rand_set,'euclidean');
            
       end
       
       
       %[h, p_ks]=kstest2(real_dist,rand_dist);
        %n=length(real_tiles);
        %p_mw=ranksum(real_dist,rand_dist);
        [n,p_mw]=ttest(real_dist, mean(rand_dist));
        %n=1;


end
plot_flag=0;

%test_bin1=bins((bins(:,1)==1) & (bins(:,3) >= 26608055) & (bins(:,2) <= 26608055))
%test_bin2=bins((bins(:,1)==1) & (bins(:,3) >= 23) & (bins(:,2) <= 23))
if (pos1(1)==12 & pos2(1)==12 & pos1(2)<69373924 & pos1(3)>69373924 & pos2(2)<85635019 & pos2(3)>85635019)  
    x=1;
%end
if plot_flag>0
    figure
    hold on

    xlabel(['chr' num2str(bins(tile1(1),1)), ':', num2str(pos1(2)), '-', num2str(pos1(3))])
    ylabel(['chr' num2str(bins(tile2(1),1)), ':', num2str(pos2(2)), '-', num2str(pos2(3))])
    
    scatter(real_tiles(:,1),real_tiles(:,2),'ko');
    axis([0 pos1(3)-pos1(2) 0 pos2(3)-pos2(2)])
    %scatter(rand_tiles(:,1),rand_tiles(:,2),'rx');
    %plot(round([dr_xy dr_xy]/x_res),[1 grid0],'w');
    %plot(grid0-round([dr_xy dr_xy]/x_res),[1 grid0],'w');

    title(['p-value: ',num2str(p_ks), ', n=', num2str(length(real_tiles), '%.3g') ])
end


end

