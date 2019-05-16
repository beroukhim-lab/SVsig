function density = tile_density( tile, bins, events , sij1dx, e_density, plot_flag )
% function returns the denisty of events in a SINGLE tile (support)

% uncertainty parameters
dr_xy=4e5; % default 4e5
h_prctile=0; %default 0
n_std=0.70; % default 0.7
dr_res=10000; % resilution of histogramp in bp

d_std=dr_xy/2;
dist_res=round(dr_xy/dr_res); 
%dist_fw=1/4; % default 1/2
min_density=4e8;

% list of unique events
breakpoint=[events(:,1:2);events(:,4:5)];
tile_pos_x=[]; tile_pos_y=[]; marg_pos_x=[]; marg_pos_y=[];
for c1=1:length(tile(:,1))
    tile1(c1)=min(tile(c1,:));
    tile2(c1)=max(tile(c1,:));
    pos1=bins(tile1(c1),1:3);
    pos2=bins(tile2(c1),1:3);
    tile_events=list_events(events(:,1:6),bins(tile1(c1),1:3),bins(tile2(c1),1:3),0,0,[],[],[],[]);
    if bins(tile(c1,1),1)==bins(tile(c1,2),1) % same chromosome
       tile_pos_x=[tile_pos_x;tile_events(:,2)];
       tile_pos_y=[tile_pos_y;tile_events(:,5)];
    else
       tile_pos_x=[tile_pos_x;tile_events(:,2).*(tile_events(:,1)<tile_events(:,4))+tile_events(:,5).*(tile_events(:,1)>tile_events(:,4))];
       tile_pos_y=[tile_pos_y;tile_events(:,5).*(tile_events(:,1)<tile_events(:,4))+tile_events(:,2).*(tile_events(:,1)>tile_events(:,4))];
    end
    marg_pos_x=[marg_pos_x;breakpoint(breakpoint(:,1)==pos1(:,1)&breakpoint(:,2)>=pos1(:,2)-d_std&breakpoint(:,2)<pos1(:,3)+d_std,2)];
    marg_pos_y=[marg_pos_y;breakpoint(breakpoint(:,1)==pos2(:,1)&breakpoint(:,2)>=pos2(:,2)-d_std&breakpoint(:,2)<pos2(:,3)+d_std,2)];  
    
    % density factor = empirical density/expected density;
    exp_density=interp1(sij1dx,e_density,mean(abs(tile_events(:,5)-tile_events(:,2))));
    density_factor(c1)=length(tile_events)*1e12/(bins(tile1(c1),3)-bins(tile1(c1),2))/(bins(tile2(c1),3)-bins(tile2(c1),2))/exp_density;    
end
density_factor=mean(density_factor);
if density_factor>1000
    density_factor=1000;
end
if density_factor<0.001
    density_factor=0.001;
end

num_events=length(tile_pos_x);



% remove background of marginal distribution by removing the 10% most
% sparse breakpoints
% s_marg_pos_x=sort(marg_pos_x);
% pad_s_marg_pos_x=[s_marg_pos_x(2);s_marg_pos_x;s_marg_pos_x(end-1)];
% d_marg_pos_x=abs(pad_s_marg_pos_x(2:end-1)-pad_s_marg_pos_x(1:end-2))+abs(pad_s_marg_pos_x(2:end-1)-pad_s_marg_pos_x(3:end));
% s_marg_pos_x(d_marg_pos_x>prctile(d_marg_pos_x,50))=[];

%tile area
tile_x=max(bins(tile1,3))-min(bins(tile1,2));
tile_y=max(bins(tile2,3))-min(bins(tile2,2));
tile_area=(bins(tile1,3)-bins(tile1,2))'*(bins(tile2,3)-bins(tile2,2));

tile_x=tile_x+2*dr_xy;
tile_y=tile_y+2*dr_xy;

%density
grid0=4000;
x_res=round(tile_x/grid0);
y_res=round(tile_y/grid0);
tile_mat=false(grid0,grid0);

convhull_c=[];
for c1=1:num_events,    
    
    ir_x=abs(marg_pos_x-tile_pos_x(c1))<d_std;
    ir_y=abs(marg_pos_y-tile_pos_y(c1))<d_std;
    
    [marg_xh,marg_xc]=hist((marg_pos_x(ir_x)-tile_pos_x(c1)),dist_res);     
    [marg_yh,marg_yc]=hist((marg_pos_y(ir_y)-tile_pos_y(c1)),dist_res);
    marg_xh=marg_xh-prctile(marg_xh,h_prctile);
    marg_xh(marg_xh<0)=0;
    marg_yh=marg_yh-prctile(marg_yh,h_prctile);
    marg_yh(marg_yh<0)=0;
     
    %Kiran added: control for the case where the sum of the marg_xh and
    %marg_yh is zero
    
    if sum(marg_xh) ~= 0 && sum(marg_yh) ~= 0 
    marg_xh=marg_xh/sum(marg_xh);
    marg_yh=marg_yh/sum(marg_yh);

    elseif sum(marg_xh)==0
    
    marg_xh = 0;
    
    
    
    end
   
    if sum(marg_yh)==0 
    
    marg_yh = 0;
    
    end
    % end of addition
    
    
    
    
    
    
    
    
    std_x = sqrt(sum(marg_xc.^2.*marg_xh))*n_std;
    std_y = sqrt(sum(marg_yc.^2.*marg_yh))*n_std;

    if std_x>d_std, 
        std_x=d_std;
    end
    if std_y>d_std, 
        std_y=d_std;
    end
%    [std_x std_y]    
%     std_x=std(marg_pos_x(ir_x)-tile_pos_x(c1));
%     std_y=std(marg_pos_y(ir_y)-tile_pos_y(c1));    
    
%     [marg_xh,marg_xc]=hist(abs(marg_pos_x(ir_x)-tile_pos_x(c1)),dist_res);     
%     [marg_yh,marg_yc]=hist(abs(marg_pos_y(ir_y)-tile_pos_y(c1)),dist_res);
%     s_marg_xh=sort(marg_xh);
%     s_marg_yh=sort(marg_yh);
%     for c2=1:length(marg_xh)-1,
%         if s_marg_xh(c2+1)<=s_marg_xh(c2),
%             s_marg_xh(c2+1)=s_marg_xh(c2)+1e-10;
%         end
%         if s_marg_yh(c2+1)<=s_marg_yh(c2),
%             s_marg_yh(c2+1)=s_marg_yh(c2)+1e-10;
%         end
%     end
%     
%     if max(marg_xh(:))*dist_fw>min(marg_xh(:)),
%         std_x=interp1(s_marg_xh',[dist_res-1:-1:0]*(marg_xc(2)-marg_xc(1))',max(marg_xh(:))*dist_fw);
%     else
%         std_x=d_std;
%     end
%     if max(marg_yh(:))*dist_fw>min(marg_yh(:)),
%         std_y=interp1(s_marg_yh',[dist_res-1:-1:0]*(marg_yc(2)-marg_yc(1))',max(marg_yh(:))*dist_fw);
%     else
%         std_y=d_std;
%     end
%     
    stdx_t=sqrt(d_std^2/(sum(ir_x)+1)+std_x^2);
    stdy_t=sqrt(d_std^2/(sum(ir_y)+1)+std_y^2);
    
    nnxy=sum(abs(tile_pos_x-tile_pos_x(c1))<stdx_t&abs(tile_pos_y-tile_pos_y(c1))<stdy_t);
    if nnxy==0,
        nnxy=1;
    end
    e_x=stdx_t/sqrt(nnxy);
    e_y=stdy_t/sqrt(nnxy);

    if plot_flag==2,
        sprintf('%u ',[std_x std_y stdx_t stdy_t nnxy e_x e_y])
    end
    
    x_coor(c1)=tile_pos_x(c1)-min(bins(tile1,2))+dr_xy;
    y_coor(c1)=tile_pos_y(c1)-min(bins(tile2,2))+dr_xy;
    
    x_c1=round((x_coor(c1)-e_x)/x_res)+1; x_c2=round((x_coor(c1)+e_x)/x_res);
    y_c1=round((y_coor(c1)-e_y)/y_res)+1 ; y_c2=round((y_coor(c1)+e_y)/y_res);
    x_range=x_c1:x_c2;
    y_range=y_c1:y_c2;
    
    convhull_c=[convhull_c; x_c1 y_c1; x_c1 y_c2; x_c2 y_c1; x_c2 y_c2];
       
    tile_mat(y_range,x_range)=true;
%    [sum(ir_x) sum(ir_y) std_x std_y nnxy]
end

[K,cluster_area] = convhull(convhull_c);
%density=[sum(tile_mat(:))/((grid0-floor(2*dr_xy/x_res))*(grid0-floor(2*dr_xy/y_res)))*tile_area num_events];
%cluster_area=sum(tile_mat(:))*x_res*y_res;
cluster_area=cluster_area*x_res*y_res;
if cluster_area>min_density,
    density=[cluster_area num_events];
else
    density=[min_density num_events];
end

if plot_flag>0,
    figure
    imagesc(tile_mat)
    hold on
    p_h=patch(convhull_c(K,1),convhull_c(K,2),[0.8 0.8 0.8]);
    set(p_h,'FaceAlpha',0.8)
    xlabel(['chr' num2str(bins(tile1(1),1))])
    ylabel(['chr' num2str(bins(tile2(1),1))])
    scatter(round(x_coor/x_res),round(y_coor/y_res),'kx');
    plot(round([dr_xy dr_xy]/x_res),[1 grid0],'w');
    plot(grid0-round([dr_xy dr_xy]/x_res),[1 grid0],'w');
    plot([1 grid0],round([dr_xy dr_xy]/y_res),'w');
    plot([1 grid0],grid0-round([dr_xy dr_xy]/y_res),'w');
    xtick_pos=get(gca,'XTick');
    for c1=1:length(xtick_pos),
        xtick_label{c1}=sprintf('%2.2fM',((xtick_pos(c1)-1)*x_res+(min(bins(tile1,2))-dr_xy))/1e6);
    end
    set(gca,'XTickLabel',xtick_label)
    ytick_pos=get(gca,'YTick');
        for c1=1:length(xtick_pos),
        ytick_label{c1}=sprintf('%2.2fM',((ytick_pos(c1)-1)*y_res+(min(bins(tile2,2))-dr_xy))/1e6);
    end
    set(gca,'YTickLabel',ytick_label)
    colormap(summer(256))
end

end

