function sij1dx = length_dist_1d_bins(events,chsize,num_1d_bins)


%diffe is the distance between the two breakpoints for each event, sorted
diffe=sort(abs(events(events(:,1)==events(:,4),5)-events(events(:,1)==events(:,4),2)));

nume=length(diffe);
bins_loc=round(linspace(1,nume,num_1d_bins));

sij1dx = diffe(bins_loc)';
sij1dx(end)=max(chsize);

end