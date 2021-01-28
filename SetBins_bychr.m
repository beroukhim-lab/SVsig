
function [bins, numbins] = SetBins_bychr(events,chsize,CHR)


for c1 = CHR
    firstbreaks = find(events(:,1)==c1);
    secondbreaks = find(events(:,4)==c1);


    bins(c1,1) = c1;
    bins(c1,2) = 1;   
    bins(c1,3) = chsize(c1);
    bins(c1,4) = length(firstbreaks) + length(secondbreaks);
    
end

numbins=length(bins);
      

