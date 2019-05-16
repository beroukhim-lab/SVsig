function len_factor = model_length_dist( events,bins, CHR, sij1dx)
%@param CHR is chromosome
%@param sij1dx is one dimensional sij on the x-axis
%@param bins is a list of the start and end positions of the bins
%@param events is the list of events (i.e. junctions)

% fraction of events per length interval 
sij1dy0 = EventLengthDist_G(sij1dx,events,0);
sij1dy=sum(sij1dy0,2);

% since sij1dx was choosen such that there are equal number of events I assume here equal distribution  
d_sij1dx=diff(sij1dx)';
sd_sij1dx(1)=d_sij1dx(1)/2;
for c1=2:length(d_sij1dx),
    sd_sij1dx(1,c1)=sd_sij1dx(c1-1)+sum(d_sij1dx(c1-1:c1))/2;
end
sij1_area=sij1dy(1:end-1).*d_sij1dx;

for c1=CHR,
    numbins_chr(c1)=sum(bins(:,1)==c1);
end

numbins=length(bins);
len_factor=sparse(sum(numbins_chr.*(numbins_chr+1))/2,length(d_sij1dx));
ct=1;
for c1=CHR,
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    for c2=firstbin:lastbin,
        diag_bin = sum(bins(c2,2:3),2)/2;        
        diag_bin_size = (bins(c2,3)-bins(c2,2))/2;
        sij1dx_b=sij1dx<diag_bin_size;
        last_diag = find(sij1dx_b,1,'last');
        len_factor(ct,:) = sij1dx_b(2:end);
        len_factor(ct,last_diag) = (diag_bin_size - sij1dx(last_diag))/(sij1dx(last_diag+1)-sij1dx(last_diag));
        len_factor(ct,:) = len_factor(ct,:)/sum(len_factor(ct,:));
        ct=ct+1;
        for c3=c2+1:lastbin,
            upper_diag_bins = bins(c3,2:3)-diag_bin;
            sij1dx_b = sij1dx<upper_diag_bins(2)&sij1dx>=upper_diag_bins(1);
            if sum(sij1dx_b)==0, % both of upper_diag_bins falls within one sij1dx bin
                sx_bin=find(sij1dx>upper_diag_bins(1),1);
%                len_factor(ct,sx_bin-1)=diff(upper_diag_bins)/d_sij1dx(sx_bin-1);'
                len_factor(ct,sx_bin-1)=1;
            else % ... cover more than one sij1dx bins
                first_sx=find(sij1dx_b,1);
                last_sx=find(sij1dx_b,1,'last');
                len_factor(ct,first_sx-1) = (sij1dx(first_sx)-upper_diag_bins(1))/d_sij1dx(first_sx-1);
                len_factor(ct,last_sx) = (upper_diag_bins(2)-sij1dx(last_sx))/d_sij1dx(last_sx);
                len_factor(ct,first_sx:last_sx-1) = 1;
%                len_factor(ct,first_sx:last_sx-1) = len_factor(ct,first_sx:last_sx-1)/sum(len_factor(ct,first_sx:last_sx-1));
                len_factor(ct,:) = len_factor(ct,:)/sum(len_factor(ct,:));
            end
            ct=ct+1;
        end
    end
end
        
        
return