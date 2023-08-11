% length distribution per direction
function sij1dy = EventLengthDist_G_testlength(sij1dx,events,plot_flag, length_factor)

% returns genome-wide avergae number of events per bp by length for each SV-type

sij1dypp=zeros(length(sij1dx),1);
sij1dy=zeros(length(sij1dx)-1,4);

DIRECTpp = find(events(:,1)==events(:,4) & events(:,3)==1 & events(:,6)==1);
DIRECTpm = find(events(:,1)==events(:,4) & events(:,3)==1 & events(:,6)==2);
DIRECTmp = find(events(:,1)==events(:,4) & events(:,3)==2 & events(:,6)==1);
DIRECTmm = find(events(:,1)==events(:,4) & events(:,3)==2 & events(:,6)==2);
    
if ~isempty(DIRECTpp)
    sij1dypp(:,1) = histc(abs(events(DIRECTpp,2)-events(DIRECTpp,5)),sij1dx);
%    sij1dy(:,1) = sij1dypp(1:end-1)./diff(sij1dx)'/sum(sij1dypp(1:end-1));
    sij1dy(:,1) = sij1dypp(1:end-1)./diff(sij1dx)';

end
if ~isempty(DIRECTpm)
    sij1dypm(:,1) = histc(abs(events(DIRECTpm,2)-events(DIRECTpm,5)),sij1dx);
%    sij1dy(:,2) = sij1dypm(1:end-1)./diff(sij1dx)'/sum(sij1dypm(1:end-1));
    sij1dy(:,2) = sij1dypm(1:end-1)./diff(sij1dx)';

end
if ~isempty(DIRECTmp)
    sij1dymp(:,1) = histc(abs(events(DIRECTmp,2)-events(DIRECTmp,5)),sij1dx);
%    sij1dy(:,3) = sij1dymp(1:end-1)./diff(sij1dx)'/sum(sij1dymp(1:end-1));
    sij1dy(:,3) = sij1dymp(1:end-1)./diff(sij1dx)';

end
if ~isempty(DIRECTmm)
    sij1dymm(:,1) = histc(abs(events(DIRECTmm,2)-events(DIRECTmm,5)),sij1dx);
%    sij1dy(:,4) = sij1dymm(1:end-1)./diff(sij1dx)'/sum(sij1dymm(1:end-1));
    sij1dy(:,4) = sij1dymm(1:end-1)./diff(sij1dx)';

end

sij1dy(end+1,:)=0;



%%%code to test robustness of length distribution
sij1dy=sij1dy.^(length_factor);

if plot_flag
   sij1dy_0=sum(sij1dy,2);
   sij1dy_n=sij1dy_0./sum(sij1dy_0(1:end-1).*diff(sij1dx'));
   loglog(sij1dx,sij1dy_n)
end