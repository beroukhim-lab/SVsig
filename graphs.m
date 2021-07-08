

%sij1dy1=interp1(sij1dx,sij1dy(:,1),sij1dx,'pchip');
%plot(log(sij1dx),log(sij1dy(:,1)),'o',log(sij1dx'),log(sij1dy1),':.');

sij1dy1000 = sum(sij1dy,3);     
sij1dy1000 = sij1dy1000/sum(sij1dy1000).*(100/10);
sij1dy1000interp=interp1(sij1dx,sij1dy1000,sij1dx,'pchip');
plot(log(sij1dx),log(sij1dy1000),'*k',log(sij1dx(1:99)),log(sij1dy1000interp(1:99)),'-k');
%plot((sij1dx),(sij1dy1000),'*m',(sij1dx(1:99)),(sij1dy1000interp(1:99)),'-m');
hold on

sij1dy50=times50.sij1dy(:,1,:);
for ca=1:4
    sij1dy50(:,:,ca) = sij1dy50(:,:,ca) * annot_frac(ca);
end
sij1dy50 = sum(sij1dy50,3);     % summing across types of events (4)
sij1dy50 = (sij1dy50/sum(sij1dy50)).*(50/5);
sij1dy50interp=interp1(times50.sij1dx,sij1dy50,times50.sij1dx,'pchip');
plot(log(times50.sij1dx),log(sij1dy50),'*r',log(times50.sij1dx(1,1:49)),log(sij1dy50interp(1,1:49)),'-r');
hold on

hold on

times100=load('/Users/shu/220210121_eventratios_withmarg_a0sijdx.mat', 'sij1dx', 'sij1dy');
sij1dy100=times100.sij1dy(:,1,:);
sij1dy100 = sum(sij1dy100,3);     % summing across types of events (4)
sij1dy100 = sij1dy100/sum(sij1dy100);
sij1dy100interp=interp1(times100.sij1dx,sij1dy100,times100.sij1dx,'pchip');
plot(log(times100.sij1dx),log(sij1dy100),'*b',log(times100.sij1dx),log(sij1dy100interp),'-b');

hold on

xline(log(1e6))

sij1dyten=tentimes.sij1dy(:,1,:);
events=events00;
nume=length(events);
annot_frac(1) = sum(events(:,3)==1&events(:,6)==1)/nume;
annot_frac(2) = sum(events(:,3)==1&events(:,6)==2)/nume;
annot_frac(3) = sum(events(:,3)==2&events(:,6)==1)/nume;
annot_frac(4) = sum(events(:,3)==2&events(:,6)==2)/nume;
for ca=1:4
    sij1dyten(:,:, ca) = sij1dyten(:,:, ca) * annot_frac(ca);
end
sij1dyten;
sij1dyten = sum(sij1dyten,3);     % summing across types of events (4)
sij1dyten = sij1dyten/sum(sij1dyten) 

sij1dy1teninterpavg=interp1(tentimes.sij1dx,sij1dyten,tentimes.sij1dx,'pchip');
plot(log(tentimes.sij1dx),log(sij1dyten),'*b',log(tentimes.sij1dx),log(sij1dy1teninterpavg),'-b');



hist1=histogram(diffe);
hi=hist1.BinEdges;
bye=hist1.Values;
plot(log(hi(1:238)),log(bye))




hold on

% annot_frac(1) = sum(events00(:,3)==1&events00(:,6)==1)/289095;
% annot_frac(2) = sum(events00(:,3)==1&events00(:,6)==2)/289095;
% annot_frac(3) = sum(events00(:,3)==2&events00(:,6)==1)/289095;
% annot_frac(4) = sum(events00(:,3)==2&events00(:,6)==2)/289095;
annot_frac=[0.181729189366817,0.337961569726215,0.300223110050330,0.180086130856639]


%comparing different sij1dx binning 
times1000=load('/Users/shu/sij1dx1000.mat', 'sij1dx', 'sij1dy');
%sij1dy1000=times1000.sij1dy(:,1,:);
sij1dy1000=sij1dy(:,1,:);
for ca=1:4
    sij1dy1000(:,:,ca) = sij1dy1000(:,:,ca) * annot_frac(ca);
end
sij1dy1000 = sum(sij1dy1000,3);     % summing across types of events (4)
sij1dy1000 = (sij1dy1000/sum(sij1dy1000)).*(1);
sij1dy1000interp=interp1(sij1dx,sij1dy1000,sij1dx,'pchip');
plot(log(sij1dx),log(sij1dy1000),'*k',log(sij1dx(1,1:99)),log(sij1dy1000interp(1,1:99)),'-k');
hold on

times100=load('/Users/shu/sij1dx100.mat', 'sij1dx', 'sij1dy');
sij1dy100=times100.sij1dy(:,1,:);
for ca=1:4
    sij1dy100(:,:,ca) = sij1dy100(:,:,ca) * annot_frac(ca);
end
sij1dy100 = sum(sij1dy100,3);     % summing across types of events (4)
sij1dy100 = (sij1dy100/sum(sij1dy100)).*(100/100);
sij1dy100interp=interp1(times100.sij1dx,sij1dy100,times100.sij1dx,'pchip');
plot(log(times100.sij1dx),log(sij1dy100),'*c',log(times100.sij1dx(1,1:99)),log(sij1dy100interp(1,1:99)),'-c');
hold on

times50=load('/Users/shu/sij1dx50.mat', 'sij1dx', 'sij1dy');
sij1dy50=times50.sij1dy(:,1,:);
for ca=1:4
    sij1dy50(:,:,ca) = sij1dy50(:,:,ca) * annot_frac(ca);
end
sij1dy50 = sum(sij1dy50,3);     % summing across types of events (4)
sij1dy50 = (sij1dy50/sum(sij1dy50)).*(50/10);
sij1dy50interp=interp1(times50.sij1dx,sij1dy50,times50.sij1dx,'pchip');
plot(log(times50.sij1dx),log(sij1dy50),'*r',log(times50.sij1dx(1,1:49)),log(sij1dy50interp(1,1:49)),'-r');
hold on

times25=load('/Users/shu/sij1dx25.mat', 'sij1dx', 'sij1dy');
sij1dy25=times25.sij1dy(:,1,:);
for ca=1:4
    sij1dy25(:,:,ca) = sij1dy25(:,:,ca) * annot_frac(ca);
end
sij1dy25 = sum(sij1dy25,3);     % summing across types of events (4)
sij1dy25 = (sij1dy25/sum(sij1dy25)).*(25/10);
sij1dy25interp=interp1(times25.sij1dx,sij1dy25,times25.sij1dx,'pchip');
plot(log(times25.sij1dx),log(sij1dy25),'*g',log(times25.sij1dx(1,1:24)),log(sij1dy25interp(1,1:24)),'-g');
hold on

times10=load('/Users/shu/sijdxsijdy.mat', 'sij1dx', 'sij1dy');
sij1dy10=times10.sij1dy(:,1,:);
for ca=1:4
    sij1dy10(:,:,ca) = sij1dy10(:,:,ca) * annot_frac(ca);
end

sij1dy10 = sum(sij1dy10,3);     % summing across types of events (4)
sij1dy10 = (sij1dy10/sum(sij1dy10)).*(1);
sij1dy10interp=interp1(times10.sij1dx,sij1dy10,times10.sij1dx,'pchip');
plot(log(times10.sij1dx),log(sij1dy10),'*c',log(times10.sij1dx),log(sij1dy10interp),'-c');
hold on

xline(log(1e6))






