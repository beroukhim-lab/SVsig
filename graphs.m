

%sij1dy1=interp1(sij1dx,sij1dy(:,1),sij1dx,'pchip');
%plot(log(sij1dx),log(sij1dy(:,1)),'o',log(sij1dx'),log(sij1dy1),':.');

sij1dy1000 = sum(sij1dy,3);     
sij1dy1000 = sij1dy1000/sum(sij1dy1000);
sij1dy1000interp=interp1(sij1dx,sij1dy1000,sij1dx,'pchip');
plot(log(sij1dx),log(sij1dy1000),'*r',log(sij1dx),log(sij1dy1000interp),'-r');

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




%comparing different sij1dx binning 
times1000=load('/Users/shu/sij1dx1000.mat', 'sij1dx', 'sij1dy');
sij1dy1000=times1000.sij1dy(:,1,:);
sij1dy1000 = sum(sij1dy1000,3);     % summing across types of events (4)
sij1dy1000 = sij1dy1000/sum(sij1dy1000);
sij1dy1000interp=interp1(times1000.sij1dx,sij1dy1000,times1000.sij1dx,'pchip');
plot(log(times1000.sij1dx),log(sij1dy1000),'*y',log(times1000.sij1dx),log(sij1dy1000interp),'-y');
hold on

times100=load('/Users/shu/sij1dx100.mat', 'sij1dx', 'sij1dy');
sij1dy100=times100.sij1dy(:,1,:);
sij1dy100 = sum(sij1dy100,3);     % summing across types of events (4)
sij1dy100 = sij1dy100/sum(sij1dy100);
sij1dy100interp=interp1(times100.sij1dx,sij1dy100,times100.sij1dx,'pchip');
plot(log(times100.sij1dx),log(sij1dy100),'*c',log(times100.sij1dx),log(sij1dy100interp),'-c');
hold on

times50=load('/Users/shu/sij1dx50.mat', 'sij1dx', 'sij1dy');
sij1dy50=times50.sij1dy(:,1,:);
sij1dy50 = sum(sij1dy50,3);     % summing across types of events (4)
sij1dy50 = sij1dy50/sum(sij1dy50);
sij1dy50interp=interp1(times50.sij1dx,sij1dy50,times50.sij1dx,'pchip');
plot(log(times50.sij1dx),log(sij1dy50),'*r',log(times50.sij1dx),log(sij1dy50interp),'-r');
hold on

times25=load('/Users/shu/sij1dx25.mat', 'sij1dx', 'sij1dy');
sij1dy25=times25.sij1dy(:,1,:);
sij1dy25 = sum(sij1dy25,3);     % summing across types of events (4)
sij1dy25 = sij1dy25/sum(sij1dy25);
sij1dy25interp=interp1(times25.sij1dx,sij1dy25,times25.sij1dx,'pchip');
plot(log(times25.sij1dx),log(sij1dy25),'*g',log(times25.sij1dx),log(sij1dy25interp),'-g');
hold on

times10=load('/Users/shu/sijdxsijdy.mat', 'sij1dx', 'sij1dy');
sij1dy10=times10.sij1dy(:,1,:);
sij1dy10 = sum(sij1dy10,3);     % summing across types of events (4)
sij1dy10 = sij1dy10/sum(sij1dy10);
sij1dy10interp=interp1(times10.sij1dx,sij1dy10,times10.sij1dx,'pchip');
plot(log(times10.sij1dx),log(sij1dy10),'*b',log(times10.sij1dx),log(sij1dy10interp),'-b');
hold on





