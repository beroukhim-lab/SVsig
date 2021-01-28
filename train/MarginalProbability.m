function R = MarginalProbability(bins_event_tble,events,numbins)

disp ('calculating marginal probabilities...');

%bins_event_tble: col1, col2 are bin indices, col3 is tile index
%R=zeros(numbins,num_annot);
R=zeros(numbins,1);
for c1=1:length(events),
%     R(bins_event_tble(c1,1),(events(c1,3)-1)*2+events(c1,6))=R(bins_event_tble(c1,1),(events(c1,3)-1)*2+events(c1,6))+1;
%     R(bins_event_tble(c1,2),(events(c1,3)-1)*2+events(c1,6))=R(bins_event_tble(c1,2),(events(c1,3)-1)*2+events(c1,6))+1;
    R(bins_event_tble(c1,1))=R(bins_event_tble(c1,1))+1;
    R(bins_event_tble(c1,2))=R(bins_event_tble(c1,2))+1;
end

R = bsxfun(@rdivide,2*R,sum(R,1));


%so this function essentially counts how many breakpoints in each bin,
%divies by total number of events --> breakpoint density 
%so this function shouldn't be called marginal probability right?.. pls
%help
%marginal probabiiliyt should be only distribution of either i or j
%breakpoint densities right


