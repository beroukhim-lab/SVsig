
%mix model
global CHR

% generate background rates for the break invasion model
break_invasion_model

% calculates length factor function
%len_factor = model_length_dist(events00,bins, CHR, sij1dx);

% calculates multiplicative model
double_break_join_model

% train mix model
mfull00=mfull{1}+mfull{2}+mfull{3}+mfull{4};


%find the optimal alpha for combination of background models
[mix_model,opt_alpha] = mix_model_param( mfull00, p, p_mult, events00, bins, CHR);


disp(strcat('optimal alpha', num2str(opt_alpha)));
save('20240719_testmodel.mat')


