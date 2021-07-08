
%mix model
global CHR

% generate background rates for the break invasion model
%break_invasion_model
break_invasion_copy

% calculates length factor function
%len_factor = model_length_dist(events00,bins, CHR, sij1dx);

% calculates multiplicative model
double_break_join_model

% train mix model
mfull00=mfull{1}+mfull{2}+mfull{3}+mfull{4};


%p= pcawg_background.p;
%p_mult= pcawg_background.p_mult;
%find the optimal alpha for combination of break invasion and double break
%join models
[mix_model,opt_alpha] = mix_model_param( mfull00, p, p_mult, events00, bins, CHR );

opt_alpha

%save('model1model220191011.mat')


