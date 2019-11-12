
%mix model


% generate background rates for the break invasion model
break_invasion_model

% calculates length factor function
len_factor = model_length_dist(events00,bins, CHR, sij1dx);

% calculates multiplicative model
double_break_join_model

% train mix model
mfull00=mfull{1}+mfull{2}+mfull{3}+mfull{4};

%find the optimal alpha for combination of break invasion and double break
%join models
[mix_model,opt_alpha] = mix_model_param( mfull00, p, p_mult, events00, bins, CHR );

opt_alpha

%save('model1model220191011.mat')
% p-values for mix model
%[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull00, mix_model, [], [], 1);
%[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull00, mix_model, bins, events, sij1dx, chsize, CHR, [], [], 1);

%do this 20 times to check for stability and BIC value% 



