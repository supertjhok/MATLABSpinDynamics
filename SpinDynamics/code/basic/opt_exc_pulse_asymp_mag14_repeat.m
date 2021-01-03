% Repeated optimization runs, store current maximum
% Gradually transform refocusing pulse (tref,aref) from rectangular(180) to
% rectangular (135)
% Use Colm's spin dynamics code for greater speed

function opt_exc_pulse_asymp_mag14_repeat(tref,pref,aref,delt)

delta=0.01; % Shortening of reversed-phase segments at every step
pulse_num=28; % refocusing pulse = rectangular(180)

% Load initial excitation pulse
tmp=load('dat_files\results_mag13.mat');
results_sort=tmp.results;
out.texc=results_sort{pulse_num,1}; out.pexc=results_sort{pulse_num,2};
out.tref=tref; out.pref=pref; out.aref=aref;

count=1;
params.delt=delt;
params.delta=delta;

results={};
while(1)
    params.texc=out.texc;
    params.pexc=out.pexc;
    params.tref=out.tref;
    params.pref=out.pref;
    params.aref=out.aref;
    
    [out]=opt_exc_pulse_asymp_mag14(params);
        
    results{count,1}=out.texc;
    results{count,2}=out.pexc;
    results{count,3}=out.tref;
    results{count,4}=out.pref;
    results{count,5}=out.echo_pk;
    results{count,6}=out.echo_rms;
    save results_mag14.mat results
    
    % Shorten excitation pulse and repeat optimization
    count=count+1;
    
    disp(out.echo_pk)
    disp(out.echo_rms)
    disp(out.tref)
end