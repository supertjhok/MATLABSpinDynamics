% Optimize CPMG excitation pulse, precalculate refocusing matrix
% Soumyajit Mandal, 09/21/10
% Allow arbitrary pulse amplitudes, 02/25/11

function [texc,pexc,echo_pk,echo_rms]=opt_exc_pulse_asymp_mag2(nseg,neff,del_w,len_acq,start)

T_90=pi/2; % normalized

delt=0.01;
tvect=-len_acq/2:delt:len_acq/2;
echo=zeros(1,length(tvect));

% Excitation pulse definition
% Use nonlinear function minimization - all segment times must be positive
if isempty(start)
    start=10*T_90*rand(1,nseg); % Random initial condition
end
lb=zeros(1,nseg); % Lower bound
ub=10*T_90*ones(1,nseg); % Upper bound

% trust-region-reflective algorithm (fmincon default) will not work for
% this problem because of the constraints, so use interior-point or
% active-set algorithms instead
options=optimset('Algorithm','interior-point','Display','iter','TolFun',1e-4,'MaxFunEvals',5000);
%options=optimset('Algorithm','active-set','Display','iter','TolFun',1e-4,'MaxFunEvals',5000);

% Segments have alternating 90/270 phase and constant amplitude
aexc=ones(1,nseg);
pexc=(pi/2)*ones(1,nseg);
for i=1:2:nseg
    pexc(i)=3*pi/2;
end

texc=fmincon(@(params)fit_function(params,pexc,aexc,neff,del_w,len_acq,tvect),start,[],[],[],[],lb,ub,[],options);
[masy]=cpmg_van_spin_dynamics_asymp_mag(texc,pexc,aexc,neff,del_w,len_acq);

for i=1:length(tvect)
    echo(i)=sum(masy.*exp(-1i*del_w*tvect(i)));
end

echo_pk=max(abs(echo));
echo_rms=sqrt(trapz(del_w,abs(masy).^2));

function val=fit_function(texc,pexc,aexc,neff,del_w,len_acq,tvect)

echo=zeros(1,length(tvect));
[masy]=cpmg_van_spin_dynamics_asymp_mag(texc,pexc,aexc,neff,del_w,len_acq);

for i=1:length(tvect)
    echo(i)=sum(masy.*exp(-1i*del_w*tvect(i)));
end

%val=-0.01*max(abs(echo)); % Optimize echo peak
val=-sqrt(trapz(del_w,abs(masy).^2))-0.01*max(abs(echo)); % Optimize echo RMS + peak