% Optimize CPMG excitation pulse, precalculate refocusing matrix
% Soumyajit Mandal, 09/21/10
% --------------------------------------------------------------
% Allow arbitrary pulse amplitudes, 02/25/11
% Use VAN_EXC as a starting point, 08/22/11
% Gradually make excitation pulse shorter -> fewer segments (09/12/11),
% shorter segments (09/13/11)
% Gradually make excitation pulse longer -> more segments (09/16/11)
% Use Colm's optimized code for speed -> 09/19/11.
% Gradually make excitation pulse longer -> stretch segments (03/29/12)
% Gradually make excitation pulse shorter -> reduce number of segments
% (04/06/12)
% --------------------------------------------------------------
% Using Colm's code only ~17% of the time is being spent in calculating the
% fitness function, i.e., the spin dynamics (for 100-segment long
% excitation pulses and 10 echoes). The rest is fmincon overhead, so
% further optimization of the spin dynamics code won't speed up the
% optimization process. Net speedup relative to my Matlab code is ~2.5x,
% which is far less than the 8x - 10x speedup of the spin dynamics itself.
% In order to get more speedup the optimization routine has to be written
% in C++ as well!
% --------------------------------------------------------------

function [out]=opt_exc_pulse_asymp_mag14(params)

texc=params.texc;
pexc=params.pexc;
tref=params.tref;
pref=params.pref;
aref=params.aref;
len_acq=params.delt;
delta=params.delta;

T_90=pi/2; % normalized
delt=0.01*T_90;
tvect=-len_acq/2:delt:len_acq/2;
echo=zeros(1,length(tvect));

tref(1)=tref(1)+delta*T_90; % Lengthen free precession periods to keep TE constant
tref(3)=tref(3)+delta*T_90;
tref(2)=tref(2)-2*delta*T_90; % Shorten refocusing pulse
[neff,del_w]=calc_rot_axis_arba(tref,pref,aref); % Recalculate refocusing axis

start=pexc;
nexc=length(texc);
aexc=ones(1,nexc); % Segments have arbitrary phase and constant amplitude

% Excitation pulse definition
% Use nonlinear function minimization - all segment times must be positive
lb=zeros(1,nexc); % Lower bound
ub=2*pi*ones(1,nexc); % Upper bound

% trust-region-reflective algorithm (fmincon default) will not work for
% this problem because of the constraints, so use interior-point or
% active-set algorithms instead
options=optimset('Algorithm','interior-point','Display','iter','TolFun',1e-4,'MaxFunEvals',10000);
%options=optimset('Algorithm','active-set','Display','iter','TolFun',1e-4,'MaxFunEvals',5000);

pexc=fmincon(@(params)fit_function(params,texc,aexc,neff,del_w,len_acq),start,[],[],[],[],lb,ub,[],options);
[masy]=cpmg_van_spin_dynamics_asymp_mag2(texc,pexc,aexc,neff,del_w,len_acq);

for i=1:length(tvect)
    echo(i)=sum(masy.*exp(-1i*del_w*tvect(i)));
end

out.texc=texc;
out.pexc=pexc;
out.tref=tref;
out.pref=pref;
out.aref=aref;
out.echo_pk=max(abs(echo));
out.echo_rms=sqrt(trapz(tvect,abs(echo).^2));

function val=fit_function(pexc,texc,aexc,neff,del_w,len_acq)

[masy]=cpmg_van_spin_dynamics_asymp_mag2(texc,pexc,aexc,neff,del_w,len_acq);

% Calculate time-domain echo
%echo=zeros(1,length(tvect));
%for i=1:length(tvect)
%    echo(i)=sum(masy.*exp(-1i*del_w*tvect(i)));
%end

% Phase inversion leaves behind only the symmetric part of the spectrum,
% i.e., the real component of the time-domain echo
%echo = real(echo);

% Optimize echo peak
%val=-0.01*max(abs(echo));

% Optimize echo RMS + echo peak
%val=-0.01*(sqrt(trapz(tvect,abs(echo).^2))+max(abs(echo)));

masy=real(masy);
val=-sqrt(trapz(del_w,masy.^2))-0.33*trapz(del_w,masy);