function [echo_int_all] = sim_cpmg_tuned_probe_img(params)

% Read in parameters
NE=params.NE;
TE=params.TE;
Tgrad=params.Tgrad;
rho=params.rho;
T1map=params.T1map;
T2map=params.T2map;
pxz=params.pxz;
FOV=params.FOV;

% Define pulse system parameters
[sp, pp]=set_params_tuned_JMR;
T_90=pp.T_90; % Nominal T_90 pulse length
T_180=2*T_90;

% Check pulse sequence timing
if pp.tacq>(TE-T_180)
    pp.tacq=TE-T_180; % Maximum acquisition period = TE-T_180
end

% Set size of simulation domain
sp.numptsy=400; 
sp.maxoffs=5; % y-axis
siz=size(rho);
sp.nx=siz(1);
sp.nz=siz(2); % (x,z) plane
sp.rho=rho;
sp.T1map=T1map;
sp.T2map=T2map; % Set sample properties

% Set size of image
px=pxz(1); 
pz=pxz(2);

% Predicted maximum number of echoes before simulation instability occurs
% (in the absence of a gradient - gradients will make matters worse)
dw=2*sp.maxoffs/sp.numptsy; % Normalized to w1
NE_max=2*T_180/(dw*TE);
display(['NE_max = ' num2str(NE_max)]);
if NE>NE_max
    disp('Warning: Attempt to simulate beyond NE_max!');
    disp('Accuracy of results is not guaranteed');
end

% Set plotting parameters
sp.plt_tx = 0;
sp.plt_rx = 0;
sp.plt_sequence = 0; % Plots on/off
sp.plt_axis = 0;
sp.plt_mn = 0;
sp.plt_echo = 0;
sp.plt_output = 1;
sp.plt_fields = 0;

% Create fields and sample parameters
sp=create_fields_single_sided(sp);

% Create receiver sensitivity map
sp.w_1r=sp.w_1; % , assuming single transmit/receive coil

% Set acquisition parameters
tacq=(pi/2)*pp.tacq/T_90; % Normalized acquisition window length
tdw=(pi/2)*pp.tdw/T_90; % Normalized receiver dwell time
nacq=round(tacq/tdw)+1; % Number of acquired time domain points
tvect=linspace(-tacq/2,tacq/2,nacq)'; % Acquisition vector, size: [nacq,1]
isoc=exp(1i*tvect*sp.del_w); % Isochromats (without additional gradients)


% Create pulse sequence (in normalized time)
% --------------------------------------------
% Pre-calculate all pulses for speed
Rtot={};
pp_in.tp=pp.T_90;
pp_in.tdel=2*pp.T_90;


pp_in.phi=pi/2;
pp_in.amp=1;
pp_out=calc_pulse_shape(sp,pp,pp_in); % Exc pulse y (phase = pi/2)
Rtot{1}=calc_rotation_matrix(sp,pp_out);

pp_in.phi=3*pi/2;
pp_out=calc_pulse_shape(sp,pp,pp_in); % Exc pulse -y (phase = 3*pi/2)
Rtot{2}=calc_rotation_matrix(sp,pp_out);

pp_in.tp=pp.T_180; 
pp_in.phi=0;
pp_out=calc_pulse_shape(sp,pp,pp_in); % Ref pulse x (phase = 0)
Rtot{3}=calc_rotation_matrix(sp,pp_out);

pp_in.tp=pp.T_180; pp_in.phi=pi/2;
pp_out=calc_pulse_shape(sp,pp,pp_in); % Ref pulse y (phase = pi/2)
Rtot{4}=calc_rotation_matrix(sp,pp_out);

% Excitation pulses 1 and 2, including timing correction
texc=[pi/2 -1]; aexc=[1 0];
pexc1=[1 0]; pexc2=[2 0];  % pexc is pulse type
acq_exc=[0 0]; gexc=[0 0]; % gexc is the gradient

% Encoding periods 1 and 2
tenc=[(pi/2)*Tgrad/T_90 pi (pi/2)*Tgrad/T_90]; aenc=[0 1 0];
penc1=[0 3 0]; penc2=[0 4 0];
acq_enc=[0 0 0]; genc=[1 0 0]; % Note: Gradient is enabled

% Refocusing cycles
nref=3; % Segments in refocusing cycle
tref=zeros(1,nref*NE); pref1=tref; pref2=tref;
aref=tref; acq_ref=tref; gref=tref;
tfp=(pi/2)*(TE-pp.T_180)/(2*T_90); % Free precession period (normalized)
for i=1:NE
    tref((i-1)*nref+1:i*nref)=[tfp pi tfp];
    pref1((i-1)*nref+1:i*nref)=[0 3 0]; % Pulse type
    pref2((i-1)*nref+1:i*nref)=[0 4 0]; % Pulse type
    aref((i-1)*nref+1:i*nref)=[0 1 0];
    acq_ref((i-1)*nref+1:i*nref)=[0 0 1];
    gref((i-1)*nref+1:i*nref)=[0 0 0]; % Gradient
end

% Create complete pulse sequences
% Fixed terms
pp.tp=[texc tenc tref]; pp.amp=[aexc aenc aref];
pp.acq=[acq_exc acq_enc acq_ref]; pp.grad=[gexc genc gref];
pp.Rtot=Rtot;

% Variable terms
pul1=[pexc1 penc1 pref1]; % Sequence 1: (y,x)
pul2=[pexc2 penc1 pref1]; % Sequence 2: (-y,x)
pul3=[pexc1 penc2 pref2]; % Sequence 3: (y,y)
pul4=[pexc2 penc2 pref2]; % Sequence 4: (-y,y)

% Calculate gradient strength vectors
Tgradn=(pi/2)*Tgrad/T_90; % Normalized gradient length
% Estimate maximum gradient frequency offsets using specified FOV
wxmax=pi*px^2/(2*FOV(1)*Tgradn); wzmax=pi*pz^2/(2*FOV(2)*Tgradn);
gradx=wxmax*linspace(-1,1,px); % x-gradient steps
gradz=wzmax*linspace(-1,1,pz); % z-gradient steps

%Get receiver transfer function
[sp.tf]=tuned_probe_rx_tf(sp,pp); % Filtering by tuned receiver


% Generate image
echo_int_all=zeros(px,pz,NE);
for i=1:px % Parallelize for speed
    spc=sp; ppc=pp; % Create local variables
    gradxc=gradx; gradzc=gradz;
    
    for j=1:pz
        % Create gradient field
        spc.del_wg=gradxc(i)*spc.del_wx+gradzc(j)*spc.del_wz;
        
        % Run PAP phase cycle: (1-2)
        ppc.pul=pul1; [macq1,mrx1]=calc_macq_tuned_probe_relax4(spc,ppc);
        ppc.pul=pul2; [macq2,mrx2]=calc_macq_tuned_probe_relax4(spc,ppc);
        %mrx_x=(mrx1-mrx2); % Apply phase cycle
        mrx_x = macq1;
        % Run PAP phase cycle: (3-4)
        ppc.pul=pul3; [macq3,mrx3]=calc_macq_tuned_probe_relax4(spc,ppc);
        ppc.pul=pul4; [macq4,mrx4]=calc_macq_tuned_probe_relax4(spc,ppc);
       % mrx_y=(mrx3-mrx4); % Apply phase cycle
        mrx_y = macq3;
        % Calculate time-domain echoes
        echo_rx_x=isoc*mrx_x'; % Size: [nacq,NE]
        echo_rx_y=isoc*mrx_y'; % Size: [nacq,NE]
        echo_rx_xy=imag(echo_rx_x)-1i*real(echo_rx_y); % Raw echo shapes
        echo_int_all(i,j,:)=trapz(tvect,echo_rx_xy)'; % Estimate echo integrals, size: [NE,1]
        
        disp(['Running: ' num2str(round(100*j/pz)) '% of row ' num2str(i)])
    end
    disp(['------Row ' num2str(i) ' (of ' num2str(px) ') completed.------'])
end

if sp.plt_output % Plot selected outputs
    eplt=2; % Echo number to plot
    echo_int_eplt=echo_int_all(:,:,eplt); % Select data for chosen echo
    
    figure; % Plot k-space of selected echo
    subplot(1,3,1); imagesc(real(echo_int_eplt));
    colorbar; title('k-space (Real)')
    subplot(1,3,2); imagesc(imag(echo_int_eplt));
    colorbar; title('k-space (Imag)')
    subplot(1,3,3); imagesc(abs(echo_int_eplt));
    colorbar; title('k-space (Mag)')
    
    figure; % Plot spin density image
    img_eplt=ifftshift(ifft2(echo_int_eplt)); % Estimate image
    imagesc(abs(img_eplt));
    colorbar; title('Spin density image (Mag)');
    end
end
function [pp_out]=calc_pulse_shape(sp,pp,pp_in)

T_90=pp.T_90;
tdeln=(pi/2)*pp_in.tdel/T_90; % Normalized delay
amp_zero=pp.amp_zero; % Minimum amplitude for calculations

% Add delay to RF pulse to account for ring down, create structure
pp_curr=pp;
pp_curr.tp = [pp_in.tp pp_in.tdel];
pp_curr.phi = [pp_in.phi 0];
pp_curr.amp = [pp_in.amp 0];

% Calculate RF pulse
sp.plt_rx=0; % Turn off plotting
[tvect2, Icr2,tvect,Icr] = tuned_probe_lp(sp,pp_curr);
% pp_out.tf1=tf1; pp_out.tf2=tf2;

delt=(pi/2)*(tvect(2)-tvect(1))/T_90; % Convert to normalized time
texc=delt*ones(1,length(tvect));
pexc=atan2(imag(Icr),real(Icr));
aexc=abs(Icr);
aexc(aexc<amp_zero)=0; % Threshold amplitude
aexc = (aexc - min(aexc)) / ( max(aexc) - min(aexc));
% aexc(aexc>1)=1;
aexc(aexc<amp_zero)=0; % Threshold amplitude


% Remove added delay from RF pulse
pp_out.tp=[texc -tdeln]; pp_out.phi=[pexc 0]; pp_out.amp=[aexc 0];
pp_out.acq=zeros(1,length(texc)+1);

end
