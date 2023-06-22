%% George Kenefati HW4B
close all;clear;clc;format short;format compact;

%% Parameter Controls

% Sinusoidal stimulus
f = input('What is the frequency in Hz?'); %Hz
omega=2*pi*f; %radians 
x=linspace(0,3.5,1000)'; %distance vector
delx=3.5/1000; %x-axis step size

% Constants
Pzo=2; %dyne/m^2
rho=1; %g/cm^3
m=0.01; %g/cm^2
s=2e8*exp(-1.5.*x); %dyne/cm^3
r=5*exp(2.25.*x); %dyne-s/cm^3

% -------------------------------------------------------------------------
%% Q4
% Uncomment a block to overwrite impedance parameters
% % Q4a
% m=0; %
% s=2e8*exp(-1.5.*x);
% r=5*exp(2.25.*x);

% % Q4b
% m=0.01;
% s=2e8*exp(-1.5.*x);
% r=(5*exp(2.25.*x)) / 10; %

% % Q4c
% m=0.01;
% s=2e8*exp(-1.5.*x) * 10; %
% r=5*exp(2.25.*x);

%--------------------------------------------------------------------------
%% Q5
% Uncomment a block to overwrite impedance parameters
% % Q5a
% m=0.015;

% % Q5b
% m=0.015;
% s=2e8*exp(-1.5.*x);
% r=5*exp(2.25.*x);

% % Q5c
% m=0.015;
% s=2e8*exp(-1.5.*x);
% r=(5*exp(2.25.*x)) / 10; %

%% Response Calulations

% Mechanical Impedance
Z = s./(1i*omega) + m.*1i*omega + r;
Zr = real(Z);
Zi = imag(Z);

% Admittance
Y = 1./Z;
Yr = real(Y);
Yi = imag(Y);
Ya = abs(Y);
Yp = angle(Y);

% Pressure at Basilar Membrane (BM)
integralY = cumsum(Y)*delx;
exp_arg = -2*omega*rho*integralY;
P = Pzo*exp(exp_arg);
Pr = real(P);

% Amplitude and Phase of Pressure
Pa = Pzo*exp(real(exp_arg));
Pp = (imag(exp_arg))/(2*pi); % unwrap phase into cycles

% Velocity of BM
V = 2.*Pr.*Yr;
Va = 2.*Pa.*Ya;
Vp = Pp + Yp;
% ---Uncomment to query point of resonance for Q5a---%
% [q5max_val, q5max_ind]=max(Va);
% q5pt_res=x(q5max_ind); % point of resonance

% Displacement of BM
Xa = (omega*Va) * 1e7; %cm -> nm
X = Xa.*cos(omega*x) * 1e7;
Xp = Vp - pi/2;

%--------------------------------------------------------------------------
% Before starting Q6, what is the original point of resonance along x?
[q6max_val_og, q6max_ind_og]=max(Xa);
q6pt_res_og=x(q6max_ind_og); % point of resonance
q6max_val_og; % amplitude at pt_res

%% Q6

% Before running this section, comment Q4 and Q5 parameters, and re-run
% cell 2 to reset original parameters.

% % Group the indices 2 before, at, and 2 after the point of resonance
% indx=[q6max_ind-2,q6max_ind-1,q6max_ind,q6max_ind+1,q6max_ind+2];
% 
% % Gradual increase in scale of negative resistance
% r_mag=[2,4,8,4,2];
% 
% % Decrease impedance at these indices (across a span of ~19 mm)
% for i=1:length(indx)
%     Z(indx) = s(indx)./(1i*omega) + m.*1i*omega - r_mag(i)*r(indx);
% end
% 
% Zr = real(Z);
% Zi = imag(Z);
% 
% %%% RECALCULATE Y, P, V, and X
% % Admittance
% Y = 1./Z;
% Yr = real(Y);
% Yi = imag(Y);
% Ya = abs(Y);
% Yp = angle(Y);
% 
% % Pressure at Basilar Membrane (BM)
% integralY = cumsum(Y)*delx;
% P = Pzo*exp(-2*omega*rho*integralY);
% Pr = real(P);
% 
% % Amplitude and Phase of Pressure
% Pa = Pzo*exp(real(-2*omega*rho*integralY));
% Pp = (imag(-2*omega*rho*integralY))/(2*pi);
% 
% % Velocity of BM
% V = 2.*Pr.*Yr;
% Va = 2.*Pa.*Ya;
% Vp = Pp + Yp;
% 
% % Recalculate displacement
% Xa = (omega*Va) * 1e7;
% X = Xa.*cos(omega*t) * 1e7;
% Xp = Vp - pi/2;
% 
% % Obtain new max impedance values
% [q6max_val_new, q6max_ind_new]=max(Xa);
% q6pt_res_new=x(q6max_ind_new); % point of resonance
% q6max_val_new; % amplitude at pt_res

% ------------------------------------------------------------------------
%% Q7
% Before running this section, comment Q4, Q5, and Q6 sections, and re-run
% cell 2 to reset original parameters.

% % Gradual increase in scale of positive stiffness
% s_mag=[2,4,8]; % first 10.5mm of cochlea base
% 
% % First, query original impedance values at cochlea base
% q7_ogZ = zeros(1,numel(s_mag)); %preallocate array
% for i=1:length(s_mag)
%     q7_ogZ(i)=Z(i);
% end
% 
% % Increase effective stiffness by the above scalars
% for i=1:length(s_mag)
%     Z(i) = ( s_mag(i)*s(i) )./(1i*omega) + m.*1i*omega + r(i);
% end
% 
% % Recalculate Zr and Zi
% Zr = real(Z);
% Zi = imag(Z);
% 
% %%% RECALCULATE Y, P, V, and X
% % Admittance
% Y = 1./Z;
% Yr = real(Y);
% Yi = imag(Y);
% Ya = abs(Y);
% Yp = angle(Y);
% 
% % Pressure at Basilar Membrane (BM)
% integralY = cumsum(Y)*delx;
% P = Pzo*exp(-2*omega*rho*integralY);
% Pr = real(P);
% 
% % Amplitude and Phase of Pressure
% Pa = Pzo*exp(real(-2*omega*rho*integralY));
% Pp = (imag(-2*omega*rho*integralY))/(2*pi);
% 
% % Velocity of BM
% V = 2.*Pr.*Yr;
% Va = 2.*Pa.*Ya;
% Vp = Pp + Yp;
% 
% % Recalculate displacement
% Xa = (omega*Va) * 1e7;
% X = Xa.*cos(omega*t) * 1e7;
% Xp = Vp - pi/2;
% 
% % Obtain new max impedance values
% [q7max_val, q7max_ind]=max(Xa);
% q6pt_res=x(q7max_ind); % point of resonance
% q7max_val; % amplitude at pt_res
% 
% % Query new impedance values at cochlea base
% q7_newZ = zeros(1,numel(s_mag)); %preallocate array
% for i=1:length(s_mag)
%     q7_newZ(i)=Z(i);
% end
% 
% % What is the mean difference in Z caused by the cochlear implant?
% q7_Zdiff = q7_newZ - q7_ogZ;
% q7_Zdiff_mean = mean(q7_Zdiff)

%% Plots

% Pressure
figure(1);
subplot(2,2,1);plot(x,Pa);title('Pressure Amplitude over the BM');
xlabel('Distance (cm)');ylabel('Pressure Amplitude (dynes/(cm^2))');

subplot(2,2,2);plot(x,Pp);title('Pressure Phase over the BM');
xlabel('Distance (cm)');ylabel('Pressure Phase (cycles)'); grid on;

subplot(2,2,3);plot(x,Pr);title('Pressure Wave over the BM');
xlabel('Distance (cm)');ylabel('Pressure Wave (dynes/(cm^2))'); grid on;

%% Velocity
figure(2);
subplot(2,2,1);plot(x,Va);title('Velocity Amplitude over the BM');
xlabel('Distance (cm)'); ylabel('Velocity Amplitude (cm/s)');

subplot(2,2,2);plot(x,Vp);title('Velocity Phase over the BM');
xlabel('Distance (cm)'); ylabel('Velocity Phase (cycles)'); grid on;

subplot(2,2,3);plot(x,V);title('Velocity Wave over the BM');
xlabel('Distance (cm)'); ylabel('Velocity Wave (cm/s)'); grid on;

%% Displacement
figure(3);
subplot(2,2,1);plot(x,Xa);title('Displacement Amplitude over the BM');
xlabel('Distance (nm)'); ylabel('Velocity Amplitude (nm/s)');

subplot(2,2,2);plot(x,Xp);title('Displacement Phase over the BM');
xlabel('Distance (nm)'); ylabel('Displacement Phase (cycles)'); grid on;

subplot(2,2,3);plot(x,X);title('Displacement Wave over the BM');
xlabel('Distance (nm)'); ylabel('Displacement Wave (nm/s)'); grid on;

%% Mechanical Impedance
figure(4);

% find where Zi crosses the x-axis
for i=1:length(Zi)
    if Zi(i) > -0.7 && Zi(i) < 0.7
        zcross=i*delx;
    end
end

subplot(2,1,1);plot(x,Zr); hold on; plot(x,Zi); hold on; 
xline(zcross);title('Mechanical Impedance (Real and Imaginary)'); 
txtlbl=sprintf('%f',zcross);
text(zcross, -1, txtlbl, 'Color', 'k', 'FontSize',14);
xlabel('Distance (cm)'); ylabel('Mechanical Impedance g/(s-cm^2)'); 
grid on; legend('Real','Imaginary','Location', 'Northwest');
title(legend,'Parts');

% Mechanical Admittance
subplot(2,1,2);plot(x,Yr); hold on; plot(x,Yi);
title('Mechanical Admittance (Real and Imaginary)');
xlabel('Distance (cm)'); ylabel('Mechanical Admittance (s-cm^2)/g'); 
grid on; legend('Real','Imaginary','Location', 'Northwest');
title(legend,'Parts');