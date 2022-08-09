clc;
% Import Earthquake ground acceleration data from .xls file
[EW] = xlsread('katanap.xls');

t = EW(:,1);  % time history of earthquake
Ag = EW(:,2);  % ground acceleration
dt = t(2)-t(1);  % unit time = 0.02

Tn = 0.01:0.01:3;  % natural period of vibration
zi = 0.025:0.025:0.1;  % Damping ratio  

m = 1;           % taking unit lumped mass of SDOF
wn = 2*pi./Tn;     %natural circular frequency of vibration
k = m*wn.^2;      % stiffness 

colorstring = 'kbgr';
sgtitle('Response Spectra','Fontweight',"bold");
for i = 1:4
    Sd = zeros(1,length(Tn));
    for j = 1:length(Tn)
        [u,v] = find_disp(Ag,zi(i),m,dt,wn(j),k(j));
        Sd(j) = max(abs(u));  
    end
   
    Sv = wn.*Sd;
    Sa_g = wn.*Sv/9.80665;
    
    set(gcf,'position',[50,50,100,500]);
    subplot(3,1,1)
    plot(Tn,Sd, 'Color', colorstring(i));
    title('Displacement Response Spectrum','Fontweight',"bold")
    xlabel('Period (sec)'),ylabel('Spectral Displacemnet (m)')
    legend('\zeta = 0.025','\zeta = 0.05','\zeta = 0.075','\zeta = 0.1', 'Location','Best');
    hold on;
    
    subplot(3,1,2)
    plot(Tn,Sv, 'Color', colorstring(i));
    title('Velocity Response Spectrum','Fontweight',"bold")
    xlabel('Period (sec)'),ylabel('Spectral Velocity (m/s)')
    legend('\zeta = 0.025','\zeta = 0.05','\zeta = 0.075','\zeta = 0.1', 'Location','Best');
    hold on;
    
    subplot(3,1,3)
    plot(Tn,Sa_g , 'Color', colorstring(i));
    title('Acceleration Response Spectrum','Fontweight',"bold")
    xlabel('Period (sec)'),ylabel('Spectral Acceleration (g)')
    legend('\zeta = 0.025','\zeta = 0.05','\zeta = 0.075','\zeta = 0.1', 'Location','Best');
    hold on;
end




%%Function  
% for calculation of response due to ground acc. of unit time step at time j
% Method based on interpolation excitation

function [u,v] = find_disp(Ag,zi,m,dt,wn,k)

sq = (1-zi^2)^(1/2);
ep = exp(-zi*wn*dt);
wd = wn*sq; % natural frequency of damped vibration
si = sin(wd*dt);
co = cos(wd*dt);

% Coefficients in reccurence formulas (?<1)
A = ep*(zi/sq*si+co);  % related to free vibration with damping
B = ep*(1/wd*si);  % related to free vibration with damping
C = 1/k*(2*zi/(wn*dt)+ep*(((1-2*zi^2)/(wd*dt)-zi/sq)*si-co*(1+2*zi/(wn*dt))));  % related to response to step force
D = 1/k*(1-2*zi/(wn*dt)+ep*(si*(2*zi^2-1)/(wd*dt)+2*zi/(wn*dt)*co));  % related to response to ramp force
 
Av = -ep*(wn/sq*si);
Bv = ep*(co-zi/sq*si);
Cv = 1/k*(-1/dt+ep*((wn/sq+zi/(dt*sq))*si+co/dt));
Dv = 1/(k*dt)*(1-ep*(zi/sq*si+co));

u = zeros(length(Ag),1);
v = zeros(length(Ag),1);

% setting intial coditions
u(1)=0;
v(1)=0;

for i=1:length(Ag)-1
    u(i+1)=A*u(i)+B*v(i)+C*(-m*Ag(i)*9.80665)+D*(-m*Ag(i+1)*9.80665);
    v(i+1)=Av*u(i)+Bv*v(i)+Cv*(-m*Ag(i)*9.80665)+Dv*(-m*Ag(i+1)*9.80665);
end
end
