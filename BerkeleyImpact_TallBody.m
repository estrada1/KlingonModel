function [t,VAR,Output] = BerkeleyImpact_TallBody
%===========================================================================
% File: BerkeleyImpact_TallBody.m created on Wed Oct  7 2015 by MotionGenesis 5.7.
% Advanced Student Licensee: Matt Estrada (until December 2015).
% Portions copyright (c) 2009-2015 Motion Genesis LLC.  Rights reserved.
% Paid-up MotionGenesis Advanced Student licensees are granted the right
% right to distribute this code for legal student-academic (non-professional) purposes only,
% provided this copyright notice appears in all copies and distributions.
%===========================================================================
% The software is provided "as is", without warranty of any kind, express or    
% implied, including but not limited to the warranties of merchantability or    
% fitness for a particular purpose. In no event shall the authors, contributors,
% or copyright holders be liable for any claim, damages or other liability,     
% whether in an action of contract, tort, or otherwise, arising from, out of, or
% in connection with the software or the use or other dealings in the software. 
%===========================================================================
%  FIX THIS RIGHT HERE
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
Ffoam_bottom=0; Ffoam_top=0;
mewBottom=0; mewTop=0; gammapp=0; phipp=0; xpp=0; ypp=0; Fx_tail=0; Fy_tail=0; T_tail=0; AttachPt_xp=0; AttachPt_yp=0; AttachPt_x_WorldFrame=0; AttachPt_y_WorldFrame=0; E_rebound=0; E_tail=0; FoamContactBottom=0; FoamContactTop=0; FootAttached=0;
 FootContact=0; Fx_contact=0; Fx_rebound=0; Fy_contact=0; Fy_fricBottom=0; Fy_fricTop=0; Fy_rebound=0; F_hardstopBottom=0; F_hardstopTop=0; GravityPotentialEnergy=0; HardstopBottomCompress=0; HardstopContactBottom=0; HardstopContactTop=0;
 HardstopTopCompress=0; H_body=0; H_tail=0; KineticEnergy=0; KineticEnergy_body=0; KineticEnergy_tail=0; Lx_body=0; Lx_tail=0; Ly_body=0; Ly_tail=0; rebound_mag=0; rebound_x=0; TailCompress=0; TailContact=0; vBottom=0; vTop=0; x_Acm=0;
 x_bottom=0; x_Ccm=0; x_FoamBottom=0; x_FoamTop=0; x_hardstop=0; x_tail=0; x_top=0; y_Acm=0; y_bottom=0; y_Ccm=0; y_FoamBottom=0; y_FoamTop=0; y_tail=0; y_top=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
Awidth                          =  3;                      % cm                  Constant
b_hardstop                      =  25;                     % N/m/s               Constant
b_tail                          =  0.01;                   % N/s                 Constant
b_top                           =  1;                      % N/m/s               Constant
d_body                          =  10;                     % cm                  Constant
d_foam                          =  3;                      % cm                  Constant
Fpre_rebound                    =  4;                      % N                   Constant
g                               =  9.80665;                % m/sec^2             Constant
IAzz                            =  2.90e-4;                % kg*m^2              Constant
ICzz                            =  0.000002;               % kg*m^2              Constant
k_hardstop                      =  10000;                  % N/m                 Constant
k_rebound                       =  14.28571428571428;      % N/m                 Constant
k_tail                          =  0.1;                    % N                   Constant
k_top                           =  100;                    % N/m                 Constant
l_foam                          =  0.02;                   % m                   Constant
l_tail                          =  0.23;                   % m                   Constant
mA                              =  .2;                     % kg                  Constant
mC                              =  .00001;                 % kg                  Constant
mQ                              =  .02;                    % kg                  Constant
phin                            =  90;                     % deg                 Constant

gamma                           = -35;                     % deg                 Initial Value
phi                             =  90;                     % deg                 Initial Value
x                               = -.1;                     % m                   Initial Value
y                               =  4;                      % m                   Initial Value
gammap                          =  0;                      % rad/sec             Initial Value
phip                            =  0;                      % rad/sec             Initial Value
xp                              =  1;                      % m/sec               Initial Value
yp                              =  .25;                    % m/sec               Initial Value
AttachPt_x                      =  0.0;                    % UNITS               Initial Value
AttachPt_y                      =  0.0;                    % UNITS               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  0.5;                    % sec                 Final Time
tStep                           =  0.005;                  % sec                 Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-07;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions.  UnitSystem: kilogram, meter, second
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
Awidth = Awidth * 0.01;
d_body = d_body * 0.01;
d_foam = d_foam * 0.01;
phin = phin * DEGtoRAD;
gamma = gamma * DEGtoRAD;
phi = phi * DEGtoRAD;

% Evaluate constants
mewTop = -0.001;
mewBottom = -0.1;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files



%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
FootContact = ceil(0.5*sign(x));
HardstopContactTop = ceil(0.5*sign(x-l_foam*cos(gamma)-(d_body+d_foam)*sin(gamma)));
HardstopContactBottom = ceil(0.5*sign(x+d_foam*sin(gamma)-l_foam*cos(gamma)));
TailContact = ceil(0.5*sign(x-l_tail*cos(gamma+phi)-(Awidth+l_foam)*cos(gamma)));
FootAttached = ceil(0.5*sign(sqrt(AttachPt_x^2+AttachPt_y^2))+0.5*FootContact);
HardstopTopCompress = ceil(0.5*sign(xp+l_foam*sin(gamma)*gammap-(d_body+d_foam)*cos(gamma)*gammap));
F_hardstopTop = -HardstopContactTop*(k_top*(x-l_foam*cos(gamma)-(d_body+d_foam)*sin(gamma))+b_top*(xp+l_foam*sin(gamma)*gammap-(d_body+d_foam)*cos(gamma)*gammap)*HardstopTopCompress);
HardstopBottomCompress = ceil(0.5*sign(xp+d_foam*cos(gamma)*gammap+l_foam*sin(gamma)*gammap));
F_hardstopBottom = -HardstopContactBottom*(k_hardstop*(x+d_foam*sin(gamma)-l_foam*cos(gamma))+b_hardstop*(xp+d_foam*cos(gamma)*gammap+l_foam*sin(gamma)*gammap)*HardstopBottomCompress);
vTop = yp - l_foam*cos(gamma)*gammap - (d_body+d_foam)*sin(gamma)*gammap;
vBottom = yp + d_foam*sin(gamma)*gammap - l_foam*cos(gamma)*gammap;
Fy_fricTop = -mewTop*F_hardstopTop*vTop/(1.0E-6+abs(vTop));
TailCompress = ceil(0.5*sign(xp+(Awidth+l_foam)*sin(gamma)*gammap+l_tail*sin(gamma+phi)*(gammap+phip)));
rebound_x = AttachPt_x*(1-FootContact);
rebound_mag = sqrt(AttachPt_y^2+rebound_x^2);
Fx_rebound = FootAttached*rebound_x*(k_rebound+Fpre_rebound/(1.0E-8+rebound_mag));
Fy_rebound = AttachPt_y*FootAttached*(k_rebound+Fpre_rebound/(1.0E-8+rebound_mag));

% Quantities that were specified
T_tail = k_tail*(phin-phi) - b_tail*phip;
Fx_tail = -TailContact*(k_hardstop*(x-l_tail*cos(gamma+phi)-(Awidth+l_foam)*cos(gamma))+b_hardstop*(xp+(Awidth+l_foam)*sin(gamma)*gammap+l_tail*sin(gamma+phi)*(gammap+phip))*TailCompress);
Fy_tail = -mewBottom*Fx_tail*sign(yp-(Awidth+l_foam)*cos(gamma)*gammap-l_tail*cos(gamma+phi)*(gammap+phip));
AttachPt_xp = -FootAttached*xp;
AttachPt_yp = -FootAttached*yp;


% Quantities to be specified
Ffoam_bottom = 0;
Ffoam_top = 0;
Fy_fricBottom = -mewBottom*vBottom*(Ffoam_bottom+F_hardstopBottom)/(1.0E-6+abs(vBottom));
xpp = ((l_tail*mQ*sin(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))-l_tail*mQ*(Awidth+l_foam)*sin(phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+  ...
l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))-(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(  ...
gamma))))*(T_tail+l_tail*(Fx_tail*sin(gamma+phi)+(g*mQ-Fy_tail)*cos(gamma+phi))-l_tail*mQ*(Awidth+l_foam)*sin(phi)*gammap^2)-(l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-l_tail^2*mQ^2*(Awidth+l_foam)*sin(  ...
phi)*cos(gamma+phi)-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))*(l_tail*Fx_tail*sin(gamma+phi)+(Awidth+l_foam)*Fx_tail*sin(gamma)+g*mA*(Awidth+  ...
l_foam)*cos(gamma)+g*mC*(Awidth+l_foam)*cos(gamma)+l_tail*(g*mQ-Fy_tail)*cos(gamma+phi)+(Awidth+l_foam)*(g*mQ-Fy_tail)*cos(gamma)+d_foam*sin(gamma)*Fy_fricBottom+d_foam*cos(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopBottom+  ...
l_foam*sin(gamma)*F_hardstopTop-l_tail*mQ*(Awidth+l_foam)*sin(phi)*(gammap^2-(gammap+phip)^2)-l_foam*cos(gamma)*Fy_fricBottom-l_foam*cos(gamma)*Fy_fricTop-(d_body+d_foam)*sin(gamma)*Fy_fricTop-(d_body+d_foam)*cos(gamma)*F_hardstopTop))/(  ...
l_tail^2*mQ^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(mA+mC+mQ)*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(  ...
Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC+mQ)*sin(phi)^2-(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(  ...
gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+  ...
l_foam)*sin(gamma))))-(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(  ...
mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma))))) - ((l_tail^2*mQ^2*cos(gamma+phi)^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+  ...
2*l_tail*(Awidth+l_foam)*cos(phi)))+(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-2*l_tail*mQ*cos(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(  ...
gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma))))+(ICzz+mQ*l_tail^2)*((mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+  ...
l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*(Awidth+l_foam)*cos(gamma)*gammap^2+mQ*((Awidth+l_foam)*cos(gamma)*gammap^2+l_tail*cos(  ...
gamma+phi)*(gammap+phip)^2)-Fx_tail-Fx_rebound-F_hardstopBottom-F_hardstopTop)-(l_tail*mQ*(Awidth+l_foam)*(mA+mC+mQ)*cos(gamma)*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(  ...
gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))-l_tail*mQ*cos(gamma+phi)*(  ...
l_tail*mQ*sin(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))-(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*(mA*(Awidth+l_foam)*sin(gamma)+mC*(  ...
Awidth+l_foam)*sin(gamma)+mQ*((Awidth+l_foam)*sin(gamma)+2*l_tail*sin(gamma+phi)))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*(Awidth+l_foam)*sin(gamma)*gammap^2+mQ*((Awidth+l_foam)*sin(gamma)*gammap^2+l_tail*sin(  ...
gamma+phi)*(gammap+phip)^2)-Fy_tail-Fy_rebound-Fy_fricBottom-Fy_fricTop))/((mA+mC+mQ)*(l_tail^2*mQ^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(mA+mC+mQ)*((ICzz+  ...
l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC+mQ)*sin(phi)^  ...
2-(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*  ...
sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))-(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*l_tail*mQ*cos(  ...
gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma))))));
ypp = ((l_tail*mQ*(Awidth+l_foam)*(mA+mC+mQ)*sin(gamma)*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(  ...
Awidth+l_foam)*sin(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))-l_tail*mQ*sin(gamma+phi)*(l_tail*mQ*cos(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(  ...
Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))-(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*((Awidth+l_foam)*cos(gamma)+2*  ...
l_tail*cos(gamma+phi)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*(Awidth+l_foam)*cos(gamma)*gammap^2+mQ*((Awidth+l_foam)*cos(gamma)*gammap^2+l_tail*cos(gamma+phi)*(gammap+phip)^2)-Fx_tail-Fx_rebound-F_hardstopBottom-F_hardstopTop)-(  ...
l_tail^2*mQ^2*sin(gamma+phi)^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*((mA+mC+mQ)*(ICzz+l_tail*mQ*(  ...
l_tail+(Awidth+l_foam)*cos(phi)))-2*l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))+(ICzz+mQ*l_tail^2)*((mA*(Awidth+l_foam)*sin(gamma)+mC*(  ...
Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))))*(g*mA+g*mC+g*mQ+mA*(  ...
Awidth+l_foam)*sin(gamma)*gammap^2+mC*(Awidth+l_foam)*sin(gamma)*gammap^2+mQ*((Awidth+l_foam)*sin(gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Fy_tail-Fy_rebound-Fy_fricBottom-Fy_fricTop))/((mA+mC+mQ)*(l_tail^2*mQ^2*(IAzz+ICzz+  ...
mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(mA+mC+mQ)*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(  ...
Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC+mQ)*sin(phi)^2-(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+  ...
l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))-(  ...
mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(  ...
gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))))) - ((l_tail*mQ*(Awidth+l_foam)*sin(phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+  ...
l_foam)*sin(gamma)))+l_tail*mQ*cos(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))-(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*(mA*(Awidth+  ...
l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma))))*(T_tail+l_tail*(Fx_tail*sin(gamma+phi)+(g*mQ-Fy_tail)*cos(gamma+phi))-l_tail*mQ*(Awidth+l_foam)*sin(phi)*gammap^2)-(l_tail^2*mQ^2*(  ...
Awidth+l_foam)*sin(phi)*sin(gamma+phi)+l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+  ...
l_foam)*cos(gamma))))*(l_tail*Fx_tail*sin(gamma+phi)+(Awidth+l_foam)*Fx_tail*sin(gamma)+g*mA*(Awidth+l_foam)*cos(gamma)+g*mC*(Awidth+l_foam)*cos(gamma)+l_tail*(g*mQ-Fy_tail)*cos(gamma+phi)+(Awidth+l_foam)*(g*mQ-Fy_tail)*cos(gamma)+  ...
d_foam*sin(gamma)*Fy_fricBottom+d_foam*cos(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopTop-l_tail*mQ*(Awidth+l_foam)*sin(phi)*(gammap^2-(gammap+phip)^2)-l_foam*cos(gamma)*Fy_fricBottom-l_foam*  ...
cos(gamma)*Fy_fricTop-(d_body+d_foam)*sin(gamma)*Fy_fricTop-(d_body+d_foam)*cos(gamma)*F_hardstopTop))/(l_tail^2*mQ^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(  ...
mA+mC+mQ)*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(  ...
mA+mC+mQ)*sin(phi)^2-(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(  ...
mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))-(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*  ...
l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))));
gammapp = ((l_tail^2*mQ^2-(mA+mC+mQ)*(ICzz+mQ*l_tail^2))*((Awidth+l_foam)*Fx_tail*sin(gamma)+g*mA*(Awidth+l_foam)*cos(gamma)+g*mC*(Awidth+l_foam)*cos(gamma)+(Awidth+l_foam)*(g*mQ-Fy_tail)*cos(gamma)+l_tail*mQ*(Awidth+l_foam)*sin(phi)*(  ...
gammap+phip)^2+d_foam*sin(gamma)*Fy_fricBottom+d_foam*cos(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopTop-l_foam*cos(gamma)*Fy_fricBottom-l_foam*cos(gamma)*Fy_fricTop-(d_body+d_foam)*sin(  ...
gamma)*Fy_fricTop-(d_body+d_foam)*cos(gamma)*F_hardstopTop)+(l_tail^2*mQ^2*(Awidth+l_foam)*sin(gamma)*cos(gamma+phi)^2+l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi))-l_tail*mQ*(Awidth+l_foam)*cos(gamma)*cos(  ...
gamma+phi))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*(Awidth+l_foam)*cos(gamma)*gammap^2+mQ*((  ...
Awidth+l_foam)*cos(gamma)*gammap^2+l_tail*cos(gamma+phi)*(gammap+phip)^2)-Fx_tail-Fx_rebound-F_hardstopBottom-F_hardstopTop)-(l_tail^2*mQ^2-(mA+mC+mQ)*(ICzz+mQ*l_tail^2))*T_tail-(l_tail^2*mQ^2*(Awidth+l_foam)*cos(gamma)*sin(gamma+phi)^  ...
2+l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi))-l_tail*mQ*(Awidth+l_foam)*sin(gamma)*sin(gamma+phi))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(  ...
Awidth+l_foam)*cos(gamma))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*(Awidth+l_foam)*sin(gamma)*gammap^2+mQ*((Awidth+l_foam)*sin(gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Fy_tail-Fy_rebound-Fy_fricBottom-  ...
Fy_fricTop))/(l_tail^2*mQ^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(mA+mC+mQ)*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+  ...
ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC+mQ)*sin(phi)^2-(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(  ...
l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(  ...
Awidth+l_foam)*sin(gamma))))-(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*  ...
l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))));
phipp = (((mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))^2+(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(  ...
gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))*(T_tail+l_tail*(Fx_tail*sin(gamma+phi)+(g*mQ-Fy_tail)*cos(gamma+phi))-l_tail*mQ*(Awidth+  ...
l_foam)*sin(phi)*gammap^2)-(l_tail^2*mQ^2-(mA+mC+mQ)*(ICzz+mQ*l_tail^2))*(l_tail*Fx_tail*sin(gamma+phi)+(Awidth+l_foam)*Fx_tail*sin(gamma)+g*mA*(Awidth+l_foam)*cos(gamma)+g*mC*(Awidth+l_foam)*cos(gamma)+l_tail*(g*mQ-Fy_tail)*cos(gamma+  ...
phi)+(Awidth+l_foam)*(g*mQ-Fy_tail)*cos(gamma)+d_foam*sin(gamma)*Fy_fricBottom+d_foam*cos(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopBottom+l_foam*sin(gamma)*F_hardstopTop-l_tail*mQ*(Awidth+l_foam)*sin(phi)*(gammap^2-(gammap+phip)^  ...
2)-l_foam*cos(gamma)*Fy_fricBottom-l_foam*cos(gamma)*Fy_fricTop-(d_body+d_foam)*sin(gamma)*Fy_fricTop-(d_body+d_foam)*cos(gamma)*F_hardstopTop))/(l_tail^2*mQ^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+  ...
l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(mA+mC+mQ)*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+  ...
l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC+mQ)*sin(phi)^2-(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(  ...
l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))-(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(  ...
l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(  ...
Awidth+l_foam)*cos(gamma))))) - ((l_tail*mQ*cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+  ...
mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))-(mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))-  ...
l_tail*mQ*sin(gamma+phi)*((mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+  ...
l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*(Awidth+l_foam)*cos(gamma)*gammap^2+mQ*((Awidth+l_foam)*cos(gamma)*gammap^2+l_tail*cos(gamma+phi)*(gammap+phip)^2)-Fx_tail-Fx_rebound-  ...
F_hardstopBottom-F_hardstopTop)-(l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+  ...
mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))-(mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))-  ...
l_tail*mQ*cos(gamma+phi)*((mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+  ...
l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*(Awidth+l_foam)*sin(gamma)*gammap^2+mQ*((Awidth+l_foam)*sin(gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Fy_tail-  ...
Fy_rebound-Fy_fricBottom-Fy_fricTop))/((mA+mC+mQ)*(l_tail^2*mQ^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)))+(mA+mC+mQ)*((ICzz+l_tail*mQ*(l_tail+(Awidth+  ...
l_foam)*cos(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2+mQ*(l_tail^2+(Awidth+l_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi))))-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC+mQ)*sin(phi)^2-(mA*(Awidth+l_foam)*  ...
sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*(  ...
Awidth+l_foam)*sin(gamma)+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma))))-(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+  ...
l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma))))));

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros(1,10);
VAR(1) = gamma;
VAR(2) = phi;
VAR(3) = x;
VAR(4) = y;
VAR(5) = gammap;
VAR(6) = phip;
VAR(7) = xp;
VAR(8) = yp;
VAR(9) = AttachPt_x;
VAR(10) = AttachPt_y;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
gamma = VAR(1);
phi = VAR(2);
x = VAR(3);
y = VAR(4);
gammap = VAR(5);
phip = VAR(6);
xp = VAR(7);
yp = VAR(8);
AttachPt_x = VAR(9);
AttachPt_y = VAR(10);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros(1,10);
VARp(1) = gammap;
VARp(2) = phip;
VARp(3) = xp;
VARp(4) = yp;
VARp(5) = gammapp;
VARp(6) = phipp;
VARp(7) = xpp;
VARp(8) = ypp;
VARp(9) = AttachPt_xp;
VARp(10) = AttachPt_yp;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
FoamContactTop = ceil(0.5*sign(x-d_foam*sin(gamma)));
FoamContactBottom = ceil(0.5*sign(x+d_foam*sin(gamma)));
x_tail = x - l_tail*cos(gamma+phi) - (Awidth+l_foam)*cos(gamma);
y_tail = y - l_tail*sin(gamma+phi) - (Awidth+l_foam)*sin(gamma);
x_Acm = x - (Awidth+l_foam)*cos(gamma);
y_Acm = y - (Awidth+l_foam)*sin(gamma);
x_FoamTop = x - d_foam*sin(gamma);
y_FoamTop = y + d_foam*cos(gamma);
x_FoamBottom = x + d_foam*sin(gamma);
y_FoamBottom = y - d_foam*cos(gamma);
x_hardstop = x - l_foam*cos(gamma) - (d_body+d_foam)*sin(gamma);
x_Ccm = x - (Awidth+l_foam)*cos(gamma);
y_Ccm = y - (Awidth+l_foam)*sin(gamma);
x_top = x - l_foam*cos(gamma) - (d_body+d_foam)*sin(gamma);
y_top = y + (d_body+d_foam)*cos(gamma) - l_foam*sin(gamma);
x_bottom = x + d_foam*sin(gamma) - l_foam*cos(gamma);
y_bottom = y - d_foam*cos(gamma) - l_foam*sin(gamma);
AttachPt_x_WorldFrame = AttachPt_x + x;
AttachPt_y_WorldFrame = AttachPt_y + y;
Fx_contact = Ffoam_bottom + Ffoam_top + Fx_tail + Fx_rebound + F_hardstopBottom + F_hardstopTop;
Fy_contact = Fy_tail + Fy_rebound + Fy_fricBottom + Fy_fricTop;
KineticEnergy = 0.5*IAzz*gammap^2 + 0.5*ICzz*(gammap+phip)^2 + 0.5*mA*(xp^2+yp^2+(Awidth+l_foam)^2*gammap^2+2*(Awidth+l_foam)*sin(gamma)*gammap*xp-2*(Awidth+l_foam)*cos(gamma)*gammap*yp) + 0.5*mC*(xp^2+yp^2+(Awidth+l_foam)^2*gammap^2+2*(  ...
Awidth+l_foam)*sin(gamma)*gammap*xp-2*(Awidth+l_foam)*cos(gamma)*gammap*yp) - 0.5*mQ*(2*(Awidth+l_foam)*cos(gamma)*gammap*yp+2*l_tail*cos(gamma+phi)*yp*(gammap+phip)-xp^2-yp^2-(Awidth+l_foam)^2*gammap^2-l_tail^2*(gammap+phip)^2-2*(  ...
Awidth+l_foam)*sin(gamma)*gammap*xp-2*l_tail*sin(gamma+phi)*xp*(gammap+phip)-2*l_tail*(Awidth+l_foam)*cos(phi)*gammap*(gammap+phip));
GravityPotentialEnergy = g*(mA*(y-(Awidth+l_foam)*sin(gamma))+mC*(y-(Awidth+l_foam)*sin(gamma))+mQ*(y-l_tail*sin(gamma+phi)-(Awidth+l_foam)*sin(gamma)));
KineticEnergy_body = 0.5*IAzz*gammap^2 + 0.5*mA*(xp^2+yp^2+(Awidth+l_foam)^2*gammap^2+2*(Awidth+l_foam)*sin(gamma)*gammap*xp-2*(Awidth+l_foam)*cos(gamma)*gammap*yp);
KineticEnergy_tail = -0.5*mQ*(2*(Awidth+l_foam)*cos(gamma)*gammap*yp+2*l_tail*cos(gamma+phi)*yp*(gammap+phip)-xp^2-yp^2-(Awidth+l_foam)^2*gammap^2-l_tail^2*(gammap+phip)^2-2*(Awidth+l_foam)*sin(gamma)*gammap*xp-2*l_tail*sin(gamma+phi)*  ...
xp*(gammap+phip)-2*l_tail*(Awidth+l_foam)*cos(phi)*gammap*(gammap+phip));
E_rebound = 0.5*rebound_mag*(2*Fpre_rebound+k_rebound*rebound_mag);
E_tail = 0.5*k_tail*(phin-phi)^2;
Lx_body = mA*(xp+(Awidth+l_foam)*sin(gamma)*gammap);
Ly_body = mA*(yp-(Awidth+l_foam)*cos(gamma)*gammap);
H_body = IAzz*gammap - mA*(y*xp+(Awidth+l_foam)*cos(gamma)*yp-x*yp-(Awidth+l_foam)*sin(gamma)*xp-(Awidth+l_foam)*(Awidth+l_foam-x*cos(gamma)-y*sin(gamma))*gammap);
Lx_tail = mQ*(xp+(Awidth+l_foam)*sin(gamma)*gammap+l_tail*sin(gamma+phi)*(gammap+phip));
Ly_tail = mQ*(yp-(Awidth+l_foam)*cos(gamma)*gammap-l_tail*cos(gamma+phi)*(gammap+phip));
H_tail = -mQ*(y*xp+l_tail*cos(gamma+phi)*yp+(Awidth+l_foam)*cos(gamma)*yp+l_tail*x*cos(gamma+phi)*(gammap+phip)+l_tail*y*sin(gamma+phi)*(gammap+phip)-x*yp-l_tail^2*(gammap+phip)-l_tail*sin(gamma+phi)*xp-(Awidth+l_foam)*sin(gamma)*xp-  ...
l_tail*(Awidth+l_foam)*cos(phi)*(phip+2*gammap)-(Awidth+l_foam)*(Awidth+l_foam-x*cos(gamma)-y*sin(gamma))*gammap);

Output = zeros(1,61);
Output(1) = t;
Output(2) = x;
Output(3) = y;
Output(4) = x_Acm;
Output(5) = y_Acm;
Output(6) = x_tail;
Output(7) = y_tail;
Output(8) = gamma;
Output(9) = phi;

Output(10) = x_FoamTop;
Output(11) = y_FoamTop;
Output(12) = x_FoamBottom;
Output(13) = y_FoamBottom;
Output(14) = AttachPt_x_WorldFrame;
Output(15) = AttachPt_y_WorldFrame;
Output(16) = x_hardstop;
Output(17) = x_Ccm;
Output(18) = y_Ccm;
Output(19) = x_top;
Output(20) = y_top;
Output(21) = x_bottom;
Output(22) = y_bottom;

Output(23) = t;
Output(24) = Ffoam_top;
Output(25) = Ffoam_bottom;
Output(26) = Fx_tail;
Output(27) = Fy_tail;
Output(28) = T_tail;
Output(29) = Fx_rebound;
Output(30) = Fy_rebound;
Output(31) = Fx_contact;
Output(32) = Fy_contact;
Output(33) = Fy_fricTop;
Output(34) = Fy_fricBottom;
Output(35) = F_hardstopTop;
Output(36) = F_hardstopBottom;

Output(37) = t;
Output(38) = FoamContactTop;
Output(39) = FoamContactBottom;
Output(40) = TailContact;
Output(41) = FootAttached;

Output(42) = KineticEnergy;
Output(43) = KineticEnergy_body;
Output(44) = KineticEnergy_tail;
Output(45) = E_rebound;
Output(46) = E_tail;
Output(47) = GravityPotentialEnergy;

Output(48) = xp;
Output(49) = yp;
Output(50) = gammap;
Output(51) = phip;
Output(52) = xpp;
Output(53) = ypp;
Output(54) = gammapp;
Output(55) = phipp;

Output(56) = Lx_body;
Output(57) = Ly_body;
Output(58) = H_body;
Output(59) = Lx_tail;
Output(60) = Ly_tail;
Output(61) = H_tail;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 7 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files BerkeleyImpact_TallBody.i  (i=1, ..., 7)\n' );
      fprintf( 1, '\n Note: Plots are automatically generated by issuing the OutputPlot command in MotionGenesis\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''BerkeleyImpact_TallBody.1'' );\n' );
      fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t              x              y            x_Acm          y_Acm         x_tail         y_tail          gamma           phi\n' );
      fprintf( 1,                '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)           (rad)          (rad)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros(1,7);
      FileIdentifier(1) = fopen('BerkeleyImpact_TallBody.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.1'); end
      fprintf(FileIdentifier(1), '%% FILE: BerkeleyImpact_TallBody.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              x              y            x_Acm          y_Acm         x_tail         y_tail          gamma           phi\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)           (rad)          (rad)\n\n' );
      FileIdentifier(2) = fopen('BerkeleyImpact_TallBody.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.2'); end
      fprintf(FileIdentifier(2), '%% FILE: BerkeleyImpact_TallBody.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%   x_FoamTop      y_FoamTop    x_FoamBottom   y_FoamBottom   AttachPt_x_WorldFrame AttachPt_y_WorldFrame  x_hardstop        x_Ccm          y_Ccm          x_top          y_top        x_bottom       y_bottom\n' );
      fprintf(FileIdentifier(2), '%%      (m)            (m)            (m)            (m)                (m)                   (m)               (m)            (m)            (m)            (m)            (m)            (m)            (m)\n\n' );
      FileIdentifier(3) = fopen('BerkeleyImpact_TallBody.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.3'); end
      fprintf(FileIdentifier(3), '%% FILE: BerkeleyImpact_TallBody.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t          Ffoam_top    Ffoam_bottom      Fx_tail        Fy_tail        T_tail       Fx_rebound     Fy_rebound     Fx_contact     Fy_contact     Fy_fricTop    Fy_fricBottom  F_hardstopTop F_hardstopBottom\n' );
      fprintf(FileIdentifier(3), '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)         (UNITS)\n\n' );
      FileIdentifier(4) = fopen('BerkeleyImpact_TallBody.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.4'); end
      fprintf(FileIdentifier(4), '%% FILE: BerkeleyImpact_TallBody.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t       FoamContactTop  FoamContactBottom  TailContact   FootAttached\n' );
      fprintf(FileIdentifier(4), '%%   (second)        (UNITS)          (UNITS)         (UNITS)        (UNITS)\n\n' );
      FileIdentifier(5) = fopen('BerkeleyImpact_TallBody.5', 'wt');   if( FileIdentifier(5) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.5'); end
      fprintf(FileIdentifier(5), '%% FILE: BerkeleyImpact_TallBody.5\n%%\n' );
      fprintf(FileIdentifier(5), '%% KineticEnergy KineticEnergy_body KineticEnergy_tail    E_rebound       E_tail     GravityPotentialEnergy\n' );
      fprintf(FileIdentifier(5), '%%      (J)              (J)                (J)              (J)            (J)                (J)\n\n' );
      FileIdentifier(6) = fopen('BerkeleyImpact_TallBody.6', 'wt');   if( FileIdentifier(6) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.6'); end
      fprintf(FileIdentifier(6), '%% FILE: BerkeleyImpact_TallBody.6\n%%\n' );
      fprintf(FileIdentifier(6), '%%      x''             y''           gamma''          phi''            x''''            y''''          gamma''''         phi''''\n' );
      fprintf(FileIdentifier(6), '%%     (m/s)          (m/s)        (rad/sec)      (rad/sec)        (m/s)          (m/s)        (rad/sec)      (rad/sec)\n\n' );
      FileIdentifier(7) = fopen('BerkeleyImpact_TallBody.7', 'wt');   if( FileIdentifier(7) == -1 ), error('Error: unable to open file BerkeleyImpact_TallBody.7'); end
      fprintf(FileIdentifier(7), '%% FILE: BerkeleyImpact_TallBody.7\n%%\n' );
      fprintf(FileIdentifier(7), '%%    Lx_body        Ly_body        H_body         Lx_tail        Ly_tail        H_tail\n' );
      fprintf(FileIdentifier(7), '%%    (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(10:22) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(23:36) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(37:41) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(5), Output(42:47) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(6), Output(48:55) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(7), Output(56:61) );  end
end


%===========================================================================
function WriteNumericalData( fileIdentifier, Output )
%===========================================================================
numberOfOutputQuantities = length( Output );
if( numberOfOutputQuantities > 0 ),
   for( i = 1 : numberOfOutputQuantities ),
      fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
   end
   fprintf( fileIdentifier, '\n' );
end
end



%===========================================================================
function [functionsToEvaluateForEvent, eventTerminatesIntegration1Otherwise0ToContinue, eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1] = EventDetection( t, VAR, uSimulink )
%===========================================================================
% Detects when designated functions are zero or cross zero with positive or negative slope.
% Step 1: Uncomment call to mdlDerivatives and mdlOutputs.
% Step 2: Change functionsToEvaluateForEvent,                      e.g., change  []  to  [t - 5.67]  to stop at t = 5.67.
% Step 3: Change eventTerminatesIntegration1Otherwise0ToContinue,  e.g., change  []  to  [1]  to stop integrating.
% Step 4: Change eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1,  e.g., change  []  to  [1].
% Step 5: Possibly modify function EventDetectedByIntegrator (if eventTerminatesIntegration1Otherwise0ToContinue is 0).
%---------------------------------------------------------------------------
% mdlDerivatives( t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
% mdlOutputs(     t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
functionsToEvaluateForEvent = [];
eventTerminatesIntegration1Otherwise0ToContinue = [];
eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1 = [];
eventDetectedByIntegratorTerminate1OrContinue0 = eventTerminatesIntegration1Otherwise0ToContinue;
end


%===========================================================================
function [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents )
%===========================================================================
isIntegrationFinished = eventDetectedByIntegratorTerminate1OrContinue0( nIndexOfEvents );
if( ~isIntegrationFinished ),
   SetNamedQuantitiesFromMatrix( VAR );
%  Put code here to modify how integration continues.
   VAR = SetMatrixFromNamedQuantities;
end
end



%===========================================================================
function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile )
%===========================================================================
OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',tStep, 'Events',@EventDetection );
t = tInitial;                 epsilonT = 0.001*tStep;                   tFinalMinusEpsilonT = tFinal - epsilonT;
printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + tStep;
printCounterFile   = 0;       isIntegrationFinished = 0;
mdlDerivatives( t, VAR, 0 );
while 1,
   if( (integrateForward && t >= tFinalMinusEpsilonT) || (~integrateForward && t <= tFinalMinusEpsilonT) ), isIntegrationFinished = 1;  end
   shouldPrintToScreen = printIntScreen && ( isIntegrationFinished || printCounterScreen <= 0.01 );
   shouldPrintToFile   = printIntFile   && ( isIntegrationFinished || printCounterFile   <= 0.01 );
   if( isIntegrationFinished || shouldPrintToScreen || shouldPrintToFile ),
      Output = mdlOutputs( t, VAR, 0 );
      OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );
      if( isIntegrationFinished ), break;  end
      if( shouldPrintToScreen ), printCounterScreen = printIntScreen;  end
      if( shouldPrintToFile ),   printCounterFile   = printIntFile;    end
   end
   [TimeOdeArray, VarOdeArray, timeEventOccurredInIntegrationStep, nStatesArraysAtEvent, nIndexOfEvents] = ode45( @mdlDerivatives, [t tAtEndOfIntegrationStep], VAR, OdeMatlabOptions, 0 );
   if( isempty(timeEventOccurredInIntegrationStep) ),
      t = TimeOdeArray( length(TimeOdeArray) );
      VAR = VarOdeArray( length(TimeOdeArray), : );
      printCounterScreen = printCounterScreen - 1;
      printCounterFile   = printCounterFile   - 1;
      if( abs(tAtEndOfIntegrationStep - t) >= abs(epsilonT) ), warning('numerical integration failed'); break;  end
      tAtEndOfIntegrationStep = t + tStep;
   else
      t = timeEventOccurredInIntegrationStep;
      VAR = nStatesArraysAtEvent;
      printCounterScreen = 0;
      printCounterFile   = 0;
      [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents );
   end
end
end

%================================================
end    % End of function BerkeleyImpact_TallBody
%================================================
