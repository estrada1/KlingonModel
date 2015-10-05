function [t,VAR,Output] = BerkeleyImpact_LongBody
%===========================================================================
% File: KlingonModel/BerkeleyImpact_LongBody.m created on Wed Sep 23 2015 by MotionGenesis 5.7.
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
eventDetectedByIntegratorTerminate1OrContinue0 = [];
mew=0; gammapp=0; phipp=0; xpp=0; ypp=0; Ffoam_bottom=0; Ffoam_top=0; Fx_tail=0; T_tail=0; AttachPt_xp=0; AttachPt_yp=0; AttachPt_x_WorldFrame=0; AttachPt_y_WorldFrame=0; FoamContactBottom=0; FoamContactTop=0; FootAttached=0; FootContact=0;
 Fx_contact=0; Fx_rebound=0; Fy_contact=0; Fy_fricBottom=0; Fy_fricTop=0; Fy_rebound=0; F_hardstopBottom=0; F_hardstopTop=0; GravityPotentialEnergy=0; HardstopContactBottom=0; HardstopContactTop=0; KineticEnergy=0; MechanicalEnergy=0; rebound_x=0;
 rebound_y=0; TailContact=0; vBottom=0; vTop=0; x_Acm=0; x_Ccm=0; x_FoamBottom=0; x_FoamTop=0; x_hardstop=0; x_tail=0; y_Acm=0; y_Ccm=0; y_FoamBottom=0; y_FoamTop=0; y_tail=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
Awidth                          =  3;                      % cm                  Constant
b_foam                          =  2.5;                    % N*sec/m             Constant
b_hardstop                      =  25;                     % N/m/s               Constant
b_tail                          =  0;%0.1;                    % N/s                 Constant
d_body                          =  5;                      % cm                  Constant
d_foam                          =  3;                      % cm                  Constant
Fpre_rebound                    =  2;                      % N                   Constant
g                               =  9.80665;                % m/sec^2             Constant
IAzz                            =  0.0029;                 % kg*m^2              Constant
ICzz                            =  0.000002;               % kg*m^2              Constant
k_foam                          =  85;                     % N/m                 Constant
k_hardstop                      =  10000;                  % N/m                 Constant
k_rebound                       =  14.28571428571428;      % N/m                 Constant
k_tail                          =  0;%1;                    % N                   Constant
l_foam                          =  0.02;                   % m                   Constant
l_tail                          =  0.23;                   % m                   Constant
mA                              =  .2;                     % kg                  Constant
mC                              =  .01;                    % kg                  Constant
mQ                              =  .02;                    % kg                  Constant
phin                            = -50;                     % deg                 Constant

gamma                           = -15;                     % deg                 Initial Value
phi                             = -30;                     % deg                 Initial Value
x                               = -.1;                     % m                   Initial Value
y                               =  4;                      % m                   Initial Value
gammap                          =  0;                      % rad/sec             Initial Value
phip                            =  0;                      % rad/sec             Initial Value
xp                              =  1;                      % m/sec               Initial Value
yp                              =  .5;                     % m/sec               Initial Value
AttachPt_x                      =  0;                      % UNITS               Initial Value
AttachPt_y                      =  0;                      % UNITS               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  1;                      % sec                 Final Time
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
mew = -0.1;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files



%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
FootContact = ceil(0.5*sign(x));
FoamContactTop = ceil(0.5*sign(x-d_foam*sin(gamma)));
FoamContactBottom = ceil(0.5*sign(x+d_foam*sin(gamma)));
HardstopContactTop = ceil(0.5*sign(x-d_foam*sin(gamma)-l_foam*cos(gamma)));
HardstopContactBottom = ceil(0.5*sign(x+(d_body+d_foam)*sin(gamma)-l_foam*cos(gamma)));
TailContact = ceil(0.5*sign(x+(d_body+d_foam)*sin(gamma)-l_tail*cos(gamma+phi)-(Awidth+l_foam)*cos(gamma)));
FootAttached = ceil(0.5*sign(sqrt(AttachPt_x^2+AttachPt_y^2))+0.5*FootContact);
F_hardstopTop = -HardstopContactTop*(k_hardstop*(x-d_foam*sin(gamma)-l_foam*cos(gamma))-b_hardstop*(d_foam*cos(gamma)*gammap-xp-l_foam*sin(gamma)*gammap));
F_hardstopBottom = HardstopContactBottom*(k_hardstop*(l_foam*cos(gamma)-x-(d_body+d_foam)*sin(gamma))-b_hardstop*(xp+l_foam*sin(gamma)*gammap+(d_body+d_foam)*cos(gamma)*gammap));
vTop = yp - d_foam*sin(gamma)*gammap - l_foam*cos(gamma)*gammap;
vBottom = yp + (d_body+d_foam)*sin(gamma)*gammap - l_foam*cos(gamma)*gammap;
rebound_x = AttachPt_x + (Awidth+l_foam)*cos(gamma);
rebound_y = AttachPt_y + (Awidth+l_foam)*sin(gamma);
Fx_rebound = (k_rebound+Fpre_rebound/(1.0E-8+sqrt((Awidth+l_foam)^2+AttachPt_x^2+AttachPt_y^2+2*(Awidth+l_foam)*AttachPt_x*cos(gamma)+2*(Awidth+l_foam)*AttachPt_y*sin(gamma))))*FootAttached*rebound_x;
Fy_rebound = (k_rebound+Fpre_rebound/(1.0E-8+sqrt((Awidth+l_foam)^2+AttachPt_x^2+AttachPt_y^2+2*(Awidth+l_foam)*AttachPt_x*cos(gamma)+2*(Awidth+l_foam)*AttachPt_y*sin(gamma))))*FootAttached*rebound_y;

% Quantities that were specified
T_tail = k_tail*(phin-phi) - b_tail*phip;
Ffoam_top = -FoamContactTop*(k_foam*(x-d_foam*sin(gamma))+b_foam*(xp-d_foam*cos(gamma)*gammap));
Ffoam_bottom = -FoamContactBottom*(k_foam*(x+d_foam*sin(gamma))+b_foam*(xp+d_foam*cos(gamma)*gammap));
Fx_tail = k_hardstop*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-x-(d_body+d_foam)*sin(gamma))*TailContact;
AttachPt_xp = -FootAttached*xp;
AttachPt_yp = -FootAttached*yp;

Fy_fricTop = -mew*vTop*(Ffoam_top+F_hardstopTop)/(1.0E-6+abs(vTop));
Fy_fricBottom = -mew*vBottom*(Ffoam_bottom+F_hardstopBottom)/(1.0E-6+abs(vBottom));
xpp = -(((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(  ...
gamma)+(d_body+d_foam)*cos(gamma)))-l_tail*mQ*(mA+mC+mQ)*sin(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*  ...
l_tail*(d_body+d_foam)*sin(phi)))-l_tail*mQ*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))*(cos(gamma+phi)*(mA*(  ...
Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+  ...
d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))))*(T_tail+l_tail*(Fx_tail*sin(gamma+phi)+g*mQ*cos(gamma+phi))-l_tail*mQ*((Awidth+l_foam)*sin(phi)-(d_body+d_foam)*cos(phi))*gammap^2)+(l_tail^2*mQ^2*cos(  ...
gamma+phi)^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))+(ICzz+l_tail*mQ*(l_tail+(  ...
Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-2*l_tail*mQ*cos(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(  ...
d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))+(ICzz+mQ*l_tail^2)*((mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(  ...
gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+  ...
2*l_tail*(d_body+d_foam)*sin(phi)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))*gammap^2-Fx_tail-Ffoam_bottom*cos(gamma)-Ffoam_top*cos(gamma)-Fx_rebound-mQ*((d_body+d_foam)*sin(  ...
gamma)*gammap^2-(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*cos(gamma+phi)*(gammap+phip)^2)-F_hardstopBottom-F_hardstopTop)+((mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(  ...
gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-l_tail*mQ*(mA+mC+mQ)*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-l_tail^2*mQ^2*cos(gamma+phi)*(cos(  ...
gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(  ...
gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))))*(d_foam*Ffoam_top+g*mQ*(d_body+d_foam)*sin(gamma)+d_foam*sin(gamma)*Fy_fricTop+d_foam*cos(gamma)*F_hardstopTop+l_foam*cos(gamma)*  ...
Fy_fricBottom+l_foam*cos(gamma)*Fy_fricTop-d_foam*Ffoam_bottom-l_tail*Fx_tail*sin(gamma+phi)-(Awidth+l_foam)*Fx_tail*sin(gamma)-(d_body+d_foam)*Fx_tail*cos(gamma)-g*l_tail*mQ*cos(gamma+phi)-g*mQ*(Awidth+l_foam)*cos(gamma)-g*mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))-(Awidth+l_foam)*(sin(gamma)*Fx_rebound+cos(gamma)*(g*mA-Fy_rebound))-l_tail*mQ*((d_body+d_foam)*cos(phi)*gammap^2+(Awidth+l_foam)*sin(phi)*(gammap+phip)^2-(Awidth+l_foam)*sin(phi)*  ...
gammap^2-(d_body+d_foam)*cos(phi)*(gammap+phip)^2)-l_foam*sin(gamma)*F_hardstopBottom-l_foam*sin(gamma)*F_hardstopTop-(d_body+d_foam)*sin(gamma)*Fy_fricBottom-(d_body+d_foam)*cos(gamma)*F_hardstopBottom)-(l_tail*mQ*sin(gamma+phi)*(ICzz+  ...
l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-(ICzz+mQ*  ...
l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(  ...
gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-l_tail*mQ*cos(gamma+phi)*(l_tail*mQ*sin(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+  ...
d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))-(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+l_foam)*sin(  ...
gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)+2*l_tail*sin(gamma+phi)))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*((Awidth+l_foam)*sin(  ...
gamma)+(d_body+d_foam)*cos(gamma))*gammap^2+mQ*((Awidth+l_foam)*sin(gamma)*gammap^2+(d_body+d_foam)*cos(gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Ffoam_bottom*sin(gamma)-Ffoam_top*sin(gamma)-Fy_rebound-Fy_fricBottom-  ...
Fy_fricTop))/(l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))+(mA+  ...
mC+mQ)^2*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*  ...
l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi))))-l_tail^2*mQ^2*(cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*  ...
cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))^2-(mA+mC+mQ)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+  ...
l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(ICzz+  ...
mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))))-(mA+mC+mQ)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(  ...
ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))));
ypp = ((l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(  ...
d_body+d_foam)*cos(gamma)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(mA*(Awidth+  ...
l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-l_tail*mQ*sin(gamma+phi)*(l_tail*mQ*cos(gamma+phi)*(IAzz+ICzz+mA*(Awidth+  ...
l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))-(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*  ...
sin(phi)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)+2*l_tail*cos(gamma+phi)-(d_body+d_foam)*sin(gamma)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))*gammap^2-Fx_tail-Ffoam_bottom*cos(gamma)-Ffoam_top*cos(gamma)-Fx_rebound-mQ*((d_body+d_foam)*sin(gamma)*gammap^2-(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*cos(gamma+phi)*(gammap+phip)^  ...
2)-F_hardstopBottom-F_hardstopTop)-(l_tail*mQ*(mA+mC+mQ)*cos(gamma+phi)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*  ...
l_tail*(d_body+d_foam)*sin(phi)))-(mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-l_tail*mQ*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(  ...
cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(  ...
gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))))*(T_tail+l_tail*(Fx_tail*sin(gamma+phi)+g*mQ*cos(gamma+phi))-l_tail*mQ*((Awidth+l_foam)*sin(phi)-(d_body+d_foam)*cos(phi))*gammap^2)-(  ...
l_tail^2*mQ^2*sin(gamma+phi)^2*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))+(ICzz+  ...
l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-2*l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*  ...
sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))))+(ICzz+mQ*l_tail^2)*((mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(  ...
l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+  ...
l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))*gammap^2+mQ*((Awidth+l_foam)*sin(gamma)*gammap^2+(d_body+d_foam)*cos(  ...
gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Ffoam_bottom*sin(gamma)-Ffoam_top*sin(gamma)-Fy_rebound-Fy_fricBottom-Fy_fricTop)-(l_tail*mQ*(mA+mC+mQ)*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+  ...
d_foam)*sin(phi)))-(mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-l_tail^2*mQ^2*  ...
sin(gamma+phi)*(cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))))*(d_foam*Ffoam_top+g*mQ*(d_body+d_foam)*sin(gamma)+d_foam*sin(gamma)*Fy_fricTop+d_foam*cos(gamma)*F_hardstopTop+l_foam*cos(  ...
gamma)*Fy_fricBottom+l_foam*cos(gamma)*Fy_fricTop-d_foam*Ffoam_bottom-l_tail*Fx_tail*sin(gamma+phi)-(Awidth+l_foam)*Fx_tail*sin(gamma)-(d_body+d_foam)*Fx_tail*cos(gamma)-g*l_tail*mQ*cos(gamma+phi)-g*mQ*(Awidth+l_foam)*cos(gamma)-g*mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))-(Awidth+l_foam)*(sin(gamma)*Fx_rebound+cos(gamma)*(g*mA-Fy_rebound))-l_tail*mQ*((d_body+d_foam)*cos(phi)*gammap^2+(Awidth+l_foam)*sin(phi)*(gammap+phip)^2-(Awidth+l_foam)*sin(phi)*  ...
gammap^2-(d_body+d_foam)*cos(phi)*(gammap+phip)^2)-l_foam*sin(gamma)*F_hardstopBottom-l_foam*sin(gamma)*F_hardstopTop-(d_body+d_foam)*sin(gamma)*Fy_fricBottom-(d_body+d_foam)*cos(gamma)*F_hardstopBottom))/(l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+  ...
ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))+(mA+mC+mQ)^2*((ICzz+l_tail*mQ*(l_tail+(Awidth+  ...
l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(  ...
d_body+d_foam)*sin(phi))))-l_tail^2*mQ^2*(cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+  ...
l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))^2-(mA+mC+mQ)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(  ...
gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(  ...
gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))))-(mA+mC+mQ)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+  ...
d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(ICzz+mQ*l_tail^2)*(mA*(  ...
Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))));
gammapp = -((mA+mC+mQ)*(l_tail^2*mQ^2-(mA+mC+mQ)*(ICzz+mQ*l_tail^2))*(d_foam*Ffoam_top+g*mQ*(d_body+d_foam)*sin(gamma)+d_foam*sin(gamma)*Fy_fricTop+d_foam*cos(gamma)*F_hardstopTop+l_foam*cos(gamma)*Fy_fricBottom+l_foam*cos(gamma)*  ...
Fy_fricTop-d_foam*Ffoam_bottom-l_tail*Fx_tail*sin(gamma+phi)-(Awidth+l_foam)*Fx_tail*sin(gamma)-(d_body+d_foam)*Fx_tail*cos(gamma)-g*l_tail*mQ*cos(gamma+phi)-g*mQ*(Awidth+l_foam)*cos(gamma)-g*mC*((Awidth+l_foam)*cos(gamma)-(d_body+  ...
d_foam)*sin(gamma))-(Awidth+l_foam)*(sin(gamma)*Fx_rebound+cos(gamma)*(g*mA-Fy_rebound))-l_tail*mQ*((d_body+d_foam)*cos(phi)*gammap^2+(Awidth+l_foam)*sin(phi)*(gammap+phip)^2-(Awidth+l_foam)*sin(phi)*gammap^2-(d_body+d_foam)*cos(phi)*(  ...
gammap+phip)^2)-l_foam*sin(gamma)*F_hardstopBottom-l_foam*sin(gamma)*F_hardstopTop-(d_body+d_foam)*sin(gamma)*Fy_fricBottom-(d_body+d_foam)*cos(gamma)*F_hardstopBottom)+((mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((  ...
Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-l_tail^2*mQ^2*cos(gamma+phi)^2*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(  ...
d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-l_tail*mQ*sin(gamma+phi)*((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-l_tail*mQ*cos(gamma+phi)*(mA*(  ...
Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(  ...
gamma))*gammap^2-Fx_tail-Ffoam_bottom*cos(gamma)-Ffoam_top*cos(gamma)-Fx_rebound-mQ*((d_body+d_foam)*sin(gamma)*gammap^2-(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*cos(gamma+phi)*(gammap+phip)^2)-F_hardstopBottom-F_hardstopTop)-(mA+mC+  ...
mQ)*((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-l_tail^2*mQ^2-l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+  ...
l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-l_tail*mQ*cos(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))*(T_tail+  ...
l_tail*(Fx_tail*sin(gamma+phi)+g*mQ*cos(gamma+phi))-l_tail*mQ*((Awidth+l_foam)*sin(phi)-(d_body+d_foam)*cos(phi))*gammap^2)-((mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(  ...
gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-l_tail^2*mQ^2*sin(gamma+phi)^2*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(  ...
gamma)-(d_body+d_foam)*sin(gamma)))-l_tail*mQ*cos(gamma+phi)*((mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(  ...
d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))*gammap^2+mQ*((Awidth+l_foam)*sin(  ...
gamma)*gammap^2+(d_body+d_foam)*cos(gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Ffoam_bottom*sin(gamma)-Ffoam_top*sin(gamma)-Fy_rebound-Fy_fricBottom-Fy_fricTop))/(l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((  ...
Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))+(mA+mC+mQ)^2*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(  ...
phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi))))-l_tail^2*mQ^2*(  ...
cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(  ...
gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))^2-(mA+mC+mQ)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+  ...
l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(  ...
d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))))-(mA+mC+mQ)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(  ...
gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))));
phipp = ((mA+mC+mQ)*((mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))^2+(mA*(Awidth+l_foam)*cos(gamma)+mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(  ...
Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi))))*(T_tail+l_tail*(Fx_tail*sin(gamma+phi)+g*mQ*cos(gamma+phi))-l_tail*mQ*((Awidth+l_foam)*sin(phi)-(d_body+d_foam)*cos(phi))*gammap^  ...
2)+(l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-(mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+  ...
l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-l_tail*mQ*cos(gamma+phi)*((mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+  ...
l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(  ...
Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))))*(g*mA+g*mC+g*mQ+mA*(Awidth+l_foam)*sin(gamma)*gammap^2+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))*gammap^2+mQ*((  ...
Awidth+l_foam)*sin(gamma)*gammap^2+(d_body+d_foam)*cos(gamma)*gammap^2+l_tail*sin(gamma+phi)*(gammap+phip)^2)-Ffoam_bottom*sin(gamma)-Ffoam_top*sin(gamma)-Fy_rebound-Fy_fricBottom-Fy_fricTop)-(mA+mC+mQ)*((mA+mC+mQ)*(ICzz+l_tail*mQ*(  ...
l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-l_tail^2*mQ^2-l_tail*mQ*sin(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*  ...
cos(gamma)))-l_tail*mQ*cos(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))*(d_foam*Ffoam_top+g*mQ*(d_body+d_foam)*sin(  ...
gamma)+d_foam*sin(gamma)*Fy_fricTop+d_foam*cos(gamma)*F_hardstopTop+l_foam*cos(gamma)*Fy_fricBottom+l_foam*cos(gamma)*Fy_fricTop-d_foam*Ffoam_bottom-l_tail*Fx_tail*sin(gamma+phi)-(Awidth+l_foam)*Fx_tail*sin(gamma)-(d_body+d_foam)*  ...
Fx_tail*cos(gamma)-g*l_tail*mQ*cos(gamma+phi)-g*mQ*(Awidth+l_foam)*cos(gamma)-g*mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))-(Awidth+l_foam)*(sin(gamma)*Fx_rebound+cos(gamma)*(g*mA-Fy_rebound))-l_tail*mQ*((d_body+d_foam)*  ...
cos(phi)*gammap^2+(Awidth+l_foam)*sin(phi)*(gammap+phip)^2-(Awidth+l_foam)*sin(phi)*gammap^2-(d_body+d_foam)*cos(phi)*(gammap+phip)^2)-l_foam*sin(gamma)*F_hardstopBottom-l_foam*sin(gamma)*F_hardstopTop-(d_body+d_foam)*sin(gamma)*  ...
Fy_fricBottom-(d_body+d_foam)*cos(gamma)*F_hardstopBottom)-(l_tail*mQ*cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(  ...
d_body+d_foam)*cos(gamma)))*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))-(mA+mC+mQ)*(ICzz+l_tail*mQ*(l_tail+(  ...
Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))-l_tail*mQ*sin(  ...
gamma+phi)*((mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))^2-(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((  ...
Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))))*(mA*(Awidth+l_foam)*cos(gamma)*gammap^2+mC*((Awidth+l_foam)*cos(gamma)-(  ...
d_body+d_foam)*sin(gamma))*gammap^2-Fx_tail-Ffoam_bottom*cos(gamma)-Ffoam_top*cos(gamma)-Fx_rebound-mQ*((d_body+d_foam)*sin(gamma)*gammap^2-(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*cos(gamma+phi)*(gammap+phip)^2)-F_hardstopBottom-  ...
F_hardstopTop))/(l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi)))+(  ...
mA+mC+mQ)^2*((ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))^2-(ICzz+mQ*l_tail^2)*(IAzz+ICzz+mA*(Awidth+l_foam)^2+mC*((Awidth+l_foam)^2+(d_body+d_foam)^2)+mQ*(l_tail^2+(Awidth+l_foam)^2+(d_body+d_foam)^2+2*  ...
l_tail*(Awidth+l_foam)*cos(phi)+2*l_tail*(d_body+d_foam)*sin(phi))))-l_tail^2*mQ^2*(cos(gamma+phi)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*  ...
cos(gamma)))-sin(gamma+phi)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))))^2-(mA+mC+mQ)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+  ...
l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma)))*(2*l_tail*mQ*sin(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(ICzz+  ...
mQ*l_tail^2)*(mA*(Awidth+l_foam)*sin(gamma)+mC*((Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))+mQ*(l_tail*sin(gamma+phi)+(Awidth+l_foam)*sin(gamma)+(d_body+d_foam)*cos(gamma))))-(mA+mC+mQ)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((  ...
Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))*(2*l_tail*mQ*cos(gamma+phi)*(ICzz+l_tail*mQ*(l_tail+(Awidth+l_foam)*cos(phi)+(d_body+d_foam)*sin(phi)))-(  ...
ICzz+mQ*l_tail^2)*(mA*(Awidth+l_foam)*cos(gamma)+mC*((Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma))+mQ*(l_tail*cos(gamma+phi)+(Awidth+l_foam)*cos(gamma)-(d_body+d_foam)*sin(gamma)))));

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
x_tail = x + (d_body+d_foam)*sin(gamma) - l_tail*cos(gamma+phi) - (Awidth+l_foam)*cos(gamma);
y_tail = y - l_tail*sin(gamma+phi) - (Awidth+l_foam)*sin(gamma) - (d_body+d_foam)*cos(gamma);
x_Acm = x - (Awidth+l_foam)*cos(gamma);
y_Acm = y - (Awidth+l_foam)*sin(gamma);
x_FoamTop = x - d_foam*sin(gamma);
y_FoamTop = y + d_foam*cos(gamma);
x_FoamBottom = x + d_foam*sin(gamma);
y_FoamBottom = y - d_foam*cos(gamma);
x_hardstop = x - d_foam*sin(gamma) - l_foam*cos(gamma);
x_Ccm = x + (d_body+d_foam)*sin(gamma) - (Awidth+l_foam)*cos(gamma);
y_Ccm = y - (Awidth+l_foam)*sin(gamma) - (d_body+d_foam)*cos(gamma);
AttachPt_x_WorldFrame = AttachPt_x + x;
AttachPt_y_WorldFrame = AttachPt_y + y;
Fx_contact = Ffoam_bottom + Ffoam_top + Fx_tail + Fx_rebound + F_hardstopBottom + F_hardstopTop;
Fy_contact = Fy_rebound + Fy_fricBottom + Fy_fricTop;
KineticEnergy = 0.5*IAzz*gammap^2 + 0.5*ICzz*(gammap+phip)^2 + 0.5*mA*(xp^2+yp^2+(Awidth+l_foam)^2*gammap^2+2*(Awidth+l_foam)*sin(gamma)*gammap*xp-2*(Awidth+l_foam)*cos(gamma)*gammap*yp) - 0.5*mC*(2*(Awidth+l_foam)*cos(gamma)*gammap*yp-  ...
xp^2-yp^2-(Awidth+l_foam)^2*gammap^2-(d_body+d_foam)^2*gammap^2-2*(Awidth+l_foam)*sin(gamma)*gammap*xp-2*(d_body+d_foam)*sin(gamma)*gammap*yp-2*(d_body+d_foam)*cos(gamma)*gammap*xp) - 0.5*mQ*(2*(Awidth+l_foam)*cos(gamma)*gammap*yp+2*  ...
l_tail*cos(gamma+phi)*yp*(gammap+phip)-xp^2-yp^2-(Awidth+l_foam)^2*gammap^2-(d_body+d_foam)^2*gammap^2-l_tail^2*(gammap+phip)^2-2*(Awidth+l_foam)*sin(gamma)*gammap*xp-2*(d_body+d_foam)*sin(gamma)*gammap*yp-2*(d_body+d_foam)*cos(gamma)*  ...
gammap*xp-2*l_tail*sin(gamma+phi)*xp*(gammap+phip)-2*l_tail*(Awidth+l_foam)*cos(phi)*gammap*(gammap+phip)-2*l_tail*(d_body+d_foam)*sin(phi)*gammap*(gammap+phip));
GravityPotentialEnergy = g*(mA*(y-(Awidth+l_foam)*sin(gamma))+mC*(y-(Awidth+l_foam)*sin(gamma)-(d_body+d_foam)*cos(gamma))+mQ*(y-l_tail*sin(gamma+phi)-(Awidth+l_foam)*sin(gamma)-(d_body+d_foam)*cos(gamma)));
MechanicalEnergy = GravityPotentialEnergy + KineticEnergy;

Output = zeros(1,38);
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

Output(19) = t;
Output(20) = Ffoam_top;
Output(21) = Ffoam_bottom;
Output(22) = Fx_tail;
Output(23) = T_tail;
Output(24) = Fx_rebound;
Output(25) = Fy_rebound;
Output(26) = Fx_contact;
Output(27) = Fy_contact;
Output(28) = Fy_fricTop;
Output(29) = Fy_fricBottom;
Output(30) = F_hardstopTop;
Output(31) = F_hardstopBottom;

Output(32) = t;
Output(33) = FoamContactTop;
Output(34) = FoamContactBottom;
Output(35) = TailContact;
Output(36) = FootAttached;

Output(37) = KineticEnergy;
Output(38) = MechanicalEnergy;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 5 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files KlingonModel/BerkeleyImpact_LongBody.i  (i=1, ..., 5)\n' );
      fprintf( 1, '\n Note: Plots are automatically generated by issuing the OutputPlot command in MotionGenesis\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''KlingonModel/BerkeleyImpact_LongBody.1'' );\n' );
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
      FileIdentifier = zeros(1,5);
      FileIdentifier(1) = fopen('KlingonModel/BerkeleyImpact_LongBody.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file KlingonModel/BerkeleyImpact_LongBody.1'); end
      fprintf(FileIdentifier(1), '%% FILE: KlingonModel/BerkeleyImpact_LongBody.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              x              y            x_Acm          y_Acm         x_tail         y_tail          gamma           phi\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)           (rad)          (rad)\n\n' );
      FileIdentifier(2) = fopen('KlingonModel/BerkeleyImpact_LongBody.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file KlingonModel/BerkeleyImpact_LongBody.2'); end
      fprintf(FileIdentifier(2), '%% FILE: KlingonModel/BerkeleyImpact_LongBody.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%   x_FoamTop      y_FoamTop    x_FoamBottom   y_FoamBottom   AttachPt_x_WorldFrame AttachPt_y_WorldFrame  x_hardstop        x_Ccm          y_Ccm\n' );
      fprintf(FileIdentifier(2), '%%      (m)            (m)            (m)            (m)                (m)                   (m)               (m)            (m)            (m)\n\n' );
      FileIdentifier(3) = fopen('KlingonModel/BerkeleyImpact_LongBody.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file KlingonModel/BerkeleyImpact_LongBody.3'); end
      fprintf(FileIdentifier(3), '%% FILE: KlingonModel/BerkeleyImpact_LongBody.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t          Ffoam_top    Ffoam_bottom      Fx_tail        T_tail       Fx_rebound     Fy_rebound     Fx_contact     Fy_contact     Fy_fricTop    Fy_fricBottom  F_hardstopTop F_hardstopBottom\n' );
      fprintf(FileIdentifier(3), '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)         (UNITS)\n\n' );
      FileIdentifier(4) = fopen('KlingonModel/BerkeleyImpact_LongBody.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file KlingonModel/BerkeleyImpact_LongBody.4'); end
      fprintf(FileIdentifier(4), '%% FILE: KlingonModel/BerkeleyImpact_LongBody.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t       FoamContactTop  FoamContactBottom  TailContact   FootAttached\n' );
      fprintf(FileIdentifier(4), '%%   (second)        (UNITS)          (UNITS)         (UNITS)        (UNITS)\n\n' );
      FileIdentifier(5) = fopen('KlingonModel/BerkeleyImpact_LongBody.5', 'wt');   if( FileIdentifier(5) == -1 ), error('Error: unable to open file KlingonModel/BerkeleyImpact_LongBody.5'); end
      fprintf(FileIdentifier(5), '%% FILE: KlingonModel/BerkeleyImpact_LongBody.5\n%%\n' );
      fprintf(FileIdentifier(5), '%% KineticEnergy MechanicalEnergy\n' );
      fprintf(FileIdentifier(5), '%%      (J)             (J)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(10:18) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(19:31) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(32:36) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(5), Output(37:38) );  end
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



%=============================================================
end    % End of function KlingonModel/BerkeleyImpact_LongBody
%=============================================================
