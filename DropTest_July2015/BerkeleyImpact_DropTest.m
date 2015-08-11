function [t,VAR,Output] = BerkeleyImpact_DropTest
%===========================================================================
% File: BerkeleyImpact_DropTest.m created on Sat Aug  1 2015 by MotionGenesis 5.4.
% Advanced Student Licensee: Matt Estrada (until December 2015).
% Portions copyright (c) 1988-2012 Paul Mitiguy and 2009-2012 Motion Genesis.
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
gammapp=0; phipp=0; xpp=0; ypp=0; Ffoam_bottom=0; Ffoam_top=0; Fx_tail=0; F_hardstop=0; T_tail=0; AttachPt_xp=0; AttachPt_yp=0; AttachPt_x_WorldFrame=0; AttachPt_y_WorldFrame=0; Fcontact=0; Fmag_rebound=0; FoamContactBottom=0; FoamContactTop=0;
 FootAttached=0; FootContact=0; Fx_rebound=0; Fy_rebound=0; HardstopContact=0; TailContact=0; x_Acm=0; x_FoamBottom=0; x_FoamTop=0; x_hardstop=0; x_tail=0; y_Acm=0; y_FoamBottom=0; y_FoamTop=0; y_tail=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
Awidth                          =  1;                      % cm                  Constant
b_foam                          =  2.5;                    % N*sec/m             Constant
b_hardstop                      =  25;                     % N/m/s               Constant
d_foam                          =  3;                      % cm                  Constant
Fpre_rebound                    =  2;                      % N                   Constant
g                               =  9.80665;                % m/sec^2             Constant
IAzz                            =  0.0029;                 % kg*m^2              Constant
ICzz                            =  0.000002;               % kg*m^2              Constant
k_foam                          =  85;                     % N/m                 Constant
k_hardstop                      =  10000;                  % N/m                 Constant
k_rebound                       =  14.28571428571428;      % N/m                 Constant
k_tail                          =  0;                      % N                   Constant
l_foam                          =  0.02;                   % m                   Constant
l_tail                          =  0.23;                   % m                   Constant
mA                              =  .2;                     % kg                  Constant
mC                              =  .01;                    % kg                  Constant
mQ                              =  .02;                    % kg                  Constant

gamma                           =  0;                      % deg                 Initial Value
phi                             =  0;                      % deg                 Initial Value
x                               = -.5;                     % m                   Initial Value
y                               =  4;                      % m                   Initial Value
gammap                          =  0;                      % rad/sec             Initial Value
phip                            =  0;                      % rad/sec             Initial Value
xp                              =  0;                      % m/sec               Initial Value
yp                              =  0;                      % m/sec               Initial Value
AttachPt_x                      =  0;                      % UNITS               Initial Value
AttachPt_y                      =  0;                      % UNITS               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  1;                      % sec                 Final Time
integStp                        =  0.005;                  % sec                 Integration Step
printIntScreen                  =  1;                      % 0 or +integer       Print-Integer
printIntFile                    =  1;                      % 0 or +integer       Print-Integer
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-07;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions.  UnitSystem: kilogram, meter, second
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
Awidth = Awidth * 0.01;
d_foam = d_foam * 0.01;
gamma = gamma * DEGtoRAD;
phi = phi * DEGtoRAD;

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, integStp, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files



%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
FootContact = ceil(0.5*sign(x));
FoamContactTop = ceil(0.5*sign(x-d_foam*sin(gamma)));
FoamContactBottom = ceil(0.5*sign(x+d_foam*sin(gamma)));
HardstopContact = ceil(0.5*sign(x-l_foam*cos(gamma)));
TailContact = ceil(0.5*sign(x-l_tail*cos(gamma+phi)));
FootAttached = ceil(0.5*sign(sqrt(AttachPt_x^2+AttachPt_y^2))+0.5*FootContact);
Fx_rebound = (k_rebound*AttachPt_x+Fpre_rebound*cos(gamma))*FootAttached;
Fy_rebound = (k_rebound*AttachPt_y+Fpre_rebound*sin(gamma))*FootAttached;


% Quantities that were specified
T_tail = -k_tail*phi;
Ffoam_top = -FoamContactTop*(k_foam*(x-d_foam*sin(gamma))+b_foam*(xp-d_foam*cos(gamma)*gammap));
Ffoam_bottom = -FoamContactBottom*(k_foam*(x+d_foam*sin(gamma))+b_foam*(xp+d_foam*cos(gamma)*gammap));
Fx_tail = -k_hardstop*(x-l_tail*cos(gamma+phi))*TailContact;
AttachPt_xp = -FootAttached*xp;
AttachPt_yp = -FootAttached*yp;
F_hardstop = -HardstopContact*(k_hardstop*(x-l_foam*cos(gamma))+b_hardstop*(xp+l_foam*sin(gamma)*gammap));

xpp = -((Awidth+l_foam)*(mA+mC)*(l_tail^2*mQ^2*sin(phi)*cos(gamma+phi)+(mA+mC+mQ)*(ICzz+mQ*l_tail^2)*sin(gamma))*(d_foam*Ffoam_top-d_foam*Ffoam_bottom-g*l_tail*mQ*sin(gamma+phi)-g*mA*(Awidth+l_foam)*sin(gamma)-g*mC*(Awidth+l_foam)*sin(  ...
gamma))-(T_tail+g*l_tail*mQ*sin(gamma+phi))*(l_tail*mQ*(mA+mC+mQ)*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)*sin(gamma+phi)-(mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(  ...
Awidth+l_foam)*sin(gamma))-l_tail*mQ*(Awidth+l_foam)*(mA+mC)*sin(phi)*(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)))-((Awidth+l_foam)*(mA+mC)*(ICzz+mQ*l_tail^2)*sin(gamma)*(l_tail*mQ*cos(gamma+phi)+  ...
mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))+l_tail*mQ*cos(gamma+phi)*(l_tail*mQ*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)*sin(gamma+phi)-(ICzz+mQ*l_tail^2)*(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+  ...
l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))))*(Ffoam_bottom*sin(gamma)+Ffoam_top*sin(gamma)+F_hardstop*sin(gamma)+Fy_rebound-mA*(Awidth+l_foam)*sin(gamma)*gammap^2-mC*(Awidth+l_foam)*sin(gamma)*gammap^2-l_tail*mQ*sin(gamma+phi)*(  ...
gammap+phip)^2)-(l_tail^2*mQ^2*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)*cos(gamma+phi)^2-2*l_tail*mQ*(ICzz+mQ*l_tail^2)*cos(gamma+phi)*(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*  ...
cos(gamma))-(ICzz+mQ*l_tail^2)*((mA+mC+mQ)*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))^2))*(g*mA+g*mC+g*mQ+Ffoam_bottom*cos(gamma)+Ffoam_top*cos(  ...
gamma)+F_hardstop*cos(gamma)+Fx_rebound-mA*(Awidth+l_foam)*cos(gamma)*gammap^2-mC*(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*mQ*cos(gamma+phi)*(gammap+phip)^2))/(mC^2*(Awidth+l_foam)^2*(mA+mC+mQ)*(ICzz+mQ*l_tail^2)+mA*(Awidth+l_foam)^2*(  ...
mA+2*mC)*(mA+mC+mQ)*(ICzz+mQ*l_tail^2)+l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(ICzz+mQ*l_tail^2)*(mA+mC+mQ)^2*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC)^2*  ...
sin(phi)^2);
ypp = -((Awidth+l_foam)*(mA+mC)*(l_tail^2*mQ^2*sin(phi)*sin(gamma+phi)-(mA+mC+mQ)*(ICzz+mQ*l_tail^2)*cos(gamma))*(d_foam*Ffoam_top-d_foam*Ffoam_bottom-g*l_tail*mQ*sin(gamma+phi)-g*mA*(Awidth+l_foam)*sin(gamma)-g*mC*(Awidth+l_foam)*sin(  ...
gamma))+(l_tail*mQ*(ICzz+mQ*l_tail^2)*cos(gamma+phi)*(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))-l_tail^2*mQ^2*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)*sin(gamma+phi)*  ...
cos(gamma+phi)-(Awidth+l_foam)*(mA+mC)*(ICzz+mQ*l_tail^2)*sin(gamma)*(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma)))*(g*mA+g*mC+g*mQ+Ffoam_bottom*cos(gamma)+Ffoam_top*cos(gamma)+F_hardstop*cos(  ...
gamma)+Fx_rebound-mA*(Awidth+l_foam)*cos(gamma)*gammap^2-mC*(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*mQ*cos(gamma+phi)*(gammap+phip)^2)-(T_tail+g*l_tail*mQ*sin(gamma+phi))*((mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(l_tail*mQ*cos(gamma+phi)+mA*(  ...
Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))-l_tail*mQ*(mA+mC+mQ)*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)*cos(gamma+phi)-l_tail*mQ*(Awidth+l_foam)*(mA+mC)*sin(phi)*(l_tail*mQ*sin(gamma+phi)+mA*(  ...
Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma)))-(l_tail^2*mQ^2*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)*sin(gamma+phi)^2-2*l_tail*mQ*(ICzz+mQ*l_tail^2)*sin(gamma+phi)*(l_tail*mQ*sin(gamma+phi)+mA*(  ...
Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))-(ICzz+mQ*l_tail^2)*((mA+mC+mQ)*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))^2))*(  ...
Ffoam_bottom*sin(gamma)+Ffoam_top*sin(gamma)+F_hardstop*sin(gamma)+Fy_rebound-mA*(Awidth+l_foam)*sin(gamma)*gammap^2-mC*(Awidth+l_foam)*sin(gamma)*gammap^2-l_tail*mQ*sin(gamma+phi)*(gammap+phip)^2))/(mC^2*(Awidth+l_foam)^2*(mA+mC+mQ)*(  ...
ICzz+mQ*l_tail^2)+mA*(Awidth+l_foam)^2*(mA+2*mC)*(mA+mC+mQ)*(ICzz+mQ*l_tail^2)+l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(ICzz+mQ*l_tail^2)*(mA+mC+mQ)^2*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-  ...
l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC)^2*sin(phi)^2);
gammapp = ((mA+mC+mQ)*(T_tail+g*l_tail*mQ*sin(gamma+phi))*((mA+mC+mQ)*(ICzz+mQ*l_tail^2)-l_tail^2*mQ^2-l_tail*mQ*(Awidth+l_foam)*(mA+mC)*cos(phi))+(Awidth+l_foam)*(mA+mC)*(l_tail^2*mQ^2*cos(gamma)*sin(gamma+phi)^2-(mA+mC+mQ)*(ICzz+mQ*  ...
l_tail^2)*cos(gamma)-l_tail^2*mQ^2*sin(gamma)*sin(gamma+phi)*cos(gamma+phi))*(Ffoam_bottom*sin(gamma)+Ffoam_top*sin(gamma)+F_hardstop*sin(gamma)+Fy_rebound-mA*(Awidth+l_foam)*sin(gamma)*gammap^2-mC*(Awidth+l_foam)*sin(gamma)*gammap^2-  ...
l_tail*mQ*sin(gamma+phi)*(gammap+phip)^2)-(mA+mC+mQ)*(l_tail^2*mQ^2-(mA+mC+mQ)*(ICzz+mQ*l_tail^2))*(d_foam*Ffoam_top-d_foam*Ffoam_bottom-g*l_tail*mQ*sin(gamma+phi)-g*mA*(Awidth+l_foam)*sin(gamma)-g*mC*(Awidth+l_foam)*sin(gamma))-(  ...
Awidth+l_foam)*(mA+mC)*(l_tail^2*mQ^2*sin(gamma)*cos(gamma+phi)^2-(mA+mC+mQ)*(ICzz+mQ*l_tail^2)*sin(gamma)-l_tail^2*mQ^2*cos(gamma)*sin(gamma+phi)*cos(gamma+phi))*(g*mA+g*mC+g*mQ+Ffoam_bottom*cos(gamma)+Ffoam_top*cos(gamma)+F_hardstop*  ...
cos(gamma)+Fx_rebound-mA*(Awidth+l_foam)*cos(gamma)*gammap^2-mC*(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*mQ*cos(gamma+phi)*(gammap+phip)^2))/(mC^2*(Awidth+l_foam)^2*(mA+mC+mQ)*(ICzz+mQ*l_tail^2)+mA*(Awidth+l_foam)^2*(mA+2*mC)*(mA+mC+  ...
mQ)*(ICzz+mQ*l_tail^2)+l_tail^2*mQ^2*(mA+mC+mQ)*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(ICzz+mQ*l_tail^2)*(mA+mC+mQ)^2*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC)^2*sin(phi)^2);
phipp = (((mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))-l_tail*mQ*sin(gamma+phi)*(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(  ...
gamma))*(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))-l_tail*mQ*cos(gamma+phi)*((mA+mC+mQ)*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+  ...
l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))^2))*(Ffoam_bottom*sin(gamma)+Ffoam_top*sin(gamma)+F_hardstop*sin(gamma)+Fy_rebound-mA*(Awidth+l_foam)*sin(gamma)*gammap^2-mC*(Awidth+l_foam)*sin(gamma)*gammap^2-l_tail*mQ*sin(gamma+phi)*(  ...
gammap+phip)^2)-(mA+mC+mQ)*((mA+mC+mQ)*(ICzz+mQ*l_tail^2)-l_tail^2*mQ^2-l_tail*mQ*(Awidth+l_foam)*(mA+mC)*cos(phi))*(d_foam*Ffoam_top-d_foam*Ffoam_bottom-g*l_tail*mQ*sin(gamma+phi)-g*mA*(Awidth+l_foam)*sin(gamma)-g*mC*(Awidth+l_foam)*sin(  ...
gamma))-(mA+mC+mQ)*(T_tail+g*l_tail*mQ*sin(gamma+phi))*((mA+mC+mQ)*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))^2-(l_tail*mQ*cos(  ...
gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))^2)-((mA+mC+mQ)*(ICzz+mQ*l_tail^2)*(l_tail*mQ*sin(gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))-l_tail*mQ*cos(gamma+phi)*(l_tail*mQ*sin(  ...
gamma+phi)+mA*(Awidth+l_foam)*sin(gamma)+mC*(Awidth+l_foam)*sin(gamma))*(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))-l_tail*mQ*sin(gamma+phi)*((mA+mC+mQ)*(IAzz+ICzz+mQ*l_tail^2+mA*(Awidth+  ...
l_foam)^2+mC*(Awidth+l_foam)^2)-(l_tail*mQ*cos(gamma+phi)+mA*(Awidth+l_foam)*cos(gamma)+mC*(Awidth+l_foam)*cos(gamma))^2))*(g*mA+g*mC+g*mQ+Ffoam_bottom*cos(gamma)+Ffoam_top*cos(gamma)+F_hardstop*cos(gamma)+Fx_rebound-mA*(Awidth+l_foam)*  ...
cos(gamma)*gammap^2-mC*(Awidth+l_foam)*cos(gamma)*gammap^2-l_tail*mQ*cos(gamma+phi)*(gammap+phip)^2))/(mC^2*(Awidth+l_foam)^2*(mA+mC+mQ)*(ICzz+mQ*l_tail^2)+mA*(Awidth+l_foam)^2*(mA+2*mC)*(mA+mC+mQ)*(ICzz+mQ*l_tail^2)+l_tail^2*mQ^2*(mA+  ...
mC+mQ)*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-(ICzz+mQ*l_tail^2)*(mA+mC+mQ)^2*(IAzz+mA*(Awidth+l_foam)^2+mC*(Awidth+l_foam)^2)-l_tail^2*mQ^2*(Awidth+l_foam)^2*(mA+mC)^2*sin(phi)^2);

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
x_tail = x - l_tail*cos(gamma+phi);
y_tail = y - l_tail*sin(gamma+phi);
x_Acm = x - (Awidth+l_foam)*cos(gamma);
y_Acm = y - (Awidth+l_foam)*sin(gamma);
x_FoamTop = x - d_foam*sin(gamma);
y_FoamTop = y + d_foam*cos(gamma);
x_FoamBottom = x + d_foam*sin(gamma);
y_FoamBottom = y - d_foam*cos(gamma);
x_hardstop = x - l_foam*cos(gamma);
AttachPt_x_WorldFrame = AttachPt_x + x;
AttachPt_y_WorldFrame = AttachPt_y + y;
Fmag_rebound = sqrt(Fx_rebound^2+Fy_rebound^2);
Fcontact = Ffoam_bottom + Ffoam_top + Fx_tail + F_hardstop + Fx_rebound;

Output = zeros(1,28);
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

Output(17) = t;
Output(18) = Ffoam_top;
Output(19) = Ffoam_bottom;
Output(20) = Fx_tail;
Output(21) = T_tail;
Output(22) = Fmag_rebound;
Output(23) = Fcontact;

Output(24) = t;
Output(25) = FoamContactTop;
Output(26) = FoamContactBottom;
Output(27) = TailContact;
Output(28) = FootAttached;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 4 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files BerkeleyImpact_DropTest.i  (i=1,2,3,4)\n' );
      fprintf( 1, '\n Note: Plots are automatically generated by issuing the OutputPlot command in MotionGenesis\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''BerkeleyImpact_DropTest.1'' );\n' );
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
      FileIdentifier = zeros(1,4);
      FileIdentifier(1) = fopen('BerkeleyImpact_DropTest.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file BerkeleyImpact_DropTest.1'); end
      fprintf(FileIdentifier(1), '%% FILE: BerkeleyImpact_DropTest.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              x              y            x_Acm          y_Acm         x_tail         y_tail          gamma           phi\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)           (rad)          (rad)\n\n' );
      FileIdentifier(2) = fopen('BerkeleyImpact_DropTest.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file BerkeleyImpact_DropTest.2'); end
      fprintf(FileIdentifier(2), '%% FILE: BerkeleyImpact_DropTest.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%   x_FoamTop      y_FoamTop    x_FoamBottom   y_FoamBottom   AttachPt_x_WorldFrame AttachPt_y_WorldFrame  x_hardstop\n' );
      fprintf(FileIdentifier(2), '%%      (m)            (m)            (m)            (m)                (m)                   (m)               (m)\n\n' );
      FileIdentifier(3) = fopen('BerkeleyImpact_DropTest.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file BerkeleyImpact_DropTest.3'); end
      fprintf(FileIdentifier(3), '%% FILE: BerkeleyImpact_DropTest.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t          Ffoam_top    Ffoam_bottom      Fx_tail        T_tail      Fmag_rebound     Fcontact\n' );
      fprintf(FileIdentifier(3), '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
      FileIdentifier(4) = fopen('BerkeleyImpact_DropTest.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file BerkeleyImpact_DropTest.4'); end
      fprintf(FileIdentifier(4), '%% FILE: BerkeleyImpact_DropTest.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t       FoamContactTop  FoamContactBottom  TailContact   FootAttached\n' );
      fprintf(FileIdentifier(4), '%%   (second)        (UNITS)          (UNITS)         (UNITS)        (UNITS)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(10:16) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(17:23) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(24:28) );  end
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
% Step 1: Uncomment call to mdlDerivatives.
% Step 2: Change functionsToEvaluateForEvent,                      e.g., change  []  to  [t - 5.67]  to stop at t = 5.67.
% Step 3: Change eventTerminatesIntegration1Otherwise0ToContinue,  e.g., change  []  to  [1]  to stop integrating.
% Step 4: Change eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1,  e.g., change  []  to  [1].
%----------------------------------------------------------------------
% mdlDerivatives( t, VAR, uSimulink );   % UNCOMMENT FOR EVENT HANDLING
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
function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, integStp, absError, relError, VAR, printIntScreen, printIntFile )
%===========================================================================
OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',integStp, 'Events',@EventDetection );
t = tInitial;                 epsilonT = 0.001*integStp;                tFinalMinusEpsilonT = tFinal - epsilonT;
printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + integStp;
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
      tAtEndOfIntegrationStep = t + integStp;
   else
      t = timeEventOccurredInIntegrationStep;
      VAR = nStatesArraysAtEvent;
      printCounterScreen = 0;
      printCounterFile   = 0;
      [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents );
   end
end
end


%=======================================================
end   % End of embedded function BerkeleyImpact_DropTest
%=======================================================
