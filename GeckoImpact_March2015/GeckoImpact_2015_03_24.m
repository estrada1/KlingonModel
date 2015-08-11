function [t,VAR,Output] = GeckoImpact
%===========================================================================
% File: GeckoImpact.m created on Thu Mar 26 2015 by MotionGenesis 5.4.
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
gammapp=0; phipp=0; thetapp=0; xpp=0; ypp=0; Fx_arm=0; Fx_foot=0; Fx_tail=0; T_arm=0; T_tail=0; AttachPt_xp=0; AttachPt_yp=0; ArmContact=0; ArmPotentialEnergy=0; FootAttached=0; FootContact=0; FootPotentialEnergy=0; Fx_rebound=0; Fy_rebound=0;
 GravityPotentialEnergy=0; KineticEnergy=0; MechanicalEnergy=0; ReboundPotentialEnergy=0; SpringPotentialEnergy=0; TailContact=0; TailPotentialEnergy=0; WallPotentialEnergy=0; x_arm=0; x_tail=0; y_arm=0; y_tail=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.80665;                % m/sec^2             Constant
IAzz                            =  0.4;                    % kg*m^2              Constant
IBzz                            =  0.000002;               % kg*m^2              Constant
ICzz                            =  0.000002;               % kg*m^2              Constant
k_arm                           =  100;                    % N                   Constant
k_foot                          =  10000;                  % N/m                 Constant
k_rebound                       =  1000;                   % N/m                 Constant
k_tail                          =  100;                    % N                   Constant
k_wall                          =  1000;                   % N                   Constant
l_arm                           =  0.5;                    % m                   Constant
l_tail                          =  0.5;                    % m                   Constant
mA                              =  2;                      % kg                  Constant
mB                              =  .00002;                 % kg                  Constant
mC                              =  .00002;                 % kg                  Constant

gamma                           =  45;                     % deg                 Initial Value
phi                             =  0;                      % deg                 Initial Value
theta                           =  0;                      % deg                 Initial Value
x                               = -2;                      % m                   Initial Value
y                               =  4;                      % m                   Initial Value
gammap                          =  0;                      % rad/sec             Initial Value
phip                            =  0;                      % rad/sec             Initial Value
thetap                          =  0;                      % rad/sec             Initial Value
xp                              =  5;                      % m/sec               Initial Value
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

% Unit conversions
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
gamma = gamma * DEGtoRAD;
phi = phi * DEGtoRAD;
theta = theta * DEGtoRAD;

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, integStp, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files
PlotOutputFiles;


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
FootContact = ceil(0.5*sign(x));
FootAttached = ceil(0.5*sign(sqrt(AttachPt_x^2+AttachPt_y^2))+0.5*FootContact);
Fy_rebound = k_rebound*AttachPt_y*FootAttached;
ypp = -(g*mA+g*mB+g*mC-Fy_rebound)/(mA+mB+mC);
ArmContact = ceil(0.5*sign(x+l_arm*sin(gamma-theta)));
TailContact = ceil(0.5*sign(x-l_tail*sin(gamma+phi)));
Fx_rebound = k_rebound*AttachPt_x*FootAttached;


% Quantities that were specified
T_arm = -k_arm*theta;
T_tail = k_tail*phi;
Fx_foot = -k_foot*x*FootContact;
Fx_arm = -k_wall*(x+l_arm*sin(gamma-theta))*ArmContact;
Fx_tail = -k_wall*(x-l_tail*sin(gamma+phi))*TailContact;
AttachPt_xp = -FootAttached*xp;
AttachPt_yp = -FootAttached*yp;

xpp = (Fx_arm+Fx_foot+Fx_tail+Fx_rebound)/(mA+mB+mC);
gammapp = (T_arm+T_tail)/IAzz;
thetapp = (T_tail+l_arm*Fx_arm*cos(gamma-theta))/IAzz + (IAzz+IBzz)*(T_arm-l_arm*Fx_arm*cos(gamma-theta))/(IAzz*IBzz);
phipp = -(T_arm-l_tail*Fx_tail*cos(gamma+phi))/IAzz - (IAzz+ICzz)*(T_tail+l_tail*Fx_tail*cos(gamma+phi))/(IAzz*ICzz);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros(1,12);
VAR(1) = gamma;
VAR(2) = phi;
VAR(3) = theta;
VAR(4) = x;
VAR(5) = y;
VAR(6) = gammap;
VAR(7) = phip;
VAR(8) = thetap;
VAR(9) = xp;
VAR(10) = yp;
VAR(11) = AttachPt_x;
VAR(12) = AttachPt_y;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
gamma = VAR(1);
phi = VAR(2);
theta = VAR(3);
x = VAR(4);
y = VAR(5);
gammap = VAR(6);
phip = VAR(7);
thetap = VAR(8);
xp = VAR(9);
yp = VAR(10);
AttachPt_x = VAR(11);
AttachPt_y = VAR(12);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros(1,12);
VARp(1) = gammap;
VARp(2) = phip;
VARp(3) = thetap;
VARp(4) = xp;
VARp(5) = yp;
VARp(6) = gammapp;
VARp(7) = phipp;
VARp(8) = thetapp;
VARp(9) = xpp;
VARp(10) = ypp;
VARp(11) = AttachPt_xp;
VARp(12) = AttachPt_yp;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
x_arm = x + l_arm*sin(gamma-theta);
y_arm = y + l_arm*cos(gamma-theta);
x_tail = x - l_tail*sin(gamma+phi);
y_tail = y - l_tail*cos(gamma+phi);
KineticEnergy = 0.5*IAzz*gammap^2 + 0.5*ICzz*(gammap+phip)^2 + 0.5*IBzz*(gammap-thetap)^2 + 0.5*mA*(xp^2+yp^2) + 0.5*mB*(xp^2+yp^2) + 0.5*mC*(xp^2+yp^2);
GravityPotentialEnergy = g*(mA+mB+mC)*y;
WallPotentialEnergy = 0.5*k_wall*(ArmContact*x_arm^2+TailContact*x_tail^2);
ArmPotentialEnergy = 0.5*k_arm*theta^2;
FootPotentialEnergy = 0.5*k_foot*x^2*FootContact;
TailPotentialEnergy = 0.5*k_tail*gamma^2;
ReboundPotentialEnergy = 0.5*k_rebound*(AttachPt_x^2+AttachPt_y^2);
SpringPotentialEnergy = ReboundPotentialEnergy + ArmPotentialEnergy + FootPotentialEnergy + TailPotentialEnergy + WallPotentialEnergy;
MechanicalEnergy = GravityPotentialEnergy + SpringPotentialEnergy + KineticEnergy;

Output = zeros(1,36);
Output(1) = t;
Output(2) = x;
Output(3) = y;
Output(4) = x_arm;
Output(5) = y_arm;
Output(6) = x_tail;
Output(7) = y_tail;
Output(8) = gamma;
Output(9) = theta;
Output(10) = phi;

Output(11) = x;
Output(12) = y;

Output(13) = t;
Output(14) = x;
Output(15) = y;
Output(16) = x_arm;
Output(17) = y_arm;

Output(18) = t;
Output(19) = FootContact;
Output(20) = ArmContact;
Output(21) = TailContact;
Output(22) = FootAttached;

Output(23) = t;
Output(24) = theta*RADtoDEG;
Output(25) = gamma*RADtoDEG;
Output(26) = phi*RADtoDEG;

Output(27) = t;
Output(28) = Fx_foot;
Output(29) = T_arm;
Output(30) = Fx_arm;
Output(31) = Fx_tail;

Output(32) = t;
Output(33) = MechanicalEnergy;
Output(34) = KineticEnergy;
Output(35) = GravityPotentialEnergy;
Output(36) = SpringPotentialEnergy;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 7 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files GeckoImpact.i  (i=1, ..., 7)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t              x              y            x_arm          y_arm         x_tail         y_tail          gamma          theta           phi\n' );
      fprintf( 1,                '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)           (rad)          (rad)          (rad)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros(1,7);
      FileIdentifier(1) = fopen('GeckoImpact.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file GeckoImpact.1'); end
      fprintf(FileIdentifier(1), '%% FILE: GeckoImpact.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              x              y            x_arm          y_arm         x_tail         y_tail          gamma          theta           phi\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)            (m)            (m)            (m)            (m)           (rad)          (rad)          (rad)\n\n' );
      FileIdentifier(2) = fopen('GeckoImpact.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file GeckoImpact.2'); end
      fprintf(FileIdentifier(2), '%% FILE: GeckoImpact.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       x              y\n' );
      fprintf(FileIdentifier(2), '%%      (m)            (m)\n\n' );
      FileIdentifier(3) = fopen('GeckoImpact.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file GeckoImpact.3'); end
      fprintf(FileIdentifier(3), '%% FILE: GeckoImpact.3\n%%\n' );
      fprintf(FileIdentifier(3), '%%       t              x              y            x_arm          y_arm\n' );
      fprintf(FileIdentifier(3), '%%   (second)          (m)            (m)          (UNITS)        (UNITS)\n\n' );
      FileIdentifier(4) = fopen('GeckoImpact.4', 'wt');   if( FileIdentifier(4) == -1 ), error('Error: unable to open file GeckoImpact.4'); end
      fprintf(FileIdentifier(4), '%% FILE: GeckoImpact.4\n%%\n' );
      fprintf(FileIdentifier(4), '%%       t         FootContact    ArmContact     TailContact   FootAttached\n' );
      fprintf(FileIdentifier(4), '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
      FileIdentifier(5) = fopen('GeckoImpact.5', 'wt');   if( FileIdentifier(5) == -1 ), error('Error: unable to open file GeckoImpact.5'); end
      fprintf(FileIdentifier(5), '%% FILE: GeckoImpact.5\n%%\n' );
      fprintf(FileIdentifier(5), '%%       t            theta          gamma           phi\n' );
      fprintf(FileIdentifier(5), '%%   (second)         (deg)          (deg)          (deg)\n\n' );
      FileIdentifier(6) = fopen('GeckoImpact.6', 'wt');   if( FileIdentifier(6) == -1 ), error('Error: unable to open file GeckoImpact.6'); end
      fprintf(FileIdentifier(6), '%% FILE: GeckoImpact.6\n%%\n' );
      fprintf(FileIdentifier(6), '%%       t           Fx_foot         T_arm         Fx_arm         Fx_tail\n' );
      fprintf(FileIdentifier(6), '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
      FileIdentifier(7) = fopen('GeckoImpact.7', 'wt');   if( FileIdentifier(7) == -1 ), error('Error: unable to open file GeckoImpact.7'); end
      fprintf(FileIdentifier(7), '%% FILE: GeckoImpact.7\n%%\n' );
      fprintf(FileIdentifier(7), '%%       t       MechanicalEnergy  KineticEnergy GravityPotentialEnergy  SpringPotentialEnergy\n' );
      fprintf(FileIdentifier(7), '%%   (second)         (UNITS)         (UNITS)            (UNITS)                (UNITS)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:10) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:10) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(11:12) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(13:17) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(4), Output(18:22) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(5), Output(23:26) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(6), Output(27:31) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(7), Output(32:36) );  end
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
function PlotOutputFiles
%===========================================================================
if( printIntFile == 0 ),  return;  end

figure;
data = load( 'GeckoImpact.2' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'y (m)' );
xlabel('x (m)');   ylabel('y (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'GeckoImpact.3' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', 'LineWidth',3 );
legend( 'x (m)', 'y (m)', 'x_arm', 'y_arm' );
xlabel('t (second)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'GeckoImpact.4' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', 'LineWidth',3 );
legend( 'FootContact', 'ArmContact', 'TailContact', 'FootAttached' );
xlabel('t (second)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'GeckoImpact.5' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'theta (deg)', 'gamma (deg)', 'phi (deg)' );
xlabel('t (second)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'GeckoImpact.6' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', 'LineWidth',3 );
legend( 'Fx_foot', 'T_arm', 'Fx_arm', 'Fx_tail' );
xlabel('t (second)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;

figure;
data = load( 'GeckoImpact.7' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', 'LineWidth',3 );
legend( 'MechanicalEnergy', 'KineticEnergy', 'GravityPotentialEnergy', 'SpringPotentialEnergy' );
xlabel('t (second)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;
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


%===========================================
end   % End of embedded function GeckoImpact
%===========================================
