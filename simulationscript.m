%% Simlulation script
% Short commands for all the big chunks I'm working with

filename = 'BerkeleyImpact_TallBody';
trialconditions = 'LongTail'; 

%% Run simulation 
eval(filename);
%%
MGworkspace = importMG(filename); 

%% Animate
animateImpact(MGworkspace, [filename '  ' trialconditions]); 

%% Plot interesting values 

load(MGworkspace); 

figure(2)
plot(t,Fx_contact, t, Fy_contact, t, Fx_rebound, t, Fy_rebound)
legend('Fx Normal', 'Fy Tangential')
xlabel('time [sec]')
ylabel('Force [N]')
title([trialconditions ' Contact'])
%savefig([trialconditions '_Forces'])

figure(3)
plot(t,KineticEnergy, t, KineticEnergy_body, t, KineticEnergy_tail)
title([trialconditions ' Energy'])
legend('Kinetic Energy', 'Body K', 'Tail K')
ylabel('Energy [J]')
xlabel('time [s]')
%savefig([trialconditions '_Energy'])

figure(4)
L = sqrt((Lx_body + Lx_tail).^2 + (Ly_body + Ly_tail).^2); 
H = H_body + H_tail; 
figure(4) 
plot(t, L, t, H)
legend('Linear Momentum', 'Angular Momentum') 

figure(5) 
total = KineticEnergy +  E_rebound + E_tail  +GravityPotentialEnergy-min(GravityPotentialEnergy); 
plot(t,KineticEnergy, t, E_rebound, t, E_tail, t, GravityPotentialEnergy -min(GravityPotentialEnergy), t, total)
title([trialconditions ' Energy'])
h_legend = legend('Kinetic Energy', 'Rebound Potential Energy', 'Tail Potential Energy', 'Gravity Potential Energy', 'total');
set(h_legend,'FontSize',14);
ylabel('Energy [J]')
xlabel('time [s]')


