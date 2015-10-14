%% Simlulation script
% Short commands for all the big chunks I'm working with

filename = 'BerkeleyImpact_TallBody';
batch = 'Oct 8 Skype'; 
trial = 'No tail';
time = char(datetime); 
mkdir(batch,time)

%% Run simulation 
eval(filename);

%%
exportdirectory = [batch '/' time]; % Save video/images for the trial
exportpath = [exportdirectory '/' trial];
MGworkspace = importMG(filename, exportpath); 

%% Animate
animateImpact(MGworkspace, exportpath); 

%% Performance Spec's 

load(MGworkspace); 

ReboundMag = sqrt(Fx_rebound.^2 + Fy_rebound.^2);
MaxCompression = min(Fx_contact);
tMaxCompression = t(find(Fx_contact == MaxCompression,1));
MaxRebound = max(ReboundMag);
tMaxRebound = t(find(ReboundMag == MaxRebound,1));

% Plot interesting values 


figure(2)
plot(t,Fx_contact, t, Fy_contact, t, Fx_rebound, t, Fy_rebound, t, ReboundMag, t, T_tail)
legend('Fx Normal', 'Fy Tangential', 'Fx Rebound', 'Fy Rebound', 'ReboundMag', 'Ttail')
xlabel('time [sec]')
ylabel('Force [N]')
title([trial ' Contact'])
savefig([exportpath ' Forces'])
saveas(gcf,[exportpath ' Forces'],'png')


figure(3)
plot(t,KineticEnergy, t, KineticEnergy_body, t, KineticEnergy_tail)
title([trial ' Energy'])
legend('Kinetic Energy', 'Body K', 'Tail K')
ylabel('Energy [J]')
xlabel('time [s]')

figure(4)
L = sqrt((Lx_body + Lx_tail).^2 + (Ly_body + Ly_tail).^2); 
H = H_body + H_tail; 
figure(4) 
plot(t, L, t, H)
legend('Linear Momentum', 'Angular Momentum') 

figure(5) 
GravOffset = 8.4; 
total = KineticEnergy +  E_rebound + E_tail  +GravityPotentialEnergy- GravOffset; 
plot(t,KineticEnergy, t, E_rebound, t, E_tail, t, GravityPotentialEnergy - GravOffset, t, total)
h_legend = legend('Kinetic Energy', 'Rebound Potential Energy', 'Tail Potential Energy', 'Gravity Potential Energy', 'total');
set(h_legend,'FontSize',14);
ylabel('Energy [J]')
xlabel('time [s]')
title([trial ' Energy'])
savefig([exportpath ' Energy'])
saveas(gcf,[exportpath ' Energy'],'png')