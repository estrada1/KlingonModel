% Importing all variables output from Motion Genesis 
% Saves a workspace in the directory
% Then, just open the .mat workspace in your script and have all the
% variables to work with
% Matt Estrada
% August 7 2015

function MGworkspace = importMG(filename,exportLocation)

[t,x,y,x_Acm,y_Acm,x_tail,y_tail,gam,phi] = BerkeleyImpact_import1_kinematics([filename '.1']);
[x_FoamTop,y_FoamTop,x_FoamBottom,y_FoamBottom,AttachPt_x_WorldFrame,AttachPt_y_WorldFrame, x_hardstop, x_Ccm, y_Ccm, x_top, y_top, x_bottom, y_bottom] = BerkeleyImpact_import2_kinematics([filename '.2']);
[t,Ffoam_top,Ffoam_bottom,Fx_tail, Fy_tail, T_tail,Fx_rebound, Fy_rebound, Fx_contact, Fy_contact, Fy_fricTop, Fy_fricBottom,F_hardstopTop, F_hardstopBottom] = BerkeleyImpact_import3_forces([filename '.3']);
[t,FoamContactTop,FoamContactBottom,TailContact,FootAttached] = BerkeleyImpact_import4_contact([filename '.4']);
[KineticEnergy,KineticEnergy_body,KineticEnergy_tail,E_rebound,E_tail, GravityPotentialEnergy] = BerkeleyImpact_import5_energy([filename '.5']);
[vx,vy,vgamma,vphi,ax,ay,agamma,aphi] = BerkeleyImpact_import6_vel_acc([filename '.6']);
[Lx_body,Ly_body,H_body,Lx_tail,Ly_tail,H_tail] = BerkeleyImpact_import7_Momentum([filename '.7']); 

save(exportLocation);
MGworkspace = exportLocation; 