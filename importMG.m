% Importing all variables output from Motion Genesis 
% Saves a workspace in the directory
% Then, just open the .mat workspace in your script and have all the
% variables to work with
% Matt Estrada
% August 7 2015

function MGworkspace = importMG(file)

[t,x,y,x_Acm,y_Acm,x_tail,y_tail,gam,phi] = BerkeleyImpact_import1_kinematics([file '.1']);
[x_FoamTop,y_FoamTop,x_FoamBottom,y_FoamBottom,AttachPt_x_WorldFrame,AttachPt_y_WorldFrame, x_hardstop] = BerkeleyImpact_import2_kinematics([file '.2']);
[t,Ffoam_top,Ffoam_bottom,Fx_tail,T_tail,Fx_rebound, Fy_rebound, Fx_contact, Fy_contact, Fy_fricTop, Fy_fricBottom,F_hardstopTop, F_hardstopBottom] = BerkeleyImpact_import3_forces([file '.3']);
[t,FoamContactTop,FoamContactBottom,TailContact,FootAttached] = BerkeleyImpact_import4_contact([file '.4']);

save(file);
MGworkspace = file; 