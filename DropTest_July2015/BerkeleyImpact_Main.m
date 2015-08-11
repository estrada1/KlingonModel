clear all; close all; clc; 
file = 'BerkeleyImpact_DropTest';
eval(file);
[t,x,y,x_Acm,y_Acm,x_tail,y_tail,gamma,phi] = BerkeleyImpact_import1_kinematics([file '.1']);
[x_FoamTop,y_FoamTop,x_FoamBottom,y_FoamBottom,AttachPt_x_WorldFrame,AttachPt_y_WorldFrame] = BerkeleyImpact_import2_kinematics([file '.2']);
[t,Ffoam_top,Ffoam_bottom,Fx_tail,T_tail,Fmag_rebound] = BerkeleyImpact_import3_forces([file '.3']);
[t,FoamContactTop,FoamContactBottom,TailContact,FootAttached] = BerkeleyImpact_import4_contact([file '.4']);

figure
plot(x_Acm,y_Acm, x_tail, y_tail)


%% Quick-and-Dirty animation for impact

figure
axis_bounds = [-.5, .2, 3.5, 4.5];
axis(axis_bounds)

% writerObj = VideoWriter('Berkeley.avi');
% open(writerObj);
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
for k= 1:length(t); 
    axis(axis_bounds)
    hold on 
    plot([0;0],[-10,10],'k','LineWidth',3)          % Wall
    plot(x_Acm(k),y_Acm(k),'r*');                           % COM Position
    plot([x_Acm(k)+0.25*sin(gamma(k)); x_Acm(k)-0.25*sin(gamma(k))],[y_Acm(k)-0.25*cos(gamma(k));y_Acm(k)+0.25*cos(gamma(k))], 'r--')
    %plot(x_arm(k),y_arm(k),'bo');                   % Arm Position 
    %plot([x(k);x_arm(k)],[y(k);y_arm(k)],'b');
    %plot(x_tail(k),y_tail(k),'go');                 % Tail Position 
    plot([x_Acm(k);x_tail(k)],[y_Acm(k);y_tail(k)],'g');
    plot(x_FoamTop(k), y_FoamTop(k), 'b*');               % Foam
    plot(x_FoamBottom(k), y_FoamBottom(k), 'b*');
    plot(x(k), y(k), 'k*');                               % Foot attachment
    hold off
    % frame = getframe;
    % writeVideo(writerObj,frame);
    pause(.05)
    clf
end
%close(writerObj);

% 
%% Plot forces nicely, to export 
% figure
% subplot(3,1,1)
% plot(t,Fx_foot, t,Fx_arm, t,Fx_tail)
% title('Compressive Forces')
% ylabel('Force [N]')
% legend('Fx foot','Fx arm','Fx tail')
% 
% subplot(3,1,2)
% plot( t,Fmag_rebound)
% title('Adhesion Required')
% ylabel('Force [N]')
% legend('Fmag rebound')
% 
% subplot(3,1,3)
% plot( t,T_arm, t,T_tail)
% title('Internal Torques')
% ylabel('Torque [Nm]')
% legend('T arm','T tail')
% xlabel('time [s]')

figure
plot(t,FoamContactTop,t,FoamContactBottom)
legend('FoamContactTop','FoamContactBottom')
