function animateImpact(MGworkspace,fileloc)

% Load in variables
load(MGworkspace); 

% Set up figure 
figure
axis_bounds = [-.5, .7, 3.5, 4.5];
axis(axis_bounds)
A_width = .03; % TODO: change this from being hardcoded
sq_aspect = 1.5; % taller than wide

% Save file if a filename is specified
if nargin>1
    writerObj = VideoWriter([fileloc '.avi'])
    open(writerObj);
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
end

% Plot a frame at each timestep
for k= 1:length(t); 
    hold on 
    plot([0;0],[-10,10],'k','LineWidth',3)          % Wall
    
    % Plot square for rigid body
    sqdx = [A_width, A_width, -A_width, -A_width]; 
    sqdy = [-A_width*sq_aspect, A_width*3*sq_aspect,A_width*3*sq_aspect,-A_width*sq_aspect];
    sqx = [x_Acm(k); x_Acm(k); x_Acm(k); x_Acm(k)];
    sqy = [y_Acm(k); y_Acm(k); y_Acm(k); y_Acm(k)];
    rot = [cos(gam(k)) -sin(gam(k)) 0; sin(gam(k)) cos(gam(k)) 0; 0 0 1];
    sqdxy = rot*[sqdx; sqdy; ones(1,4)];
    sqx = sqx' +sqdxy(1,:); sqy = sqy' +sqdxy(2,:);
    fill(sqx,sqy,'r'); 
    
    plot(x_Acm(k),y_Acm(k),'r*');                           % COM Position
    plot([x_Acm(k)+0.25*sin(gam(k)); x_Acm(k)-0.25*sin(gam(k))],[y_Acm(k)-0.25*cos(gam(k));y_Acm(k)+0.25*cos(gam(k))], 'r--')
    plot(x_top(k),y_top(k),'bo');                   % Top Contact
    plot(x_bottom(k),y_bottom(k),'bo');             % Bottom Contact
    plot(x_tail(k),y_tail(k),'go');                 % Tail Position 
    plot([x_Ccm(k);x_tail(k)],[y_Ccm(k);y_tail(k)],'g');
    %plot(x_FoamTop(k), y_FoamTop(k), 'b*');               % Foam
    %plot(x_FoamBottom(k), y_FoamBottom(k), 'b*');
    plot(x(k), y(k), 'k*');                               % Foot attachment
    plot([x(k) AttachPt_x_WorldFrame(k)],[y(k) AttachPt_y_WorldFrame(k)], 'p-','LineWidth',3) % Rebound spring
    
    % Visualize forces
    scaleForce = 1/20; 
    scaleTorque = 50; 
    
    % Tail force
    plot([x_tail(k) (x_tail(k)-Fx_tail(k)*scaleForce)],[y_tail(k) (y_tail(k)-Fy_tail(k)*scaleForce)],'c', 'LineWidth',3) 
    % Lower hardstop
    plot([x_bottom(k) (x_bottom(k)-F_hardstopBottom(k)*scaleForce)],[y_bottom(k) (y_bottom(k)-Fy_fricBottom(k)*scaleForce)],'c', 'LineWidth',3) 
    % Lower hardstop
    plot([x_top(k) (x_top(k)-F_hardstopTop(k)*scaleForce)],[y_top(k) (y_top(k)-Fy_fricTop(k)*scaleForce)],'c', 'LineWidth',3) 
    % Tail Torque
    plot(x_Acm(k), y_Acm(k), 'oc', 'MarkerSize', abs(T_tail(k)*scaleTorque)+.000001) ;
    
    % Denote Time
    dim = [0.2 0.5 0.3 0.3];
    str = {'Time:',num2str(t(k))};
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'Fontsize',14);
    
    
    hold off
    axis(axis_bounds)

    %axis equal
        

    if nargin>1
        frame = getframe;
        writeVideo(writerObj,frame);
    end

    pause(.05)
    clf
end

if nargin>1
    close(writerObj)
end


