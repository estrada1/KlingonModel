function animateImpact(MGworkspace,filename)

% Load in variables
load(MGworkspace); 

% Set up figure 
figure
axis_bounds = [-.5, .2, 3.5, 4.5];
axis(axis_bounds)
A_width = .03; % TODO: change this from being hardcoded
sq_aspect = 1.5; % taller than wide

% Save file if a filename is specified
if nargin>1
    writerObj = VideoWriter([filename '.avi']);
    open(writerObj);
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
end

% Plot a frame at each timestep
for k= 1:length(t); 
    axis(axis_bounds)
    hold on 
    plot([0;0],[-10,10],'k','LineWidth',3)          % Wall
    
    % Plot square for rigid bodt
    sqdx = [A_width, A_width, -A_width, -A_width]; 
    sqdy = [-A_width*sq_aspect, A_width*sq_aspect,A_width*sq_aspect,-A_width*sq_aspect];
    sqx = [x_Acm(k); x_Acm(k); x_Acm(k); x_Acm(k)];
    sqy = [y_Acm(k); y_Acm(k); y_Acm(k); y_Acm(k)];
    rot = [cos(gam(k)) -sin(gam(k)) 0; sin(gam(k)) cos(gam(k)) 0; 0 0 1];
    sqdxy = rot*[sqdx; sqdy; ones(1,4)];
    sqx = sqx' +sqdxy(1,:); sqy = sqy' +sqdxy(2,:);
    fill(sqx,sqy,'r'); 
    
    plot(x_Acm(k),y_Acm(k),'r*');                           % COM Position
    plot([x_Acm(k)+0.25*sin(gam(k)); x_Acm(k)-0.25*sin(gam(k))],[y_Acm(k)-0.25*cos(gam(k));y_Acm(k)+0.25*cos(gam(k))], 'r--')
    %plot(x_arm(k),y_arm(k),'bo');                   % Arm Position 
    %plot([x(k);x_arm(k)],[y(k);y_arm(k)],'b');
    plot(x_tail(k),y_tail(k),'go');                 % Tail Position 
    plot([x_Acm(k);x_tail(k)],[y_Acm(k);y_tail(k)],'g');
    plot(x_FoamTop(k), y_FoamTop(k), 'b*');               % Foam
    plot(x_FoamBottom(k), y_FoamBottom(k), 'b*');
    plot(x(k), y(k), 'k*');                               % Foot attachment
    plot([x_Acm(k) AttachPt_x_WorldFrame(k)],[y_Acm(k) AttachPt_y_WorldFrame(k)], 'p-','LineWidth',3) % Rebound spring
    hold off
    axis equal
    
    if nargin>1
        frame = getframe;
        writeVideo(writerObj,frame);
    end

    pause(.05)
    clf
end

if nargin>1
    close(writerObj);
end


