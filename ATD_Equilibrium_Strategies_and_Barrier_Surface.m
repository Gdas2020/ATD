%% All co-ordinates are in moving frame of reference
clear all
close all

% Target Co-ordinate in moving frame
xT0 = 0.0;     
xT = 1.0;

xD0 = xT0+0.6;  % Defender initial position  

vD = 1.0;

vA = 0.7;   % Intruder Speed in global frame
vT = 0.2;   % Target Speed

dx = 0.1;
xD = [linspace(xT0,xD0,ceil((xD0-xT0)/dx)+1),linspace(xD0,xT,ceil((xT-xD0)/dx)+1)];

xAsave1 = [];
xAsave2 = [];
xAsave3 = [];
xAsave4 = [];

xA = [];
dist_A = [];

figure(); hold on; grid on;
axis equal;

t1 = 0.1;
t2 = 0.1;

%% Linear Region

for i = 1:length(xD) 
    % Defender strategy
    if xD(1,i) > xD0
        lambda = -1;
    elseif xD(1,i) < xD0
        lambda = 1;
    else
        lambda = 0;
    end
    
    % optimal Attacker strategy
    % when y<0 
    eta = sqrt((1-lambda*vT)^2-vA^2)/vA;

    cosA = -lambda/sqrt(eta^2+1); 
    sinA = eta/sqrt(eta^2+1);
    
    if lambda ==-1
        Aopt = atan(sinA/cosA);
    elseif lambda == 1 
        Bopt = atan(sinA/cosA);
    else
        lambda = 0;
    end
    
    % Attacker speen in moving frame
    vAx = vA*cosA - vT;
    vAy = vA*sinA;
    
    time = abs(xD(1,i)-xD0)/vD;

    dist_A(1,1) = xD(1,i)-vAx*time;
    dist_A(2,1) = -vAy*time;
    
    xA = [dist_A(1,1); dist_A(2,1)];
    xAsave1 = [xAsave1, xA];
    
    plot([xD(1,i) xA(1,1)] , [0 xA(2,1)],'k-')  
    drawnow();
    pause(t1);

    % when y>0 
    eta = -sqrt((1-lambda*vT)^2-vA^2)/vA;

    cosA = -lambda/sqrt(eta^2+1);   % A* = Optimal Attacker Strategy
    sinA = eta/sqrt(eta^2+1);
    
    vAx = vA*cosA -vT;
    vAy = vA*sinA;
    
    time = abs(xD(1,i)-xD0)/vD;    

    dist_A(1,1) = xD(1,i)-vAx*time;
    dist_A(2,1) = -vAy*time;
    
    xA = [dist_A(1,1); dist_A(2,1)];
    xAsave2 = [xAsave2, xA];
    
    plot([xD(1,i) xA(1,1)] , [0 xA(2,1)],'k-')
    drawnow();
    pause(t1);

end
%% For Intruder ahead of Defender

dA = 0.3;
A = [linspace(Aopt, pi, ceil((pi-Aopt)/dA)+1), -linspace(pi, Aopt, ceil((pi-Aopt)/dA)+1)];
for i = 1:length(A)
        % when y<0 
        
    cosA = cos(A(i));
    sinA = sin(A(i));

    vAx = vA*cosA - vT;
    vAy = vA*sinA;
    
    time = abs(xT-xD0);

    dist_A(1,1) = -vAx*time + time;
    dist_A(2,1) = -vAy*time;
    
    xA = [dist_A(1,1); dist_A(2,1)];
    xAsave3 = [xAsave3, [xA(1,1)+xD0; xA(2,1)]];
    
    plot([time+xD0 xA(1,1)+xD0] , [0 xA(2,1)],'b-')
    drawnow();
    pause(t1);

end

%% For Defender ahead of Intruder

dB = 0.3;
B = [linspace(pi+Bopt,0, ceil(abs(pi+Bopt))/dB+1), -linspace(0, pi+Bopt, ceil(abs(pi+Bopt))/dB+1)];

count = 0;
for i = 1:length(B)
        % when y<0 
        
    cosA = cos(B(i));
    sinA = sin(B(i));

    vAx = vA*cosA - vT;
    vAy = vA*sinA;
    
    time = abs(xD0-xT0);

    dist_A(1,1) = -vAx*time;
    dist_A(2,1) = -vAy*time;
    
    if abs(dist_A(1,1)) > xD0
        dist_A(1,1) = xD0;
    end
    

    xA = [dist_A(1,1); dist_A(2,1)];
    xAsave4 = [xAsave4, xA];
    
    plot([0 xA(1,1)] , [time-xD0 xA(2,1)],'b-')
    drawnow();
    pause(t1);
end

%% Plots

plot(xAsave1(1,:),xAsave1(2,:),'r-', 'linewidth', 1.5)
    drawnow();
    pause(t2);

plot(xAsave2(1,:),xAsave2(2,:),'r-', 'linewidth', 1.5)
    drawnow();
    pause(t2);
  
plot(xAsave3(1,:),xAsave3(2,:),'r-', 'linewidth', 1.5)
    drawnow();
    pause(t2);
   
plot(xAsave4(1,:),xAsave4(2,:),'r-', 'linewidth', 1.5)
    drawnow();
    pause(t2);

title('Trajectories for optimal strategy')
xlabel('x')
ylabel('y')
title("Optimal trajectories with v_{A} = " + vA + " and v_{T} = " + vT)

% adding axis lines
axh = gca; % use current axes
line(get(axh,'XLim'), [0 0], 'Color', 'k', 'LineStyle', '-','linewidth',1);
line([0 0], get(axh,'YLim'), 'Color', 'k', 'LineStyle', '-','linewidth',1);