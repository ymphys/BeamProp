function beam_path(R,psi,gamma,B)
% Plots the path of an electron beam in the presence of a magnetic field
% Sam Collopy, Sem-sol, February 2021

% fundamental constants
c = 2.998e8; % speed of light [m/s]
qe = 1.609e-19; % fundamental charge [C]
me = 9.11e-31; %electron mass [kg]

v0 = c*sqrt(1-1/gamma^2); % beam velocity
wB = qe*B/(gamma*me); % cyclontron frequency

% Find limits to find first zero of test function
theta_test = linspace(0,pi/2, 100000);
fn_test = (R*sin(psi))^2 - ((v0*sin(theta_test)/wB).^2).*(2-2*cos(wB*R*cos(psi)./(v0*cos(theta_test))));
% find first sign flip
index_test = find(fn_test < 0, 1, 'first');

r0 = R*sin(psi)
r0max = 2*v0/wB

if(r0 > r0max)
    disp('Error: Target out of range.')
end

if(1)

    % Plot optimization function as a function of theta
    figure()
    plot(theta_test,fn_test)
    grid on
end

% Find solution for theta
rootfn = @(theta) rootfn(R,psi,v0,wB,theta);
theta = fzero(rootfn,[theta_test(index_test-1), theta_test(index_test)]);

tmax = R*cos(psi)/(v0*cos(theta)); % time at target impact
a0 = abs(v0*sin(theta)/wB) % 2D cyclotron radius

phi = atan2((1-cos(wB*tmax)),sin(wB*tmax));
%r0 = R*sin(psi)

t = linspace(0,tmax);

x = a0*(cos(wB*t)-1);
y = a0*sin(wB*t);
z = v0*cos(theta)*t;

xp = cos(phi)*x + sin(phi)*y;
yp = cos(phi)*y - sin(phi)*x;

xT = 0;
yT = R*sin(psi);
zT = R*cos(psi);

figure()

plot3(xp,yp,z);
grid on;

xlabel('x')
ylabel('y')
zlabel('z')
hold on
plot3([0,xT],[0,yT],[0,zT])
axis equal
%plot3(xT,yT,zT)
% >> xlim([0 20])
% >> ylim([0 200])
% >> zlim([0 1000])
% >> set(gca, 'DataAspectRatio', [1 10 50])
% >> set(gca, 'fontsize', 20)

% Translate angles to angles from target
%thetaT = theta - psi;
%phiT = acos(sin(psi)^2*cos(phi)+cos(psi)^2);

xp = sin(theta)*cos(phi)*cos(psi) - cos(theta)*sin(psi);
yp = -sin(theta)*sin(phi);
zp = sin(theta)*cos(phi)*sin(psi) + cos(theta)*cos(psi);

thetaT = atan(xp/zp);
phiT = atan(yp/zp);

disp(['Targeting angles are', num2str(thetaT*180/pi), ' deg. polar, ', num2str(phiT*180/pi), ' deg. elevation.' ])
end

function fn_val = root_fn(R,psi,v0,wB,theta)

    %fn_val = R*sin(psi) -(v0*sin(theta)/wB)*sqrt(2-2*cos(wB*R*cos(psi)/(v0*cos(theta))));
    fn_val = (R*sin(psi))^2 -(v0*sin(theta)/wB)^2*(2-2*cos(wB*R*cos(psi)/(v0*cos(theta))));
endfunction