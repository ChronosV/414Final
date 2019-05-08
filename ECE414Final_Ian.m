%Get a TF
s = tf('s');

%Create constants
Jg = 6.2e-6;    %Gearbox inertia
Js = 1.4e-7;    %Angular sensor inertia
Gv = 5;         %Voltage amplifier gain
Ks = 10;    %touch sensor, in V/m
Jc = 5.5271e-6;     %inertia of channel
Rb = 1e-2;      %Radius of ball, in meters
Rh = 6e-3;      %Height difference between center of ball and top of rail, in meters
Mb = 1.1305e-2;     %Mass of ball
Jb = 0.4*Mb*(Rb)^2; %Ball inertia
g = 9.81;       %Acceleration of gravity
N = 10;         %Chosen Gear Ratio
Jx = Jg + (Js/(N^2));   %Constant inertia, neglecting motor, w/o ball and beam inertia
% Jx = Jg +(Js + Jb)/(N^2);   %Constant inertia, neglecting motor, w/o beam inertia
% Jx = Jg + (Js + Jc + Jb)/(N^2);   %Constant inertia, neglecting motor


%% Find the Nominal plants
% Need to do this for all for motors


% for n = 1:4;
%     motor = motornumber(n);
%     motor.Jeff = motor.Jm + Jx;
%     
%     Kv = 60/(2*pi*motor.Kt(3));
%     Ke = motor.Kt(3);
%     numerator = (-Gv.*motor.Kt(3).*Kv*g)/(motor.L.*motor.Jeff*N*Mb.*(s^5));
%     delta = 1 - ((-motor.R/(motor.L.*s)) - (motor.Bm/(motor.Jeff.*s)) - ((Ke.*motor.Kt(3))/(motor.L.*motor.Jeff.*(s^2)))) + ((-motor.R/(motor.L.*s)).*(-motor.Bm/(motor.Jeff.*s)));
%     
%     motor.G = numerator/delta;
%     
%     Plant(n) = motor;
%     
%     figure(n);
%     rlocus(Plant(n).G);
% end

motor = motornumber(4);
motor.Jeff = motor.Jm + Jx;

Kv = 60/(2*pi*motor.Kt(3));
Ke = motor.Kt(3);
numerator = (-Gv.*motor.Kt(3).*Kv*g)/(motor.L.*motor.Jeff*N*Mb.*(s^5));
delta = 1 - ((-motor.R/(motor.L.*s)) - (motor.Bm/(motor.Jeff.*s)) - ((Ke.*motor.Kt(3))/(motor.L.*motor.Jeff.*(s^2)))) + ((-motor.R/(motor.L.*s)).*(-motor.Bm/(motor.Jeff.*s)));

motor.G = numerator/delta;
Gnom = motor.G;

Gans = minreal(Gnom);
rltool(Gans);


% for m = 1:4
%     
%     motor = motornumber(m); % Choose one of the four available motors
%     Inertia calculations
%     motor.Jsys = motor.Jm + d.Jx;   % Calculate total system inertia
%     
%     whole plant
%     Numerator = (d.Gv .* motor.Kt * d.Ks * d.Rp)./(motor.L .* motor.Jsys);
%     sOneTerm = (motor.R/motor.L)+(motor.Bm./motor.Jsys);
%     sZeroTerm = (((motor.Kt).^2)+ motor.R * motor.Bm)./(motor.L .* motor.Jsys);
%     
%     Velocity
%     Znum =(motor.Kt)./(motor.L.*motor.Jsys);
%     Zden =(s^2 + (motor.R/motor.L)*s+(motor.Bm./motor.Jsys).*s+((motor.Kt.^2 + motor.R*motor.Bm)./(motor.L.*motor.Jsys)));
%     
%     Current
%     Curr_plant_num = (d.Gv/motor.L)*(s+(motor.Bm./motor.Jsys));
%     Curr_plant_den = s^2+((motor.R/motor.L)+(motor.Bm./motor.Jsys))*s + ((((motor.Kt).^2)+motor.R*motor.Bm)./(motor.L.*motor.Jsys));
%     
%     Calculate open loop transfer functions
%     for j = 1:3
%         motor.G(j) = Numerator(j)/(s*(s^2 + sOneTerm(j) * s + sZeroTerm(j))); % Calculate plant transfer function
%         Z(j) = Znum(j)/Zden(j);
%         motor.Vel(j) = d.Ks*Z(j)*d.Gv*d.Rp;
%         motor.Cur(j) = Curr_plant_num(j)/Curr_plant_den(j);
%         motor.Volt(j) = d.Gv;
%     end
%     
%     Save the data so we have it for later
%     Plant(m) = motor;
%     
%     Plot Root Locus of plant transfer function
%     figure(m+10); clf;
% 	rlocus(Plant(m).G);
% end










