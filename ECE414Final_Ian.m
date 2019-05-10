%Get a TF
s = tf('s');

%Create constants
Jg = 6.2e-6;    %Gearbox inertia
Js = 1.4e-7;    %Angular sensor inertia
Gv = 5;         %Voltage amplifier gain
Ks = 10;    %touch sensor, in V/m
%Jc = 5.5271e-6;     %inertia of channel, which can be neglected
Rb = 1e-2;      %Radius of ball, in meters
Rh = 6e-3;      %Height difference between center of ball and top of rail, in meters
% Mb = 1.1305e-2;     %Mass of ball
Mb = 2.1111;   %NOT the mass of the ball
% Jb = 0.4*Mb*(Rb)^2; %Ball inertia, which can be neglected
g = 9.81;       %Acceleration of gravity
N = 10;         %Chosen Gear Ratio
Jt = 1.8;       %Track inertia, calculate by treating it as a rod
Jx = Jg + (Jt/(N^2));   %Constant inertia, neglecting motor
% Jx = Jg +(Js + Jb)/(N^2);   %Constant inertia, neglecting motor, w/o beam inertia
% Jx = Jg + (Js + Jc + Jb)/(N^2);   %Constant inertia, neglecting motor


%% Find the Nominal plants
% Need to do this for all for motors


for n = 1:4;
    motor = motornumber(n);
    Jeff = motor.Jm + Jx;

    Kv = 10;
    Ke = motor.Kt(3);
    Ke = Kt(3);
    L = motor.L;
    R = motor.R;
    Bm = motor.Bm;

    num = (-Gv.*Ke.*Kv.*g)/(L.*Jeff.*N.*Mb.*(s^5));
    den = 1 + (((R/L) + (Bm/Jeff))/s) + ((Ke^2 + R.*Bm)/(L.*Jeff))/(s^2);
    
    G(n) = num/den;
    
    hold on;
    figure(2);
    pzmap(G(n));
    grid on;
    
    hold on;
    figure(3);
    pzmap(G(n));
    grid on;

%     figure(1);
%     rlocus(Plant(n).G);
end

figure(2);
title('Pole-Zero Map of Nominal Plants');
legend('Motor 1', 'Motor 2', 'Motor 3', 'Motor 4');

figure(3);
title('Pole-Zero Map of Nominal Plants');
axis([-0.4 0.1 -0.006 0.006]);
legend('Motor 1', 'Motor 2', 'Motor 3', 'Motor 4');

% figure(1);
% subplot(2, 2, 1);
% rlocus(G(1));
% % pzmap(G(1));
% title('Root Locus for Motor 1');
% % grid;
% 
% subplot(2, 2, 2);
% rlocus(G(2));
% % pzmap(G(2));
% title('Root Locus for Motor 2');
% % grid;
% 
% subplot(2, 2, 3);
% rlocus(G(3));
% % pzmap(G(3));
% title('Root Locus for Motor 3');
% % grid;
% 
% subplot(2, 2, 4);
% rlocus(G(4));
% % pzmap(G(4));
% title('Root Locus for Motor 4');
% % grid;

% motor = motornumber(3);
% Jeff = motor.Jm + Jx;
% 
% % Kv = 60/(2*pi*motor.Kt(3));
% Kv = 10;
% Kt = motor.Kt;
% Ke = Kt(3);
% L = motor.L;
% R = motor.R;
% Bm = motor.Bm;
% % numerator = (-Gv.*motor.Kt(3).*Kv.*g)/(motor.L.*motor.Jeff.*N.*Mb.*(s^5));
% % delta = 1 - ((-motor.R/(motor.L.*s)) - (motor.Bm/(motor.Jeff.*s)) - ((Ke.*motor.Kt(3))/(motor.L.*motor.Jeff.*(s^2)))) + ((-motor.R/(motor.L.*s)).*(-motor.Bm/(motor.Jeff.*s)));
% 
% num = (-Gv.*Ke.*Kv.*g)/(L.*Jeff.*N.*Mb.*(s^5));
% den = 1 + (((R/L) + (Bm/Jeff))/s) + ((Ke^2 + R.*Bm)/(L.*Jeff))/(s^2);
% 
% % num = (-Gv*motor.Kt(3)*Kv*g)/(motor.L*motor.Jeff*N*Mb);
% % den = [1, ((motor.R/motor.L) + (motor.Bm/motor.Jeff)), ((Ke^2 + motor.R*motor.Bm)/(motor.L*motor.Jeff)), 0, 0, 0];
% % motor.G = tf(num, den);
% % Gnom = motor.G;
% 
% % motor.G = numerator/delta;
% G = num/den;
% 
% Gans = minreal(G);
% % rltool(Gans);



