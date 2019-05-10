%% ECE 414 Project
% Nick Fortin and Quang Luong
% Date: 4/23/19


%% Given Constants
close all

% System Parts
Jg = 6.2e-6;            % gearbox inertia (kg-m^2)
Js = 1.4e-7;            % angular sensor inertia (kg-m^2)
Gv = 5;                 % voltage amplifier gain (5V/V)
K = 1;                          % pre amp gain
% K = .94;

% Characteristics that are the same for motors 1, 2, 3, and 4
Jm = 5e-5;              % motor inertia (kg-m^2)
Bm = 3e-6;              % motor viscous friction (Nm-s/rad)

N = linspace(10,50,100);        % turn ratio

Bro = 1e-2;             % ball outer radius (1 cm)

% move to calculated section
hr = 5.697e-3;          % height from top of track to ball
M = (1 + 2/5 * (Bro/hr) ^ 2);

% Motor 1 Characteristics
Kt_1 = [0.2025 0.225 0.2475];   % motor 1 torque constant
R_1 = 8;                        % motor 1 resistance
L_1 = 25e-3;                    % motor 1 inductance

% Motor 2 Characteristics
Kt_2 = [0.1575 0.175 0.1925];   % motor 2 torque constant
R_2 = 6;                        % motor 2 resistance
L_2 = 16e-3;                    % motor 2 inductance

% Motor 3 Characteristics
Kt_3 = [0.1125 0.125 0.1375];   % motor 3 torque constant
R_3 = 4;                        % motor 3 resistance
L_3 = 7.5e-3;                   % motor 3 inductance

% Motor 4 Characteristics
Kt_4 = [0.2475 0.275 0.3025];   % motor 4 torque constant
R_4 = 12;                       % motor 4 resistance
L_4 = 32e-3;                    % motor 4 inductance

Kv = 0.1*100;                   % sensor V/m constant
g = 9.81;                        % acceleration due to gravity

% Inertias
Jb_steel = 1.3404e-6;           % inertia of steel ball, kgm^2
Jb_aluminum = 4.523e-7;         % inertia of aluminum ball, kgm^2
% J_track = 3.08e-3;              % inertia of track, kgm^2
J_track = 1.8;              % moment of inertia of track
%% Calculated Values

% Effective Inertias 
% J_eff_steel = Jm + Js + (Jg + J_track + Jb_steel)./(N.^2);    
% J_eff_aluminum = Jm + Js + (Jg + J_track + Jb_aluminum)./(N.^2);
% J_eff = J_eff_steel;

J_eff_steel = Jm + Js + Jg + J_track./(N.^2);    
J_eff_aluminum = Jm + Js + Jg + J_track./(N.^2);
J_eff = J_eff_steel;


% This plot shows the inertia for both the aluminum and steel balls
figure(21)
plot(N, J_eff_steel, N, J_eff_aluminum);
grid on
xlabel('Gear Ratio, N');
ylabel('Inertia');
title('Gear Ratio vs Inertia for both types of ball');
legend('Steel Ball', 'Aluminum Ball', 'Location', 'northwest');

% This plot shows the differnce in the effective inertias
figure(2)
plot(N, J_eff_steel-J_eff_aluminum);
grid on
xlabel('Gear Ratio, N');
ylabel('Inertia');
title('Gear Ratio vs Inertia for both types of ball differnce');
legend('Difference', 'Location', 'northwest');


%% Find the variation of Kt with N constant
close all
% To understand the effect of variance on Kt, motor 1 was used as the
% effects on the other motors would have similar effects.
N = 10;     % Decleration of turn ratio at 10 as constant
J_eff_10 = Jm + Js + Jg + J_track./N^2;  % Find Jeff using N = 10
Kt_vary = linspace(0.2025, 0.2475, 100);
for i = 1:length(Kt_vary)
    G_num = [((-Gv*Kt_vary(i).*Kv*g)/(L_1*J_eff_10*N*M))];
    G_denom = [1, (R_1/L_1 + Bm/J_eff_10), ((Kt_vary(i).*Kt_vary(i))/(L_1*J_eff_10) + (R_1*Bm)/(L_1*J_eff_10)), 0, 0, 0 ];
    G_1_vary_kt(i) = tf(G_num, G_denom);
    G_1_vary_kt_zpk(i) = zpk(G_1_vary_kt(i));
end

% This plots the pole zero map for the varying Kt to determine the effect
% of the varying kt on pole and zero locations
figure(3)
pzmap(G_1_vary_kt);
grid on;

% This plots the variance in Kt versus the gain of the plant
figure(4)
plot(Kt_vary, [G_1_vary_kt_zpk.K]);
grid on
xlabel('Kt');
ylabel('Gain, K');
title('Varianace in Kt vs Gain');

% It is determined that the effect of varying Kt has little effect on the
% pole and zero locations. For the gain, the values differ linearly.


%% Find the variation of N with constant Kt
close all
% To understand the effect of variance on N, motor 1 was used as the
% effects on the other motors would have similar effects.
Kt = 0.225;                     % Keep kt constant at average value 
N = linspace(10,50,100);        % Set up range for turn ratios
for i = 1:length(N)
    G_num = [((-Gv*Kt*Kv*g)./(L_1.*J_eff(i).*N(i)*M))];
    G_denom = [1, (R_1/L_1 + Bm./J_eff(i)), ((Kt*Kt)./(L_1*J_eff(i)) + (R_1*Bm)./(L_1.*J_eff(i))), 0, 0, 0 ];
    G_1_vary_N(i) = tf(G_num, G_denom);
    G_1_vary_N_zpk(i) = zpk(G_1_vary_N(i));
end

% This plots the pole zero map for the varying N to determine the effect
% of the varying N on pole and zero locations
figure(5)
pzmap(G_1_vary_N);
grid on;

% This plots the variance N versus the gain of the plant
figure(6)
plot(N, [G_1_vary_N_zpk.K]);
grid on
xlabel('N');
ylabel('Gain, K');
title('Varianace in N vs Gain');

% It is determined that the effects of varying N is more significant for
% the gain values as they are related to a negative inverse relationship.
% For the pole and zero locations, the effect of N is very slightly more
% significant then the effect of varying Kt. With the farther away closer
% pole being father left.


%% Variance for all motors
close all
% This section focuses on the effects on pole locations by varying N while
% making Kt be the extremes and average values.
% Motor 1
N = linspace(10,50,100);
for j = 1:length(Kt_1)  % Iterates through 3 Kt values for Motor 1 
for i = 1:length(N)
    G1_num = [((-Gv*Kt_1(j)*Kv*g)./(L_1.*J_eff(i).*N(i)*M))];
    G1_denom = [1, (R_1/L_1 + Bm./J_eff(i)), ((Kt_1(j)*Kt_1(j))./(L_1*J_eff(i)) + (R_1*Bm)./(L_1.*J_eff(i))), 0, 0, 0 ];
    G_1(i) = tf(G1_num, G1_denom);
    G_1_zpk(i) = zpk(G_1(i));
end
hold on
figure(7)
subplot(3,1,j)          % This plots the pzmap for all 3 values of Kt
pzmap(G_1);
grid on
end

% Motor 2
N = linspace(10,50,100);
for j = 1:length(Kt_2)  % Iterates through 3 Kt values for Motor 2
for i = 1:length(N)
    G2_num = [((-Gv*Kt_2(j)*Kv*g)./(L_2.*J_eff(i).*N(i)*M))];
    G2_denom = [1, (R_2/L_2 + Bm./J_eff(i)), ((Kt_2(j)*Kt_2(j))./(L_2*J_eff(i)) + (R_2*Bm)./(L_2.*J_eff(i))), 0, 0, 0 ];
    G_2(i) = tf(G2_num, G2_denom);
    G_2_zpk(i) = zpk(G_2(i));
end
hold on
figure(8)
subplot(3,1,j)          % This plots the pzmap for all 3 values of Kt
pzmap(G_2);
grid on
end

% Motor 3
N = linspace(10,50,100);
for j = 1:length(Kt_3)  % Iterates through 3 Kt values for Motor 3
for i = 1:length(N)
    G3_num = [((-Gv*Kt_3(j)*Kv*g)./(L_3.*J_eff(i).*N(i)*M))];
    G3_denom = [1, (R_3/L_3 + Bm./J_eff(i)), ((Kt_3(j)*Kt_3(j))./(L_3*J_eff(i)) + (R_3*Bm)./(L_3.*J_eff(i))), 0, 0, 0 ];
    G_3(i) = tf(G3_num, G3_denom);
    G_3_zpk(i) = zpk(G_3(i));
end
hold on
figure(9)
subplot(3,1,j)          % This plots the pzmap for all 3 values of Kt
pzmap(G_3);
grid on
end

% Motor 4
N = linspace(10,50,100);
for j = 1:length(Kt_4)  % Iterates through 3 Kt values for Motor 4
for i = 1:length(N)
    G4_num = [((-Gv*Kt_4(j)*Kv*g)./(L_4.*J_eff(i).*N(i)*M))];
    G4_denom = [1, (R_4/L_4 + Bm./J_eff(i)), ((Kt_4(j)*Kt_4(j))./(L_4*J_eff(i)) + (R_4*Bm)./(L_4.*J_eff(i))), 0, 0, 0 ];
    G_4(i) = tf(G4_num, G4_denom);
    G_4_zpk(i) = zpk(G_4(i));
end
hold on
figure(10)
subplot(3,1,j)          % This plots the pzmap for all 3 values of Kt
pzmap(G_4);
grid on
end
%i suc pepe
% It is determined that a turn ratio of 10 makes the pole locations farther
% to the left. This helps with the design process, therefore for the later
% sections N = 10. It was also not clear which motor was the best therefore
% further investigation is necessary.


%% Single Motor Plots to determine pole locations for N = 10
close all
N = 10;
G1_num_Single = [((-Gv*0.225*Kv*g)./(L_1.*J_eff_10.*N*M))];
G1_denom_Single = [1, (R_1/L_1 + Bm./J_eff_10), ((0.225*0.225)./(L_1*J_eff_10) + (R_1*Bm)./(L_1.*J_eff_10)), 0, 0, 0 ];
G_1_Single = tf(G1_num_Single, G1_denom_Single);
G_1_zpk = zpk(G_1_Single);

G2_num_Single = [((-Gv*0.175*Kv*g)./(L_2.*J_eff_10.*N*M))];
G2_denom_Single = [1, (R_2/L_2 + Bm./J_eff_10), ((0.175*0.175)./(L_2*J_eff_10) + (R_2*Bm)./(L_2.*J_eff_10)), 0, 0, 0 ];
G_2_Single = tf(G2_num_Single, G2_denom_Single);
G_2_zpk = zpk(G_2_Single);

G3_num_Single = [((-Gv*0.125*Kv*g)./(L_3.*J_eff_10.*N*M))];
G3_denom_Single = [1, (R_3/L_3 + Bm./J_eff_10), ((0.125*0.125)./(L_3*J_eff_10) + (R_3*Bm)./(L_3.*J_eff_10)), 0, 0, 0 ];
G_3_Single = tf(G3_num_Single, G3_denom_Single);
G_3_zpk = zpk(G_3_Single);

G4_num_Single = [((-Gv*0.275*Kv*g)./(L_4.*J_eff_10.*N*M))];
G4_denom_Single = [1, (R_4/L_4 + Bm./J_eff_10), ((0.275*0.275)./(L_4*J_eff_10) + (R_4*Bm)./(L_4.*J_eff_10)), 0, 0, 0 ];
G_4_Single = tf(G4_num_Single, G4_denom_Single);
G_4_zpk = zpk(G_4_Single);

% This plots the pzmap for all single motors
figure(11)
subplot(4,1,1)
pzmap(G_1_Single);
grid on
subplot(4,1,2)
pzmap(G_2_Single);
grid on
subplot(4,1,3)
pzmap(G_3_Single);
grid on
subplot(4,1,4)
pzmap(G_4_Single);
grid on

% From the above plots, looking at the slower of the two poles to the left
% of the imaginary axis, motor 1 and 4 can be considered as the pole
% locations are very similar. The next sections will look at these motors
% to determine the best motor to use.


%% Controller design for Motors 1 and 4
close all
% Motor 4
% figure(11)
% rlocus(G_1_Single)
C_1 = C;
C = tf(C);
% rltool(G_1_Single)
L_num1 = conv2(G_1_Single.Numerator{1,:},C.Numerator{1,:});
L_denom1 = conv2(G_1_Single.Denominator{1,:}, C.Denominator{1,:});
T_num1 = L_num1;
T_denom1 = L_num1 + L_denom1;
Tf1 = tf(T_num1, T_denom1);

% % Motor 4
% % figure(11)
% % rlocus(G_4_Single)
% C_4 = tf(C_4);
% % rltool(G_4_Single)
% L_num4 = conv2(G_4_Single.Numerator{1,:},C_4.Numerator{1,:});
% L_denom4 = conv2(G_4_Single.Denominator{1,:}, C_4.Denominator{1,:});
% T_num4 = L_num4;
% T_denom4 = L_num4 + L_denom4;
% Tf4 = tf(T_num4, T_denom4);

% Create an input pulse signal
% [Vel1, Acc1, Pos1] = spulse([0 2.5 2.5 5], 0.1962, 1);
[Vel1, Acc1, Pos1] = spulse([0 0.25 0.75 1] * 7, 0.0936, 1);
t = linspace(0, 120, 1000);

% Motor 1
figure(21)
lsim(Tf1,Vel1(t),t)
title('Velocity vs Time for Closed Loop Tf for Motor 1');
figure(22)
lsim(Tf1,Acc1(t),t)
title('Acceleration vs Time for Closed Loop Tf for Motor 1');
figure(23)
lsim(Tf1,Pos1(t),t)
title('Position vs Time for Closed Loop Tf for Motor 1');

% Motor 4

%% Control Effort Transfer function characteristics for Motor 1
Tu1 = Tf1/G_1_Single;         % Find the control effort from T/G
Tu_Simplify_1 = minreal(Tu1)  % Simplify for use in report
figure(13)
lsim(Tu1,Pos1(t),t)                 % Find the response to the desired input
title('Control Effort for Motor 1');


%% Input to Current Transfer Function Characteristics for Motor 1
% Find the input to current path for I(s)/R(s)
C_1_zpk = zpk(C);
N = 10;                     % Set N = 10
for i = 1:length(Kt_1)         
    s4 = tf([1 0 0 0 0],[1]);   % make a s^4 term
    Is = ((K*C.*Gv)/L_1)*s4;  % Numerator for I(s)/R(s)
    Rs1 = tf([1, (R_1/L_1 + Bm./J_eff_10), ((Kt_1(i)*Kt_1(i))./(L_1*J_eff_10) + (R_1*Bm)./(L_1.*J_eff_10)), 0, 0, 0],[1]);
    Rs2 = tf([((C.*Gv*Kt_1(i)*g*Kv)/(L_1*J_eff_10*N*M))],[1]);
    Rs = Rs1 - Rs2;             % Denominator for I(s)/R(s)
    IRS(i) = Is/Rs;                % Find I(s)/R(s)
    IRS_Simplify1 = minreal(IRS) % Simplify for use in report
    figure(14)
    hold on
    lsim(IRS(i),Pos1(t),t)               % Find the response to the desired input
    title('Input to Current Response Motor 1');
    hold off
end
legend('Kt = 0.2025', 'Kt = 0.225','Kt = 0.2475', 'Location', 'northeast');


%% Input to phase Transfer Function Characteristics for Motor 1
% Find the input to phase path for O(s)/R(s)
for i = 1:length(Kt_1)    
    s2 = tf([1 0 0],[1]);       % make a s^2 term
    Os = ((K*C*Gv*Kt_1)/(L_1*J_eff_10*N))*s2;     % Numerator for O(s)/R(s)   
    ORS = Os/Rs;                % Find O(s)/R(s)
    ORS_Simplify1 = minreal(ORS) % Simplify fo use in report
    figure(15)
    hold on
    lsim(ORS(i),Pos1(t),t)               % Find the response to the desired input
    title('Input to Phase Response Motor 1');
    hold off
end
legend('Kt = 0.2025', 'Kt = 0.225','Kt = 0.2475', 'Location', 'northeast');


%% Input to Output Transfer Function Characteristics for Motor 1
% Find the input to phase path for Y(s)/R(s)
for i = 1:length(Kt_1)    
    Ys = ((K*C*Gv*Kt_1*-g*Kv)/(L_1*J_eff_10*N*M));     % Numerator for Y(s)/R(s)   
    YRS = Ys/Rs;                % Find Y(s)/R(s)
    YRS_Simplify1 = minreal(YRS) % Simplify fo use in report
    figure(16)
    hold on
    lsim(YRS(i),Pos1(t),t)               % Find the response to the desired input
    title('Input to Output Response Motor 1');
    hold off
end
legend('Kt = 0.2025', 'Kt = 0.225','Kt = 0.2475', 'Location', 'northeast');


%% Control Effort Transfer function characteristics for Motor 4
Tu4 = Tf4/G_4_Single;         % Find the control effort from T/G
Tu_Simplify_4 = minreal(Tu4)   % Simplify for use in report
figure(17)
lsim(Tu4,Pos1(t),t)                % Find the response to the desired input
title('Control Effort for Motor 4');


%% Input to Current Transfer Function Characteristics for Motor 4
% Find the input to current path for I(s)/R(s)
C_4_zpk = zpk(C_4);
N = 10;                     % Set N = 10
for i = 1:length(Kt_4)         
    s4 = tf([1 0 0 0 0],[1]);   % make a s^4 term
    Is = ((K*C_4.*Gv)/L_4)*s4;  % Numerator for I(s)/R(s)
    Rs1 = tf([1, (R_4/L_4 + Bm./J_eff_10), ((Kt_4(i)*Kt_4(i))./(L_4*J_eff_10) + (R_4*Bm)./(L_4.*J_eff_10)), 0, 0, 0],[1]);
    Rs2 = tf([((C_4.*Gv*Kt_4(i)*g*Kv)/(L_4*J_eff_10*N*M))],[1]);
    Rs = Rs1 - Rs2;             % Denominator for I(s)/R(s)
    IRS(i) = Is/Rs;                % Find I(s)/R(s)
    IRS_Simplify4 = minreal(IRS) % Simplify for use in report
    figure(18)
    hold on
    lsim(IRS(i),Pos1(t),t)               % Find the response to the desired input
    title('Input to Current Response Motor 4');
    hold off
end
legend('Kt = 0.1125', 'Kt = 0.125','Kt = 0.1375', 'Location', 'northeast');


%% Input to phase Transfer Function Characteristics for Motor 4
% Find the input to phase path for O(s)/R(s)
for i = 1:length(Kt_4)    
    s2 = tf([1 0 0],[1]);       % make a s^2 term
    Os = ((K*C_4*Gv*Kt_4)/(L_4*J_eff_10*N))*s2;     % Numerator for O(s)/R(s)   
    ORS = Os/Rs;                % Find O(s)/R(s)
    ORS_Simplify4 = minreal(ORS) % Simplify fo use in report
    figure(19)
    hold on
    lsim(ORS(i),Pos1(t),t)               % Find the response to the desired input
    title('Input to Phase Response Motor 4');
    hold off
end
legend('Kt = 0.1125', 'Kt = 0.125','Kt = 0.1375', 'Location', 'northeast');


%% Input to Output Transfer Function Characteristics for Motor 4
% Find the input to phase path for Y(s)/R(s)
for i = 1:length(Kt_4)    
    Ys = ((K*C_4*Gv*Kt_4*-g*Kv)/(L_4*J_eff_10*N*M));     % Numerator for Y(s)/R(s)   
    YRS = Ys/Rs;                % Find Y(s)/R(s)
    YRS_Simplify4 = minreal(YRS) % Simplify fo use in report
    figure(20)
    hold on
    lsim(YRS(i),Pos1(t),t)               % Find the response to the desired input
    title('Input to Output Response Motor 4');
    hold off
end
legend('Kt = 0.1125', 'Kt = 0.125','Kt = 0.1375', 'Location', 'northeast');

