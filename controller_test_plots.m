%% Controller design for Motors 1 and 4
close all
% Motor 4
% figure(11)
% rlocus(G_1_Single)
% C_1 = C;
C = tf(C);
% rltool(G_1_Single)
L_num = conv2(G.num{1, :}, C.num{1, :});
L_den = conv2(G.den{1, :}, C.den{1, :});
T_num = L_num;
T_den = L_num + L_den;
Tf = tf(T_num, T_den);

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
[Vel, Acc, Pos] = spulse([0 0.25 0.75 1] * 7, 0.0936, 1);
t = linspace(0, 120, 1000);

% Motor 1
figure(21);
lsim(Tf, Vel(t), t);
title('Velocity vs Time for Closed Loop Tf for Motor');
figure(22);
lsim(Tf, Acc(t), t);
title('Acceleration vs Time for Closed Loop Tf for Motor');
figure(23);
lsim(Tf, Pos(t), t);
title('Position vs Time for Closed Loop Tf for Motor');

% Motor 4

%% Control Effort Transfer function characteristics for Motor 1
Tu = Tf/G;         % Find the control effort from T/G
Tu_Simplify = minreal(Tu);    % Simplify for use in report
figure(13);
lsim(Tu, Pos(t), t);                 % Find the response to the desired input
title('Control Effort for Motor');


%% Input to Current Transfer Function Characteristics for Motor 1
% Find the input to current path for I(s)/R(s)
C_zpk = zpk(C);
% N = 10;          % Set N = 10
K = 1;           % Pre-amp gain.

s4 = tf([1 0 0 0 0], [1]);   % make a s^4 term
Is = ((K*C.*Gv)/L)*s4;       % Numerator for I(s)/R(s)
Rs1 = tf([1, (R/L + Bm./Jeff), ((Kt*Kt)./(L*Jeff) + (R*Bm)./(L.*Jeff)), 0, 0, 0],[1]);
Rs2 = tf([((C.*Gv*Kt*g*Kv)/(L*Jeff*N*Mb))], [1]);
Rs = Rs1 - Rs2;              % Denominator for I(s)/R(s)
IRS = Is/Rs;                 % Find I(s)/R(s)
IRS_Simplify = minreal(IRS); % Simplify for use in report
figure(14);
hold on;
lsim(IRS, Pos(t), t);               % Find the response to the desired input
title('Input to Current Response Motor');
hold off;

% legend('Kt = 0.2025', 'Kt = 0.225', 'Kt = 0.2475', 'Location', 'northeast');


%% Input to phase Transfer Function Characteristics for Motor 1
% Find the input to phase path for O(s)/R(s)
  
s2 = tf([1 0 0], [1]);                % make a s^2 term
Os = ((K*C*Gv*Kt)/(L*Jeff*N))*s2;     % Numerator for O(s)/R(s)   
ORS = Os/Rs;                          % Find O(s)/R(s)
ORS_Simplify = minreal(ORS);          % Simplify fo use in report
figure(15);
hold on;
lsim(ORS, Pos(t), t);                 % Find the response to the desired input
title('Input to Phase Response Motor');
hold off;

% legend('Kt = 0.2025', 'Kt = 0.225', 'Kt = 0.2475', 'Location', 'northeast');


%% Input to Output Transfer Function Characteristics for Motor 1
% Find the input to phase path for Y(s)/R(s)
  
Ys = ((K*C*Gv*Kt*-g*Kv)/(L*Jeff*N*Mb));     % Numerator for Y(s)/R(s)   
YRS = Ys/Rs;                                % Find Y(s)/R(s)
YRS_Simplify = minreal(YRS);                % Simplify fo use in report
figure(16);
hold on;
lsim(YRS, Pos(t), t);                       % Find the response to the desired input
title('Input to Output Response Motor');
hold off;

% legend('Kt = 0.2025', 'Kt = 0.225', 'Kt = 0.2475', 'Location', 'northeast');


%% Control Effort Transfer function characteristics for Motor 4
%Tu4 = Tf4/G_4_Single;         % Find the control effort from T/G
%Tu_Simplify_4 = minreal(Tu4)   % Simplify for use in report
%figure(17)
%lsim(Tu4,Pos1(t),t)                % Find the response to the desired input
%title('Control Effort for Motor 4');


%% Input to Current Transfer Function Characteristics for Motor 4
% Find the input to current path for I(s)/R(s)
%C_4_zpk = zpk(C_4);
%N = 10;                     % Set N = 10
%for i = 1:length(Kt_4)         
%    s4 = tf([1 0 0 0 0],[1]);   % make a s^4 term
%    Is = ((K*C_4.*Gv)/L_4)*s4;  % Numerator for I(s)/R(s)
%    Rs1 = tf([1, (R_4/L_4 + Bm./J_eff_10), ((Kt_4(i)*Kt_4(i))./(L_4*J_eff_10) + (R_4*Bm)./(L_4.*J_eff_10)), 0, 0, 0],[1]);
%    Rs2 = tf([((C_4.*Gv*Kt_4(i)*g*Kv)/(L_4*J_eff_10*N*M))],[1]);
%    Rs = Rs1 - Rs2;             % Denominator for I(s)/R(s)
%    IRS(i) = Is/Rs;                % Find I(s)/R(s)
%    IRS_Simplify4 = minreal(IRS) % Simplify for use in report
%    figure(18)
%    hold on
%    lsim(IRS(i),Pos1(t),t)               % Find the response to the desired input
%    title('Input to Current Response Motor 4');
%    hold off
%end
%legend('Kt = 0.1125', 'Kt = 0.125','Kt = 0.1375', 'Location', 'northeast');


%% Input to phase Transfer Function Characteristics for Motor 4
% Find the input to phase path for O(s)/R(s)
%for i = 1:length(Kt_4)    
%    s2 = tf([1 0 0],[1]);       % make a s^2 term
%    Os = ((K*C_4*Gv*Kt_4)/(L_4*J_eff_10*N))*s2;     % Numerator for O(s)/R(s)   
%    ORS = Os/Rs;                % Find O(s)/R(s)
%    ORS_Simplify4 = minreal(ORS) % Simplify fo use in report
%    figure(19)
%    hold on
%    lsim(ORS(i),Pos1(t),t)               % Find the response to the desired input
%    title('Input to Phase Response Motor 4');
%    hold off
%end
%legend('Kt = 0.1125', 'Kt = 0.125','Kt = 0.1375', 'Location', 'northeast');
%

%% Input to Output Transfer Function Characteristics for Motor 4
% Find the input to phase path for Y(s)/R(s)
%for i = 1:length(Kt_4)    
%    Ys = ((K*C_4*Gv*Kt_4*-g*Kv)/(L_4*J_eff_10*N*M));     % Numerator for Y(s)/R(s)   
%    YRS = Ys/Rs;                % Find Y(s)/R(s)
%    YRS_Simplify4 = minreal(YRS) % Simplify fo use in report
%    figure(20)
%    hold on
%    lsim(YRS(i),Pos1(t),t)               % Find the response to the desired input
%    title('Input to Output Response Motor 4');
%    hold off
%end
%legend('Kt = 0.1125', 'Kt = 0.125','Kt = 0.1375', 'Location', 'northeast');

