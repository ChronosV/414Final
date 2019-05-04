function Motor=MotorData(num)

%returns the motor parameters

Motor = struct();

if num==1
    Motor.Jm = 5.0e-5;
    Motor.Bm = 3.0e-6;
    Motor.Kt(1) = 0.225;  %nominal
    Motor.Kt(2) = 0.2025; %min
    Motor.Kt(3) = 0.2475; %max
    Motor.R = 8;
    Motor.L = 25e-3;

elseif num==2
    Motor.Jm = 5.0e-5;
    Motor.Bm = 3.0e-6;
    Motor.Kt(1) = 0.175;  %nominal
    Motor.Kt(2) = 0.1575; %min
    Motor.Kt(3) = 0.1925; %max
    Motor.R = 6;
    Motor.L = 16e-3;

elseif num==3
    Motor.Jm = 5.0e-5;
    Motor.Bm = 3.0e-6;
    Motor.Kt(1) = 0.125;  %nominal
    Motor.Kt(2) = 0.1125; %min
    Motor.Kt(3) = 0.1375; %max
    Motor.R = 4;
    Motor.L = 7.5e-3;

elseif num==2
    Motor.Jm = 5.0e-5;
    Motor.Bm = 3.0e-6;
    Motor.Kt(1) = 0.275;  %nominal
    Motor.Kt(2) = 0.2475; %min
    Motor.Kt(3) = 0.3025; %max
    Motor.R = 12;
    Motor.L = 32e-3;
end
end