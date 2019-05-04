function Motor=MotorData(num)

%returns the motor parameters

Motor = struct();

if num==1
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt(1) = 0.225;  %nominal
    Kt(2) = 0.2025; %min
    Kt(3) = 0.2475; %max
    Re = 8;
    L = 25e-3;

elseif num==2
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt(1) = 0.175;  %nominal
    Kt(2) = 0.1575; %min
    Kt(3) = 0.1925; %max
    Re = 6;
    L = 16e-3;

elseif num==3
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt(1) = 0.125;  %nominal
    Kt(2) = 0.1125; %min
    Kt(3) = 0.1375; %max
    Re = 4;
    L = 7.5e-3;

elseif num==2
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt(1) = 0.275;  %nominal
    Kt(2) = 0.2475; %min
    Kt(3) = 0.3025; %max
    Re = 12;
    L = 32e-3;
end
end