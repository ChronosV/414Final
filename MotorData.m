function Motor=MotorData(num)

%returns the motor parameters

Motor = struct();

if num==1
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt = 0.225;
    Re = 8;
    L = 25e-3;

elseif num==2
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt = 0.175;
    Re = 6;
    L = 16e-3;

elseif num==3
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt = 0.125;
    Re = 4;
    L = 7.5e-3;

elseif num==2
    Jm = 5.0e-5;
    Bm = 3.0e-6;
    Kt = 0.275;
    Re = 12;
    L = 32e-3;
end
end