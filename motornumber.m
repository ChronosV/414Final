function Motor=motornumber(num)
%%
% Function to return motor parameters
%
% Idea taken from Dalton Binette
%
% Written 5/1/2018 
%
%%

Motor = struct();

if num==1
    Motor.Jm = 5e-5;
    Motor.Bm = 3e-6;
    Motor.Kt(1) = .225 * 1.1;   % max
    Motor.Kt(2) = .225 * 0.9;   % min
    Motor.Kt(3) = .225;         % nominal
    Motor.R = 8;
    Motor.L = .025;

elseif num==2
    Motor.Jm = 5e-5;
    Motor.Bm = 3e-6;
    Motor.Kt(1) = .175 * 1.1;   % max
    Motor.Kt(2) = .175 * 0.9;   % min
    Motor.Kt(3) = .175;         % nominal 
    Motor.R = 6;
    Motor.L = .016;

elseif num==3
    Motor.Jm = 5e-5;
    Motor.Bm = 3e-6;
    Motor.Kt(1) = .125 * 1.1;   % max
    Motor.Kt(2) = .125 * 0.9;   % min
    Motor.Kt(3) = .125;         % nominal  
    Motor.R = 4;
    Motor.L = .0075;

elseif num==4
    Motor.Jm = 5e-5;
    Motor.Bm = 3e-6;
    Motor.Kt(1) = .275 * 1.1;   % max
    Motor.Kt(2) = .275 * 0.9;   % min
    Motor.Kt(3) = .275;         % nominal  
    Motor.R = 12;
    Motor.L = .032;
end
end