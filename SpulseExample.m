% spulse.m usage examples

t1 = 1; % when pulse starts rising

t2 = 2; % when pulse reaches its final amplitude

t3 = 3; % when pulse starts falling

t4 = 4; % when pulse returns to zero

T = [t1 t2 t3 t4]; % pack times into a vector

F = 5;  % pulse amplitude

% consider each pulse transition type
Type = true; % one half sine wave transition
[fun1,dfun1,ifun1] = spulse(T,F,Type);

Type = false; % one quarter sine wave transition
[fun2,dfun2,ifun2] = spulse(T,F,Type);


t = linspace(0,5,300); % time points where to evaluate pulses

% the outputs of spulse are functions! so simply evaluate them at the
% desired data points

subplot(2,1,1)
plot(t,fun1(t),t,fun2(t))
xlabel('Time')
ylabel('Ball Velocity')
legend('half sine','quarter sine')


subplot(2,1,2)
plot(t,ifun1(t),t,ifun2(t))
xlabel('Time')
ylabel('Ball Position')
legend('half sine','quarter sine')


% so mess with T and F so that the ball position starts at zero and ends at
% five or a little less.
% if you plot the derivative dfun(t), you get the ball acceleration. The
% larger the acceleration, the more control effort is needed. You know this
% from driving your car!!!

% lsim simulates linear systems (Duh!) for arbitrary inputs. The input you
% give it are data points in time t, and the system input u(t) at each of
% the data points. These data points must be close to together to make the
% simulation behave properly since lsim does not interpolate u(t) to get
% inputs between the data points. You should not need more than 1000 data
% points to make it work well. In fact it might work well with one half
% that many data points.

