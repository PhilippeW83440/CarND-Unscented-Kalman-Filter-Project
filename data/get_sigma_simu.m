
L=load('laser_simu.txt');
[n_row, n_col]=size(L)

dt=0.05; % second
v_prev = 0;
accelerations = [];

%yawrate_prev = 0;
%yaw_accelerations = [];

for i=1:n_row
    vx = L(i, 1);
    vy = L(i, 2);
    v = sqrt(vx*vx+vy*vy);
    %yawrate = L(i, 9);
    if (i > 1)
        acc = (v - v_prev) / dt;
        accelerations = [accelerations; acc];
        
        %yaw_acc = (yawrate - yawrate_prev) / dt;
        %yaw_accelerations = [yaw_accelerations; yaw_acc];
    end
    v_prev = v;
    %yawrate_prev = yawrate;
end

mu_a = mean(accelerations);
std_a = std(accelerations);
fprintf('laser:: mu_a=%f std_a=%f\n', mu_a, std_a)

% mu_yawdd = mean(yaw_accelerations);
% std_yawdd = std(yaw_accelerations);
% fprintf('laser: mu_yawdd=%f std_yawdd=%f\n', mu_yawdd, std_yawdd)

figure(1);
plot(accelerations);

%figure(2);
%plot(yaw_accelerations);
