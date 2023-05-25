alltime = 500;
omega = 1;
A = 0.08;
nu = 0.1;
% pose1 - начальный тип статической остойчивости коробля - всегда первый непорежденный
% pose2 - конечный тип статической остойчивости коробля
pose1 = 1;
pose2 = 5;

% a0, a1, a3, a5 - коэффициенты полинома пятого порядка, задающего функцию восстанавливающего момента коробля.
% Порядок типов остойчивости: 

% 	1-ый неповрежденный     1-ый		2-ой		3-ий		4-ый		5-ый

a0=[		0               0           -0.2		0       	-0.07 		0.07]; 
a1=[		0.64            0.25		0.64		-0.64		-0.64		-0.64]; 
a3=[		-0.1            -0.1		-0.1 		2.5 		2.5 		2.5]; 
a5=[		-0.07           -0.05		-0.07		-1.3		-1.3		-1.3];

f1=figure;
f2=figure;

for time = 0:alltime
    proc = time/alltime;
    recovery = [(a5(pose1)+(a5(pose2)-a5(pose1))*proc) 0 (a3(pose1)+(a3(pose2)-a3(pose1))*proc) 0 (a1(pose1)+(a1(pose2)-a1(pose1))*proc) a0(pose1)+(a0(pose2)-a0(pose1))*proc];
    [time_dif,theta_dif] = ode45(@(t,y) DoDt(t,y, omega, A, nu, recovery),[0 alltime],[0 0]);
    figure(f1);
    plot(theta_dif(:,1),theta_dif(:,2), "b", 'LineWidth', 1.5)
    axis([-1 1 -1 1])
    drawnow  
    figure(f2);
    plot(time_dif,theta_dif(:,1), "b", 'LineWidth', 1.5)
    drawnow  
end

function DthetaDtime = DoDt(ti, th, omega, A, nu, recovery)
    DthetaDtime = [th(2);A*cos(omega*ti) - nu*th(2) - (recovery(6)+recovery(5)*th(1)+recovery(3)*th(1).^3+recovery(1)*th(1).^5)];
end