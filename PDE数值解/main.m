%%
% a(r,theta) = 1 + r^2
% f(r, theta) = 8 * r^2 + 4,
% g(r, theta) = 0
% u(r, theta) = r^2 - 1
clc;clear;close all;

N = 160;              
M = 2*N;             
dr = 1/N;
dtheta = 2*pi/M;

r = linspace(0,1,N+1)';           
theta = linspace(0,2*pi,M+1);     
[TH,R] = meshgrid(theta,r);

u_exact = R.^2 - 1;
f = 8*R.^2 + 4;
a = 1 + R.^2;

A = sparse((N+1)*M,(N+1)*M);
F = zeros((N+1)*M,1);

for j = 1:M
    for i = 1:N+1
        idx = i + (j-1)*(N+1);

        if i==1
            idx_r1 = idx+1;

            A(idx,idx) = -4*a(i,j)/dr^2;
            A(idx,idx_r1) = 4*a(i,j)/dr^2;

            F(idx) = f(i,j);
            
        elseif i==N+1 
            A(idx,idx)=1;
            F(idx)=0;
            
        else 
            rp = r(i)+dr/2;
            rm = r(i)-dr/2;
            a_p = (a(i,j)+a(i+1,j))/2;
            a_m = (a(i,j)+a(i-1,j))/2;

            idx_ip1 = idx+1;
            idx_im1 = idx-1;

           
            Ar_plus = a_p/(dr^2) + a_p/(2*r(i)*dr);
            Ar_minus = a_m/(dr^2) - a_m/(2*r(i)*dr);
            A(idx,idx_ip1)=Ar_plus;
            A(idx,idx_im1)=Ar_minus;
            A(idx,idx)=-(Ar_plus+Ar_minus)-2*a(i,j)/(r(i)^2*dtheta^2);

            jp = mod(j,M)+1;
            jm = mod(j-2,M)+1;
            idx_jp = i+(jp-1)*(N+1);
            idx_jm = i+(jm-1)*(N+1);

            A(idx,idx_jp)=a(i,j)/(r(i)^2*dtheta^2);
            A(idx,idx_jm)=a(i,j)/(r(i)^2*dtheta^2);

            F(idx)=f(i,j);
        end
    end
end

u_numeric = A\F;
u_numeric = reshape(u_numeric,N+1,M);
u_numeric = [u_numeric u_numeric(:,1)]; 


error = abs(u_numeric-u_exact);
max_error = max(error(:));
fprintf('Max error = %e\n',max_error);

[X,Y] = pol2cart(TH,R);
mesh(X,Y,u_numeric);
xlabel('x');ylabel('y');zlabel('u');
title('Finite Difference Solution');
%% 
% a(x,y) = 1 + x^2 + y^2, a(r,theta) = 1 + r^2
% f(r, theta) = 2 * e^{r^2}[4*r^2+8*r^4+2]
% g(r, theta) = 0
% u(r, theta) = (r^2 - 1) * e^{r^2}
clc;clear;close all;

N = 200;              
M = 2*N;             
dr = 1/N;
dtheta = 2*pi/M;

r = linspace(0,1,N+1)';           
theta = linspace(0,2*pi,M+1);     
[TH,R] = meshgrid(theta,r);

u_exact = (R .^ 2 - 1).*exp(R.^2);
f = 2 .* exp(R.^2).*(4.*R.^2+8.*R.^4+2);
a = 1 + R.^2;

A = sparse((N+1)*M,(N+1)*M);
F = zeros((N+1)*M,1);

for j = 1:M
    for i = 1:N+1
        idx = i + (j-1)*(N+1);

        if i==1
            idx_r1 = idx+1;

            A(idx,idx) = -4*a(i,j)/dr^2;
            A(idx,idx_r1) = 4*a(i,j)/dr^2;

            F(idx) = f(i,j);
            
        elseif i==N+1 
            A(idx,idx)=1;
            F(idx)=0;
            
        else 
            rp = r(i)+dr/2;
            rm = r(i)-dr/2;
            a_p = (a(i,j)+a(i+1,j))/2;
            a_m = (a(i,j)+a(i-1,j))/2;

            idx_ip1 = idx+1;
            idx_im1 = idx-1;

           
            Ar_plus = a_p/(dr^2) + a_p/(2*r(i)*dr);
            Ar_minus = a_m/(dr^2) - a_m/(2*r(i)*dr);
            A(idx,idx_ip1)=Ar_plus;
            A(idx,idx_im1)=Ar_minus;
            A(idx,idx)=-(Ar_plus+Ar_minus)-2*a(i,j)/(r(i)^2*dtheta^2);

            jp = mod(j,M)+1;
            jm = mod(j-2,M)+1;
            idx_jp = i+(jp-1)*(N+1);
            idx_jm = i+(jm-1)*(N+1);

            A(idx,idx_jp)=a(i,j)/(r(i)^2*dtheta^2);
            A(idx,idx_jm)=a(i,j)/(r(i)^2*dtheta^2);

            F(idx)=f(i,j);
        end
    end
end

u_numeric = A\F;
u_numeric = reshape(u_numeric,N+1,M);
u_numeric = [u_numeric u_numeric(:,1)]; 


error = abs(u_numeric-u_exact);
max_error = max(error(:));
fprintf('Max error = %e\n',max_error);

[X,Y] = pol2cart(TH,R);
mesh(X,Y,u_numeric);
xlabel('x');ylabel('y');zlabel('u');
title('Finite Difference Solution');
%% 计算收敛阶，以第一个算例为例
clc;clear;close all
e = [5.048852e-03 1.288651e-03 3.241234e-04 8.117250e-05 2.030317e-05];
p = log(e(1:end-1) ./ e(2:end)) ./ log(2);
disp('计算得到的收敛阶p值为:')
disp(p)
mean_p = mean(p);
fprintf('平均收敛阶为 %.4f\n',mean_p)

