clc
clear;
close all;
tic
% Problem:
% sovling the wave equation, grad^2 u = alpha*d^2u/dt2, in a rectangular range
% [x0,x1,y0,y1] with the following boundary ocnditions:
%   for x=x0:  u=a;    % Dirichlet boundary condition 
%          or  dudx=a; % Neumann boundary condition
%   for x=x1:  u=b;   
%          or  dudx=b;
%   for y=y0:  u=c; 
%          or  dudx=c;    
%   for y=y1:  u=d; 
%          or  dudx=d;
% There should also be an initial condition, that is u0(x,y) at t=0;
% 
x0=0; x1=5; y0=0; y1=5;
b_types =[0,0,1,1]; % boundary types: '0' for Dirichlet boundary; '1' for Neumann
b_values=[0,0,0,0]; % boundary values
alpha=0.04;

%%%
% define the simulation area
dx=0.1;  % grid szie
dy=dx;
x=x0:dx:x1;
y=y0:dy:y1;
u0=zeros(numel(x),numel(y));


x_ini=x>1&x<2;
y_ini=y>1&y<2;
u0(x_ini,y_ini)=0.4*exp(-(x(x_ini)-1.5).^2/0.5^2)'*exp(-(y(y_ini)-1.5).^2/0.5^2);

%%% scale the simulation to Cellular Network
Vmax=0.4;  %[V] the max read voltage used for reading the memristor 
% scale the boundary values
bv4Mat=b_values;
% for Neuuman boundaries, scale the values according to dimesion change (from dx or dy to 1)
bv4Mat(b_types==1) = b_values(b_types==1)*dx; 
% The capacitor value should be defined by the heat transfer coefficient
% TODO should be also scaled by other values
C = alpha*(dx*dy); 

MatrixA = [ 0 , 1 , 0 ; 1, -3 ,1 ; 0, 1, 0];
MatrixB = [ 0 , 0 , 0; 0 , 0 , 0; 0, 0, 0];
Input_temp = zeros(numel(x),numel(y));

I = 0;
R_x = 1;
dt = 1e-2;
T = 2;
VxMatInt = u0;
MatrixU = u0;
[VxMatHist, VxStable, VyMatHist, VyStable] = simulate(VxMatInt,MatrixU,T,C,R_x,I,dt,MatrixA,MatrixB,Vmax,b_types,bv4Mat); 

u=VyStable; 
figure()
surf(x,y,u');
xlabel('x')
ylabel('y')
zlim([-0.4,0.4])
zlabel('Temp')
title('Scaled')
[caz,cel] = view(37,31);

toc
%% creating a GIF
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testNew_wave.gif';
for n = 1:10:size(VyMatHist,3)
    % Draw plot for y = x.^n
    u=VyMatHist(:,:,n); 
    surf(x,y,u','edgecolor','none');
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis([x0 x1 y0 y1 -0.4 0.4])
    [caz,cel] = view(37,31);
    %shading interp
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end