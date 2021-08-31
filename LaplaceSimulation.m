clc
clear;
close all;

% Problem:
% sovling the Laplace equation, grad^2 u = 0, in a rectangular range
% [x0,x1,y0,y1] with the following boundary ocnditions:
%   for x=x0:  u=a;    % Dirichlet boundary condition 
%          or  dudx=a; % Neumann boundary condition
%   for x=x1:  u=b;   
%          or  dudx=b;
%   for y=y0:  u=c; 
%          or  dudx=c;    
%   for y=y1:  u=d; 
%          or  dudx=d;
x0=0; x1=1; y0=0; y1=2;
b_types =[0,0,1,1]; % boundary types: '0' for Dirichlet boundary; '1' for Neumann
b_values=[0.4,0,1,0]; % boundary values

%%%
% define the simulation area
dx=0.2;  % grid szie
dy=dx;
x=x0:dx:x1;
y=y0:dy:y1;


%%% scale the simulation to Cellular Network
Vmax=0.4;  %[V] the max read voltage used for reading the memristor 
% scale the boundary values
bv4Mat=b_values;
% for Neuuman boundaries, scale the values according to dimesion change (from dx or dy to 1)
bv4Mat(b_types==1) = b_values(b_types==1)*dx; 

MatrixA = [ 0 , 1 , 0 ; 1, -3 ,1 ; 0, 1, 0];
MatrixB = [ 0 , 0 , 0; 0 , 0 , 0; 0, 0, 0];
Input_temp = zeros(numel(x),numel(y));

I = 0;
C = 0.4;  % this value could be as small as possible
R_x = 1;
dt = 1e-2;
T = 10;
N = T/dt;
t_Vec = linspace(0,T,N+2);
VxMatInt = Input_temp;
MatrixU = Input_temp;
[VxMatHist, VxStable, VyMatHist, VyStable] = simulate(VxMatInt,MatrixU,T,C,R_x,I,dt,MatrixA,MatrixB,Vmax,b_types,bv4Mat); 

% scale the value back
u=VyStable; 
figure()
surf(x,y,u');
xlabel('x')
ylabel('y')
zlabel('Potential')
title('Scaled')
[caz,cel] = view(37,31);
%% creating a GIF
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testNew_Laplace2.gif';
for n = 1:18:size(VyMatHist,3)
    % Draw plot for y = x.^n
    u=VyMatHist(:,:,n); 
    surf(x,y,u');
    xlabel('x')
    ylabel('y')
    zlabel('Potential')
     zlim([-0.3,0.4])
    %axis([1 20 1 20 -0.4 0.4])
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