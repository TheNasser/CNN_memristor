clc
clear;
close all;

% Problem:
% sovling the heat equation, grad^2 u = alpha*du/dt, in a rectangular range
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
x0=0; x1=5; y0=0; y1=10;
b_types =[0,0,0,0]; % boundary types: '0' for Dirichlet boundary; '1' for Neumann
b_values=[0,0,0,0]; % boundary values
alpha=4;

%%%
% define the simulation area
dx=0.5;  % grid szie
dy=dx;
x=x0:dx:x1;
y=y0:dy:y1;
u0=zeros(numel(x),numel(y));
u0(x>1&x<4,y>1&y<9)=0.4;
u0(x>2&x<4,y>4&y<6)=0;

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
T = 10;
N = T/dt;
t_Vec = linspace(0,T,N+2);
VxMatInt = u0;
MatrixU = u0;
[VxMatHist, VxStable, VyMatHist, VyStable] = simulate(VxMatInt,MatrixU,T,C,R_x,I,dt,MatrixA,MatrixB,Vmax,b_types,bv4Mat); 
heatmap(VyStable,'Colormap',jet);

u=VyStable; 
figure()
surf(x,y,u');
xlabel('x')
ylabel('y')
zlim([0,0.4])
zlabel('Temp')
title('Scaled')
[caz,cel] = view(37,31);
%%
PlotOutput(t_Vec ,VxMatHist ,VyMatHist,N,3,11,T);
PlotOutput(t_Vec ,VxMatHist ,VyMatHist,N,6,3,T);
%% creating a GIF
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testNew_heat-dissipation2.gif';
for n = 1:18:size(VyMatHist,3)
    % Draw plot for y = x.^n
    u=VyMatHist(:,:,n); 
    surf(x,y,u');
    xlabel('x')
    ylabel('y')
    zlabel('Temp')
    axis([x0 x1 y0 y1 0 0.4])
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