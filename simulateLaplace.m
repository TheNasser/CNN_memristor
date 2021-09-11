% runs a simulation of the CNN , with the given Input Matrix of each cell 

% (VxMatInt: the inital Vx voltage for all cells ,VuMat the initail vilatge for all cells 
%T: simulation time ,C : capacitor for every cell,R_x : resistor for every cell ,I : bias current in every cell ,
%dt is the time intervals , MatA is matrix A,MatB is Matrix B ,Vmax,si the max voltage,
%b_types is the boundary types 1 for neuman 0 for direchlet,b_values the value for the boundaries (both should be of size 4 for the 4 boundaries 
function [VxMatHist, VxMat, VyMatHist, VyMat] = simulate(VxMatInt,VuMat,T,C,R_x,I,dt,MatA,MatB,Vmax,b_types,b_values)
    [m,n]=size(VxMatInt);
    VxMatHist=zeros(m,n,numel(T));
    VyMatHist=zeros(m,n,numel(T));
    
    VxMat = VxMatInt;
    VyMat = 0.5 * ( abs(VxMat + Vmax) - abs(VxMat -Vmax)); 
    
    dVxdt=zeros(m,n);
    
    count = 1;
    t = 0;
    while t <= T
        % boundary conditions
        if b_types(1)==0
                VxMat(1,:)=b_values(1);
        else
                VxMat(1,:)=VxMat(2,:)-b_values(1);
        end
        if b_types(2)==0
                VxMat(end,:)=b_values(2);
        else
                VxMat(end,:)=VxMat(end-1,:)+b_values(2);
        end
        if b_types(3)==0
                VxMat(:,1)=b_values(3);
        else
                VxMat(:,1)=VxMat(:,2)-b_values(3);
        end
        if b_types(4)==0
                VxMat(:,end)=b_values(4);
        else
                VxMat(:,end)=VxMat(:,end-1)+b_values(4);
        end
        
        % within the matrix (non-boundary elements)
        for r = 2:m-1
            for c = 2:n-1
                SumA = sum(MatA.*VyMat(r-1:r+1,c-1:c+1),'all');
                SumB = sum(MatB.*VuMat(r-1:r+1,c-1:c+1),'all');
                ddVx_dtdt = (-(1/R_x)*VxMat(r,c) + SumA + SumB + I)/C;
                % to implement wave function, revision on the circuit is
                % needed, to implement I=d^2V/dt^2
                dVxdt(r,c) = dVxdt(r,c) + ddVx_dtdt*dt;
            end
        end
        VxMat = VxMat + dVxdt*dt;
    
        VyMat = 0.5 * ( abs(VxMat + Vmax) - abs(VxMat -Vmax)); 
        VxMatHist(:,:,count) = VxMat;
        VyMatHist(:,:,count) = VyMat;
        count = count + 1;
        t = t + dt;
    end
end