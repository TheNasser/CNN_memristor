% runs a simulation of the CNN , with the given Input Matrix of each cell 
function [VxMatHist, VxMat, VyMatHist, VyMat] = simulate(VxMatInt,VuMat,T,C,R_x,I,dt,MatA,MatB,Vmax,b_types,b_values)
    [m,n]=size(VxMatInt);
    VxMatHist=zeros(m,n,numel(T));
    VyMatHist=zeros(m,n,numel(T));
    
    VxMat = VxMatInt;
    VyMat = 0.5 * ( abs(VxMat + Vmax) - abs(VxMat -Vmax)); 
    
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
                dV_x = (-(1/R_x)*VxMat(r,c) + SumA + SumB + I)/C;
                VxMat(r,c) = VxMat(r,c) + dV_x * dt;
            end
        end
    
    VyMat = 0.5 * ( abs(VxMat + Vmax) - abs(VxMat -Vmax)); 
    
    VxMatHist(:,:,count) = VxMat;
    VyMatHist(:,:,count) = VyMat;
    count = count + 1;
    t = t + dt;
    end
end