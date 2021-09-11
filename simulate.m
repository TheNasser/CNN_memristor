% runs a simulation of the CNN , with the given Input Matrix of each cell 

% (VxMatInt: the inital Vx voltage for all cells ,VuMat the initail vilatge for all cells 
%T: simulation time ,C : capacitor for every cell,R_x : resistor for every cell ,I : bias current in every cell ,
%dt is the time intervals , MatA is matrix A,MatB is Matrix B.
function [VxMatHist, VxStable, VyMatHist, VyStable] = simulate(V_xMatInt,V_uMat,m,n,T,C,R_x,I,dt,MatA,MatB)
    V_yMat = zeros(m,n);
    V_yMat = computeVy(V_yMat,V_xMatInt,m,n);
   % V_yMat = V_xMatInt; % should be computeVy  in regular 
    V_xMatNew = V_xMatInt;
    V_xMatHist = V_xMatInt;
    V_yMatHist = V_yMat;
    count = 1;
    t = 0;
    while t <= T
        for r = 1:m
            for c = 1:n
             
          
                
                V_xMatNew(r,c) = V_xMatNew(r,c) + computeVx(C,R_x,I,dt,MatA,MatB,V_xMatNew(r,c),V_yMat,V_uMat,m,n,r,c);
          
                
                
            end
        end
    V_xMatHist(:,:,count) = V_xMatNew;
    
    % Matrix-wise operation would be fine, do not need to loop in the matrix
     V_yMat = V_xMatNew;
   % V_yMat = 0.5 * ( abs(V_xMatNew + Vmax) - abs(V_xMatNew -Vmax)); //  don't need Vy 
    
    V_yMat = computeVy(V_yMat,V_xMatNew,m,n);
    V_yMatHist(:,:,count) = V_yMat;
    count = count + 1;
    t = t + dt;
    end

    VxMatHist = V_xMatHist;
    VxStable = V_xMatNew;
    VyMatHist = V_yMatHist;
    VyStable = V_yMat;
end