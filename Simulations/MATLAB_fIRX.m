function IRX = fIRX(n,d,r)

    if d == 2
% creating IRX matrix
        IRX = zeros(n*d,3);
        
% creating d columns as tiled identity matrix w/ n rows
        IRX(:,1:2) = repmat(eye(d),n,1);
        
% specifying components of submatrix (to add to IRX matrix)
% s = submatrix in both 2D and 3D
% s matrix should specify r vectors wrt origin
        for i=1:n
            s = zeros(2,1);

% r already specifies distances wrt to origin --> just need to plug them in
            s(1,1) = -r(2,i);
            s(2,1) = r(1,i);
            IRX(2*(i-1)+1:2*i,3) = s;
        end  
    end

    if d == 3

% creating IRX matrix
        IRX = zeros(n*d,6); 

% creating d columns as tiled identity matrix
        IRX(:,1:3) = repmat(eye(d),n,1);

% defining resistance matrix
        for i = 1:n        
            s = zeros(3);
            
% creating top triangle of 3x3 matrix to add to IRX
            s(1,2) = r(3,i);
            s(1,3) = -r(2,i);
            s(2,3) = r(1,i);

% flipping matrix over center diagonal
            s = s - s';

% tacking on new matrix to IRX matrix
            IRX(3*(i-1)+1:3*i,4:6) = s;
        end
    end       
end
