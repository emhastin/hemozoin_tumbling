classdef MATLAB_sim_hemozoin
   
    properties
      radius
      MM 
      mu
      lab_vel_tot
      brick
      position_matrix_rotated
      pos_tot
      rotation_matrix_cell
      bd
      rot_mat
      p_rot
      p_r
      Rdt
      body_trans_vel_tot
      body_rot_tot_up
      lab_ang_tot
   end
   
   methods
      
       
       
       function obj=MATLAB_sim_hemozoin(e, d, mu, p, X, Y, Z)
          
N = size(p,1);

cm = [mean(X) mean(Y) mean(Z)];  % center of mass = this should be the origin [0 0 0]

e = e*10^(-9);                   % converting deltaS from nm to m
epsvec = e*ones(1,N);
                                 % r matrix = list of original coordinates in m
r = p';                          % p' = correct orientation for r matrix
IRX = fIRX(N,d,r);               % IRX = matrix of kinematic velocities of each small sphere (node)
G = fmakeGregstv2(r,mu,epsvec);  % G = matrix of hydrodynamic velocities of each small sphere (node) using Regularized Stokeslet eqn
RM = IRX'*inv(G)*IRX;            % RM = resistance matrix
obj.MM = inv(RM);                % MM = mobility matrix, accounting for forces and torques on small spheres

       end

       
       % all calculations below are done from center of each brick (not
       % individual spheres)
       
       function obj=rotation_translation(obj, mu, p, radius, dt, t)
           T = 293.15;                                     % temperature of plasmodium food vacuole in K (~98.38F = avg human temp)
           kb = 1.38*10^(-23);                             % Boltzmann constant
           obj.radius=radius;                              % change in each time step (sec?)
           avg_vel = 0;                                    % average velocity of hemozoin crystal (due to Brownian motion)
           Rt = eye(3);                                    % rotational matrix to convert body basis to lab basis 
           stdev_trans = sqrt(2*abs(obj.MM)*kb*T*dt);      % standard deviation of displacement of hemozoin crystal
           stdev_ang = sqrt(2*abs(obj.MM)*kb*T*dt);
           obj.mu = mu;

           
%           for i = 1:4*length(t)  
%     rot_x = -2.5;
%     rot_y = 0.6;
%     rot_z = -0.6;
%     
    rot_x = normrnd(avg_vel,stdev_ang(4,4))/dt;     % generates random # for omega (angular velocity in radians/sec) 
    rot_y = normrnd(avg_vel,stdev_ang(5,5))/dt;     % diagonals of MM are forces in x,y,z followed by torques in x,y,z
    rot_z = normrnd(avg_vel,stdev_ang(6,6))/dt;     % stdev = root mean square, only nonzero on diagonal
    
    body_rot_tot(:,1) = [rot_x; rot_y; rot_z];
    
    obj.body_rot_tot_up = body_rot_tot;
    
% converting from body basis into lab basis with rotations
    w = [rot_x rot_y rot_z];                 % making vector of omega in order to calculate magnitude
    % could add LJ velocity here or add 76-91 in other loop
    
    
    omegas = (1/norm(w)^2)*[ rot_x^2,  rot_x* rot_y,  rot_x* rot_z;... % calculating sub components of R matrix 
                             rot_y* rot_x,  rot_y^2,  rot_y* rot_z;...
                             rot_z* rot_x,  rot_z* rot_y,  rot_z^2];
   
    omegax = (1/norm(w))*[0, - rot_z,  rot_y;...              % making omega cross (wx)
                           rot_z, 0, - rot_x;...
                         - rot_y,  rot_x, 0];
    
    obj.Rdt = omegas+cos(norm(w)*dt)*(eye(3)-omegas)+sin(norm(w)*dt)*omegax;    % infinitesimal rotation based on dt
    Rt = Rt*obj.Rdt;
    
%     if abs(det(Rt)-1)>1e-6
%         det(Rt)
%         error('hi')
%     
%     end
%    
    obj.rotation_matrix_cell = Rt;     % rotational matrix changes for each time step according to angular velocity

    trans_vel_x = normrnd(avg_vel,stdev_trans(1,1))/dt;    % generates random translational velocity from Gaussian distribution
    trans_vel_y = normrnd(avg_vel,stdev_trans(2,2))/dt;    % in x y and z orientations for displacements
    trans_vel_z = normrnd(avg_vel,stdev_trans(3,3))/dt;
    
    obj.body_trans_vel_tot(:,1) = [trans_vel_x; trans_vel_y; trans_vel_z];
    
% body basis - orientation with respect to each hemozoin movement, lab basis - overall movement relative to the system of reference
%     obj.lab_vel_tot(:,:) = Rt*obj.body_trans_vel_tot(:,:);      % adjusting velocity to include translation relative to the reference system
%     obj.lab_ang_tot(:,:) = Rt*obj.body_rot_tot_up(:,:);
    
% end
% initializing displacement at 0
% obj.posx = zeros(length(p),1);
% obj.posy = zeros(length(p),2);
% obj.posz = zeros(length(p),3);
% 
% % randomizing starting positions
% obj.posx(1) = radius/2*(-1+(2).*rand);
% obj.posy(1) = radius/2*(-1+(2).*rand);
% obj.posz(1) = radius/2*(-1+(2).*rand);
%%
%-------------------
% counter = 2;
% for j = 2:length(t)+1
% %    distance(j)=sqrt(posx(j-1)^2+posy(j-1)^2+posz(j-1)^2);
% %     if (distance(j)>radius-0.0000005)
% % %         obj.radius-0.0000004
% %     posx(j) = posx(j-1)+lab_vel_tot(1,j-1)*dt-posx(j-1)*0.000001*mu/distance(j)*dt;  % adding distance from initial position
% %     posy(j) = posy(j-1)+lab_vel_tot(2,j-1)*dt-posy(j-1)*0.000001*mu/distance(j)*dt;  % finding distance by using lab velocity*time step
% %     posz(j) = posz(j-1)+lab_vel_tot(3,j-1)*dt-posz(j-1)*0.000001*mu/distance(j)*dt;
% %     else
%     %-------------------
%     while true
%         obj.posx(j) = obj.posx(j-1)+obj.lab_vel_tot(1,counter-1)*dt;
%         obj.posy(j) = obj.posy(j-1)+obj.lab_vel_tot(2,counter-1)*dt;
%         obj.posz(j) = obj.posz(j-1)+obj.lab_vel_tot(3,counter-1)*dt;
%         distance(j)=sqrt(obj.posx(j)^2+obj.posy(j)^2+obj.posz(j)^2);
%         if ~ (distance(j) > radius)
%             counter = counter+1;
%             break
%         else
%             counter = counter+1;            
%     %-------------------
%    obj.posx = obj.posx+obj.lab_vel_tot(1,1)*dt;  % adding distance from initial position
%    obj.posy = obj.posy+obj.lab_vel_tot(2,1)*dt;  % finding distance by using lab velocity*time step
%    obj.posz = obj.posz+obj.lab_vel_tot(3,1)*dt;
%        
% %     end
% 
%         end
%     end
% end
%%
% obj.p_rot = p';
% obj.position_matrix_rotated(:,:) = obj.p_rot';


    
%     obj.rot_mat = obj.rotation_matrix_cell;
%     obj.p_rot = obj.rot_mat*obj.p_rot;
%     
%     obj.p_r = [obj.p_rot(1,:);obj.p_rot(2,:);obj.p_rot(3,:)]';
%     obj.position_matrix_rotated(:,:) = obj.p_r;


%     obj.pos_tot = [posx; posy; posz]';      % total movement in x,y,z for each time step

       end
       
   end      
     
       
   end
    
    
