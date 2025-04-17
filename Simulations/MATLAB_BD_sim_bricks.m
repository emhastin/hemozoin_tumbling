clc
clear all
close all

deltaS_xy = 85;            % distance between nodes in x and y directions (nm)
deltaS_z = 97.5;           % distance between nodes in z direction (nm)
e = deltaS_xy/3;           % smaller deltaS = more Ns (# of small spheres), more accurate and converging to actual #, e = radius of each node
mu = 0.01;                 % dynamic viscosity of FV (10 cP) in N*sec/m^2 = Pa*sec
d = 3;                     % moving in 3 dimensions

% defining x y z ranges for hemozoin brick in nm (170x170x585 nm (10^-9 m))
X = -85:deltaS_xy:85;     
Y = -85:deltaS_xy:85;
Z = -292.5:deltaS_z:292.5;  

% finding point locations for all vertices in grid
[xx yy zz] = meshgrid(X, Y, Z);
coords = [xx(:) yy(:) zz(:)];      % coords = list of every vertex in 3D grid

% making grid lines on x y z shape
[X1 Y1 Z1] = meshgrid(xx([1 end]),yy,zz);
X1 = permute(X1,[2 1 3]); 
Y1 = permute(Y1,[2 1 3]); 
Z1 = permute(Z1,[2 1 3]);
X1(end+1,:,:) = NaN; 
Y1(end+1,:,:) = NaN; 
Z1(end+1,:,:) = NaN;
[X2 Y2 Z2] = meshgrid(xx,yy([1 end]),zz);
X2(end+1,:,:) = NaN; 
Y2(end+1,:,:) = NaN; 
Z2(end+1,:,:) = NaN;
[X3 Y3 Z3] = meshgrid(xx,yy,zz([1 end]));
X3 = permute(X3,[3 1 2]); 
Y3 = permute(Y3,[3 1 2]); 
Z3 = permute(Z3,[3 1 2]);
X3(end+1,:,:) = NaN; 
Y3(end+1,:,:) = NaN; 
Z3(end+1,:,:) = NaN;

% coordinates yz face on x boundary
ind1 = coords(:,1) == -85;
yz_face1 = coords(ind1,:);

% coordinates yz face on x boundary
ind2 = coords(:,1) == 85;
yz_face2 = coords(ind2,:);

% coordinates xz face on y boundary
ind3 = coords(:,2) == -85;
xz_face1 = coords(ind3,:);

% coordinates xz face on y boundary
ind4 = coords(:,2) == 85;
xz_face2 = coords(ind4,:);

% coordinates xy face on z boundary
ind5 = coords(:,3) == -292.5; 
xy_face1 = coords(ind5,:);

% coordinates xy face on z boundary
ind6 = coords(:,3) == 292.5;  
xy_face2 = coords(ind6,:);

%%----------------------------------------------------------------------------------
% Visualization
% Plotting hemozoin
figure(1)

hold on
plot3(yz_face1(:,1),yz_face1(:,2),yz_face1(:,3), 'o')
plot3(yz_face2(:,1),yz_face2(:,2),yz_face2(:,3), 'o')
plot3(xz_face1(:,1),xz_face1(:,2),xz_face1(:,3), 'o')
plot3(xz_face2(:,1),xz_face2(:,2),xz_face2(:,3), 'o')
plot3(xy_face1(:,1),xy_face1(:,2),xy_face1(:,3), 'o')
plot3(xy_face2(:,1),xy_face2(:,2),xy_face2(:,3), 'o')

% figure('Renderer','opengl')
h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], [Z1(:);Z2(:);Z3(:)]);
set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-')

axis equal
box on
axis on
view(500,500)
grid on

title('Hemozoin Geometry')
xlabel('x (nm)') % x-axis label
ylabel('y (nm)') % y-axis label
zlabel('z (nm)') % z-axis label
pbaspect([1 1 1])
%-------------------------------------------------------------------------------

% concatenating all face points to make large list of points
p = cat(1,yz_face1,yz_face2,xz_face1,xz_face2,xy_face1,xy_face2);

% deleting repeats that overlap at each edge
p = unique(p, 'rows')*10^-9;  % points are all in m scale (converting nm to m)


%-----------------------------------------------------------------------------
% Visualization
% Confirming no overlaps

figure(2)

hold on
plot3(p(:,1),p(:,2),p(:,3), 'o')

axis equal
box on
axis on
view(500,500)
grid on

title('Hemozoin Geometry')
xlabel('x (nm)') % x-axis label
ylabel('y (nm)') % y-axis label
zlabel('z (nm)') % z-axis label
axis equal

% sets parameters for experiment 
n_hemozoin=35; % number of hemozoin crystals
dt = 0.0033;  % time step for each frame (in seconds)
t = 0:dt:5;   % time length of simulation (in seconds)
radius = 3e-6;    % FV radius (in meters)


% looping through each individual crystal to establish starting positions
for h = 1:n_hemozoin
  
%     initializing displacement at 0
    center_posx{h} = zeros(1);
    center_posy{h} = zeros(1);
    center_posz{h} = zeros(1);
    
%     randomizing starting positions
    center_posx{h}(1) = radius/2*(-1+(2).*rand);
    center_posy{h}(1) = radius/2*(-1+(2).*rand);
    center_posz{h}(1) = radius/2*(-1+(2).*rand);
    hemo_CM{1,h}(1,:) = [center_posx{h}(1) center_posy{h}(1) center_posz{h}(1)];
        
    hemo(h)=hemozoin_v11(e, d, mu, p, X,Y,Z) ;

%      starting node locations
    p_r1{1,h}(:,1) = p(:,1) + center_posx{h}(1,1);
    p_r1{1,h}(:,2) = p(:,2) + center_posy{h}(1,1);
    p_r1{1,h}(:,3) = p(:,3) + center_posz{h}(1,1);
    
%       distances from brick center of mass to node locations
    node_distance_CM{1,h}(:,:) = [p_r1{1,h}(:,1)-hemo_CM{1,h}(1,1) p_r1{1,h}(:,2)-hemo_CM{1,h}(1,2) p_r1{1,h}(:,3)-hemo_CM{1,h}(1,3)];
       
    position_matrix_rotated{1,h}(:,:) = p_r1{1,h}(:,:);

end


%   establishing starting positions of all crystals
for h=1:n_hemozoin

    
    posx{1,h}(1,:) = p_r1{1,h}(:,1);
    posy{1,h}(1,:) = p_r1{1,h}(:,2);
    posz{1,h}(1,:) = p_r1{1,h}(:,3);
    
    pos_tot{1,h}(:,:) = [posx{1,h}(1,:); posy{1,h}(1,:); posz{1,h}(1,:)];

    brick{1,h}(:,1) = p_r1{1,h}(:,1);
    brick{1,h}(:,2) = p_r1{1,h}(:,2);
    brick{1,h}(:,3) = p_r1{1,h}(:,3);
    
    tot_displacement{1,h}(:,:) = [mean(posx{1,h}(1,:)); mean(posy{1,h}(1,:)); mean(posz{1,h}(1,:))];
    tot_rotation{1,h}(:,:) = eye(3);
    
    distance{1,h}(:,1) = sqrt(posx{1,h}(1,:).^2+posy{1,h}(1,:).^2+posz{1,h}(1,:).^2);
    
    rotation_matrix_cell{1,h}(:,:) = eye(3);
    
end

%   looping through each time step and displacing crystals based on physical constraints

for i = 2:length(t)
  
  % computing Brownian displacements/rotations at each time step
  for h = 1:n_hemozoin  
    
    hemo(h)=hemo(h).rotation_translation(mu,p,  radius, dt,t);
    
    % saving rotation matrix and random velocity in x,y,z
    rotation_matrix_cell_dt{i,h}(:,:) = hemo(1,h).rotation_matrix_cell;
    
    rotation_matrix_cell{i,h}(:,:) = rotation_matrix_cell_dt{i,h}(:,:)*rotation_matrix_cell{i-1,h}(:,:);
    
    lab_vel_tot{i,h}(:,:) = rotation_matrix_cell{i,h}(:,:)*hemo(1,h).body_trans_vel_tot;
    pos_tot_brownian{i,h}(:,:) = lab_vel_tot{i,h}(:,:) * dt;
    
    posx{i,h}(1,:) = brick{i-1,h}(:,1) + (lab_vel_tot{i,h}(1,:)*dt);
    posy{i,h}(1,:) = brick{i-1,h}(:,2) + (lab_vel_tot{i,h}(2,:)*dt);
    posz{i,h}(1,:) = brick{i-1,h}(:,3) + (lab_vel_tot{i,h}(3,:)*dt);
            
    pos_tot{i,h}(:,:) = [posx{i,h}(1,:); posy{i,h}(1,:); posz{i,h}(1,:)];
    hemo_CM{i,h}(1,:) = [mean(posx{i,h}(1,:)) mean(posy{i,h}(1,:)) mean(posz{i,h}(1,:))];
   
    % saving rotation matrix used at each crystal time step for LJ potentials
    rotation_matrix_timestep{i,h}(:,:) = rotation_matrix_cell{i,h}(:,:);
  
    position_matrix_rotated{i,h}(:,:) =  rotation_matrix_cell{i,h}(:,:)*node_distance_CM{1,h}(:,:)' + hemo_CM{i,h}'; 
  
    brick_brownian{i,h}(:,1) = position_matrix_rotated{i,h}(1,:); 
    brick_brownian{i,h}(:,2) = position_matrix_rotated{i,h}(2,:); 
    brick_brownian{i,h}(:,3) = position_matrix_rotated{i,h}(3,:);
    
    tot_force{i,h}(:,:) = zeros(3,1);
    tot_torque{i,h}(:,:) = zeros(3,1);
  
  
  end
  
  
  
  % computing LJ displacements/rotations at each time step
  for h1=1:n_hemozoin-1
        for h2=h1+1:n_hemozoin % looping through all crystals

            % find center point of brick 1
            centroid_brick{i,h1}(1,:) = [mean(brick_brownian{i,h1}(:,1)) mean(brick_brownian{i,h1}(:,2)) mean(brick_brownian{i,h1}(:,3))];
                
            % find center point of brick 2
            centroid_brick{i,h2}(1,:) = [mean(brick_brownian{i,h2}(:,1)) mean(brick_brownian{i,h2}(:,2)) mean(brick_brownian{i,h2}(:,3))];
            
            crystal1_f = zeros(size(pos_tot{1,h1},2),size(pos_tot{1,h2},2),3);
            crystal2_f = zeros(size(pos_tot{1,h1},2),size(pos_tot{1,h2},2),3);
            crystal1_t = zeros(size(pos_tot{1,h1},2),size(pos_tot{1,h2},2),3);
            crystal2_t = zeros(size(pos_tot{1,h1},2),size(pos_tot{1,h2},2),3); 
            
            % looping through all nodes on both spheres
            for n1=1:size(pos_tot{1,h1},2)
                for n2=1:size(pos_tot{1,h2},2) 

                node_distance(n1,n2) = sqrt((brick_brownian{i,h1}(n1,1)-brick_brownian{i,h2}(n2,1)).^2+(brick_brownian{i,h1}(n1,2)-brick_brownian{i,h2}(n2,2)).^2+(brick_brownian{i,h1}(n1,3)-brick_brownian{i,h2}(n2,3)).^2);
            
                if node_distance(n1,n2) < 5e-7 %% if crystals come within 500 nm of each other
                    
                distance_vector = [brick_brownian{i,h1}(n1,1)-brick_brownian{i,h2}(n2,1) brick_brownian{i,h1}(n1,2)-brick_brownian{i,h2}(n2,2) brick_brownian{i,h1}(n1,3)-brick_brownian{i,h2}(n2,3)];
                
                % describing the unit vector direction that the force will go in
                n_hat(1,1) = distance_vector(1,1)/node_distance(n1,n2);
                n_hat(1,2) = distance_vector(1,2)/node_distance(n1,n2);
                n_hat(1,3) = distance_vector(1,3)/node_distance(n1,n2);
                
                % defining forces and torques on each sphere (Lennard-Jones potential)
                r = node_distance(n1,n2);
                % using radius of each node (e) in meters
                sigma = e*2*10^(-9);
                % eps = epsilon from lennard jones potential (Joules)
                eps = 2*10^-21;

                
                % forces in all three directions, sphere on crystal 1
                crystal1_f(n1,n2,1) = 48*eps*((sigma/r)^13)*n_hat(1,1); 
                crystal1_f(n1,n2,2) = 48*eps*((sigma/r)^13)*n_hat(1,2); 
                crystal1_f(n1,n2,3) = 48*eps*((sigma/r)^13)*n_hat(1,3); 
                
                % making sure LJ potentials don't kick crystals out of
                % boundary, moving back ~500 nm
                if abs(crystal1_f(n1,n2,1)) > 5e-16
                    crystal1_f(n1,n2,1) = 5e-16*n_hat(1,1);
                end 
            
                if abs(crystal1_f(n1,n2,2)) > 5e-16
                    crystal1_f(n1,n2,2) = 5e-16*n_hat(1,2);
                end
            
                if abs(crystal1_f(n1,n2,3)) > 5e-16
                    crystal1_f(n1,n2,3) = 5e-16*n_hat(1,3);
                end
               
                % opposite forces in all three directions, sphere on crystal 2
                crystal2_f(n1,n2,1) = -1*(crystal1_f(n1,n2,1));
                crystal2_f(n1,n2,2) = -1*(crystal1_f(n1,n2,2));
                crystal2_f(n1,n2,3) = -1*(crystal1_f(n1,n2,3));

                % torques in all three directions, sphere on crystal 1
                distance_cm_vector1 = [brick_brownian{i,h1}(n1,1)-centroid_brick{i,h1}(1,1) brick_brownian{i,h1}(n1,2)-centroid_brick{i,h1}(1,2) brick_brownian{i,h1}(n1,3)-centroid_brick{i,h1}(1,3)];

                crystal1_t(n1,n2,1) = distance_cm_vector1(1,2)*crystal1_f(n1,n2,3)-distance_cm_vector1(1,3)*crystal1_f(n1,n2,2);
                crystal1_t(n1,n2,2) = distance_cm_vector1(1,3)*crystal1_f(n1,n2,1)-distance_cm_vector1(1,1)*crystal1_f(n1,n2,3);
                crystal1_t(n1,n2,3) = distance_cm_vector1(1,1)*crystal1_f(n1,n2,2)-distance_cm_vector1(1,2)*crystal1_f(n1,n2,1);
       
                % torques in all three directions, sphere on crystal 2
                distance_cm_vector2 = [brick_brownian{i,h2}(n2,1)-centroid_brick{i,h2}(1,1) brick_brownian{i,h2}(n2,2)-centroid_brick{i,h2}(1,2) brick_brownian{i,h2}(n2,3)-centroid_brick{i,h2}(1,3)];
                
                crystal2_t(n1,n2,1) = distance_cm_vector2(1,2)*crystal2_f(n1,n2,3)-distance_cm_vector2(1,3)*crystal2_f(n1,n2,2); 
                crystal2_t(n1,n2,2) = distance_cm_vector2(1,3)*crystal2_f(n1,n2,1)-distance_cm_vector2(1,1)*crystal2_f(n1,n2,3);
                crystal2_t(n1,n2,3) = distance_cm_vector2(1,1)*crystal2_f(n1,n2,2)-distance_cm_vector2(1,2)*crystal2_f(n1,n2,1);
               
                % this will end with 74x74x3 matrix of individual values for forces/torques on each specific sphere in x,y,z direction
            
            
                end
                end
            end
        
        % total force/torque in vacuole basis on each crystal
        tot_brick1_force = sum(crystal1_f,[1 2]);
        tot_brick2_force = sum(crystal2_f,[1 2]);
        tot_brick1_torque = sum(crystal1_t,[1 2]);
        tot_brick2_torque = sum(crystal2_t,[1 2]);
                
        % combining forces/torques to multiply by inverse rotation matrix to convert back to crystal basis  
        tot_brick1_force = rotation_matrix_timestep{i,h1}(:,:)'*[tot_brick1_force(1,1,1); tot_brick1_force(1,1,2); tot_brick1_force(1,1,3)];
        tot_brick1_torque = rotation_matrix_timestep{i,h1}(:,:)'*[tot_brick1_torque(1,1,1); tot_brick1_torque(1,1,2); tot_brick1_torque(1,1,3)];
                
        tot_brick2_force = rotation_matrix_timestep{i,h2}(:,:)'*[tot_brick2_force(1,1,1); tot_brick2_force(1,1,2); tot_brick2_force(1,1,3)];
        tot_brick2_torque = rotation_matrix_timestep{i,h2}(:,:)'*[tot_brick2_torque(1,1,1); tot_brick2_torque(1,1,2); tot_brick2_torque(1,1,3)];
        
        tot_force{i,h1}(:,1) = tot_force{i,h1}(:,:) + tot_brick1_force;
        tot_force{i,h2}(:,1) = tot_force{i,h2}(:,:) + tot_brick2_force;
        
        tot_torque{i,h1}(:,1) = tot_torque{i,h1}(:,:) + tot_brick1_torque;
        tot_torque{i,h2}(:,1) = tot_torque{i,h2}(:,:) + tot_brick2_torque;
        end
  end   
  
  % converting total forces/torques of each brick into LJ displacements
  for h = 1:n_hemozoin
      
      tot_brick_force_mem{i,h}(:,:) = zeros(1,3);
      tot_brick_torque_mem{i,h}(:,:) = zeros(1,3);
      
        % keeping crystals inside membrane
        for n = 1:size(pos_tot{1,h},2)
        % distance of node to origin
        rad_node = sqrt(brick_brownian{i,h}(n,1).^2+brick_brownian{i,h}(n,2).^2+brick_brownian{i,h}(n,3).^2);
                
        if abs(radius - rad_node) < 5e-7 %% if crystals come within 500 nm of membrane
                    
            distance_vector_origin = [brick_brownian{i,h}(n,1) brick_brownian{i,h}(n,2) brick_brownian{i,h}(n,3)];
                    
            % describing the unit vector direction that the force will go in
            r_hat(1,1) = distance_vector_origin(1,1)/rad_node;
            r_hat(1,2) = distance_vector_origin(1,2)/rad_node;
            r_hat(1,3) = distance_vector_origin(1,3)/rad_node;
                    
            % defining forces and torques on each sphere (Lennard-Jones potential)
            distance_to_mem{i,h}(n,1) = abs(radius - rad_node);
            % using radius of each node (e) in meters
            sigma = e*2*10^(-9);
            % eps = epsilon from lennard jones potential (Joules)
            eps = 2*10^-21;
                
            % forces in all three directions, node on crystal
            crystal_f_mem{i,h}(n,1) = 48*eps*((sigma/distance_to_mem{i,h}(n,1))^13)*-r_hat(1,1); 
            crystal_f_mem{i,h}(n,2) = 48*eps*((sigma/distance_to_mem{i,h}(n,1))^13)*-r_hat(1,2); 
            crystal_f_mem{i,h}(n,3) = 48*eps*((sigma/distance_to_mem{i,h}(n,1))^13)*-r_hat(1,3); 
            
            % making sure LJ potentials don't kick crystals out of
            % boundary, moving back 500 nm
            if abs(crystal_f_mem{i,h}(n,1)) > 5e-16
                crystal_f_mem{i,h}(n,1) = 5e-16*-r_hat(1,1);
            end 
            
            if abs(crystal_f_mem{i,h}(n,2)) > 5e-16
                crystal_f_mem{i,h}(n,2) = 5e-16*-r_hat(1,2);
            end
            
            if abs(crystal_f_mem{i,h}(n,3)) > 5e-16
                crystal_f_mem{i,h}(n,3) = 5e-16*-r_hat(1,3);
            end
            
            % torques in all three directions, node on crystal
            distance_cm_vector = [brick_brownian{i,h}(n,1)-centroid_brick{i,h}(1,1) brick_brownian{i,h}(n,2)-centroid_brick{i,h}(1,2) brick_brownian{i,h}(n,3)-centroid_brick{i,h}(1,3)];

            crystal_t_mem{i,h}(n,1) = distance_cm_vector(1,2)*crystal_f_mem{i,h}(n,3)-distance_cm_vector(1,3)*crystal_f_mem{i,h}(n,2);
            crystal_t_mem{i,h}(n,2) = distance_cm_vector(1,3)*crystal_f_mem{i,h}(n,1)-distance_cm_vector(1,1)*crystal_f_mem{i,h}(n,3);
            crystal_t_mem{i,h}(n,3) = distance_cm_vector(1,1)*crystal_f_mem{i,h}(n,2)-distance_cm_vector(1,2)*crystal_f_mem{i,h}(n,1);
                 
            tot_brick_force_mem{i,h} = sum(crystal_f_mem{i,h},1);
            tot_brick_torque_mem{i,h} = sum(crystal_t_mem{i,h},1);  
            
            % combining forces/torques to multiply by inverse rotation matrix to convert back to crystal basis  

            tot_brick_force_mem2{i,h}(:,:) = rotation_matrix_timestep{i,h}(:,:)'*[tot_brick_force_mem{i,h}(1,1); tot_brick_force_mem{i,h}(1,2); tot_brick_force_mem{i,h}(1,3)];
            tot_brick_torque_mem2{i,h}(:,:) = rotation_matrix_timestep{i,h}(:,:)'*[tot_brick_torque_mem{i,h}(1,1); tot_brick_torque_mem{i,h}(1,2); tot_brick_torque_mem{i,h}(1,3)];
          
            
            tot_force{i,h}(:,1) = tot_force{i,h}(:,:) + tot_brick_force_mem2{i,h}(:,:);
            tot_torque{i,h}(:,1) = tot_torque{i,h}(:,:) + tot_brick_force_mem2{i,h}(:,:);
            
            
        end
        end
     
        % MM is defined for every brick * sum of force/torque on all spheres on brick * rotation matrix at that time point 
        repul_vel_brick(1,1) = hemo(1,h).MM(1,1)*tot_force{i,h}(1,1);
        repul_vel_brick(2,1) = hemo(1,h).MM(2,2)*tot_force{i,h}(2,1);
        repul_vel_brick(3,1) = hemo(1,h).MM(3,3)*tot_force{i,h}(3,1);
            
        repul_rot_brick(1,1) = hemo(1,h).MM(4,4)*tot_torque{i,h}(1,1);
        repul_rot_brick(2,1) = hemo(1,h).MM(5,5)*tot_torque{i,h}(2,1);
        repul_rot_brick(3,1) = hemo(1,h).MM(6,6)*tot_torque{i,h}(3,1); 
       
        % converting from body basis into vacuole basis with rotations (brick 1)
        w1 = repul_rot_brick';                 % making vector of omega in order to calculate magnitude

    
    
        omegas = (1/norm(w1)^2)* [w1(1,1)^2,  w1(1,1)* w1(1,2),  w1(1,1)* w1(1,3);... % calculating sub components of R matrix 
                    w1(1,2)* w1(1,1),  w1(1,2)^2,  w1(1,2)* w1(1,3);...
                    w1(1,3)* w1(1,1),  w1(1,3)* w1(1,2),  w1(1,3)^2];
   
        omegax = (1/norm(w1))*[0, - w1(1,3),  w1(1,2);...              % making omega cross (wx)
                    w1(1,3), 0, - w1(1,1);...
                    - w1(1,2),  w1(1,1), 0];
    
        Rdt = omegas+cos(norm(w1)*dt)*(eye(3)-omegas)+sin(norm(w1)*dt)*omegax;    % infinitesimal rotation based on dt
        
        if sum(Rdt,'all') > 0
            tot_rotation{i,h} = rotation_matrix_cell{i,h}(:,:)*Rdt;
        else 
            tot_rotation{i,h} = rotation_matrix_cell{i,h}(:,:);
        end
            

     
        lj_dis_tot_brick{i,h}(:,:) = tot_rotation{i,h} * repul_vel_brick(:,1) * dt; 
        
  end
        
% updating brick to include LJ, Brownian positions
    for h = 1:n_hemozoin
        
        tot_displacement{i,h}(:,:) = tot_displacement{i-1,h}(:,:) + lj_dis_tot_brick{i,h}(:,:) + pos_tot_brownian{i,h}(:,:);
               
         brick_updated{i,h}(:,:) = tot_rotation{i,h}(:,:)*node_distance_CM{1,h}(:,:)' + tot_displacement{i,h};


         brick{i,h}(:,1) = brick_updated{i,h}(1,:);
         brick{i,h}(:,2) = brick_updated{i,h}(2,:);
         brick{i,h}(:,3) = brick_updated{i,h}(3,:);
    
    end
end   

% writes videos
v = VideoWriter('Hemozoin_Movement.avi');
v.FrameRate = round((1/dt));
open(v);


axis on
view(500,500)
grid on

%  ----------------------visualization    

for i = 1:length(t) 
    hold off
    figure(4)
   
    for h=1:n_hemozoin
        drawnow        
        hemo(h).bd{i,1}(:,:) = boundary(brick{i,h}(:,:));
        trisurf(hemo(h).bd{i,1}(:,:),brick{i,h}(:,1),brick{i,h}(:,2),brick{i,h}(:,3),'FaceAlpha',0.1)
        hold on
   
    xlim([-3.5e-6 3.5e-6]) 
    ylim([-3.5e-6 3.5e-6])
    zlim([-3.5e-6 3.5e-6])
       % Label axes.
 xlabel('X', 'FontSize', 20);
 ylabel('Y', 'FontSize', 20);
 zlabel('Z', 'FontSize', 20);
 axis square
    end

    
% Make unit sphere
[x,y,z] = sphere;
% Scale to desired radius.
x = x * radius;
y = y * radius;
z = z * radius;
% Translate sphere to new location.
offset = 0;
% Plot as surface.
surf(x+offset,y+offset,z+offset,'FaceAlpha',0.05) 

    
    frame = getframe(gcf);
    writeVideo(v,frame);

end

close(v);



% creating tracks from brick positions (CM)
for h = 1:n_hemozoin
    for i = 1:length(t)
        brick_CM{i,h}(:,:) = [mean(brick{i,h}(:,1)) mean(brick{i,h}(:,2)) mean(brick{i,h}(:,3))];
        tracks{h}(i,1:3) = brick_CM{i,h}(1,1:3);
    end
    
    [msd,tau] = MSD1(tracks{h}(:,:), dt); 
    msd1{1,h}(:,2) = msd;
    msd1{1,h}(:,1) = tau;
    
end


figure(5) 
hold on

for h = 1:n_hemozoin
    plot(msd1{1,h}(:,1),msd1{1,h}(:,2)) 
    dlm = fitlm(msd1{1,h}(1:55,1),msd1{1,h}(1:55,2),'Intercept',false); 
    x = msd1{1,h}(1:55,1); 
    linear_coefficients(h,:) = table2array(dlm.Coefficients);
    d_app(h,1) = linear_coefficients(h,1)/1e-12/4;
end

xlabel('s [s]') 
axis([0 0.3 0 5e-13])
ylabel('MSD [m^2]') 
hold off
