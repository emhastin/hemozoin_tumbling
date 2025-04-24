%simulation of stepsize distribution for crystal with propulsion
%x is crystal long axis, y and z are crystal short axes
clear all

%below is to use for no distribution case
%onevel=70;
%propulsiondirection = [0; 1; 0]; %unit vector of propulsion direction in body frame

mu= 1e-3*(1/10^6); %absolute viscosity, kg/(um s)
kBT = 293*1.381*10^-23*10^12 %boltzmann constant * 293 K, kg um^2/s^2 


%below values can be used to simulate crystal geometry instead of sphere
%Dx =kBT*2.9977e8; %translational diffusion of crystal along long axis, um^2/s
       % using value of 2.9977e8 s/kg from MM.zip, which
%Dy= kBT*2.5676e8; %y and z are diffusion constants along short axes
%Dz= kBT*2.5676e8; % using value of 2.5676e8 s/kg from MM.zip


activeD=5.8; %translational diffusion constant from active portion of bimodal fit


%%For spherical model
%%
D = 1.08; %translational diffusion constant (measured for hemozoin in water), um^2/s
%D= 0.17 %translational diffusion constant (from theory), um^2/s
Dx=D;
Dy=D;
Dz=D;
%%
%%

avgD = (Dx+Dy+Dz)/3;
a = kBT/avgD/6/pi/mu; %sphere radius to produce translational diffusion constant above, um


fprintf('sphere radius in um: %d \n', a);

%below values can be used to simulate crystal geometry instead of sphere
%Drx = kBT*8.9930e9 ; %rotational diffusion constant about long axis, rad^2/s
             % using value of 8.9930e21 s/(kg m^2)= 8.9939e9 s/(kg um^2) from MM.zip
%Dry = kBT*2.9975e9; %rotational diffusion constant about short axes, rad^2/s
%Drz = kBT*2.9975e9; %using value of 2.9975e21 s/(kg m^2)=2.9975e9 s/(kg um^2)  from MM.zip

%% For spherical model
%%
Dr = kBT/(8*pi*mu*a^3); %rotational diffusion constant, rad^2/s
Drx=Dr;
Dry=Dr;
Drz=Dr;
%%
%%

avgDr = (Drx+Dry+Drz)/3;

tr = (pi/100)^2/avgDr ; %timescale for rotational diffusion (time to diffuse pi/100 radians), s
dt =  1e-5;%timestep, make sure it is less than rotational diffusion timescale, s
fprintf('timestep (%d) should be less than rotational diffusion timescale for rotation by pi/100 radians (%d) \n', dt, tr);


%from paper: trajectories up to 0.5 s (~100 frames)
numbertraj = 620; %3500 trajectories in expt
sample =  round(1/300/dt); %number of dt in 0.00333 s (300fps) 
sampleinterval = dt*sample; 
max_timesteps = sample * 100;

trajDisplacements = zeros(numbertraj, 2, max_timesteps/sample);
trajdata = zeros(numbertraj,max_timesteps/sample);
veldist = zeros(1,numbertraj);
omegadist = zeros(numbertraj, 3);
propulsiondirection = zeros(numbertraj,3);
%generate distribution of propulsion velocities and directions
%           veldist is for speed
%           propulsiondirection is for body frame direction of propulsion
%           omegadist is for body frame active rotation
veldiststdev= 60;  %standard deviation of each component of propulsion velocity distribution
                  %units um/s
                  %note stdev o velocity is therefore sqrt(3) veldiststdev
omegadiststdev = 2*pi*200/4; %standard deviation of each component of propulsion angular velocity
                     %units rad/s
                     

for trajcounter = 1:numbertraj     
%   These are for the velocity distribution
    vel = normrnd(0,veldiststdev,[1,3]);
    veldist(trajcounter) = norm(vel);
    propulsiondirection(trajcounter,:) = vel/norm(vel); %unit vector of propulsion direction in body frame
    omegadist(trajcounter,:) = normrnd(0,omegadiststdev,[1,3]);

%   This is for no distribution
%   vel= rand([1 3])-0.5;
%   veldist(trajcounter) = onevel;
%   propulsiondirection(trajcounter,:) = vel/norm(vel); %unit vector of propulsion direction in body frame
end

for trajcounter= 1:numbertraj
    v0= veldist(trajcounter);
    heading = zeros(3,max_timesteps); %direction of propulsion
    position = zeros(3,max_timesteps);
    position(:,1) = [0; 0; 0];
    heading(:,1) = propulsiondirection(trajcounter,:); %velocity direction in body frame
    twoDposition = zeros(2,max_timesteps/sample);
    rot=[1 0 0; 0 1 0; 0 0 1]; %rotation matrix that takes body frame to lab frame

    %rotation due to propulsion angular velocity in body frame for each
    %time step
        dthetax = dt*omegadist(trajcounter,1);
        dthetay = dt*omegadist(trajcounter,2);
        dthetaz = dt*omegadist(trajcounter,3);
        phi = sqrt(dthetax^2 + dthetay^2 + dthetaz^2);%rotation angle
        nhatx = dthetax/phi; %rotation axis
        nhaty = dthetay/phi;
        nhatz = dthetaz/phi;
        proprot = (1- cos(phi))*[nhatx*nhatx nhatx*nhaty nhatx*nhatz; nhaty*nhatx nhaty*nhaty nhaty*nhatz; nhatz*nhatx nhatz*nhaty nhatz*nhatz];
        proprot = proprot + cos(phi)*[1 0 0;0 1 0;0 0 1];
        proprot = proprot + sin(phi)*[0 -nhatz nhaty; nhatz 0 -nhatx; -nhaty nhatx 0];     


    for tn = 1:max_timesteps
        
        dposx = normrnd(0, sqrt(dt*2*Dx) );
        dposy = normrnd(0, sqrt(dt*2*Dy) );
        dposz = normrnd(0, sqrt(dt*2*Dz) );
%        position(:,tn+1) = position(:,tn) + rot*[dposx;dposy;dposz] + heading(:,tn)*v0*dt; 
        position(:,tn+1) = position(:,tn) + [dposx;dposy;dposz] + heading(:,tn)*v0*dt;     

        rot = rot*proprot; %update rotation matrix for rotation drot in body frame
           
        %Brownian rotational diffusion
        dthetax = normrnd(0, sqrt(dt*2*Drx) );
        dthetay = normrnd(0, sqrt(dt*2*Dry) );
        dthetaz = normrnd(0, sqrt(dt*2*Drz) );
        phi = sqrt(dthetax^2 + dthetay^2 + dthetaz^2);%rotation angle
        nhatx = dthetax/phi; %rotation axis
        nhaty = dthetay/phi;
        nhatz = dthetaz/phi;
        drot = (1- cos(phi))*[nhatx*nhatx nhatx*nhaty nhatx*nhatz; nhaty*nhatx nhaty*nhaty nhaty*nhatz; nhatz*nhatx nhatz*nhaty nhatz*nhatz];
        drot = drot + cos(phi)*[1 0 0;0 1 0;0 0 1];
        drot = drot + sin(phi)*[0 -nhatz nhaty; nhatz 0 -nhatx; -nhaty nhatx 0];
        %heading(:,tn+1) = drot*heading(:,tn) ; %rotational diffusion
        rot = rot*drot; %update rotation matrix for rotation drot in body frame

        heading(:,tn+1) = rot*heading(:,1);


        if mod(tn,sample) == 0
            twoDposition(:, round(tn/sample)+1) = position(1:2, tn);

        end
               
    end
        
    


    stepsize= zeros(1,max_timesteps/sample);
    nstepsize= zeros(1,max_timesteps/sample);
    for step=1:(max_timesteps/sample)
        stepsize(step) = norm( twoDposition(:,step+1) - twoDposition(:,step));
        nstepsize(step) = stepsize(step)/sqrt(sampleinterval);
  
        trajDisplacements(trajcounter,:,step) = twoDposition(:,step+1);
    end
    trajdata(trajcounter,:) = nstepsize;
end

%calculate average mean square displacement for each time lag
samplenumber = max_timesteps/sample;
trajMSD = zeros(numbertraj,samplenumber-1);
MSD = zeros(1,samplenumber-1);
for trajcounter = 1:numbertraj
    for timelag = 1:(samplenumber-1)
        for j=1:(samplenumber-timelag)
            trajMSD(trajcounter,timelag) = trajMSD(trajcounter,timelag) + (norm( trajDisplacements(trajcounter,:,j+timelag) - trajDisplacements(trajcounter,:,j) ))^2;
        end
        trajMSD(trajcounter,timelag) = trajMSD(trajcounter,timelag)/(samplenumber-timelag);
    end
    MSD =  MSD + trajMSD(trajcounter,:);
end
MSD = MSD/numbertraj;

%figure(1)
%scatter( position(1,:), position(2,:));axis equal
    

plotsteps = trajdata(:);

figure(1)
histogram(plotsteps, 'Normalization', 'pdf') 
xlim([0 10])
ylim([0 0.6])
hold on



ndx= 0:0.05:10;
dx=ndx*sqrt(sampleinterval);
sigma = 2*activeD*sampleinterval;
analyticalpdf = 1/sigma* dx.*exp(- (dx.^2)/(2*sigma)) ;
nanalyticalpdf = analyticalpdf*sqrt(sampleinterval);
plot(ndx,nanalyticalpdf)
hold off

%figure(2) 
%histogram(veldist)



figure(2)
timelags= sampleinterval*( 1:(samplenumber-1) );
plot(timelags, MSD)

numbertraj58=numbertraj;
timelags58 = timelags;
MSD58=MSD;
trajMSD58 = trajMSD;
plotsteps58= plotsteps;
save("data58.mat","timelags58","MSD58","trajMSD58","plotsteps58","D","veldiststdev","omegadiststdev","numbertraj58");