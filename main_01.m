% --SYM_1--
%Sym_1 is the main executable file for IAC2023 simulation. The simulation
%consists of a femtosatellite swarm in low earth orbit, which is collecting
%data from the world magnetic model. A perturbation is applied to the world
%magentic model, which affects the magnetic field in a scale of minutes,
%seconds, and microseconds. At the end of the simulation, the data is
%processed to return intensity maps over time of a single cell, intensity
%maps over space for multiple cells, and a field reconstruction of a cell
%over time. 

% Initialize
%Clean up 
clear 
clc
addpath("RungeKuttar4thOrderIntegrator\")
tic
%Time parameter definition-------------------------------------------------
n_orbits = 5; %Choose the number of orbits we want to analyse
t_0 = 1672531200; %mission start epoch [s] %NOTE: 00:00 1 of January, 2023
t_f = t_0+n_orbits*5.43010002e3; %mission end epoch [s] %NOTE: 00:00 2 OF January, 2023 ; 5.43010002e3 is the orbital period
dt = 0.5; %stepsize [s]
%Initial orbit design------------------------------------------------------
a = 7400; %Semi major axis [km]
e = 0.0000001; %Eccentricity
i = 10; %Inclination [rads]
Omega = 0; %Right ascention of the ascending node [rad]
w = 0; %Angle of perigee [rad]
M = 0; %Mean anomaly [rad]
Kepler = [a,e,i,Omega,w,M]; %Keplerian vector 
%Swarm definition----------------------------------------------------------
ID = 1; %identification digit
type = "MagMeasure"; %Swarm Type
%Drone manufacturing instructions
dType = "Drone_m";
dPop = 163; 
dNeed = 50; 
systems = table(dType,dPop,dNeed);
%Planetary Body parameters-------------------------------------------------
PRadius = 6371; %[km] planetary radius
rotRate = [0 0 2*pi/86400]; %[rad/s] Rotation rate relative to ECI
R_p = [0 0 0]; 
V_p = [0 0 0];
Theta_p = [0 0 0]; %NOTE: must calculate for the starting epoch
% Instantiation 
myClock = tiktok(t_0,t_f,dt); 
mySwarm = swarm(ID,type,systems,Kepler,t_0);
myEarth = celeBod(R_p,V_p,Theta_p,rotRate,PRadius,myClock); 
%Clean up
clear a dNeed dPop dt dType e i ID Kepler M Omega PRadius R_p rotRate systems t_0 t_f Theta_p type V_p w
disp("instantiation complete")
%__________________________________________________________________________
% Orbit propagation
deployBrood(mySwarm); %Apply deployment delta vee
for i=1:length(mySwarm.network{1})
    propagateOrb(mySwarm.network{1}{i},myClock);
end
disp("orbit propagation complete")
%__________________________________________________________________________
% Data collection 
%No error, yes perturbation
tic
for i=1:length(mySwarm.network{1})
    getMdata(mySwarm.network{1}{i},myEarth,myClock)
end
disp("data collection complete")
toc
%__________________________________________________________________________
%%  FIGURE 1: SWARM orbits
figure(1)
plotPath = 0; 
if plotPath == 1 
    figure(1)
    hold on
    for i=1:length(mySwarm.network{1})
        plot3(mySwarm.network{1}{i}.r(1,:),mySwarm.network{1}{i}.r(2,:),mySwarm.network{1}{i}.r(3,:))
    end
    axis equal
    grid on
    pbaspect([1 1 1])
    plot3(0,0,0,'ro')
    view([1 1 1])
    xlabel("x_{ECI}")
    ylabel("y_{ECI}")
    zlabel("z_{ECI}")
    hold off
    disp("orbit plotting complete")
end

%% FIGURE 2: SINGLE SPACECRAFT POV
%get a vector for the latitude of the satellite realative to Earth
% r_s = mySwarm.network{1}{1}.r(1:2,:);
% theta_s = zeros(length(myClock.timeline(:,1)),1);
% for i=1:length(r_s)-1
%     theta_s(i) = acos(dot(r_s(:,i),[1 0 ])/norm(r_s));
% end
% theta_se = theta_s-myEarth.Theta_t(:,3);
% disp(theta_se(1,:))
% for i = 1:length(theta_se)
%     if theta_se(i) > pi/2
%         theta_se(i) = 2*pi - theta_se(i);
%     end
% end

figure(2)
hold on
ylabel("Magnetic field intensity [nT]")
xlabel("Time since deployment [h]")
%xlim([0 ])
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{1}.D_m(:,3),'b')
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{1}.D_m(:,2),'g')
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{1}.D_m(:,1),'r')
legend('Z_{ECI}','Y_{ECI}','X_{ECI}')
% set(gca,'looseInset',get(gca,'TightInset'))
% saveas(gca,'Results/singleSatFullDataSet.jpeg')
hold off

%% FIGURE 3: CONTINOUS MONITORING DEMO 
% 1. Choose a grid size in m and create space
mySTS = stereoTaxicSpace(100000); %100 km lattice
% 2. Identify a populated cell at 1/4 period to observe a location with
% some spread
R = mySwarm.network{1}{1}.r(:,2700);
index = getIndex(R,1,1,1,mySTS);       
% 3. Monitor cell
[~,F] = monitorCel_F(myClock,mySwarm.network{1},mySTS,index);
disp("monitoring complete")

% Single voxel value over time / fixed space value over time
figure(3)
hold on 
grid on
T = myClock.timeline(:,1)/3600;
plot(T,F)
xlabel("Time since deployment [h]")
ylabel("Magnetic field intensity [nT]")
hold off

%% FIGURE 4: Distribued measurement for fixed time T = 25000 (MAGNETIC FIELD INTENSITY ONLY)
%1. Choose the time 
T = 25000; %first column is used for indexation reasons
%2. Sort the data [location, magnitude]
G_mi = zeros(length(mySwarm.network{1}),4);

for i=1:length(mySwarm.network{1})
    G_mi(i,1) = mySwarm.network{1}{i}.r(1,T);
    G_mi(i,2) = mySwarm.network{1}{i}.r(2,T);
    G_mi(i,3) = mySwarm.network{1}{i}.r(3,T);
    G_mi(i,4) = mySwarm.network{1}{i}.D_m(T,4);
end

figure(4)
scatter3(G_mi(:,1),G_mi(:,2),G_mi(:,3),10,G_mi(:,4),'filled')
colorbar
xlabel('X_{ECI} [m]')
ylabel('Y_{ECI} [m]')
zlabel('Z_{ECI} [m]')
%view([0 0 1])

%% FIGURE 5: Stereotaxic interpretation of Figure 4, for a lattice side of 10km. 
%Get the position at which each measurement was made
R = zeros(3,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    R(:,i) = mySwarm.network{1}{i}.r(:,25000);
end
%Get the magnetic field intensity for each measurement
F = zeros(1,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    F(1,i) = mySwarm.network{1}{i}.D_m(25000,4);
end
%Get index location for each measurement
mySTS2 = stereoTaxicSpace(10000);
index = zeros(3,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    index(:,i) = getIndex(R(:,i),1,1,1,mySTS2);
end
%Create cell array analogous to the final heat map. Each cell in the array
%is a cell in the stereotaxic space, and only cells within the maximum and
%minimum indexed positions are aknowledged (to avoid mapping all of space).
Map = cell(max(index(1,:))-min(index(1,:)),max(index(2,:))-min(index(2,:)));
%Assign all measurements made to a cell corresoponding its indexation
for i=1:length(index)
    if index(1,i)-min(index(1,:)) == 0 && max(index(2,:))-index(2,i) > 0
        Map{index(1,i)-min(index(1,:))+1,max(index(2,:))-index(2,i)} = horzcat(Map{index(1,i)-min(index(1,:))+1,max(index(2,:))-index(2,i)},F(i)); 
    end
    if index(1,i)-min(index(1,:)) > 0 && max(index(2,:))-index(2,i) == 0
        Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1} = horzcat(Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1},F(i));
    end
    if index(1,i)-min(index(1,:)) > 0 && max(index(2,:))-index(2,i) > 0
        Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)} = horzcat(Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)},F(i));
    end %EDIT WARNING: Although previously stable, an error started to occur on the second index of the Map variable. To correct, a +1 was removed. 
end
% Homogenize map and convert to heatmap (i.e., apply an averaging or
% interpolation rule to each cell to resolve each cell into a single value
for i=1:width(Map)
    for j = 1:height(Map)
        if Map{j,i} ~= 0
           Map{j,i} = sum(Map{j,i})/length(Map{j,i});
        else
           Map{j,i} = 0;
        end
    end
end
%Correct orientation
Map = flipud(fliplr(Map).');
%Convert cell array into a matrix array and generate heatmap
Map = cell2mat(Map);
figure(5)
heatmap(Map)
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


clear ans F G_mi i index j Map Ax n_orbits R T r_s theta_s theta_se 
runtime = toc;
%
load handel
sound(y,Fs)

clear y Fs 


