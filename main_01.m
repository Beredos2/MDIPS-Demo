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
%Time parameter definition-------------------------------------------------
n_orbits = 5; %Choose the number of orbits we want to analyse
t_0 = 1672531200; %mission start epoch [s] %NOTE: 00:00 1 of January, 2023
t_f = t_0+n_orbits*5.43010002e3; %mission end epoch [s] %NOTE: 00:00 2 OF January, 2023 ; 5.43010002e3 is the orbital period
dt = 0.5; %stepsize [s]
%Initial orbit design------------------------------------------------------
a = 6678; %Semi major axis [km]
e = 0.0000001; %Eccentricity
i = 0; %Inclination [rads]
Omega = 0; %Right ascention of the ascending node [rad]
w = 0; %Angle of perigee [rad]
M = 0; %Mean anomaly [rad]
Kepler = [a,e,i,Omega,w,M]; %Keplerian vector 
%Swarm definition----------------------------------------------------------
ID = 1; %identification digit
type = "MagMeasure"; %Swarm Type
%Drone manufacturing instructions
dType = "Drone_m";
dPop = 50; 
dNeed = 10; 
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
% Orbit propagation
deployBrood(mySwarm); %Modify v_0 for all drones
for i=1:length(mySwarm.network{1})
    propagateOrb(mySwarm.network{1}{i},myClock);
end
disp("orbit propagation complete")
% Data collection 
%No error, yes perturbation
for i=1:length(mySwarm.network{1})
    getMdata(mySwarm.network{1}{i},myEarth,myClock)
end
disp("data collection complete")
%Position knolwedge error

%Measuremetn knowledge error


% %%  Plot Orbits
% figure(2)
% hold on
% for i=1:length(mySwarm.network{1})
%     plot3(mySwarm.network{1}{i}.r(1,:),mySwarm.network{1}{i}.r(2,:),mySwarm.network{1}{i}.r(3,:))
% end
% axis equal
% grid on
% pbaspect([1 1 1])
% plot3(0,0,0,'ro')
% %hold off
% disp("orbit plotting complete")

%%  Monitor the value of a singel cell over time 
% 1. Choose a grid size and create space
mySTS = stereoTaxicSpace(100000); %The volume chosen reprsents a 100km cube. This might be too big when considering real world applications, but is a perfect excample of continous monitoring when compared against the discontinuity of the propagation modelled here
% 2. Identify a populated cell at 1/4 period
R = mySwarm.network{1}{1}.r(:,2700);
index = getIndex(R,1,1,1,mySTS); 
% 3. Monitor cell
[T,F] = monitorCel_F(myClock,mySwarm.network{1},mySTS,index);
disp("monitoring complete")


%% Single voxel value over time / fixed space value over time
figure(1)
hold on 
grid on
T = myClock.timeline(:,1)/3600;
plot(T,F)
%plot(T,F,'o')
xlabel("Time since deployment [h]")
ylabel("Magnetic field intensity [nT]")
%xlim([23217.9e4 23289.9e4])
hold off

%% Distribued measurement for fixed time T = 25000
%1. Choose the time 
T = 25000; %first column is used for indexation reasons
%2. Sort the data [location, magnitude]
G_mi = zeros(length(mySwarm.network{1}),4);

for i=1:length(mySwarm.network{1})
    G_mi(i,1) = mySwarm.network{1}{i}.r(1,T);
    G_mi(i,2) = mySwarm.network{1}{i}.r(2,T);
    G_mi(i,3) = mySwarm.network{1}{i}.r(3,T);
    G_mi(i,4) = mySwarm.network{1}{i}.D_m(T,7);
end

figure(2)
scatter3(G_mi(:,1),G_mi(:,2),G_mi(:,3),10,G_mi(:,4),'filled')
colorbar
xlabel('X_{ECI} [m]')
ylabel('Y_{ECI} [m]')
zlabel('Z_{ECI} [m]')

%% Distributed measurement for fixed time T=25000 in stereotaxic space S(cubic,10)
%Get the position and magnitude of each measurement
R = zeros(3,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    R(:,i) = mySwarm.network{1}{i}.r(:,25000);
end
F = zeros(1,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    F(1,i) = mySwarm.network{1}{i}.D_m(25000,7);
end
%Get index location for each measurement
mySTS2 = stereoTaxicSpace(500);
index = zeros(3,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    index(:,i) = getIndex(R(:,i),1,1,1,mySTS2);
end
%Create a data storage analogue of the space of interest
Map = cell(max(index(1,:))-min(index(1,:)),max(index(2,:))-min(index(2,:)));
Map2 = cell(max(index(1,:))-min(index(1,:)),max(index(2,:))-min(index(2,:)));

for i=1:length(index)
    if index(1,i)-min(index(1,:)) == 0 && max(index(2,:))-index(2,i) > 0
        Map{index(1,i)-min(index(1,:))+1,max(index(2,:))-index(2,i)} = horzcat(Map{index(1,i)-min(index(1,:))+1,max(index(2,:))-index(2,i)},F(i)); 
    end
    if index(1,i)-min(index(1,:)) > 0 && max(index(2,:))-index(2,i) == 0
        Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1} = horzcat(Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1},F(i));
    end
    if index(1,i)-min(index(1,:)) > 0 && max(index(2,:))-index(2,i) > 0
        Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1} = horzcat(Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1},F(i));
    end
end
%% Homogenize map and convert to heatmap
for i=1:width(Map)
    for j = 1:height(Map)
        if Map{j,i} ~= 0
            Map{j,i} = sum(Map{j,i})/length(Map{j,i});
            %Map2{j,i} = length(Map{j,i}); 
            %length(Map{j,i})
        else
            Map{j,i} = 0;
        end
    end
end
Map = flipud(fliplr(Map).');
Map = cell2mat(Map);
heatmap(Map)
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%% 
% figure(4)
% heatmap(Map2)



%% Plotting data of a single drone over time 
%get a vector for the latitude of the satellite realative to Earth
r_s = mySwarm.network{1}{1}.r(1:2,:);
theta_s = zeros(length(myClock.timeline(:,1)),1);
for i=1:length(r_s)-1
    theta_s(i) = acos(dot(r_s(:,i),[1 0 ])/norm(r_s));
end
theta_se = theta_s-myEarth.Theta_t(:,3);
disp(theta_se(1,:))
for i = 1:length(theta_se)
    if theta_se(i) > pi/2
        theta_se(i) = 2*pi - theta_se(i);
    end
end
%% 
figure(4)
hold on 
ylabel("Magnetic field intensity [nT]")
xlabel("Time since deployment [h]")
xlim([0 7.54])
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{1}.D_m(:,3),'b')
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{1}.D_m(:,2),'g')
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{1}.D_m(:,1),'r')
%yyaxis('right')
%ylim([-7,7])
%plot(myClock.timeline(:,1),theta_se,'k')
legend('Z_{ECI}','Y_{ECI}','X_{ECI}')
% set(gca,'looseInset',get(gca,'TightInset'))
% saveas(gca,'Results/singleSatFullDataSet.jpeg')

clear ans F G_mi i index j Map n_orbits R T
%
load handel
sound(y,Fs)
