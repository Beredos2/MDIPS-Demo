%SANDBOX
%This is where pieces of code get tested to make sure they work before
%implementation 

% %% Simulating the orbit of a single femtosatellite over time with
% Runge-Kutta propagation [OBSOLETE]
% clear
% clc
% %Create a timeline
% dT = 300; 
% T = 900000; 
% myClock = tiktok(T,dT); 
% % Keplerian Elements of IC
% a = 29599.8;
% e = 0.0001;
% i = 0.9774;
% Omega = 1.3549;
% w = 0;
% M = 0.2645;
% % Instantiate drone
% myDrone = drone(1,[a,e,i,Omega,w,M]);
% %Propagate throughout orbits
% propagate(myDrone,myDrone.r,myDrone.v,myClock.dT,myClock.nT); 
% % plot results
% hold on
% %for t=1:myClock.nT
%     plot3(myDrone.r(1,:),myDrone.r(2,:),myDrone.r(3,:))
% %end

%% Simulating an orbiting swarm
clear
clc
%Create a timeline
myClock = tiktok(1693785600,1693866200,30); 
%Define swarm
ID = 1;
type = "fun"; 
%Varaibles for table
dType = "Adam";
dPop = 10; 
dNeed = 10; 
systems = table(dType,dPop,dNeed);
%Keplerian elements 
% Keplerian Elements of IC
a = 6800.8;
e = 0.0000001;
i = 0;
Omega = 0;
w = 0;
M = 0;
Kepler = [a,e,i,Omega,w,M]; 
mySwarm = swarm(ID,type,systems,Kepler,0);
deployBrood(mySwarm); %Modify v_0 for all drones

%Propagate orbits
for i=1:length(mySwarm.network{1})
    propagateOrb(mySwarm.network{1}{i},myClock);
end

%% Plot orbits
hold on
for i=1:length(mySwarm.network{1})
    plot3(mySwarm.network{1}{i}.r(1,:),mySwarm.network{1}{i}.r(2,:),mySwarm.network{1}{i}.r(3,:))
end
figure(1)
axis equal
grid on
pbaspect([1 1 1])
plot3(0,0,0,'ro')
hold off

clear dType dPop dNeed a e i Omega w M Kepler type T ID dT

%% Calculating the angle difference between ECEF and ECI
%ECI
X = [1 0 0];
Y = [0 1 0];
Z = [0 0 1];
%Earth 
Earth = celeBod([0 0 0],[0 0 0],[0 0 0], [0 0 2*pi/86400], 6371,myClock);
%% get data 
for i=1:length(mySwarm.network{1})
    getMdata(mySwarm.network{1}{i},Earth,myClock)
end
%% Data processing 
figure(2)
hold on
subplot(2,2,1)
for id=1:length(mySwarm.network{1})
    plot3(mySwarm.network{1}{id}.r(1,60),mySwarm.network{1}{id}.r(2,60),mySwarm.network{1}{id}.r(3,60),'r*')
end
hold off
subplot(2,2,2)
hold on
for id=1:length(mySwarm.network{1})
    plot3(mySwarm.network{1}{id}.r(1,500),mySwarm.network{1}{id}.r(2,500),mySwarm.network{1}{id}.r(3,500),'r*')
end
hold off
subplot(2,2,3)
hold on
for id=1:length(mySwarm.network{1})
    plot3(mySwarm.network{1}{id}.r(1,1000),mySwarm.network{1}{id}.r(2,1000),mySwarm.network{1}{id}.r(3,1000),'r*')
end
hold off
subplot(2,2,4)
hold on
for id=1:length(mySwarm.network{1})
    plot3(mySwarm.network{1}{id}.r(1,2000),mySwarm.network{1}{id}.r(2,2000),mySwarm.network{1}{id}.r(3,2000),'r*')
end

%% signal design 
%Create the sample space: 5 minutes
t = linspace(0,60*5,2*60*5); 
f_1 = cos((2*pi/300)*t);
f_2 = 1.5*cos((pi/15)*t);
f_3 = 2*cos((6*pi)*t); 

figure(1)
hold on
xlabel('time [s]')
ylabel('amplitude')
plot(t,f_1)
plot(t,f_2)
plot(t,f_3)
legend('f_1','f_2','f_3')
hold off

figure(2)
plot(t,f_1+f_2+f_3)

%% what LEO orbit allows me to revisit the same location an exact number of times each day, so that at the start of the next day we are at the same position? 
%-Nomenclature-
% w_e: earth's rotation rate in ECI [rad/s]
% w_s: spacecraft orbit rate in ECI [rad/s]
% T_e: period of Earth's rotatoin (a day) [s]
% T_se: spacecraft revisit time [s] 

% To ensure that the spacecraft is at the same location at the start of
% each day, 
%% How to organize satellite data into an easily accessed array? We want time on the x axis (counting left 1, to reight k), and satID on the y axis, (counting top 1 to bottom n)
%Create an array of size {n,k} 
SymData = cell(length(mySwarm.network{1}),length(myClock.timeline(:,1)));
%Find each satellite
for i=1:length(mySwarm.network{1})
    %for each time step
    for j=1:length(myClock.timeline(:,1))
      SymData{i,j} = {mySwarm.network{1}{i}.r(:,j)',mySwarm.network{1}{i}.D_m(i,:)};
    end
end

%% Get a TFSV measurmeent from SymData
%For time interval 50k, plot the spatial magnetic field variation in each
%axis. 
%Define X 
X = zeros(height(SymData),1); 
for i=1:length(X)
    X(i) = SymData{i,50000}{1}(1);
end
%Define Y
Y = zeros(height(SymData),1); 
for i=1:length(Y)
    Y(i) = SymData{i,50000}{1}(2);
end
%Define Z
Z = zeros(height(SymData),1); 
for i=1:length(Z)
    Z(i) = SymData{i,50000}{1}(3);
end
%Get Magnetic field intensity 
F = zeros(height(SymData));
for i=1:length(F)
    F(i) = SymData{i,50000}{2}(7);
end

%Plot figures side-by-side
subplot(1,3,1)
plot(X,F,'*')
ylim([4e4,7e4])
subplot(1,3,2)
plot(Y,F,'*')
ylim([4e4,7e4])
subplot(1,3,3)
plot(Z,F,'*')
ylim([4e4,7e4])

%% 
scatter3(X,Y,F)

%******************ERROR***********************%
%Not plotting correctly. Check vector contents 


