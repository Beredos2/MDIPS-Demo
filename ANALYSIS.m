% MissionData: is a script that will review the MDIPS_Demo results in 
% detail to ensure that the results are relevant to the scientific
% objectives. Easy to read graphs, figures, and tables are presented that
% address the mission history and the data collected throughout the
% mission.

%Recover simulation data
load("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\MDIPS Demo\RESULTS.mat")
%Provide access to classes and function toolbox
addpath("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\/MDIPS Demo")
% 
% %% the path of an arbitrary satellite througout the mission
% r_t = mySwarm.network{1}{1}.r; 
% %t = myClock.timeline(:,1);
% figure(1)
% title("Path of an arbitrary satellite throughout the misison")
% xlabel("X_ECI")
% ylabel("Y_ECI")
% zlabel("Z_ECI")
% scatter3(r_t(1,:),r_t(2,:),r_t(3,:),10,'red','filled','o',)
% view([1,1,1])

%% Plotting relative distance between satellites over time
%This script calculates the spread of the satellites along each of the ECEI
%axis by building an adjacency matrix of the relative positions of the
%satellites, calculating the smallest distance on each file. Finally, the
%script calcualtes the average of the smallest distances, and the spread,
%recording both for each time frame. 

%Prepare figure
% figure(10)
% %title("Max distance in each axis between drones")
% xlabel("time [hours]")
% ylabel("distance [km]")
% hold on 
% 
% figure(11)
% %title("Swarm density along each dimension")
% xlabel("time [hours]")
% ylabel("density [# drones/km]")
% hold on

figure(12)
%title("Density of the swarm over time, measured over maximum calculated volume")
xlabel("time [h]")
ylabel("Density [satellites/km^3]")

%Main loop, repeats as many times as there are time frames
for t=1:10:length(myClock.timeline)
    %Build the ajdacency matrix for the first frame
    ADJ_X = zeros(length(mySwarm.network{1}));
    ADJ_Y = zeros(length(mySwarm.network{1}));
    ADJ_Z = zeros(length(mySwarm.network{1}));
    for i=1:length(mySwarm.network{1}) % 1:n
        for j=i:length(mySwarm.network{1}) %i:n
            %Fill it in!
            ADJ_X(i,j) = mySwarm.network{1}{i}.r(1,t) - mySwarm.network{1}{j}.r(1,t); %X_ECI
            ADJ_Y(i,j) = mySwarm.network{1}{i}.r(2,t) - mySwarm.network{1}{j}.r(2,t); %Y_ECI
            ADJ_Z(i,j) = mySwarm.network{1}{i}.r(3,t) - mySwarm.network{1}{j}.r(3,t); %Z_ECI
            %Fill the other half
            ADJ_X(j,i) = -ADJ_X(i,j);
            ADJ_Y(j,i) = -ADJ_Y(i,j);
            ADJ_Z(j,i) = -ADJ_Z(i,j);
        end
    end

    %ADJ = sqrt(ADJ_X.^2+ADJ_Y.^2+ADJ_Z.^2);
    %Get the max value in each dimension 
    xmax = max(max(ADJ_X));
    ymax = max(max(ADJ_Y)); 
    zmax = max(max(ADJ_Z)); 
    %Draw the greatest distance between satellites on each axis [km]
    % figure(10)
    % plot(t/3600,xmax/1000,'b.')
    % plot(t/3600,ymax/1000,'r.')
    % plot(t/3600,zmax/1000,'g.')
    % drawnow
    % pause(0.1)
    % Draw the satellite density [#/km]
    % figure(11)
    % plot(t/3600,length(mySwarm.network{1})/(xmax/1000),'b.')
    % plot(t/3600,length(mySwarm.network{1})/(ymax/1000),'r.')
    % plot(t/3600,length(mySwarm.network{1})/(zmax/1000),'g.')
    % drawnow
    % pause(0.1)
    %Draw the combined density of the swarm [#/km]
    figure(12)
    hold on
    plot(t/3600,(163*10^9)/(ymax*xmax*zmax),'r.')
    % plot(t/3600,163./xmax,'g.')
    % plot(t/3600,log10(163/(zmax)),'b.')
    % plot(t/3600,log10(163/(ymax*zmax)),'c.')
    hold off


end

% figure(10)
% legend('\Delta x_{max}','\Delta y_{max}','\Delta z_{max}')
% xlim([0,15])
% 
% figure(11)
% legend('\delta x','\delta y', '\delta z')
% xlim([0,15])
% ylim([0,100])
% 
% clear i j t n_orbits

%% Volume density plotting over time
figure(12)
%title("Density of the swarm over time, measured over maximum calculated volume")
xlabel("time [h]")
ylabel("Density [satellites/km^3]")
%Create a storage variable
vmax = zeros(length(myClock.timeline),1);
zyAmax = zeros(length(myClock.timeline),1);
xLmax = zeros(length(myClock.timeline),1);
for t=1:50:length(myClock.timeline)
    %Calculate the volume
    vmax(t) = xmax(t)*ymax(t)*zmax(t);
    zyAmax(t) = ymax(t)*zmax(t);
    xLmax(t) = xmax(t); 
end


%Plot to figure 12
figure(12)
hold on
%plot(myClock.timeline(:,1),163./zyAmax,'.')
plot(myClock.timeline(:,1),log(163./xLmax),'.')
%plot(myClock.timeline(:,1),163./vmax,'.')

    %% Statistical analysis
    %A simple average of the shortest distance is an innadequate measure of the
    %spread, because it ignores the possibility of sparcely distributed dense
    %clusters. Instead, any understanding of distribution requires a common
    %reference. Any line in the ADJ matrix can be used to set any of the drones
    %as a reference. Because the difference of distnaces is absolute, choosing
    %any of the drones as the reference is equivalent. The first is taken for
    %convenience.
    
    %Measure distance of all drones relative to (1)
    dist_x = ADJ_X(1,1:163); 
    dist_y = ADJ_Y(1,1:163);
    dist_z = ADJ_Z(1,1:163); 
    dist_v = ADJ(1,1:163); 
    IdRef = linspace(1,163,163);
    %Calculate all standard deviations
    std_x = std(dist_x);
    std_y = std(dist_y);
    std_z = std(dist_z);
    std_v = std(dist_v); 
    %calculate the mean distance
    mean_x = mean(dist_x); 
    mean_y = mean(dist_y); 
    mean_z = mean(dist_z); 
    mean_v = mean(dist_v); 
    %concatenate output
    length_stats = [mean_x std_x; mean_y std_y; mean_z std_z; mean_v std_v];
    %Distance to angle separation 
    thetax = acos(1-(dist_x.^2/(2*norm(mySwarm.network{1}{1}.r(:,t))^2)));
    thetay = acos(1-(dist_y.^2/(2*norm(mySwarm.network{1}{1}.r(:,t))^2))); 
    thetaz = acos(1-(dist_z.^2/(2*norm(mySwarm.network{1}{1}.r(:,t))^2))); 
    thetav = acos(1-(dist_v.^2/(2*norm(mySwarm.network{1}{1}.r(:,t))^2))); 
    %mean angle
    mean_thetax = mean(thetax); 
    mean_thetay = mean(thetay);
    mean_thetaz = mean(thetaz);
    mean_thetav = mean(thetav);
    %standard devation 
    std_thetax = std(thetax);
    std_thetay = std(thetay); 
    std_thetaz = std(thetaz);
    std_thetav = std(thetav); 
    %concatenate antgular data
    angle_stats = [mean_thetax std_thetax; mean_thetay std_thetay; mean_thetaz std_thetaz; mean_thetav std_thetav];
  
    clear dist_x dist_y dist_z dist_v std_x std_y std_z std_v mean_x mean_y mean_z mean_v ADJ ADJ_X ADJ_Y ADJ_Z mean thetax thetay thetaz thetav mean_thetax mean_thetay mean_thetaz mean_thetav std_thetax std_thetay std_thetaz std_thetav
    
  %% visualise spread
  %Mean separation in degrees
    figure(13)
    hold on
    title('separation in x')
    plot(myClock.timeline(t),angle_stats(1,1),'b.') %mean
    plot(myClock.timeline(t),angle_stats(1,2),'r.') %std 
    ylabel('Degrees')
    xlabel('time')
    legend('mean statellite separation in degrees','standard deviation')
    % values in y
    figure(14)
    hold on
    title('separation in y')
    plot(myClock.timeline(t),angle_stats(2,1),'b.') %mean
    plot(myClock.timeline(t),angle_stats(2,2),'r.') %std 
    ylabel('Degrees')
    xlabel('time')
    legend('mean statellite separation in degrees','standard deviation')
    % values in z
    figure(15)
    hold on
    title('separation in z')
    plot(myClock.timeline(t),angle_stats(3,1),'b.') %mean
    plot(myClock.timeline(t),angle_stats(3,2),'r.') %std 
    ylabel('Degrees')
    xlabel('time')
    legend('mean statellite separation in degrees','standard deviation')
    %values in v
    figure(16)
    hold on
    plot(myClock.timeline(t),angle_stats(4,1),'b.') %mean
    plot(myClock.timeline(t),angle_stats(4,2),'r.') %std 
    ylabel('Degrees')
    xlabel('time')
    legend('mean statellite separation in degrees','standard deviation')    
   

%% Observation of satellite spread over time 
%Initialize figure
figure(20)
hold on
axis equal
%Create a storage array
G_mi = zeros(length(mySwarm.network{1}),4);
%Open plot loop
for t=1:100:length(myClock.timeline)
%Assign frame values
    for i=1:length(mySwarm.network{1})
        G_mi(i,1) = mySwarm.network{1}{i}.r(1,t);
        G_mi(i,2) = mySwarm.network{1}{i}.r(2,t);
        G_mi(i,3) = mySwarm.network{1}{i}.r(3,t);
        G_mi(i,4) = mySwarm.network{1}{i}.D_m(t,4);
    end
    %Plot to screen
    scatter3(G_mi(:,1),G_mi(:,2),G_mi(:,3),10,G_mi(:,4),'filled')
    colorbar
    xlabel('X_{ECI} [m]')
    ylabel('Y_{ECI} [m]')
    zlabel('Z_{ECI} [m]')
    title(myClock.timeline(t,1))
    view([1 1 1])
    drawnow
    pause(0.1)
end

clear t G_mi i 

%% Distance of drone 1 relative to every other drone
%Prepare figure
figure(21)
xlabel("ID")
ylabel("Distance [m]")
xlim([0,length(mySwarm.network{1})]);
hold on

%create drone axis
IRef = linspace(1,length(mySwarm.network{1}),length(mySwarm.network{1}));
%create storage vector
distance = zeros(1,length(mySwarm.network{1}));
%Plot all values over time! 
for t=1:100:length(myClock.timeline)
    for i=1:length(mySwarm.network{1})
        distance(i) = norm(mySwarm.network{:}{i}.r(3,t)-mySwarm.network{:}{1}.r(3,t));
    end
        plot(IRef,distance,'.')
        title("Distance of ID:1 to all other drones \n",t)
        drawnow
        pause(0.1)
end

%% Figures 30-39: Distributed measurements 
%The distributed measurement analysis is carried out y observing the
%performance of the swarm in its three topological states: linear
%expansion, planar expansion, and transitory. The time stamps used for
%analyisis are identified from the SwarmSizeTracking.fig file. 

% *** TRANSITORY-STATE *** %
% I would like to show three events side by side, measured once the swarm
% has developed. The middle point between the y and z intersection with the
% x function is taken as a suitable place for observing the transitory
% topology, as it presents a "middle ground". For the three events of
% interest, I have selected the first three steps of stable expansion,
% biased towards positive time (i.e., the far side of the step). 

% 1. Select time stamps 
t_1 = 3.9725*3600;  % seconds
t_2 = 6.87528*3600; % seconds
t_3 = 9.75028*3600; % seconds
% Get index for the timeline vector 
ti_1 = t_1*2; 
ti_2 = t_2*2;
ti_3 = t_3*3;
T = [ti_1,ti_2,ti_3]; 
clear t_1 t_2 t_3 ti_1 ti_2 ti_3

% Plotting
ID = ['a','b','c'];

figure(30);

for i=1:3
    %2. Sort the data [location, magnitude]
    G_mi = zeros(length(mySwarm.network{1}),4);
    for j=1:length(mySwarm.network{1})
        G_mi(j,1) = mySwarm.network{1}{j}.r(1,i);
        G_mi(j,2) = mySwarm.network{1}{j}.r(2,i);
        G_mi(j,3) = mySwarm.network{1}{j}.r(3,i);
        G_mi(j,4) = mySwarm.network{1}{j}.D_m(i,4);
        if i==1
            disp(j)
            disp(G_mi)
        end
    end
    % %generate figure number
    % fignum = strcat(num2str(Fig),num2str(i)); 
    % fignum = str2double(fignum);
    % %plot
    % figure(fignum)
    % hold on
    nexttile(i)
    scatter3(G_mi(:,1),G_mi(:,2),G_mi(:,3),5,G_mi(:,4),'filled')
    colorbar
    title(ID(i))
    xlabel('X_{ECI} [m]')
    ylabel('Y_{ECI} [m]')
    zlabel('Z_{ECI} [m]')
    %view([0 0 1])
end




