%----------------------CELESTIAL BODY-------------------------------------%
%celeBod(R,V,Theta,w) instantiates a celestial body in location R with
%respect to the frame of reference, velocity V, orientation Theta, rotation
%rate w, where R,V,Theta, and w are vectors in three dimensions, and whose
%values are measured with respect to an inertial frame of reference.   

classdef celeBod < handle
    properties
        R       %Position m
        V       %Velocity m/s
        Theta_0   %Orientation rad
        Theta_t %Orientation history (array)
        w       %Angular rate rad/s
        radius  %Planetary radius in km
    end

    methods
        %------------------------CONSTRUCTOR-------------------------------
        function CelestialBody = celeBod(R,V,Theta,w,radius,C_obj)
            CelestialBody.R = R;
            CelestialBody.V = V;
            CelestialBody.Theta_0 = Theta;
            CelestialBody.w = w;
            CelestialBody.radius = radius;
            %Calculate orientation history
            CelestialBody.Theta_t = zeros(length(C_obj.timeline),3); %Orientations stored in a nTx3 array, where the columns are Theta_X,Theta_Y,Theta_Z in the ECI frame, and states throughout time are read top to bottom.
            for t=1:length(C_obj.timeline)
                CelestialBody.Theta_t(t,:) = getOrientations(CelestialBody,C_obj.timeline(t));
            end
                
        end
        %-------------------------DYNAMICS---------------------------------
        function Theta_t = getOrientations(CelestialBody,t) %NOTE:Because of how the simulation works, it would be better for the program to calculate all orientations and store them, which can be later used as lookup.
            Theta_x = CelestialBody.Theta_0(1)+CelestialBody.w(1)*t;
            Theta_y = CelestialBody.Theta_0(2)+CelestialBody.w(2)*t;
            Theta_z = CelestialBody.Theta_0(3)+CelestialBody.w(3)*t;
            Theta_t = [Theta_x,Theta_y,Theta_z];
        end
        %------------------------MANETIC FIELD----------------------------%
        function [XYZ,H,D,I,F] = planetMagField(P_obj,height,latitude,longitude,epoch)
            % epoch to decimal year
            d = datetime(epoch,'ConvertFrom','epochtime','Format','yy/MM/dd');
            d = decyear(d);
            [XYZ,H,D,I,F] = wrldmagm(height,latitude,longitude,d);
            %Apply perturbation
            for i=1:3
                if i==1
                    XYZ(i) = XYZ(i)+XYZ(i)*0.333*cos(0.01885*epoch);
                end
                if i==2 
                    XYZ(i) = XYZ(i) + XYZ(i)*0.455*cos((pi/15)*epoch);
                end
                if i==3
                    XYZ(i) = XYZ(i)+XYZ(i)*3.66*cos(6*pi*epoch)-XYZ(i)*sin(longitude);
                end
            end
            %recalculate intensity
            F = norm(XYZ); 
        end
    end
end 