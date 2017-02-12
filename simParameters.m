function param = simParameters()
    
    % Main body segment
    m1 = 7.9026;       %[kg]         Mass of main body segment 
    J1 = 0.08;         %[kgm^2]      Rotational inertia of main body segment 
    lH = 0.138;        %[m]          Distance from main body CoM to hip center

    % Thigh segments
    m2 = 0.7887;       %[kg]         Mass of thigh segment
    J2 = 0.002207;     %[kgm^2]      Rotational inertia of thigh segment
    l2 = 0.0193399;    %[m]          Distance from hip center to thigh CoM
    lL2 = 0.2;         %[m]          Distance from hip center to knee center

    % Shank segments
    m3 = 0.510;        %[kg]         Mass of shank segment
    J3 = 0.006518;     %[kgm^2]      Rotational inertia of thigh segment
    l3 = 0.165235;     %[m]          Distance from knee center to shank CoM
    lL3 = 0.2385;      %[m]          Distance from knee center to foot center
    
    % Environmental constants
    g = 9.81;          %[m/s^2]      Accelerations due to gravity
    param = [m1 J1 lH m2 J2 l2 lL2 m3 J3 l3 lL3 g];
end