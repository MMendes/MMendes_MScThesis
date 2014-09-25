%% BRAINCARE
%  Miguel Mendes
%  Last update: 2014.03.10 at 23:24
% ------------------------------------------------------------------------
%       Placing Electrodes on the Scalp according to the 10-20 EEG system
% ------------------------------------------------------------------------

% clear all;
close all; clc;

%% ------------------------------------------------------------------------

% Scalp
Xscalp = sc.n(:,1);
Yscalp = sc.n(:,2);
Zscalp = sc.n(:,3);

%% ------------------------------------------------------------------------
% ELECTRODES Referencial
load elec_pos_new_ver2;

% Electrodes Coordinates and Labels
XYZ1020 = test_1020;
labels = test_chan_lab; % labels
labels{20,1} = 'C5' ;
labels{21,1} = 'C6' ;
Ref_XYZ1020 = zeros(19,3);

% Fp1(1) and Fp2(2)
T_Fp = makehgtform('zrotate',pi/2,'xrotate',0.4*pi/2);
Ref_XYZ1020(1:2,:) =  XYZ1020(1:2,:) * T_Fp(1:3,1:3);

% The rest of the electrodes
T_rest = makehgtform('zrotate',pi/2,'xrotate',0.2*pi/2);
Ref_XYZ1020(3:19,:) =  XYZ1020(3:19,:) * T_rest(1:3,1:3);

%% ------------------------------------------------------------------------
% Moving electrodes referential to (0,0,0)

x_ele = Ref_XYZ1020(:,1);
y_ele = Ref_XYZ1020(:,2);
z_ele = Ref_XYZ1020(:,3);

% Electrodes centroid
ex = mean(x_ele) ;
ey = mean(y_ele) ;
ez = mean(z_ele) ;

% Move Electrodes centroid to (0,0,0)
x_ele2 = x_ele - ex ;
y_ele2 = y_ele - ey ;
z_ele2 = z_ele - ez ;

Ref_XYZ1020 = [x_ele2, y_ele2, z_ele2];
%% ------------------------------------------------------------------------
% ROTATIONS: changing electrodes positions

% Final Coordinates: initialize
Rot_XYZ1020 = zeros(21,3);

%   Fp1(1) and Fp2(2)
T_Fp1 = makehgtform('xrotate',-0.005*pi/2);
Rot_XYZ1020(1,:) = Ref_XYZ1020(1,:) * T_Fp1(1:3,1:3);
Rot_XYZ1020(2,:) = Ref_XYZ1020(2,:);

%   F3(3) and F4(4)
T_F3 = makehgtform('zrotate',0.00*pi/2, 'yrotate',0.225*pi/2,  'xrotate',0.05*pi/2);
Rot_XYZ1020(3,:) =  Ref_XYZ1020(3,:) * T_F3(1:3,1:3);
T_F4 = makehgtform('zrotate',-0.00*pi/2, 'yrotate',-0.070*pi/2, 'xrotate',0.075*pi/2);
Rot_XYZ1020(4,:) =  Ref_XYZ1020(4,:) * T_F4(1:3,1:3);

%   C3(5) and C4(6)
T_C3 = makehgtform('yrotate',0.310*pi/2);
T_C4 = makehgtform('yrotate',-0.295*pi/2);
Rot_XYZ1020(5,:) = Ref_XYZ1020(5,:) * T_C3(1:3,1:3);  % C3
Rot_XYZ1020(6,:) = Ref_XYZ1020(6,:) * T_C4(1:3,1:3);  % C4

%   P3(7) and P4(8)
T_P3 = makehgtform('yrotate',0.1*pi/2,'xrotate',-0.1*pi/2);
Rot_XYZ1020(7,:) =  Ref_XYZ1020(7,:) * T_P3(1:3,1:3);
T_P4 = makehgtform('zrotate',0.005*pi/2, 'yrotate',-0.060*pi/2,'xrotate',-0.1*pi/2);
Rot_XYZ1020(8,:) =  Ref_XYZ1020(8,:) * T_P4(1:3,1:3);

%   O1(9) and O2(10)
T_O1 = makehgtform('xrotate',-0.2*pi/2);
Rot_XYZ1020(9,:) = Ref_XYZ1020(9,:) * T_O1(1:3,1:3);
T_O2 = makehgtform('xrotate',-0.210*pi/2);
Rot_XYZ1020(10,:) = Ref_XYZ1020(10,:) * T_O2(1:3,1:3);

%   F7(11) and F8(12)
T_F7 = makehgtform('zrotate',-0.015*pi/2,'yrotate',0.1*pi/2);
T_F8 = makehgtform('zrotate',-0.080*pi/2,'yrotate',-0.160*pi/2);
Rot_XYZ1020(11,:) =  Ref_XYZ1020(11,:) * T_F7(1:3,1:3);
Rot_XYZ1020(12,:) =  Ref_XYZ1020(12,:) * T_F8(1:3,1:3);

%   T3(13) and T4(14)
T_T3 = makehgtform('zrotate',-0.005*pi/2,'yrotate',0.10*pi/2);
Rot_XYZ1020(13,:) =  Ref_XYZ1020(13,:) * T_T3(1:3,1:3);
T_T4 = makehgtform('yrotate',-0.005*pi/2);
Rot_XYZ1020(14,:) =  Ref_XYZ1020(14,:) * T_T4(1:3,1:3);

%   T5(15) and T6(16)
T_T5 = makehgtform('zrotate',-0.105*pi/2,'yrotate',0.250*pi/2);
Rot_XYZ1020(15,:) =  Ref_XYZ1020(15,:) * T_T5(1:3,1:3);
T_T6 = makehgtform('yrotate',-0.1*pi/2);
Rot_XYZ1020(16,:) =  Ref_XYZ1020(16,:) * T_T6(1:3,1:3);

%   Fz(17) and Cz(18)
T_Fz = makehgtform('xrotate',0.08*pi/2);
Rot_XYZ1020(17,:) = Ref_XYZ1020(17,:) * T_Fz(1:3,1:3);
T_Cz = makehgtform('yrotate',-0.002*pi/2,'xrotate',-0.0625*pi/2);
Rot_XYZ1020(18,:) = Ref_XYZ1020(18,:) * T_Cz(1:3,1:3);

%   Pz(19)
T_Pz = makehgtform('xrotate',-0.12*pi/2);
Rot_XYZ1020(19,:) = Ref_XYZ1020(19,:) * T_Pz(1:3,1:3);

%   C5(20) and C6(21)
T_C5 = makehgtform('yrotate',-0.085*pi/2);
T_C6 = makehgtform('yrotate',0.1*pi/2);
Rot_XYZ1020(20,:) = Ref_XYZ1020(5,:) * T_C5(1:3,1:3);  % similar to C3
Rot_XYZ1020(21,:) = Ref_XYZ1020(6,:) * T_C6(1:3,1:3);  % similar to C4


% Final Coordinates: ending
x_ele3 = Rot_XYZ1020(:,1);
y_ele3 = Rot_XYZ1020(:,2);
z_ele3 = Rot_XYZ1020(:,3);

%% ------------------------------------------------------------------------
% SPHERICAL PROJECTION

% Cartesian to Spherical coordinates
[th_elec,phi_elec,r_elec] = cart2sph(x_ele3,y_ele3,z_ele3);
[th_model,phi_model,r_model] = cart2sph(Xscalp,Yscalp,Zscalp);

% Matrix with the 21 electrodes positions
Scalp_elec = zeros(21,3);

% Setting tolerance
tol_max = 10*pi/180; % tolerance
tol_min = 0.1*pi/180; % tolerance
tol_step = 0.1*pi/180;


for tol = tol_min : tol_step : tol_max          % incrementing tolerance
    
    for i = 1:21    % 21 electrodes
        
        if sum(Scalp_elec(i,:)) == 0
            
            diff_th = abs (th_model - th_elec(i));
            diff_phi = abs (phi_model - phi_elec(i));
            
            
            for j = 1:size(sc.n,1)
                
                if (diff_th(j) < tol && diff_phi(j) < tol)
                    
                    Scalp_elec(i,1) = Xscalp(j);
                    Scalp_elec(i,2) = Yscalp(j);
                    Scalp_elec(i,3) = Zscalp(j);
                    
                    node_id(i,1) = j;
                    
                end
            end
        end
    end
end

%% Plotting Scalp Electrodes
% figure(2)
% hold on;
% plot3(Scalp_elec(:,1),Scalp_elec(:,2),Scalp_elec(:,3),'bo');
% text(Scalp_elec(:,1),Scalp_elec(:,2),Scalp_elec(:,3),labels);
% xlabel('X axis'); ylabel('y axis'); zlabel('z axis');
% plot3(0,0,0,'ro');

%% ---------------------------------------------------------------------
% Create spherical electrodes and place them in the scalp

EEG = struct;   % EEG contains all the 21 electrodes information (n,f,e)
centers = zeros(21,3);

figure(3);
hold on;
for i = 1:21
   
    center = Scalp_elec(i,:);
    centers(i,:) = center;
    
    [Enode,Eface,Eelem]=meshasphere(center,6,20,10); % electrodes diameter: 12mm
    
    Eface(:,4) = i+1; Eelem(:,5) = i+1;
    
    EEG.(sprintf('elec_%d',i)).node = Enode;
    EEG.(sprintf('elec_%d',i)).face = Eface;
    %EEG.(sprintf('elec_%d',i)).elem = Eelem;
    
    plotmesh(Enode, Eface,'facecolor','b');

end

plotmesh(sc.n,sc.f,'facecolor','y');
xlabel('X axis'); ylabel('y axis'); zlabel('z axis');
view([1 0 0 ]);

save('EEG1020','EEG');


