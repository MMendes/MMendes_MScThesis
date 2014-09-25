%% BRAINCARE
%  Miguel Mendes
%  Last update: 2014.04.29

% ------------------------------------------------------------------------
%       Create3DHeadModel.m was created in order to create a 5-layered
%       realistic model of the human head containing: white and grey
%       matter, CSF, scalp and skull. Additionally, 21 surface electrodes
%       were added to the model according to the 10-20 EEG system.
%
%       The *.dfs files containing the surfaces information resulted from
%  the MRI BrainSuite segmentation.
%       The output *.mphtxt file is then ready to be imported into COMSOL
%  Multiphysics as a meshed volume.
% ------------------------------------------------------------------------
clear all; close all; clc; tic;

%% Load Head info
load('HeadLayers.mat');

wm  = Head5.wm ;
gm  = Head5.gm ;
csf = Head5.csf ;
sk = Head5.sk ;
sc = Head5.sc ;

%% -----------------------------------------------------------------------

% Run Electrodes1020: this will display the electrodes on the scalp surface
% according to the 10-20 system, through spherical projection.
Electrodes1020;

% Import Electrodes mesh (nodes and faces)
e1_n = EEG.elec_1.node;       e1_f = EEG.elec_1.face;
e2_n = EEG.elec_2.node;       e2_f = EEG.elec_2.face;
e3_n = EEG.elec_3.node;       e3_f = EEG.elec_3.face;
e4_n = EEG.elec_4.node;       e4_f = EEG.elec_4.face;
e5_n = EEG.elec_5.node;       e5_f = EEG.elec_5.face;
e6_n = EEG.elec_6.node;       e6_f = EEG.elec_6.face;
e7_n = EEG.elec_7.node;       e7_f = EEG.elec_7.face;
e8_n = EEG.elec_8.node;       e8_f = EEG.elec_8.face;
e9_n = EEG.elec_9.node;       e9_f = EEG.elec_9.face;
e10_n = EEG.elec_10.node;     e10_f = EEG.elec_10.face;
e11_n = EEG.elec_11.node;     e11_f = EEG.elec_11.face;
e12_n = EEG.elec_12.node;     e12_f = EEG.elec_12.face;
e13_n = EEG.elec_13.node;     e13_f = EEG.elec_13.face;
e14_n = EEG.elec_14.node;     e14_f = EEG.elec_14.face;
e15_n = EEG.elec_15.node;     e15_f = EEG.elec_15.face;
e16_n = EEG.elec_16.node;     e16_f = EEG.elec_16.face;
e17_n = EEG.elec_17.node;     e17_f = EEG.elec_17.face;
e18_n = EEG.elec_18.node;     e18_f = EEG.elec_18.face;
e19_n = EEG.elec_19.node;     e19_f = EEG.elec_19.face;
e20_n = EEG.elec_20.node;     e20_f = EEG.elec_20.face;
e21_n = EEG.elec_21.node;     e21_f = EEG.elec_21.face;

%       Note: there should be a way to the this in a smart-code way, but I
%       still don't know how to create the variable named "ei_n", where
%       i = 1, 2, 3, ... 21.


% Mergemesh 21 Electrodes
[n_21elec,f_21elec] = mergemesh( e1_n,e1_f , e2_n,e2_f   ,  e3_n,e3_f   ,  e4_n,e4_f   , e5_n,e5_f , ...
                                 e6_n,e6_f   ,  e7_n,e7_f  , e8_n,e8_f , e9_n,e9_f , e10_n,e10_f , e11_n,e11_f , e12_n,e12_f , ...
                                 e13_n,e13_f , e14_n,e14_f , e15_n,e15_f ,  ...
                                 e16_n,e16_f  , e17_n,e17_f  , e18_n,e18_f  , e19_n,e19_f  , e20_n,e20_f , e21_n,e21_f);

%n_21elec =  e17_n;
%f_21elec =  e17_f;
          
% 6; 13; 17:18

close all;
%[n_21elec,f_21elec] = mergemesh(  e17_n,e17_f  , e18_n,e18_f);



% Defining different face domains regarding the electrode-scalp intersection
[n_all,f_all] = surfboolean(sc.n,sc.f , 'all' , n_21elec , f_21elec);

n_all(:,end) = [];

%[n_f,f_f] = meshcheckrepair(n_all,f_all , 'meshfix');

f_s = f_all(find(f_all(:,4)==1),:); % Scalp outside electrodes
f_e = f_all(find(f_all(:,4)==4),:); % Electrodes outside head
f_i = f_all(find(f_all(:,4)==3),:); % Scalp inside electrodes


f_s(:,4) = 5;
f_e(:,4) = 6;
f_i(:,4) = 7;


%%
figure(4); hold on;
plotmesh(n_all,f_s,'x>0','facecolor','b');
plotmesh(n_all,f_i,'facecolor','y');
plotmesh(n_all,f_e,'facecolor','g');
xlabel('X axis'); ylabel('y axis'); zlabel('z axis');
view([1 0 0 ]);




%% Join all surfaces
[n_end,f_end] = mergemesh(wm.n,wm.f   ,   gm.n , gm.f   ,   csf.n,csf.f   ,   sk.n,sk.f , n_all,f_s , n_all,f_i , n_all,f_e);
%[n_end,f_end] = mergemesh(wm_n2,wm_f2  , br_n2,br_f2  , is_n2,is_f2   , os_n2,os_f2,  sc_n2 , sc_f2 );

[n_end,f_end] = removeisolatednode(n_end,f_end);


% size(n_end) = [ 922055    3 ] 
% size(f_end) = [ 1833671   4 ] 
      
% Check Faces domains
f_dom = unique(f_end(:,end))

figure(5); hold on;% test line by line
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(1)),:),'x>0','facecolor','g') % scalp
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(2)),:),'x>0','facecolor','m') % half-spherical electrodes
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(3)),:),'x>0','facecolor','y') % sc-elec intersections
%hold on;
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(4)),:),'x>0','facecolor','r') % WM
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(5)),:),'x>0','facecolor','g') % GM
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(6)),:),'x>0','facecolor','b') % CSF
plotmesh(n_end,f_end(find(f_end(:,end)==f_dom(7)),:),'x>0','facecolor','m') % Skull
xlabel('X axis'); ylabel('Y axis'); zlabel('Z axis');
view([-1 0 0]);


%% Create meshed volume

p1 = [ 22  60 -20];
p2 = [ 50 -34 -10];
p3 = [ 6  -52 -40];
p4 = [ 9  -76 -52];
p5 = [ 80 -10 -63];

% Find points inside the half-spherical electrodes
[th,phi,r] = cart2sph(centers(:,1),centers(:,2),centers(:,3));  % centers came from "Electrodes1020"
[x,y,z] = sph2cart(th,phi,r+2.5);                                 % move the center of the electrode distally
innerpoints = [x y z];

some_inner = [innerpoints(1,:) ;innerpoints(2,:) ;innerpoints(3,:) ; ...
              innerpoints(4,:) ;innerpoints(6,:) ;innerpoints(7,:) ; innerpoints(9,:) ;innerpoints(16,:) ; ...
              innerpoints(17,:) ;innerpoints(18,:) ;innerpoints(19,:) ;innerpoints(21,:) ];


regions=[ p1 ; p2 ; p3 ; p4 ; p5 ];% innerpoints(1:21,:)]; % innerpoints contains the interior points of the half-spherical electrodes
p_ext = [60 -100 0];

% NOTE: maybe I just need to inner points for the all electrodes because
% there are only two face domains for all the electrodes (2 and 3).

tic;
% fprintf('\nSTEP 4:\n Volume meshing has started...\n\n');
 [node,elem,face]=surf2mesh(n_end,f_end,min(n_end(:,1:3))-1,max(n_end(:,1:3))+1,1,100,regions,[p_ext],0);
 fprintf('\nSTEP 4:\n Meshed volume has been created!\n\n');
toc_s2m = toc;

figure;
plotmesh(node,face);

%% Checking and defining domains
% At this point the unique values of both faces and elements should be
% confirmed. This means that, ideally, there are just 5 face domains and 5
% elements domains, ordered from inner to outer surfaces.

% % Faces
dfaces = unique(face(:,end))
face(find(face(:,end)==dfaces(1)),4)=1;
face(find(face(:,end)==dfaces(2)),4)=2;
face(find(face(:,end)==dfaces(3)),4)=3;
face(find(face(:,end)==dfaces(4)),4)=4;
face(find(face(:,end)==dfaces(5)),4)=5;
face(find(face(:,end)==dfaces(6)),4)=6;
face(find(face(:,end)==dfaces(7)),4)=7;
% 


figure(7); hold on;
plotmesh(node,face(find(face(:,end)==1),:),'x>0','facecolor','m') % White matter
plotmesh(node,face(find(face(:,end)==2),:),'x>0','facecolor','y') % Grey matter
plotmesh(node,face(find(face(:,end)==3),:),'x>0','facecolor','b') % CSF
plotmesh(node,face(find(face(:,end)==4),:),'x>0','facecolor','g') % Skull
plotmesh(node,face(find(face(:,end)==5),:),'x>0','facecolor','r') % Scalp
plotmesh(node,face(find(face(:,end)==6),:),'x>0','facecolor','y') % intersection
plotmesh(node,face(find(face(:,end)==7),:),'x>0','facecolor','b') % half-electrode


% 
% 
% % Elements
delem = unique(elem(:,end))

    % Head layers
elem(find(elem(:,5)==delem(1)),5)=6; % Electrodes
elem(find(elem(:,5)==delem(2)),5)=1;
elem(find(elem(:,5)==delem(3)),5)=2;
elem(find(elem(:,5)==delem(4)),5)=3;
elem(find(elem(:,5)==delem(5)),5)=4;
elem(find(elem(:,5)==delem(6)),5)=5;
% 
%     % Electrodes
% elem(find(elem(:,5)==delem(8)),5)=7; 
% elem(find(elem(:,5)==delem(9)),5)=8; 
% elem(find(elem(:,5)==delem(10)),5)=9; 
% elem(find(elem(:,5)==delem(11)),5)=10; 
% elem(find(elem(:,5)==delem(12)),5)=11; 
% elem(find(elem(:,5)==delem(13)),5)=12; 
% elem(find(elem(:,5)==delem(14)),5)=13; 
% 
% elem(find(elem(:,5)==delem(15)),5)=14; 
% elem(find(elem(:,5)==delem(16)),5)=15; 
% elem(find(elem(:,5)==delem(17)),5)=16; 
% elem(find(elem(:,5)==delem(18)),5)=17; 
% elem(find(elem(:,5)==delem(19)),5)=18; 
% elem(find(elem(:,5)==delem(20)),5)=19; 
% elem(find(elem(:,5)==delem(21)),5)=20; 
%  




figure(8); hold on;
plotmesh(node,elem(find(elem(:,end)==1),:),'facecolor','w') % White matter
plotmesh(node,elem(find(elem(:,end)==2),:),'x>0','facecolor','b') % Grey matter
plotmesh(node,elem(find(elem(:,end)==3),:),'x>5','facecolor','g') % CSF
plotmesh(node,elem(find(elem(:,end)==4),:),'x>10','facecolor','m') % Skull
plotmesh(node,elem(find(elem(:,end)==5),:),'x>15','facecolor','y') % Scalp

plotmesh(node,elem(find(elem(:,end)==6),:),'facecolor','r') % Electrodes

plotmesh(node,elem(find(elem(:,end)==2),:),'z<15','facecolor','b') % Grey matter
plotmesh(node,elem(find(elem(:,end)==3),:),'z<10','facecolor','g') % CSF
plotmesh(node,elem(find(elem(:,end)==4),:),'z<5','facecolor','m') % Skull
plotmesh(node,elem(find(elem(:,end)==5),:),'z<0','facecolor','y') % Scalp


% Domains
surface_domains = [5 4 ; 4 3 ; 3 2 ; 2 1 ; 1 0 ; 1 6 ; 6 0];
%                   wm     gm   csf    sk   sc

close all;
    

% %% Saving Head Model mesh data
HeadModel = struct();
HeadModel.nodes = node;
HeadModel.faces = face;
HeadModel.elements = elem;
save('HeadModel_data','HeadModel');
% 
% %% Saving the data in a *.mphtxt file
% fprintf('\nSTEP 3:\n Writing the *.mphtxt file...\n\n');
tic; mesh2comsol_v2(node, face, elem, surface_domains, 'HeadModel_21.mphtxt'); toc_m2c = toc;
% fprintf('\nTHE END:\n Create3DHeadModel has finished.\n\n');
% toc;

