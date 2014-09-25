%  Miguel Mendes
%  Last update: 2014.05.03

% ------------------------------------------------------------------------
%   Subermal Cz grid
% ------------------------------------------------------------------------
clear all; close all; clc; tic;


%% ------------------------------------------------------------------------
% HEAD LAYERS
load('HeadLayers.mat');

wm  = Head5.wm ;
gm  = Head5.gm ;
csf = Head5.csf ;
sk = Head5.sk ;
sc = Head5.sc ;

%% ------------------------------------------------------------------------
% ELECTRODES Centers

centers = zeros(5,5,3);

centers(1,1,:) = [20        10      91.5] ;
centers(1,2,:) = [10        10      93.5] ;
centers(1,3,:) = [0.75      9.8     95.3] ;
centers(1,4,:) = [-9.9      10.3   94.1] ; %
centers(1,5,:) = [-20       9.7     92] ;

centers(2,1,:) = [20    20      91.25];
centers(2,2,:) = [10    19.8	93.7];
centers(2,3,:) = [0.6	20      94.6] ;
centers(2,4,:) = [-9.7	20      93.5] ;
centers(2,5,:) = [-20	20      91.5] ;

centers(3,1,:) = [20		29.5       90.3] ;
centers(3,2,:) = [10        29.5       92.45] ;
centers(3,3,:) = [0.0193 	30.2713     93.4340] ;
centers(3,4,:) = [-10  29.75  92.5] ; 
centers(3,5,:) = [-19.8  29.85  90] ;

centers(4,1,:) = [20        40          87.5];
centers(4,2,:) = [10        40          89.5];
centers(4,3,:) = [-0.0289 	40.3973 	90.5530] ;
centers(4,4,:) = [-10       40          89.4];
centers(4,5,:) = [-20       40          87.2] ;

centers(5,1,:) = [20        49          83.5];
centers(5,2,:) = [10     	49          85];
centers(5,3,:) = [-0.3398 	49.7903 	85.3660] ;
centers(5,4,:) = [-10       49          85];
centers(5,5,:) = [-20       49          83];

% Compute Distances
dist = zeros(5,5);
for i = 1:5
    for j = 1:5
    diff = centers(i,j,:)-centers(3,3,:);
    d = norm(diff(:));
    dist(i,j) = d;
    end
end
dist

%% ------------------------------------------------------------------------
% ELECTRODES Meshing

Grid = struct;
centers2 = zeros(25,3);

figure; hold on;
for i=1:5
    
    for j=1:5
        center = centers(i,j,:);
        centers2(5*(i-1)+j,:) = center(:);
        [en,ef,ee]=meshasphere(center(:),2.0,1,1); % radius: 2mm
        Grid.(sprintf('e%d',5*(i-1)+j)).n = en;  % final ns
        Grid.(sprintf('e%d',5*(i-1)+j)).f = ef;  % final fs
        plotmesh(en,ef,'facecolor','b')
    end 
end
plotmesh(sk.n,sk.f,'facecolor','w')


% Import Electrodes mesh (ns and fs)
e1_n = Grid.e1.n;       e1_f = Grid.e1.f;
e2_n = Grid.e2.n;       e2_f = Grid.e2.f;
e3_n = Grid.e3.n;       e3_f = Grid.e3.f;
e4_n = Grid.e4.n;       e4_f = Grid.e4.f;
e5_n = Grid.e5.n;       e5_f = Grid.e5.f;
e6_n = Grid.e6.n;       e6_f = Grid.e6.f;
e7_n = Grid.e7.n;       e7_f = Grid.e7.f;
e8_n = Grid.e8.n;       e8_f = Grid.e8.f;
e9_n = Grid.e9.n;       e9_f = Grid.e9.f;
e10_n = Grid.e10.n;     e10_f = Grid.e10.f;
e11_n = Grid.e11.n;     e11_f = Grid.e11.f;
e12_n = Grid.e12.n;     e12_f = Grid.e12.f;
e13_n = Grid.e13.n;     e13_f = Grid.e13.f;
e14_n = Grid.e14.n;     e14_f = Grid.e14.f;
e15_n = Grid.e15.n;     e15_f = Grid.e15.f;
e16_n = Grid.e16.n;     e16_f = Grid.e16.f;
e17_n = Grid.e17.n;     e17_f = Grid.e17.f;
e18_n = Grid.e18.n;     e18_f = Grid.e18.f;
e19_n = Grid.e19.n;     e19_f = Grid.e19.f;
e20_n = Grid.e20.n;     e20_f = Grid.e20.f;
e21_n = Grid.e21.n;     e21_f = Grid.e21.f;
e22_n = Grid.e22.n;     e22_f = Grid.e22.f;
e23_n = Grid.e23.n;     e23_f = Grid.e23.f;
e24_n = Grid.e24.n;     e24_f = Grid.e24.f;
e25_n = Grid.e25.n;     e25_f = Grid.e25.f;

[n_grid,f_grid] = mergemesh( e1_n,e1_f ,e2_n,e2_f , e3_n,e3_f , e4_n,e4_f ,e5_n,e5_f  ,...
                             e6_n,e6_f, e7_n,e7_f, e8_n,e8_f ,  e9_n,e9_f, e10_n,e10_f ,...
                             e11_n,e11_f , e12_n,e12_f, e13_n,e13_f ,e14_n,e14_f , e15_n,e15_f , ...
                             e16_n,e16_f,e17_n,e17_f, e18_n,e18_f , e19_n,e19_f,e20_n,e20_f ,...
                             e21_n,e21_f , e22_n,e22_f, e23_n,e23_f , e24_n,e24_f, e25_n,e25_f );
                         


% Boolean operation
[n_all,f_all] = surfboolean(sk.n,sk.f , 'all' , n_grid , f_grid);
n_all(:,end) = [];

f_i = f_all(find(f_all(:,4)==3),:); % Skull outside electrodes
f_e = f_all(find(f_all(:,4)==4),:); % Electrodes outside head
f_s = f_all(find(f_all(:,4)==1),:); % Skull inside electrodes

f_s(:,4) = 4;
f_e(:,4) = 6;
f_i(:,4) = 7;

% figure; hold on;
% plotmesh(n_all,f_s,'x>0','facecolor','b');
% plotmesh(n_all,f_i,'facecolor','y');
% plotmesh(n_all,f_e,'facecolor','r');
% xlabel('X axis'); ylabel('y axis'); zlabel('z axis');
% view([-1 0 0 ]);


%% Join all surfaces
[n_end,f_end] = mergemesh(wm.n,wm.f   ,   gm.n , gm.f   ,   csf.n,csf.f      ,   n_all,f_s   ,  n_all,f_i   ,  n_all,f_e,   sc.n,sc.f );
[n_end,f_end] = removeisolatednode(n_end,f_end);

%% Create meshed volume
p1 = [ 22  60 -20];
p2 = [ 50 -34 -10];
p3 = [ 6  -52 -40];
p4 = [ 9  -76 -52];
p5 = [ 80 -10 -63];

regions=[ p1 ; p2 ; p3 ; p4  ; p5 ];%innerpoints(13,:)]; % innerpoints contains the interior points of the half-spherical electrodes
p_ext = [60 -100 0];

[node,elem,face]=surf2mesh(n_end,f_end,min(n_end(:,1:3))-1,max(n_end(:,1:3))+1,1,100,regions,[p_ext],0);

%% Plotting

% % Faces
dfaces = unique(face(:,end))
face(find(face(:,end)==dfaces(1)),4)=1; % WM
face(find(face(:,end)==dfaces(2)),4)=2; % GM
face(find(face(:,end)==dfaces(3)),4)=3; % CSF
face(find(face(:,end)==dfaces(4)),4)=4; % SK
face(find(face(:,end)==dfaces(5)),4)=5; % circular
face(find(face(:,end)==dfaces(6)),4)=6; % half-sph
face(find(face(:,end)==dfaces(7)),4)=7; % SC
% 

figure(7); hold on;
plotmesh(node,face(find(face(:,end)==1),:),'x>0','facecolor','m') % White matter
plotmesh(node,face(find(face(:,end)==2),:),'x>0','facecolor','y') % Grey matter
plotmesh(node,face(find(face(:,end)==3),:),'x>0','facecolor','b') % CSF
plotmesh(node,face(find(face(:,end)==4),:),'x>0','facecolor','g') % Skull
plotmesh(node,face(find(face(:,end)==5),:),'x>0','facecolor','r') % intersection
plotmesh(node,face(find(face(:,end)==6),:),'x>0','facecolor','y') % half-electrode
plotmesh(node,face(find(face(:,end)==7),:),'x>0','facecolor','b') % Scalp


% % Elements
delem = unique(elem(:,end))
elem(find(elem(:,5)==delem(1)),5)=6;	% electrodes
elem(find(elem(:,5)==delem(2)),5)=1;
elem(find(elem(:,5)==delem(3)),5)=2;
elem(find(elem(:,5)==delem(4)),5)=3;
elem(find(elem(:,5)==delem(5)),5)=4;
elem(find(elem(:,5)==delem(6)),5)=5;

figure(8); hold on;
plotmesh(node,elem(find(elem(:,end)==1),:),'facecolor','w') % WM
plotmesh(node,elem(find(elem(:,end)==2),:),'x>0','facecolor','b') % GM
plotmesh(node,elem(find(elem(:,end)==3),:),'x>5','facecolor','g') % CSF
plotmesh(node,elem(find(elem(:,end)==4),:),'x>10','facecolor','m') % Skull
plotmesh(node,elem(find(elem(:,end)==5),:),'x>15','facecolor','y') % Scalp

plotmesh(node,elem(find(elem(:,end)==2),:),'z<15','facecolor','b') % GM
plotmesh(node,elem(find(elem(:,end)==3),:),'z<10','facecolor','g') % CSF
plotmesh(node,elem(find(elem(:,end)==4),:),'z<5','facecolor','m') % Skull
plotmesh(node,elem(find(elem(:,end)==5),:),'z<0','facecolor','y') % Scalp

plotmesh(node,elem(find(elem(:,end)==6),:),'x>15','facecolor','b') % Electrodes
plotmesh(node,elem(find(elem(:,end)==6),:),'z<0','facecolor','b') % Electrodes


% Domains
surface_domains = [5 4 ; 4 3 ; 3 2 ; 2 1 ; 2 6 ; 6 1; 1 0];
%                   wm     gm   csf    sk   int  half  sc


% %% Saving the data in a *.mphtxt file
tic; mesh2comsol_v2(node, face, elem, surface_domains, 'Czgrid.mphtxt'); toc_m2c = toc;

