%% BRAINCARE
%  Miguel Mendes
%  Last update: 2014.05.26 at 11:00
% ------------------------------------------------------------------------
%       Head Model
%       Subdermal Oz Strip 5
%       FORWARD Problem using the Reciprocity Theorem
% 
% ------------------------------------------------------------------------

%% STARTUP
%
% COMSOL LiveLink with MATLAB
%   COMSOL Multiphysics 43b
%   MATALB R2013b
addpath (strcat(getenv('COMSOL_ROOT'),'/mli'))
mphstart( str2double(getenv('MYLINKPORT')))

import com.comsol.model.*
import com.comsol.model.util.*

%% MODEL
model = ModelUtil.create('Model');
model.modelPath('/home/rodrigum/BrainCare/Subdermal/Oz_grid');
model.modelNode.create('mod1');

geom1 = model.geom.create('geom1', 3);  % 3D model
mesh1 = model.mesh.create('mesh1', 'geom1');
phys = model.physics.create('ec', 'ConductiveMedia', 'geom1');
model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');
model.study('std1').feature('stat').activate('ec', true);
model.geom('geom1').lengthUnit('mm');   % dimensions in mm
model.geom('geom1').run;


%% MESH

% Import five-layered Head Model
HM_imp = mesh1.feature.create('imp1', 'Import');
HM_imp.set('filename', '/home/rodrigum/BrainCare/Subdermal/Oz_grid/Ozgrid.mphtxt');
HM_imp.importData;
HM_imp.name('HeadModel');
mesh1.run;


%% DOMAINS & BOUNDARIES

% 3D domains
elecID = (7:31); 

% 2D boundaries
bound = zeros(5,5);
bound(1,:) = [51 45 29 25 15];
bound(2,:) = [47 41 27 21  9];
bound(3,:) = [55 37 31 23 11];
bound(4,:) = [49 40 35 19  6];
bound(5,:) = [53 43 33 17 13];




%% MATERIALS
%                                                                    Units:
% Conductivity values from Katrina and Narayan's paper                [S/m]
% Permittivity values found in the website (according to conductivity)  [1]

%   White Matter
wm = model.material.create('mat1');
wm.propertyGroup('def').set('electricconductivity', {'0.14'});
wm.propertyGroup('def').set('relpermittivity', {'212'});
wm.selection.set([5]);
wm.name('wm_mat');

%   Grey Matter
gm = model.material.create('mat2');
gm.propertyGroup('def').set('electricconductivity', {'0.33'});
gm.propertyGroup('def').set('relpermittivity', {'250'});
gm.selection.set([4]);
gm.name('gm_mat');

%   CSF
csf = model.material.create('mat3');
csf.propertyGroup('def').set('electricconductivity', {'1.82'});
csf.propertyGroup('def').set('relpermittivity', {'110'});
csf.selection.set([2]);
csf.name('csf_mat');

%   Skull
sk = model.material.create('mat4');
sk.propertyGroup('def').set('electricconductivity', {'0.058'}); % Skull live
sk.propertyGroup('def').set('relpermittivity', {'17'});
sk.selection.set([3 6]); % 6 is the ID of the small domain
sk.name('skull_mat');

%   Scalp
sc = model.material.create('mat5');
sc.propertyGroup('def').set('electricconductivity', {'0.43'});
sc.propertyGroup('def').set('relpermittivity', {'110'});  % Wet skin
sc.selection.set([1]); 
sc.name('scalp_mat');

%   Platinum Electrodes
pt = model.material.create('mat6');
pt.propertyGroup('def').set('electricconductivity', {'9.44e6'});
pt.propertyGroup('def').set('relpermittivity', {'0.14'}); % Missing good REF!
pt.selection.set(elecID);
pt.name('elec_mat');

%% ELECTRIC IDENTITIES

% GROUND
ground = phys.feature.create('gnd1', 'Ground', 2);
ground.selection.all;
ground.selection.set([31]); % The ID of the half-spherical boundary Oz


%% SIMULATION

% There are five (i = 1,2,3,4,5) horinzontal electrode strips

i = 5;

for j = 1:5

% TERMINAL: Inject current
term_name = sprintf('term%d',j);
term = phys.feature.create(term_name, 'Terminal', 2);
term.selection.set(bound(i,j));
term.set('I0', 1, '1');

% INSULATION: exclude the ground and the terminal from the list
bound2 = bound(:);
bound2 = bound2(find(bound2~=31)); % remove ground
bound2 = bound2(find(bound2~=bound(i,j))); % remove terminal

insul_name = sprintf('ein%d',j+1);
insul = phys.feature.create(insul_name, 'ElectricInsulation', 2);
insul.selection.set(bound2);

% SOLUTION
model.study.remove('std1');
model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');
model.study('std1').feature('stat').activate('ec', true);

sol1 = model.sol.create('sol1');
sol1.study('std1');
sol1.feature.create('st1', 'StudyStep');
sol1.feature('st1').set('study', 'std1');
sol1.feature('st1').set('studystep', 'stat');
sol1.feature.create('v1', 'Variables');
sol1.feature('v1').set('control', 'stat');
sol1.feature.create('s1', 'Stationary');
sol1.feature('s1').feature.create('fc1', 'FullyCoupled');
sol1.feature('s1').feature.create('i1', 'Iterative');
sol1.feature('s1').feature('i1').set('linsolver', 'cg');
sol1.feature('s1').feature('fc1').set('linsolver', 'i1');
sol1.feature('s1').feature('i1').feature.create('mg1', 'Multigrid');
sol1.feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
sol1.feature('s1').feature.remove('fcDef');
%sol1.feature('s1').set('nonlin','on');
sol1.attach('std1');
sol1.runAll;


% PLOTING
pg1 = model.result.create('pg1', 'PlotGroup3D');
pg1.name('Electric Potential (ec)');
pg1.set('oldanalysistype', 'noneavailable');
pg1.set('data', 'dset1');
pg1.feature.create('mslc1', 'Multislice');
pg1.feature('mslc1').name('Multislice');
pg1.feature('mslc1').set('oldanalysistype', 'noneavailable');
pg1.feature('mslc1').set('data', 'parent');
pg1.run;



if j == 1
    
    E51 = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
    save('E51','E51');
    model.save('OzE51');
    
elseif j == 2
    
    E52 = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
    save('E52','E52');
    model.save('OzE52');
    
elseif j == 3
    
    E53 = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
    save('E53','E53');
    model.save('OzE53');
       
elseif j == 4
    
    E54 = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
    save('E54','E54');
    model.save('OzE54');
        
elseif j == 5
    
    E55 = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
    save('E55','E55');
    model.save('OzE55');
        
end

term.active(false);
insul.active(false);

end

