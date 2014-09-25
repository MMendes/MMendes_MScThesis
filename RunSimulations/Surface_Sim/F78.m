%% BRAINCARE
%  Miguel Mendes
%  Last update: 2014.05.12 at 16:18
% ------------------------------------------------------------------------
%       Head Model
%       Surface Electrodes: 10-20 EEG system
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
model.modelPath('/home/rodrigum/BrainCare/EEG1020');
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
HM_imp.set('filename', '/home/rodrigum/BrainCare/EEG1020/HeadModel_21.mphtxt');
HM_imp.importData;
HM_imp.name('HeadModel');
mesh1.run;


%% Domains and Boundaries

% Head

% Electrodes
labels    = {'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F7' 'F8' 'T3' 'T4' 'T5' 'T6' 'Fz' 'Cz' 'Pz' 'C5' 'C6'};
elec_id   = [ 20    14    22   13   23   12   21   11   19   15   25   9    27   7    10   24   17   18   16   26   8 ];
elec_b1   = [34     21    37   19   40   17   35   15   31   23   44   11   47   6    42   13   27   29   25   46   9 ]; % half-spherical
elec_b2   = [33     22    38   20   39   18   36   16   32   24   43   12   48   7    41   14   28   30   26   45   10]; % circles

% Cz is the 18 position


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
csf.selection.set([3]);
csf.name('csf_mat');

%   Skull
sk = model.material.create('mat4');
sk.propertyGroup('def').set('electricconductivity', {'0.058'}); % Skull live
sk.propertyGroup('def').set('relpermittivity', {'17'});
sk.selection.set([1 6]); % 6 is the ID of the small domain
sk.name('skull_mat');

%   Scalp
sc = model.material.create('mat5');
sc.propertyGroup('def').set('electricconductivity', {'0.43'});
sc.propertyGroup('def').set('relpermittivity', {'110'});  % Wet skin
sc.selection.set([2]); 
sc.name('scalp_mat');

%   Platinum Electrodes
pt = model.material.create('mat6');
pt.propertyGroup('def').set('electricconductivity', {'9.44e6'});
pt.propertyGroup('def').set('relpermittivity', {'0.14'}); % Missing good REF!
pt.selection.set(elec_id);
pt.name('elec_mat');

%% ELECTRIC IDENTITIES

% GROUND
ground = phys.feature.create('gnd1', 'Ground', 2);
ground.selection.all;
ground.selection.set([29]); % The ID of the half-spherical boundary of Cz is 29

%% SAVE MODEL
%model.save('EEG1020');


%% SIMULATIONS

EEG = struct;

for i = 11 : 12
    
    
    %if i ~= 18
        
        
        % TERMINAL: Inject current
        term_name = sprintf('term%d',i);
        term = phys.feature.create(term_name, 'Terminal', 2);
        term.selection.set(elec_b1(i));
        term.set('I0', 1, '1');
        
        
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
        
        
         
        if i == 11
            
            F7 = struct;
            F7.GM = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
            
            save('F7','F7')
            
            model.save('F7_model');
            
        elseif i == 12
            
            F8 = struct;
            F8.GM = mpheval(model,{'ec.normJ','ec.Jx', 'ec.Jy', 'ec.Jz','ec.sigmaxx'},'selection',4,'edim','domain');
            
            save('F8','F8')
            
            model.save('F8_model');
            
        end


        % Desactivate terminal current
        term.active(false);

end
   
% save('EEG','EEG');

