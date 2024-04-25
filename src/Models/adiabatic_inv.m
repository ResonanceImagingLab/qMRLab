classdef adiabatic_inv < AbstractModel

        properties
        MRIinputs = {'PulseParams', 'DefaultTissueParams'}; % No input data needed 
        xnames = {}; % Not sure what there are for but dont think needed 
        voxelwise = 0; % No voxel wise fitting 

        % Creates sections in protocol boxes. Fill default vals later 
        Prot = struct('PulseParams', struct('Format',{{'beta(rad/s)'  'A0' 'n' 'nSamples' 'Q' 'Trf'}} ...
            ,'Mat',[]), ...
            'DefaultTissueParams', struct('Format',{{}},'Mat',[]));

                          
        buttons = {};
        options= struct();
         
end 
methods 


function obj = adiabatic_inv()
            obj.options = button2opts(obj.buttons);
        end
end 
end 