function [mat] = matprop(mat)



if ( strcmp(mat.stressstate,'plane_strain') ) 
    % Plane Strain case
      mat.D = (mat.E/((1+mat.nu)*(1-2*mat.nu)))*[ 1-mat.nu mat.nu 0; mat.nu 1-mat.nu 0; 0 0 (1-2*mat.nu)/2];
else
    % Plane Stress case
      mat.D = (mat.E/(1-mat.nu^2))*[1 mat.nu 0; mat.nu 1 0; 0  0 (1-mat.nu)/2];
end

%Element type (deprec)
mat.eltype  = 'Q4';



