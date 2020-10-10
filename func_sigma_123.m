function value = func_sigma_123(dim,x,y,z)

PMLMeasurement.depth = 6;

% PMLMeasurement.xmin_positiveside =  25;
% PMLMeasurement.xmax_negativeside =  5;
% PMLMeasurement.ymin_positiveside =  25;
% PMLMeasurement.ymax_negativeside =  5;
% PMLMeasurement.zmin_positiveside =  25;
% PMLMeasurement.zmax_negativeside =  5;

PMLMeasurement.xmin_positiveside =  24;
PMLMeasurement.xmax_negativeside =  6;
PMLMeasurement.ymin_positiveside =  Inf;
PMLMeasurement.ymax_negativeside = -Inf;
PMLMeasurement.zmin_positiveside =  Inf;
PMLMeasurement.zmax_negativeside = -Inf;


%PMLMeasurement.xmin_positiveside =  Inf;
%PMLMeasurement.xmax_negativeside = -Inf;
%PMLMeasurement.ymin_positiveside =  Inf;
%PMLMeasurement.ymax_negativeside = -Inf;
%PMLMeasurement.zmin_positiveside =  Inf;
%PMLMeasurement.zmax_negativeside = -Inf;

m=4;
sigma_max = 40;

switch dim 
    case 1
        if x>PMLMeasurement.xmin_positiveside
            value = sigma_max*(abs(x-PMLMeasurement.xmin_positiveside)/PMLMeasurement.depth)^m;
        elseif x<PMLMeasurement.xmax_negativeside
            value = sigma_max*(abs(x-PMLMeasurement.xmax_negativeside)/PMLMeasurement.depth)^m;
        else
            value = 0;
        end
    case 2
        if y>PMLMeasurement.ymin_positiveside
            value = sigma_max*(abs(y-PMLMeasurement.ymin_positiveside)/PMLMeasurement.depth)^m;
        elseif y<PMLMeasurement.ymax_negativeside
            value = sigma_max*(abs(y-PMLMeasurement.ymax_negativeside)/PMLMeasurement.depth)^m;
        else
            value = 0;
        end
    case 3
        if z>PMLMeasurement.zmin_positiveside
            value = sigma_max*(abs(z-PMLMeasurement.zmin_positiveside)/PMLMeasurement.depth)^m;
        elseif z<PMLMeasurement.zmax_negativeside
            value = sigma_max*(abs(z-PMLMeasurement.zmax_negativeside)/PMLMeasurement.depth)^m;
        else
            value = 0;
        end
end
end