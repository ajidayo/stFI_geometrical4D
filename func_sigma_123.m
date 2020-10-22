function value = func_sigma_123(PosVec)

PMLMeasurement.depth = 5;
PMLMeasurement.xmax_negativeside =  5;
PMLMeasurement.xmin_positiveside = 25;
PMLMeasurement.ymax_negativeside =  5;
PMLMeasurement.ymin_positiveside = 15;
PMLMeasurement.zmax_negativeside =  5;
PMLMeasurement.zmin_positiveside = 15;

% PMLMeasurement.depth = 5;
% PMLMeasurement.xmax_negativeside =  5;
% PMLMeasurement.xmin_positiveside =  25;
% PMLMeasurement.ymax_negativeside =  5;
% PMLMeasurement.ymin_positiveside =  25;
% PMLMeasurement.zmax_negativeside =  5;
% PMLMeasurement.zmin_positiveside =  25;

% PMLMeasurement.depth = 6;
% PMLMeasurement.xmax_negativeside =  6;
% PMLMeasurement.xmin_positiveside =  24;
% PMLMeasurement.ymax_negativeside = -Inf;
% PMLMeasurement.ymin_positiveside =  Inf;
% PMLMeasurement.zmax_negativeside = -Inf;
% PMLMeasurement.zmin_positiveside =  Inf;

% PMLMeasurement.depth = 6;
% PMLMeasurement.xmax_negativeside =  6;
% PMLMeasurement.xmin_positiveside = 24;
% PMLMeasurement.ymax_negativeside =  6;
% PMLMeasurement.ymin_positiveside = 24;
% PMLMeasurement.zmax_negativeside =  6;
% PMLMeasurement.zmin_positiveside = 24;

% PMLMeasurement.depth = 5;
% PMLMeasurement.xmax_negativeside =  5;
% PMLMeasurement.xmin_positiveside = 35;
% PMLMeasurement.ymax_negativeside =  5;
% PMLMeasurement.ymin_positiveside = 25;
% PMLMeasurement.zmax_negativeside =  5;
% PMLMeasurement.zmin_positiveside = 25;

% PMLMeasurement.depth = 1;
% PMLMeasurement.xmin_positiveside =  Inf;
% PMLMeasurement.xmax_negativeside = -Inf;
% PMLMeasurement.ymin_positiveside =  Inf;
% PMLMeasurement.ymax_negativeside = -Inf;
% PMLMeasurement.zmin_positiveside =  Inf;
% PMLMeasurement.zmax_negativeside = -Inf;

% PMLMeasurement.depth = 6;
% PMLMeasurement.xmax_negativeside =  6;
% PMLMeasurement.xmin_positiveside = 44;
% PMLMeasurement.ymax_negativeside =  6;
% PMLMeasurement.ymin_positiveside = 24;
% PMLMeasurement.zmax_negativeside =  6;
% PMLMeasurement.zmin_positiveside = 24;


m=4;
sigma_max = 40;

value = zeros(3,1);
if PosVec(1)>PMLMeasurement.xmin_positiveside
    value(1) = sigma_max*(abs(PosVec(1)-PMLMeasurement.xmin_positiveside)/PMLMeasurement.depth)^m;
elseif PosVec(1)<PMLMeasurement.xmax_negativeside
    value(1) = sigma_max*(abs(PosVec(1)-PMLMeasurement.xmax_negativeside)/PMLMeasurement.depth)^m;
else
    value(1) = 0;
end
if PosVec(2)>PMLMeasurement.ymin_positiveside
    value(2) = sigma_max*(abs(PosVec(2)-PMLMeasurement.ymin_positiveside)/PMLMeasurement.depth)^m;
elseif PosVec(2)<PMLMeasurement.ymax_negativeside
    value(2) = sigma_max*(abs(PosVec(2)-PMLMeasurement.ymax_negativeside)/PMLMeasurement.depth)^m;
else
    value(2) = 0;
end
if PosVec(3)>PMLMeasurement.zmin_positiveside
    value(3) = sigma_max*(abs(PosVec(3)-PMLMeasurement.zmin_positiveside)/PMLMeasurement.depth)^m;
elseif PosVec(3)<PMLMeasurement.zmax_negativeside
    value(3) = sigma_max*(abs(PosVec(3)-PMLMeasurement.zmax_negativeside)/PMLMeasurement.depth)^m;
else
    value(3) = 0;
end
end