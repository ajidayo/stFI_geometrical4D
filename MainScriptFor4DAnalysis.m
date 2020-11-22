clear;
clear AdjustPrimNodePos

global SpDIM EPSILON
SpDIM   = 3;
EPSILON = 10^(-7);
%% Inputs
SelectPreset = 4; % Preset = {1,2,3,4} is available. See ParameterPreset for details for each settings.
[RefMeshPresetType,MeshMeasurements,PMLMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer,SpElemPositionIdx,EqualKappaGroupIdx] = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
NodePosPrim_M = zeros(3,Num_of_Elem.SpN);
for SpNIdx = 1:Num_of_Elem.SpN
    NodePosPrim_M(:,SpNIdx) = NodePos.Prim(SpNIdx).Vec;
end
NodePosDual_M = zeros(3,Num_of_Elem.SpV);
for SpVIdx = 1:Num_of_Elem.SpV
    NodePosDual_M(:,SpVIdx) = NodePos.Dual(SpVIdx).Vec;
end

%%
FuncHandle_Sigma123Distribut = @(PosVec) func_sigma_123(PosVec,PMLMeasurements);
[SpElemProperties,STElemProperties,Num_of_Elem,SigmaSpV] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePosDual_M,FuncHandle_Sigma123Distribut);
[D0,D1,D2,D3] = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem);

Task = GenerEmptyTask;
TaskDepGraph = digraph;
[Task,STElemProperties] ...
    = Generate_InitialCondition_Tasks(Task,Num_of_Elem,SpElemProperties,STElemProperties);
[Task,STElemProperties] ...
    = Generate_BoundaryCondition_Tasks(Task,SpElemProperties,STElemProperties);
disp('Generating PML Tasks')
[Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = Generate_PML_Tasks(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph);
% [Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
%     = GenerateST_FI_Tasks_4D_ST(SpElemProperties,STElemProperties,Task,TaskDepGraph);
disp('Generating STFI Tasks')
[Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateTask_STFI_4DST_Block(SpElemProperties,STElemProperties,Task,TaskDepGraph,sC,D2,D1);
disp('Generating SpFI Tasks')
[Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph);
TaskOrder = SortTasks(Task,TaskDepGraph,STElemProperties,D1,D2);

cdt = 0.1;
disp(['cdt = ', num2str(cdt)])

[kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePosPrim_M,NodePosDual_M,SpElemProperties,STElemProperties,Num_of_Elem);

[Sigma] = ComputePMLSigma_3D_Space(SigmaSpV,cdt,NodePosPrim_M,NodePosDual_M,sG,sC,sD,SpElemProperties);
%% 
Z0 = 1;
SphereImpedanceCenterPosVec = [...
    0.5*(MeshMeasurements.XFineFromCoord+MeshMeasurements.XFineToCoord);...
    0.5*(MeshMeasurements.YFineFromCoord+MeshMeasurements.YFineToCoord);...
    0.5*(MeshMeasurements.ZFineFromCoord+MeshMeasurements.ZFineToCoord)];
SphereImpedanceRadius = -1 ...
    +min([...
    0.5*(MeshMeasurements.XFineToCoord-MeshMeasurements.XFineFromCoord);...
    0.5*(MeshMeasurements.YFineToCoord-MeshMeasurements.YFineFromCoord);...
    0.5*(MeshMeasurements.ZFineToCoord-MeshMeasurements.ZFineFromCoord)]);
%SphereImpedanceRadius = 0;
FuncHandle_ZDistribut = @(PosVec) Func_SphereImpedance(PosVec,SphereImpedanceCenterPosVec,SphereImpedanceRadius);
Z = ComputeZ_3D_Space(FuncHandle_ZDistribut,NodePosDual_M,sC,sD,SpElemProperties,Num_of_Elem);
Kappa_over_Z = kappa./Z;
%%
% SourceCenterPos = [0.25*MeshMeasurements.XCoord;0.5*MeshMeasurements.YCoord;0.5*MeshMeasurements.ZCoord];
% Source = GenerateSource(SourceCenterPos,MeshMeasurements,SpElemProperties,ElemPer,sC,cdt);

IncWaveEDirec = [0;1;0];
IncWaveHDirec = [0;0;1];
WaveNumVec = [2*pi/20;0;0];
WaveNumber = norm(WaveNumVec);
IncWaveEFuncHandle = @(w) sin(WaveNumber*w)*Heaviside(w);
IncWaveHFuncHandle = @(w) (1/Z0)*sin(WaveNumber*w)*Heaviside(w);
Source = ...
    Generate_TFSFSource(IncWaveEFuncHandle,IncWaveHFuncHandle,IncWaveEDirec,IncWaveHDirec,WaveNumVec,SpElemProperties,sG,sC,sD,NodePosPrim_M,NodePosDual_M,cdt);

Num_of_Elem.STSource=sum([Source(1:size(Source,2)).UpdNum],2);

TMM                         = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Sigma,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem);
[TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Num_of_Elem,SpElemProperties);

%%

FieldDoFs  = [zeros(Num_of_Elem.SpP+Num_of_Elem.SpS,1);zeros(Num_of_Elem.PMLSpP+Num_of_Elem.PMLSpS,1)];
% FieldDoFs  = rand(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.PMLSpP+Num_of_Elem.PMLSpS,1);

Num_of_Steps = 8000;
disp(['Number of Steps = ', num2str(Num_of_Steps)])
Time = 0;
FieldDoFs = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,FieldDoFs,Source);
Time = Num_of_Steps*cdt;
ZConst = round(0.5*MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse);
YConst = round(0.5*MeshMeasurements.YCoord/MeshMeasurements.dyCoarse);

PlotMagneticFluxDensity_ZComponent_atZEquals(ZConst,FieldDoFs,FaceArea,SpElemPositionIdx.SpP,MeshMeasurements)
%PlotMagneticFluxDensity2D_atYEquals(FieldDoFs,YConst,FaceArea,SpElemPositionIdx.SpP,MeshMeasurements)
%PlotMagneticFluxDensity3D(FieldDoFs,FaceArea,SpElemPositionIdx.SpP,MeshMeasurements)
PlotMagneticFlux3D(FieldDoFs,FaceArea,SpElemPositionIdx.SpP,MeshMeasurements)
