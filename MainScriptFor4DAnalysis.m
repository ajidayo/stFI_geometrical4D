clear;
clear AdjustPrimNodePos

global SpDIM EPSILON
SpDIM   = 3;
EPSILON = 10^(-7);
%% Inputs
SelectPreset = 2; % Preset = {1,2} is available. See ParameterPreset for details for each settings.
[RefMeshPresetType,MeshMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer] = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
RefImpedance_SpV = ones(Num_of_Elem.SpV,1);

%%
[SpElemProperties,STElemProperties,Num_of_Elem,PrimFacePos] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos);
[D0,D1,D2,D3]               = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem);

Task                        = struct;
TaskDepGraph                = digraph;

[Task,STElemProperties] ...
    = Generate_InitialCondition_Tasks(Task,Num_of_Elem,SpElemProperties,STElemProperties);
[Task,STElemProperties] ...
    = Generate_BoundaryCondition_Tasks(Task,SpElemProperties,STElemProperties);
[Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = Generate_PML_Tasks(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph);
[Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateST_FI_Tasks_4D_ST(SpElemProperties,STElemProperties,Task,TaskDepGraph);
[Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph);

TaskOrder = SortTasks(Task,TaskDepGraph,STElemProperties,D1,D2);

cdt = 0.1;
disp(['cdt = ', num2str(cdt)])

[kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties,Num_of_Elem);

Z = ComputeImpedance_for_EachSTPs(RefImpedance_SpV,sC,sD,SpElemProperties,Num_of_Elem);
% Z = ComputeZ_3D_Space(@func_Impedance,NodePos,sG,sD,SpElemProperties)
Kappa_over_Z = kappa./Z;

[Sigma] = ComputePMLSigma_3D_Space(cdt,NodePos,sG,sD,SpElemProperties,Num_of_Elem);

%% 

Source = GenerateSource(MeshMeasurements,SpElemProperties,ElemPer,sC);
Num_of_Elem.STSource        = 0;
for SpSourceIdx = 1:size(Source,2) 
    Num_of_Elem.STSource=Num_of_Elem.STSource+Source(SpSourceIdx).UpdNum;
end

TMM                         = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Sigma,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem);
[TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Num_of_Elem,SpElemProperties);

%%

FieldDoFs  = zeros(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.PMLSpP+Num_of_Elem.PMLSpS,1);
% FieldDoFs  = rand(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.PMLSpP+Num_of_Elem.PMLSpS,1);

Num_of_Steps = 1000;
disp(['Number of Steps = ', num2str(Num_of_Steps)])
Time = 0;
FieldDoFs = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,FieldDoFs,Source);

ZConst = round(0.5*MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse);
YConst = round(0.5*MeshMeasurements.YCoord/MeshMeasurements.dyCoarse);
PlotMagneticFluxDensity_ZComponent_atZEquals(ZConst,FieldDoFs,FaceArea,PrimFacePos,MeshMeasurements)
PlotMagneticFluxDensity2D_atYEquals(FieldDoFs,YConst,FaceArea,PrimFacePos,MeshMeasurements)
PlotMagneticFluxDensity3D(FieldDoFs,FaceArea,PrimFacePos,MeshMeasurements)
