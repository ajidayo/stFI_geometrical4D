% DONE 1 Properties_of_Sp_Elements, PML part
% DONE 2 RefMesh_Subgrid, Corner
% DONE 3 RefMesh_Subgrid, Boundary Conditions
% DONE 4 RefMesh_Subgrid, Given area and lengths
% DONE 5 Kappa, NonOrthogonal
% DONE 6 Kappa, Corner
% DONE 7 Source
% DONE 8 Plotfunctions

clear;
global SpDIM EPSILON
SpDIM   = 3; % Fixed at 3, immutable parameter.
EPSILON = 10^(-7);
%% Inputs
SelectPreset = 2; % Preset = {1,2} is available. See ParameterPreset for details for each settings.
[RefMeshPresetType,MeshMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer] = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
RefImpedance_SpV = ones(Num_of_Elem.SpV,1);

disp('point1')
%%
[SpElemProperties,STElemProperties,Num_of_Elem,PrimFacePos] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos);
disp('point2')
[D0,D1,D2,D3]               = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem);
disp('point2.1')

Task                        = struct;
TaskDepGraph                = digraph;
Map_SpElem_to_FirstGlobTask = struct;

[Task,STElemProperties] ...
    = Generate_InitialCondition_Tasks(Task,Num_of_Elem,SpElemProperties,STElemProperties);
[Task,STElemProperties] ...
    = Generate_BoundaryCondition_Tasks(Task,SpElemProperties,STElemProperties);
[Task,TaskDepGraph,SpElemProperties,STElemProperties,Map_SpElem_to_FirstGlobTask] ...
    = Generate_PML_Tasks(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask);
[Task,TaskDepGraph,SpElemProperties,STElemProperties,Map_SpElem_to_FirstGlobTask] ...
    = GenerateST_FI_Tasks_4D_ST(SpElemProperties,STElemProperties,Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask);
[Task,TaskDepGraph,SpElemProperties,STElemProperties,Map_SpElem_to_FirstGlobTask] ...
    = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask);

clearvars Map_SpElem_to_FirstGlobTask

TaskOrder = SortTasks(Task,TaskDepGraph,STElemProperties,D1,D2);

disp('point3')

cdt = 0.1;
disp(['cdt = ', num2str(cdt)])

[kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties,Num_of_Elem);

disp('point4')

Z = ComputeImpedance_for_EachSTPs(RefImpedance_SpV,sC,sD,SpElemProperties,Num_of_Elem);
Kappa_over_Z = kappa./Z;

disp('point4.1')

[Sigma] = ComputePMLSigma_3D_Space(cdt,NodePos,sG,sD,SpElemProperties,Num_of_Elem);

%% 
disp('point5')
Source = SourceFDTD(MeshMeasurements,SpElemProperties,ElemPer,sC);
Num_of_Elem.STSource        = 0;
for SpSourceIdx = 1:size(Source,2) 
    Num_of_Elem.STSource=Num_of_Elem.STSource+Source(SpSourceIdx).UpdNum;
end
TMM                         = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Sigma,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem);
[TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Num_of_Elem,SpElemProperties);
disp('point5.1')
%%

Num_of_PMLDualDoF_Init = size(find(SpElemProperties.SpP.PML==true),2)+size(find(SpElemProperties.SpS.PML==true),2);

FieldDoFs  = zeros(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_PMLDualDoF_Init,1);
% FieldDoFs  = rand(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_PMLDualDoF_Init,1);

Num_of_Steps = 10;
disp(['Number of Steps = ', num2str(Num_of_Steps)])
Time = 0;
disp('point6')
FieldDoFs = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,FieldDoFs,Source);
disp('point7')
ZConst = round(0.5*MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse);
PlotMagneticFluxDensity_ZComponent_atZEquals(ZConst,FieldDoFs,FaceArea,PrimFacePos,MeshMeasurements)
YConst = round(0.5*MeshMeasurements.YCoord/MeshMeasurements.dyCoarse);
PlotMagneticFluxDensity2D_atYEquals(FieldDoFs,YConst,FaceArea,PrimFacePos,MeshMeasurements)
PlotMagneticFluxDensity3D(FieldDoFs,FaceArea,PrimFacePos,MeshMeasurements)
