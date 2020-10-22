function [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer,SpElemPositionIdx] ...
    = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum)
switch RefMeshPresetType
    case 'FDTD'
        [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer] ...
            = RefMesh_FDTD(MeshMeasurements,LocalUpdateNum);
        SpElemPositionIdx = [];
    case 'FDTDWithSubgrid'
        [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer,SpElemPositionIdx] ...
            = RefMesh_Subgrid(MeshMeasurements,LocalUpdateNum);
        
    case 'BeltLike_LocalTimeStep_NoRefinement'
   
    case 'BeltLike_LocalTimeStep_LocalRefinement'
        
    otherwise
        warning('Unexpected RefMeshPresetType. No reference mesh generated.')
end
end