function [RefMeshPresetType,MeshMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset)
switch SelectPreset
    case 1
        RefMeshPresetType = 'FDTD';
        LocalUpdateNum = 2;
        MeshMeasurements.XCoord = 30;
        MeshMeasurements.YCoord = 10;
        MeshMeasurements.ZCoord = 10;
        MeshMeasurements.dx = 1;
        MeshMeasurements.dy = 1;
        MeshMeasurements.dz = 1;
        disp(['Preset: FDTD, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Discretization Width:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dx), ...
            ', ', num2str(MeshMeasurements.dy), ...
            ', ', num2str(MeshMeasurements.dz)])
        disp(['LocalUpdateNum = ', num2str(LocalUpdateNum)])
    case 2
        RefMeshPresetType = 'FDTDWithSubgrid';
        SubgridFineness = 2;
        
        MeshMeasurements.XCoord = 30;
        MeshMeasurements.YCoord = 10;
        MeshMeasurements.ZCoord = 10;
        MeshMeasurements.XFineFromCoord = 20;
        MeshMeasurements.YFineFromCoord = 2;
        MeshMeasurements.ZFineFromCoord = 2;
        MeshMeasurements.XFineToCoord   = 20;
        MeshMeasurements.YFineToCoord   = 2;
        MeshMeasurements.ZFineToCoord   = 2;
        MeshMeasurements.dx = 1;
        MeshMeasurements.dy = 1;
        MeshMeasurements.dz = 1;
        XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
        YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
        ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
        XIdx_FineFrom = round(MeshMeasurements.XFineFromCoord/MeshMeasurements.dx)+1;
        XIdx_FineTo   = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dx);
        YIdx_FineFrom = round(MeshMeasurements.YFineFromCoord/MeshMeasurements.dy)+1;
        YIdx_FineTo   = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dy);
        ZIdx_FineFrom = round(MeshMeasurements.ZFineFromCoord/MeshMeasurements.dz)+1;
        ZIdx_FineTo   = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dz);
        
        LocalUpdateNum                      = ones(XSize,YSize,ZSize);
        MeshMeasurements.LocalGridFineness  = ones(XSize,YSize,ZSize);
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridFineness*ones(XSize,YSize,ZSize);
        MeshMeasurements.LocalGridFinenes(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness*ones(XSize,YSize,ZSize);   
        
        disp(['Preset: FDTD with Subgird, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Discretization Width:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dx), ...
            ', ', num2str(MeshMeasurements.dy), ...
            ', ', num2str(MeshMeasurements.dz)])
    otherwise 
         warning('Preset undefined.')
end
end