function [RefMeshPresetType,MeshMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset)
switch SelectPreset
    case 1
        RefMeshPresetType = 'FDTD';
        LocalUpdateNum = 2;
        MeshMeasurements.XCoord = 10;
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
        MeshMeasurements.XFineFromCoord = 17;
        MeshMeasurements.YFineFromCoord = 2;
        MeshMeasurements.ZFineFromCoord = 2;
        MeshMeasurements.XFineToCoord   = 23;
        MeshMeasurements.YFineToCoord   = 8;
        MeshMeasurements.ZFineToCoord   = 8;
        MeshMeasurements.dxCoarse = 1;
        MeshMeasurements.dyCoarse = 1;
        MeshMeasurements.dzCoarse = 1;
        XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
        YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
        ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
        XIdx_FineFrom = round(MeshMeasurements.XFineFromCoord/MeshMeasurements.dxCoarse)+1;
        XIdx_FineTo   = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse);
        YIdx_FineFrom = round(MeshMeasurements.YFineFromCoord/MeshMeasurements.dyCoarse)+1;
        YIdx_FineTo   = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse);
        ZIdx_FineFrom = round(MeshMeasurements.ZFineFromCoord/MeshMeasurements.dzCoarse)+1;
        ZIdx_FineTo   = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse);
        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        disp(['Preset: FDTD with Subgird, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Discretization Width at Coarse Grids:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dxCoarse), ...
            ', ', num2str(MeshMeasurements.dyCoarse), ...
            ', ', num2str(MeshMeasurements.dzCoarse)])
    otherwise 
         warning('Preset undefined.')
end
end