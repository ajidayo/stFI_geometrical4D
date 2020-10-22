function [RefMeshPresetType,MeshMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset)
global EPSILON
switch SelectPreset
    case 1
        RefMeshPresetType = 'FDTD';
        LocalUpdateNum = 2;
        MeshMeasurements.XCoord = 20;
        MeshMeasurements.YCoord = 20;
        MeshMeasurements.ZCoord = 20;
        MeshMeasurements.dxCoarse = 1;
        MeshMeasurements.dyCoarse = 1;
        MeshMeasurements.dzCoarse = 1;
        disp(['Preset: FDTD, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Discretization Width:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dxCoarse), ...
            ', ', num2str(MeshMeasurements.dyCoarse), ...
            ', ', num2str(MeshMeasurements.dzCoarse)])
        disp(['LocalUpdateNum = ', num2str(LocalUpdateNum)])
    case 2
        RefMeshPresetType = 'FDTDWithSubgrid';
        SubgridFineness  = 4;
        SubgridUpdateNum = 4;
        
        MeshMeasurements.XCoord = 50;
        MeshMeasurements.YCoord = 30;
        MeshMeasurements.ZCoord = 30;
        MeshMeasurements.XFineFromCoord = 30;
        MeshMeasurements.YFineFromCoord = 10;
        MeshMeasurements.ZFineFromCoord = 10;
        MeshMeasurements.XFineToCoord   = 20;
        MeshMeasurements.YFineToCoord   = 20;
        MeshMeasurements.ZFineToCoord   = 20;
        MeshMeasurements.dxCoarse = 1;
        MeshMeasurements.dyCoarse = 1;
        MeshMeasurements.dzCoarse = 1;
        XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
        YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
        ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
        XIdx_FineFrom = round(MeshMeasurements.XFineFromCoord/MeshMeasurements.dxCoarse)+1;
        if abs(MeshMeasurements.XFineToCoord-MeshMeasurements.XCoord)<EPSILON
            XIdx_FineTo = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse) +1;
        else
            XIdx_FineTo = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse);
        end
        YIdx_FineFrom = round(MeshMeasurements.YFineFromCoord/MeshMeasurements.dyCoarse)+1;
        if abs(MeshMeasurements.YFineToCoord-MeshMeasurements.YCoord)<EPSILON
            YIdx_FineTo = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse) +1;
        else
            YIdx_FineTo = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse);
        end
        ZIdx_FineFrom = round(MeshMeasurements.ZFineFromCoord/MeshMeasurements.dzCoarse)+1;
        if abs(MeshMeasurements.ZFineToCoord-MeshMeasurements.ZCoord)<EPSILON
            ZIdx_FineTo = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse) +1;
        else
            ZIdx_FineTo = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse);
        end
        
        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridUpdateNum * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        disp(['Preset: FDTD with Subgird, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Spatial Discretization Width at Coarse Grids:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dxCoarse), ...
            ', ', num2str(MeshMeasurements.dyCoarse), ...
            ', ', num2str(MeshMeasurements.dzCoarse)])
        disp(['Spatial Discretization Width at Local Grids:dx, dy, dz = ',...
            num2str(MeshMeasurements.dxCoarse/SubgridFineness), ...
            ', ', num2str(MeshMeasurements.dyCoarse/SubgridFineness), ...
            ', ', num2str(MeshMeasurements.dzCoarse/SubgridFineness)])
        disp(['FineGrid Region: x = [',num2str(MeshMeasurements.XFineFromCoord),' ',num2str(MeshMeasurements.XFineToCoord),'], ',...
            'y = [',num2str(MeshMeasurements.YFineFromCoord),' ',num2str(MeshMeasurements.YFineToCoord),'], ',...
            'z = [',num2str(MeshMeasurements.ZFineFromCoord),' ',num2str(MeshMeasurements.ZFineToCoord),']. ']);
        disp(['Relative Temporal Discritization Width in Local Grids:',num2str(1/SubgridUpdateNum)])
    case 3
        RefMeshPresetType = 'FDTDWithSubgrid';
        SubgridFineness  = 2;
        SubgridUpdateNum = 2;
        
        MeshMeasurements.XCoord = 40;
        MeshMeasurements.YCoord = 30;
        MeshMeasurements.ZCoord = 30;
        MeshMeasurements.XFineFromCoord = 20;
        MeshMeasurements.YFineFromCoord = 10;
        MeshMeasurements.ZFineFromCoord = 10;
        MeshMeasurements.XFineToCoord   = 30;
        MeshMeasurements.YFineToCoord   = 20;
        MeshMeasurements.ZFineToCoord   = 20;
        MeshMeasurements.dxCoarse = 1;
        MeshMeasurements.dyCoarse = 1;
        MeshMeasurements.dzCoarse = 1;
        XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
        YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
        ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
        XIdx_FineFrom = round(MeshMeasurements.XFineFromCoord/MeshMeasurements.dxCoarse)+1;
        if abs(MeshMeasurements.XFineToCoord-MeshMeasurements.XCoord)<EPSILON
            XIdx_FineTo = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse) +1;
        else
            XIdx_FineTo = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse);
        end
        YIdx_FineFrom = round(MeshMeasurements.YFineFromCoord/MeshMeasurements.dyCoarse)+1;
        if abs(MeshMeasurements.YFineToCoord-MeshMeasurements.YCoord)<EPSILON
            YIdx_FineTo = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse) +1;
        else
            YIdx_FineTo = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse);
        end
        ZIdx_FineFrom = round(MeshMeasurements.ZFineFromCoord/MeshMeasurements.dzCoarse)+1;
        if abs(MeshMeasurements.ZFineToCoord-MeshMeasurements.ZCoord)<EPSILON
            ZIdx_FineTo = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse) +1;
        else
            ZIdx_FineTo = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse);
        end
        
        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridUpdateNum * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        disp(['Preset: FDTD with Subgird, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Spatial Discretization Width at Coarse Grids:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dxCoarse), ...
            ', ', num2str(MeshMeasurements.dyCoarse), ...
            ', ', num2str(MeshMeasurements.dzCoarse)])
        disp(['Spatial Discretization Width at Local Grids:dx, dy, dz = ',...
            num2str(MeshMeasurements.dxCoarse/SubgridFineness), ...
            ', ', num2str(MeshMeasurements.dyCoarse/SubgridFineness), ...
            ', ', num2str(MeshMeasurements.dzCoarse/SubgridFineness)])
        disp(['FineGrid Region: x = [',num2str(MeshMeasurements.XFineFromCoord),' ',num2str(MeshMeasurements.XFineToCoord),'], ',...
            'y = [',num2str(MeshMeasurements.YFineFromCoord),' ',num2str(MeshMeasurements.YFineToCoord),'], ',...
            'z = [',num2str(MeshMeasurements.ZFineFromCoord),' ',num2str(MeshMeasurements.ZFineToCoord),']. ']);
        disp(['Relative Temporal Discritization Width in Local Grids:',num2str(1/SubgridUpdateNum)])
    case 4
            RefMeshPresetType = 'FDTDWithSubgrid';
        SubgridFineness  = 2;
        SubgridUpdateNum = 2;
        
        MeshMeasurements.XCoord = 30;
        MeshMeasurements.YCoord = 20;
        MeshMeasurements.ZCoord = 20;
        MeshMeasurements.XFineFromCoord = 18;
        MeshMeasurements.YFineFromCoord = 8;
        MeshMeasurements.ZFineFromCoord = 8;
        MeshMeasurements.XFineToCoord   = 22;
        MeshMeasurements.YFineToCoord   = 12;
        MeshMeasurements.ZFineToCoord   = 12;
        MeshMeasurements.dxCoarse = 1;
        MeshMeasurements.dyCoarse = 1;
        MeshMeasurements.dzCoarse = 1;
        XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
        YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
        ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
        XIdx_FineFrom = round(MeshMeasurements.XFineFromCoord/MeshMeasurements.dxCoarse)+1;
        if abs(MeshMeasurements.XFineToCoord-MeshMeasurements.XCoord)<EPSILON
            XIdx_FineTo = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse) +1;
        else
            XIdx_FineTo = round(MeshMeasurements.XFineToCoord  /MeshMeasurements.dxCoarse);
        end
        YIdx_FineFrom = round(MeshMeasurements.YFineFromCoord/MeshMeasurements.dyCoarse)+1;
        if abs(MeshMeasurements.YFineToCoord-MeshMeasurements.YCoord)<EPSILON
            YIdx_FineTo = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse) +1;
        else
            YIdx_FineTo = round(MeshMeasurements.YFineToCoord  /MeshMeasurements.dyCoarse);
        end
        ZIdx_FineFrom = round(MeshMeasurements.ZFineFromCoord/MeshMeasurements.dzCoarse)+1;
        if abs(MeshMeasurements.ZFineToCoord-MeshMeasurements.ZCoord)<EPSILON
            ZIdx_FineTo = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse) +1;
        else
            ZIdx_FineTo = round(MeshMeasurements.ZFineToCoord  /MeshMeasurements.dzCoarse);
        end
        
        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridUpdateNum * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        disp(['Preset: FDTD with Subgird, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Spatial Discretization Width at Coarse Grids:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dxCoarse), ...
            ', ', num2str(MeshMeasurements.dyCoarse), ...
            ', ', num2str(MeshMeasurements.dzCoarse)])
        disp(['Spatial Discretization Width at Local Grids:dx, dy, dz = ',...
            num2str(MeshMeasurements.dxCoarse/SubgridFineness), ...
            ', ', num2str(MeshMeasurements.dyCoarse/SubgridFineness), ...
            ', ', num2str(MeshMeasurements.dzCoarse/SubgridFineness)])
        disp(['FineGrid Region: x = [',num2str(MeshMeasurements.XFineFromCoord),' ',num2str(MeshMeasurements.XFineToCoord),'], ',...
            'y = [',num2str(MeshMeasurements.YFineFromCoord),' ',num2str(MeshMeasurements.YFineToCoord),'], ',...
            'z = [',num2str(MeshMeasurements.ZFineFromCoord),' ',num2str(MeshMeasurements.ZFineToCoord),']. ']);
        disp(['Relative Temporal Discritization Width in Local Grids:',num2str(1/SubgridUpdateNum)])

    otherwise
        warning('Preset undefined.')
end
end