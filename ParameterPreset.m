function [RefMeshPresetType,MeshMeasurements,PMLMeasurement,LocalUpdateNum] = ParameterPreset(SelectPreset)
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
        MeshMeasurements.XFineFromCoord = 20;
        MeshMeasurements.YFineFromCoord = 10;
        MeshMeasurements.ZFineFromCoord = 10;
        MeshMeasurements.XFineToCoord   = 30;
        MeshMeasurements.YFineToCoord   = 20;
        MeshMeasurements.ZFineToCoord   = 20;
        MeshMeasurements.XSFFromCoord = 1;
        MeshMeasurements.YSFFromCoord = 1;
        MeshMeasurements.ZSFFromCoord = 1;
        MeshMeasurements.XSFToCoord   = MeshMeasurements.XCoord;
        MeshMeasurements.YSFToCoord   = MeshMeasurements.YCoord;
        MeshMeasurements.ZSFToCoord   = MeshMeasurements.ZCoord;
        PMLMeasurement.depth = 5;
        PMLMeasurement.xmax_negativeside =  5;
        PMLMeasurement.xmin_positiveside = 45;
        PMLMeasurement.ymax_negativeside =  5;
        PMLMeasurement.ymin_positiveside = 25;
        PMLMeasurement.zmax_negativeside =  5;
        PMLMeasurement.zmin_positiveside = 25;
        
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
                XIdx_TFFrom = round(MeshMeasurements.XTFFromCoord/MeshMeasurements.dxCoarse)+1;
        if abs(MeshMeasurements.XTFToCoord-MeshMeasurements.XCoord)<EPSILON
            XIdx_TFTo = round(MeshMeasurements.XTFToCoord  /MeshMeasurements.dxCoarse) +1;
        else
            XIdx_TFTo = round(MeshMeasurements.XTFToCoord  /MeshMeasurements.dxCoarse);
        end
        YIdx_TFFrom = round(MeshMeasurements.YTFFromCoord/MeshMeasurements.dyCoarse)+1;
        if abs(MeshMeasurements.YTFToCoord-MeshMeasurements.YCoord)<EPSILON
            YIdx_TFTo = round(MeshMeasurements.YTFToCoord  /MeshMeasurements.dyCoarse) +1;
        else
            YIdx_TFTo = round(MeshMeasurements.YTFToCoord  /MeshMeasurements.dyCoarse);
        end
        ZIdx_TFFrom = round(MeshMeasurements.ZTFFromCoord/MeshMeasurements.dzCoarse)+1;
        if abs(MeshMeasurements.ZTFToCoord-MeshMeasurements.ZCoord)<EPSILON
            ZIdx_TFTo = round(MeshMeasurements.ZTFToCoord  /MeshMeasurements.dzCoarse) +1;
        else
            ZIdx_TFTo = round(MeshMeasurements.ZTFToCoord  /MeshMeasurements.dzCoarse);
        end
        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridUpdateNum * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        
        MeshMeasurements.isInSFRegion = false(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.isInSFRegion(XIdx_TFFrom:XIdx_TFTo,YIdx_TFFrom:YIdx_TFTo,ZIdx_TFFrom:ZIdx_TFTo)...
            = true;
        
        
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
        disp(['Total-Field Region: x = [',num2str(MeshMeasurements.XTFFromCoord),' ',num2str(MeshMeasurements.XTFToCoord),'], ',...
            'y = [',num2str(MeshMeasurements.YTFFromCoord),' ',num2str(MeshMeasurements.YTFToCoord),'], ',...
            'z = [',num2str(MeshMeasurements.ZTFFromCoord),' ',num2str(MeshMeasurements.ZTFToCoord),']. ']);

    case 3
        RefMeshPresetType = 'FDTDWithSubgrid';
        SubgridFineness  = 2;
        SubgridUpdateNum = 2;
        
        MeshMeasurements.XCoord = 40;
        MeshMeasurements.YCoord = 30;
        MeshMeasurements.ZCoord = 30;
        MeshMeasurements.XFineFromCoord = 10;
        MeshMeasurements.YFineFromCoord = 10;
        MeshMeasurements.ZFineFromCoord = 10;
        MeshMeasurements.XFineToCoord   = 20;
        MeshMeasurements.YFineToCoord   = 20;
        MeshMeasurements.ZFineToCoord   = 20;
        MeshMeasurements.XTFFromCoord = 7;
        MeshMeasurements.YTFFromCoord = 7;
        MeshMeasurements.ZTFFromCoord = 7;
        MeshMeasurements.XTFToCoord   = 33;
        MeshMeasurements.YTFToCoord   = 23;
        MeshMeasurements.ZTFToCoord   = 23;
        PMLMeasurement.depth = 5;
        PMLMeasurement.xmax_negativeside =  5;
        PMLMeasurement.xmin_positiveside = 35;
        PMLMeasurement.ymax_negativeside =  5;
        PMLMeasurement.ymin_positiveside = 25;
        PMLMeasurement.zmax_negativeside =  5;
        PMLMeasurement.zmin_positiveside = 25;

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
        
        XIdx_TFFrom = round(MeshMeasurements.XTFFromCoord/MeshMeasurements.dxCoarse)+1;
        if abs(MeshMeasurements.XTFToCoord-MeshMeasurements.XCoord)<EPSILON
            XIdx_TFTo = round(MeshMeasurements.XTFToCoord  /MeshMeasurements.dxCoarse) +1;
        else
            XIdx_TFTo = round(MeshMeasurements.XTFToCoord  /MeshMeasurements.dxCoarse);
        end
        YIdx_TFFrom = round(MeshMeasurements.YTFFromCoord/MeshMeasurements.dyCoarse)+1;
        if abs(MeshMeasurements.YTFToCoord-MeshMeasurements.YCoord)<EPSILON
            YIdx_TFTo = round(MeshMeasurements.YTFToCoord  /MeshMeasurements.dyCoarse) +1;
        else
            YIdx_TFTo = round(MeshMeasurements.YTFToCoord  /MeshMeasurements.dyCoarse);
        end
        ZIdx_TFFrom = round(MeshMeasurements.ZTFFromCoord/MeshMeasurements.dzCoarse)+1;
        if abs(MeshMeasurements.ZTFToCoord-MeshMeasurements.ZCoord)<EPSILON
            ZIdx_TFTo = round(MeshMeasurements.ZTFToCoord  /MeshMeasurements.dzCoarse) +1;
        else
            ZIdx_TFTo = round(MeshMeasurements.ZTFToCoord  /MeshMeasurements.dzCoarse);
        end
        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridUpdateNum * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        MeshMeasurements.isInSFRegion = true(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.isInSFRegion(XIdx_TFFrom:XIdx_TFTo,YIdx_TFFrom:YIdx_TFTo,ZIdx_TFFrom:ZIdx_TFTo)...
            = false;        
        
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
        disp(['Total-Field Region: x = [',num2str(MeshMeasurements.XTFFromCoord),' ',num2str(MeshMeasurements.XTFToCoord),'], ',...
            'y = [',num2str(MeshMeasurements.YTFFromCoord),' ',num2str(MeshMeasurements.YTFToCoord),'], ',...
            'z = [',num2str(MeshMeasurements.ZTFFromCoord),' ',num2str(MeshMeasurements.ZTFToCoord),']. ']);
    case 4
        RefMeshPresetType = 'FDTDWithSubgrid';
        SubgridFineness  = 2;
        SubgridUpdateNum = 2;
        
        MeshMeasurements.XCoord = 34;
        MeshMeasurements.YCoord = 24;
        MeshMeasurements.ZCoord = 24;
        MeshMeasurements.XFineFromCoord = 19;
        MeshMeasurements.YFineFromCoord = 9;
        MeshMeasurements.ZFineFromCoord = 9;
        MeshMeasurements.XFineToCoord   = 25;
        MeshMeasurements.YFineToCoord   = 15;
        MeshMeasurements.ZFineToCoord   = 15;
        MeshMeasurements.XTFFromCoord = 7;
        MeshMeasurements.YTFFromCoord = 7;
        MeshMeasurements.ZTFFromCoord = 7;
        MeshMeasurements.XTFToCoord   = 27;
        MeshMeasurements.YTFToCoord   = 17;
        MeshMeasurements.ZTFToCoord   = 17;
        %         MeshMeasurements.XTFToCoord   = MeshMeasurements.XCoord;
        %         MeshMeasurements.YTFToCoord   = MeshMeasurements.YCoord;
        %         MeshMeasurements.ZTFToCoord   = MeshMeasurements.ZCoord;
        PMLMeasurement.depth = 5;
        PMLMeasurement.xmax_negativeside =  5;
        PMLMeasurement.ymax_negativeside =  5;
        PMLMeasurement.zmax_negativeside =  5;
        PMLMeasurement.xmin_positiveside = 29;
        PMLMeasurement.ymin_positiveside = 19;
        PMLMeasurement.zmin_positiveside = 19;


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
        
        XIdx_TFFrom = round(MeshMeasurements.XTFFromCoord/MeshMeasurements.dxCoarse)+1;
        if abs(MeshMeasurements.XTFToCoord-MeshMeasurements.XCoord)<EPSILON
            XIdx_TFTo = round(MeshMeasurements.XTFToCoord  /MeshMeasurements.dxCoarse) +1;
        else
            XIdx_TFTo = round(MeshMeasurements.XTFToCoord  /MeshMeasurements.dxCoarse);
        end
        YIdx_TFFrom = round(MeshMeasurements.YTFFromCoord/MeshMeasurements.dyCoarse)+1;
        if abs(MeshMeasurements.YTFToCoord-MeshMeasurements.YCoord)<EPSILON
            YIdx_TFTo = round(MeshMeasurements.YTFToCoord  /MeshMeasurements.dyCoarse) +1;
        else
            YIdx_TFTo = round(MeshMeasurements.YTFToCoord  /MeshMeasurements.dyCoarse);
        end
        ZIdx_TFFrom = round(MeshMeasurements.ZTFFromCoord/MeshMeasurements.dzCoarse)+1;
        if abs(MeshMeasurements.ZTFToCoord-MeshMeasurements.ZCoord)<EPSILON
            ZIdx_TFTo = round(MeshMeasurements.ZTFToCoord  /MeshMeasurements.dzCoarse) +1;
        else
            ZIdx_TFTo = round(MeshMeasurements.ZTFToCoord  /MeshMeasurements.dzCoarse);
        end

        
        LocalUpdateNum                      = ones(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.LocalGridFineness  = ones(XSize+1,YSize+1,ZSize+1);
        
        LocalUpdateNum(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo)...
            = SubgridUpdateNum * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        MeshMeasurements.LocalGridFineness(XIdx_FineFrom:XIdx_FineTo,YIdx_FineFrom:YIdx_FineTo,ZIdx_FineFrom:ZIdx_FineTo) ...
            = SubgridFineness * ones(XIdx_FineTo-XIdx_FineFrom+1,YIdx_FineTo-YIdx_FineFrom+1,ZIdx_FineTo-ZIdx_FineFrom+1);
        
        MeshMeasurements.isInSFRegion = false(XSize+1,YSize+1,ZSize+1);
        MeshMeasurements.isInSFRegion(XIdx_TFFrom:XIdx_TFTo,YIdx_TFFrom:YIdx_TFTo,ZIdx_TFFrom:ZIdx_TFTo)...
            = true;

        
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
        disp(['Total-Field Region: x = [',num2str(MeshMeasurements.XTFFromCoord),' ',num2str(MeshMeasurements.XTFToCoord),'], ',...
            'y = [',num2str(MeshMeasurements.YTFFromCoord),' ',num2str(MeshMeasurements.YTFToCoord),'], ',...
            'z = [',num2str(MeshMeasurements.ZTFFromCoord),' ',num2str(MeshMeasurements.ZTFToCoord),']. ']);


    otherwise
        warning('Preset undefined.')
end
end