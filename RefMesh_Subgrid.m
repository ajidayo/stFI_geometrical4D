function [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer] = RefMesh_Subgrid(MeshMeasurements,LocalUpdateNum)
global EPSILON
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

disp('RefMesh_Subgrid: Counting the number of elements')
VolPerCoarseGrid = zeros(XSize,YSize,ZSize); 
VolPerXRow      = zeros(YSize,ZSize);
VolPerXYPlane   = zeros(ZSize,1);
VolNum          = 0;
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            VolPerCoarseGrid(XIdx,YIdx,ZIdx) = ...
                (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx))^(3);
            VolPerXRow(YIdx,ZIdx) = VolPerXRow(YIdx,ZIdx) + VolPerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        VolPerXYPlane(ZIdx) = VolPerXYPlane(ZIdx) + VolPerXRow(YIdx,ZIdx);
    end
    VolNum = VolNum + VolPerXYPlane(ZIdx);
end

YZFacePerCoarseGrid = zeros(XSize+1,YSize,ZSize);
YZFacePerYRow    = zeros(ZSize,XSize+1);
YZFacePerYZPlane = zeros(XSize+1,1);
YZFaceNum        = 0;
for XIdx = 1:XSize    
    for ZIdx = 1:ZSize
        for YIdx = 1:YSize
            if XIdx == 1
                AdjCoarseGridIsFiner_XMinusDirec = false; 
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjCoarseGridIsFiner_XMinusDirec = true;
                else
                    AdjCoarseGridIsFiner_XMinusDirec = false;
                end
            end
            if AdjCoarseGridIsFiner_XMinusDirec == false
                YZFacePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)  ).^(3);
            else
                YZFacePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)  ).^(2);
            end
            YZFacePerYRow(ZIdx,XIdx) ...
                = YZFacePerYRow(ZIdx,XIdx) + YZFacePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        YZFacePerYZPlane(XIdx) = YZFacePerYZPlane(XIdx) + YZFacePerYRow(ZIdx,XIdx);
    end
    YZFaceNum = YZFaceNum + YZFacePerYZPlane(XIdx);
end
for XIdx = XSize+1
    for ZIdx = 1:ZSize
        for YIdx = 1:YSize
            YZFacePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XSize,YIdx,ZIdx)  ).^(2);
            YZFacePerYRow(ZIdx,XIdx) ...
                = YZFacePerYRow(ZIdx,XIdx) + YZFacePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        YZFacePerYZPlane(XIdx) = YZFacePerYZPlane(XIdx) + YZFacePerYRow(ZIdx,XIdx);
    end
    YZFaceNum = YZFaceNum + YZFacePerYZPlane(XIdx);
end
ZXFacePerCoarseGrid = zeros(XSize,YSize+1,ZSize);
ZXFacePerZRow    = zeros(XSize,YSize+1);
ZXFacePerZXPlane = zeros(YSize+1,1);
ZXFaceNum        = 0;
for YIdx = 1:YSize
    for XIdx = 1:XSize
        for ZIdx = 1:ZSize
            if YIdx == 1
                AdjCoarseGridIsFiner_YMinusDirec = false; 
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjCoarseGridIsFiner_YMinusDirec = true;
                else
                    AdjCoarseGridIsFiner_YMinusDirec = false;
                end
            end
            if AdjCoarseGridIsFiner_YMinusDirec == false
                ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)  ).^(3);
            else
                ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)  ).^(2);
            end
            ZXFacePerZRow(XIdx,YIdx) = ZXFacePerZRow(XIdx,YIdx) + ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        ZXFacePerZXPlane(YIdx) = ZXFacePerZXPlane(YIdx) + ZXFacePerZRow(XIdx,YIdx);
    end
    ZXFaceNum = ZXFaceNum + ZXFacePerZXPlane(YIdx);
end
for YIdx = YSize+1
    for XIdx = 1:XSize
        for ZIdx = 1:ZSize
            ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx) = ...
                (MeshMeasurements.LocalGridFineness(XIdx,YSize,ZIdx)  ).^(2);
            ZXFacePerZRow(XIdx,YIdx) = ZXFacePerZRow(XIdx,YIdx) + ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        ZXFacePerZXPlane(YIdx) = ZXFacePerZXPlane(YIdx) + ZXFacePerZRow(XIdx,YIdx);
    end
    ZXFaceNum = ZXFaceNum + ZXFacePerZXPlane(YIdx);
end
XYFacePerCoarseGrid = zeros(XSize,YSize,ZSize+1);
XYFacePerXRow    = zeros(YSize,ZSize+1);
XYFacePerXYPlane = zeros(ZSize+1,1);
XYFaceNum        = 0;
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            if ZIdx == 1
                AdjCoarseGridIsFiner_ZMinusDirec = false; 
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjCoarseGridIsFiner_ZMinusDirec = true;
                else
                    AdjCoarseGridIsFiner_ZMinusDirec = false;
                end
            end
            if AdjCoarseGridIsFiner_ZMinusDirec == false
                XYFacePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )  ).^(3);
            else
                XYFacePerCoarseGrid(XIdx,YIdx,ZIdx) = ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)-1)...
                    *(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)  ).^(2)...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)  ).^(2);
            end
            XYFacePerXRow(YIdx,ZIdx) = XYFacePerXRow(YIdx,ZIdx) + XYFacePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        XYFacePerXYPlane(ZIdx) = XYFacePerXYPlane(ZIdx) + XYFacePerXRow(YIdx,ZIdx);
    end
    XYFaceNum = XYFaceNum + XYFacePerXYPlane(ZIdx);
end
for ZIdx = ZSize+1
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            XYFacePerCoarseGrid(XIdx,YIdx,ZIdx) = ...
                +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZSize)  ).^(2);
            XYFacePerXRow(YIdx,ZIdx) = XYFacePerXRow(YIdx,ZIdx) + XYFacePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        XYFacePerXYPlane(ZIdx) = XYFacePerXYPlane(ZIdx) + XYFacePerXRow(YIdx,ZIdx);
    end
    XYFaceNum = XYFaceNum + XYFacePerXYPlane(ZIdx);
end

XEdgePerCoarseGrid = zeros(XSize,YSize+1,ZSize+1);
XEdgePerYRow    = zeros(ZSize+1,XSize);
XEdgePerYZPlane = zeros(XSize,1);
XEdgeNum        = 0;
for XIdx = 1:XSize
    for ZIdx = 1:ZSize
        for YIdx = 1:YSize
            if 1 == YIdx
                AdjVolIsFiner_YMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YMinusDirec = true;
                else
                    AdjVolIsFiner_YMinusDirec = false;
                end
            end
            if 1 == ZIdx
                AdjVolIsFiner_ZMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_ZMinusDirec = true;
                else
                    AdjVolIsFiner_ZMinusDirec = false;
                end
            end
            if YIdx ~= 1 && ZIdx ~= 1
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YZMinusDirec = true;
                else
                    AdjVolIsFiner_YZMinusDirec = false;
                end
            else
                AdjVolIsFiner_YZMinusDirec = false;
            end
            if      AdjVolIsFiner_YMinusDirec &&  AdjVolIsFiner_ZMinusDirec
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    =   (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx  )  ).^(3)...
                    - 2*(MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx  )  ).^(2)...
                    +   (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx  )  ).^(2)...
                    +   (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx-1)  ).^(2);
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)>= MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                        - MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
                else
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                        - MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);
                end
            elseif  AdjVolIsFiner_YMinusDirec && ~AdjVolIsFiner_ZMinusDirec
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)  ).^(2);
            elseif ~AdjVolIsFiner_YMinusDirec &&  AdjVolIsFiner_ZMinusDirec
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )  ).^(2);
            elseif ~AdjVolIsFiner_YMinusDirec && ~AdjVolIsFiner_ZMinusDirec
                if AdjVolIsFiner_YZMinusDirec
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        =   (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx  )  ).^(3)...
                        -   (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx  )  ).^(1)...
                        +   (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)  ).^(1);
                else
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)).^(3);
                end
            end
            XEdgePerYRow(ZIdx,XIdx) = XEdgePerYRow(ZIdx,XIdx) + XEdgePerCoarseGrid(XIdx,YIdx,ZIdx);      
        end
        for YIdx = YSize+1
            if ZIdx == 1
                AdjVolIsFiner_YZMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YZMinusDirec = true;
                else
                    AdjVolIsFiner_YZMinusDirec = false;
                end
            end
            XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)-1)...
                *(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx));
            if AdjVolIsFiner_YZMinusDirec
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1);
            else
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);        
            end
            XEdgePerYRow(ZIdx,XIdx) = XEdgePerYRow(ZIdx,XIdx) + XEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        XEdgePerYZPlane(XIdx) = XEdgePerYZPlane(XIdx) + XEdgePerYRow(ZIdx,XIdx);
    end
    for ZIdx = ZSize+1
        for YIdx = 1:YSize
            if YIdx ==1
                AdjVolIsFiner_YZMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YZMinusDirec = true;
                else
                    AdjVolIsFiner_YZMinusDirec = false;
                end
            end
            XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)-1)...
                *(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1));
            if AdjVolIsFiner_YZMinusDirec
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1);
            else
                XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    XEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);        
            end
            XEdgePerYRow(ZIdx,XIdx) = XEdgePerYRow(ZIdx,XIdx) + XEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for YIdx = YSize+1
            XEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1));
            XEdgePerYRow(ZIdx,XIdx) = XEdgePerYRow(ZIdx,XIdx) + XEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        XEdgePerYZPlane(XIdx) = XEdgePerYZPlane(XIdx) + XEdgePerYRow(ZIdx,XIdx);
    end
    XEdgeNum  = XEdgeNum  + XEdgePerYZPlane(XIdx);
end
YEdgePerCoarseGrid = zeros(XSize+1,YSize,ZSize+1);
YEdgePerZRow    = zeros(XSize+1,YSize);
YEdgePerZXPlane = zeros(YSize,1);
YEdgeNum        = 0;
for YIdx = 1:YSize
    for XIdx = 1:XSize
        for ZIdx = 1:ZSize
            if 1 == ZIdx
                AdjVolIsFiner_ZMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_ZMinusDirec = true;
                else
                    AdjVolIsFiner_ZMinusDirec = false;
                end
            end
            if 1 == XIdx
                AdjVolIsFiner_XMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_XMinusDirec = true;
                else
                    AdjVolIsFiner_XMinusDirec = false;
                end
            end
            if ZIdx ~= 1 && XIdx ~= 1
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_ZXMinusDirec = true;
                else
                    AdjVolIsFiner_ZXMinusDirec = false;
                end
            else
                AdjVolIsFiner_ZXMinusDirec = false;
            end
            if      AdjVolIsFiner_ZMinusDirec &&  AdjVolIsFiner_XMinusDirec
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx  )  ).^(3)...
                    - 2*(MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx  )  ).^(2)...
                    +   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx-1)  ).^(2)...
                    +   (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx  )  ).^(2);
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)>= MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        - MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
                else
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        - MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
                end
            elseif  AdjVolIsFiner_ZMinusDirec && ~AdjVolIsFiner_XMinusDirec
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )  ).^(2);
            elseif ~AdjVolIsFiner_ZMinusDirec &&  AdjVolIsFiner_XMinusDirec
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)  ).^(2);
            elseif ~AdjVolIsFiner_ZMinusDirec && ~AdjVolIsFiner_XMinusDirec
                if AdjVolIsFiner_ZXMinusDirec
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        =   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx  )  ).^(3)...
                        -   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx  )  ).^(1)...
                        +   (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)  ).^(1);
                else
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)).^(3);
                end
            end
            YEdgePerZRow(XIdx,YIdx) = YEdgePerZRow(XIdx,YIdx) ...
                + YEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for ZIdx = ZSize+1
            if XIdx == 1
                AdjVolIsFiner_ZXMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_ZXMinusDirec = true;
                else
                    AdjVolIsFiner_ZXMinusDirec = false;
                end
            end
            YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)-1)...
                *(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1));
            if AdjVolIsFiner_ZXMinusDirec
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1);
            else
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
            end
            YEdgePerZRow(XIdx,YIdx) = YEdgePerZRow(XIdx,YIdx) ...
                + YEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        YEdgePerZXPlane(YIdx) = YEdgePerZXPlane(YIdx) + YEdgePerZRow(XIdx,YIdx);
    end
    for XIdx = XSize+1
        for ZIdx = 1:ZSize
            if ZIdx == 1
                AdjVolIsFiner_ZXMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_ZXMinusDirec = true;
                else
                    AdjVolIsFiner_ZXMinusDirec = false;
                end
            end
            YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)-1)...
                *(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx));
            if AdjVolIsFiner_ZXMinusDirec
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1);
            else
                YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    YEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);        
            end
            YEdgePerZRow(XIdx,YIdx) = YEdgePerZRow(XIdx,YIdx) ...
                + YEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for ZIdx = ZSize+1
            YEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1));
            YEdgePerZRow(XIdx,YIdx) = YEdgePerZRow(XIdx,YIdx) ...
                + YEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        YEdgePerZXPlane(YIdx) = YEdgePerZXPlane(YIdx) + YEdgePerZRow(XIdx,YIdx);
    end
    YEdgeNum  = YEdgeNum  + YEdgePerZXPlane(YIdx);
end
ZEdgePerCoarseGrid = zeros(XSize+1,YSize+1,ZSize);
ZEdgePerXRow    = zeros(YSize+1,ZSize);
ZEdgePerXYPlane = zeros(ZSize,1);
ZEdgeNum        = 0;
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            if 1 == XIdx
                AdjVolIsFiner_XMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_XMinusDirec = true;
                else
                    AdjVolIsFiner_XMinusDirec = false;
                end
            end
            if 1 == YIdx
                AdjVolIsFiner_YMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YMinusDirec = true;
                else
                    AdjVolIsFiner_YMinusDirec = false;
                end
            end
            if XIdx ~= 1 && YIdx ~= 1
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_XYMinusDirec = true;
                else
                    AdjVolIsFiner_XYMinusDirec = false;
                end
            else
                AdjVolIsFiner_XYMinusDirec = false;
            end
            if      AdjVolIsFiner_XMinusDirec &&  AdjVolIsFiner_YMinusDirec
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx  ,ZIdx)  ).^(3)...
                    - 2*(MeshMeasurements.LocalGridFineness(XIdx  ,YIdx  ,ZIdx)  ).^(2)...
                    +   (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx  ,ZIdx)  ).^(2)...
                    +   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx-1,ZIdx)  ).^(2);
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)>= MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        - MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);
                else
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        - MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
                end
            elseif  AdjVolIsFiner_XMinusDirec && ~AdjVolIsFiner_YMinusDirec
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)  ).^(2);
            elseif ~AdjVolIsFiner_XMinusDirec &&  AdjVolIsFiner_YMinusDirec
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    = (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)  ).^(2)...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)-1)...
                    * (MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)  ).^(2);
            elseif ~AdjVolIsFiner_XMinusDirec && ~AdjVolIsFiner_YMinusDirec
                if AdjVolIsFiner_XYMinusDirec
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        =   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx  ,ZIdx)  ).^(3)...
                        -   (MeshMeasurements.LocalGridFineness(XIdx  ,YIdx  ,ZIdx)  ).^(1)...
                        +   (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)  ).^(1);
                else
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        = (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)).^(3);
                end
            end
            ZEdgePerXRow(YIdx,ZIdx) = ZEdgePerXRow(YIdx,ZIdx) + ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for XIdx = XSize+1
            if YIdx ==1
                AdjVolIsFiner_XYMinusDirec = false;
            else
            if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_XYMinusDirec = true;
            else
                AdjVolIsFiner_XYMinusDirec = false;
            end
            end
            ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)-1)...
                *(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx));
            if AdjVolIsFiner_XYMinusDirec
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx);
            else
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
            end
            ZEdgePerXRow(YIdx,ZIdx) = ZEdgePerXRow(YIdx,ZIdx) + ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        ZEdgePerXYPlane(ZIdx) = ZEdgePerXYPlane(ZIdx) + ZEdgePerXRow(YIdx,ZIdx);
    end
    for YIdx = YSize+1
        for XIdx = 1:XSize
            if XIdx ==1
                AdjVolIsFiner_XYMinusDirec = false;  
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_XYMinusDirec = true;
                else
                    AdjVolIsFiner_XYMinusDirec = false;
                end
            end
            ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)-1)...
                *(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx));
            if AdjVolIsFiner_XYMinusDirec
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx);
            else
                ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)=...
                    ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx)...
                    +MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);
            end
            ZEdgePerXRow(YIdx,ZIdx) = ZEdgePerXRow(YIdx,ZIdx) + ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for XIdx = XSize+1
            ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx) = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx);
            ZEdgePerXRow(YIdx,ZIdx) = ZEdgePerXRow(YIdx,ZIdx) + ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        ZEdgePerXYPlane(ZIdx) = ZEdgePerXYPlane(ZIdx) + ZEdgePerXRow(YIdx,ZIdx);
    end
    ZEdgeNum  = ZEdgeNum  + ZEdgePerXYPlane(ZIdx);
end

NodePerCoarseGrid = zeros(XSize+1,YSize+1,ZSize+1);
NodePerXRow    = zeros(YSize+1,ZSize+1);
NodePerXYPlane = zeros(ZSize+1,1);
NodeNum        = 0;
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            if XIdx == 1
                AdjVolIsFiner_XMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_XMinusDirec = true;
                else
                    AdjVolIsFiner_XMinusDirec = false;
                end
            end
            if YIdx == 1
                AdjVolIsFiner_YMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YMinusDirec = true;
                else
                    AdjVolIsFiner_YMinusDirec = false;
                end
            end
            if ZIdx == 1
                AdjVolIsFiner_ZMinusDirec = false;
            else
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_ZMinusDirec = true;
                else
                    AdjVolIsFiner_ZMinusDirec = false;
                end
            end
            if XIdx ~= 1 && YIdx ~= 1
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_XYMinusDirec = true;
                else
                    AdjVolIsFiner_XYMinusDirec = false;
                end
            else
                AdjVolIsFiner_XYMinusDirec = false;
            end
            if YIdx ~= 1 && ZIdx ~= 1
                if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                    AdjVolIsFiner_YZMinusDirec = true;
                else
                    AdjVolIsFiner_YZMinusDirec = false;
                end
            else
                AdjVolIsFiner_YZMinusDirec = false;
            end
            if ZIdx ~= 1 && XIdx ~= 1
                if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx) 
                    AdjVolIsFiner_ZXMinusDirec = true;
                else
                    AdjVolIsFiner_ZXMinusDirec = false;
                end
            else
                AdjVolIsFiner_ZXMinusDirec = false;
            end
            
            % Compute NodePerXRow 
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = 1 + (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )-1).^(3);
            if     ~AdjVolIsFiner_YMinusDirec && ~AdjVolIsFiner_ZMinusDirec
                if  ~AdjVolIsFiner_YZMinusDirec
                    NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)-1);
                else
                    NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        +(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)-1);
                end
            elseif  AdjVolIsFiner_YMinusDirec && ~AdjVolIsFiner_ZMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx  )-1);
            elseif ~AdjVolIsFiner_YMinusDirec &&  AdjVolIsFiner_ZMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx-1)-1);
            elseif  AdjVolIsFiner_YMinusDirec &&  AdjVolIsFiner_ZMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)-1);
            end
            if     ~AdjVolIsFiner_ZMinusDirec && ~AdjVolIsFiner_XMinusDirec
                if  ~AdjVolIsFiner_ZXMinusDirec
                    NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)-1);
                else
                    NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)-1);
                end
            elseif  AdjVolIsFiner_ZMinusDirec && ~AdjVolIsFiner_XMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx-1)-1);
            elseif ~AdjVolIsFiner_ZMinusDirec &&  AdjVolIsFiner_XMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx  )-1);
            elseif  AdjVolIsFiner_ZMinusDirec &&  AdjVolIsFiner_XMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)-1);
            end
            if     ~AdjVolIsFiner_XMinusDirec && ~AdjVolIsFiner_YMinusDirec
                if  ~AdjVolIsFiner_XYMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)-1);
                else
                    NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                        +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)-1);
                end
            elseif  AdjVolIsFiner_XMinusDirec && ~AdjVolIsFiner_YMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx  ,ZIdx)-1);
            elseif ~AdjVolIsFiner_XMinusDirec &&  AdjVolIsFiner_YMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx  ,YIdx-1,ZIdx)-1);
            elseif  AdjVolIsFiner_XMinusDirec &&  AdjVolIsFiner_YMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)-1);
            end
            if AdjVolIsFiner_XMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)-1).^(2);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx)-1).^(2);
            end
            if AdjVolIsFiner_YMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)-1).^(2);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx)-1).^(2);
            end
            if AdjVolIsFiner_ZMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)-1).^(2);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) = NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    +(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx  )-1).^(2);
            end
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for XIdx = XSize+1
            if MeshMeasurements.LocalGridFineness(XIdx-1,max([YIdx-1 1]),ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_XYMinusDirec = true;
            else
                AdjVolIsFiner_XYMinusDirec = false;
            end
            if MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,max([ZIdx-1 1]))>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_ZXMinusDirec = true;
            else
                AdjVolIsFiner_ZXMinusDirec = false;
            end
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = 1+(MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)-1)^2;
            if AdjVolIsFiner_XYMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)-1);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)-1);
            end
            if AdjVolIsFiner_ZXMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)-1);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx)-1);
            end
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        NodePerXYPlane(ZIdx) = NodePerXYPlane(ZIdx) + NodePerXRow(YIdx,ZIdx);
    end
    for YIdx = YSize+1
        for XIdx = 1:XSize
            if MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,max([ZIdx-1 1]))>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_YZMinusDirec = true;
            else
                AdjVolIsFiner_YZMinusDirec = false;
            end            
            if MeshMeasurements.LocalGridFineness(max([XIdx-1 1]),YIdx-1,ZIdx)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_XYMinusDirec = true;
            else
                AdjVolIsFiner_XYMinusDirec = false;
            end
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = 1+(MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)-1)^2;
            if AdjVolIsFiner_YZMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)-1);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     + (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)-1);
            end
            if AdjVolIsFiner_XYMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx)-1);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                     + (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx)-1);
            end
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for XIdx = XSize+1
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx);
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        NodePerXYPlane(ZIdx) = NodePerXYPlane(ZIdx) + NodePerXRow(YIdx,ZIdx);
    end
    NodeNum = NodeNum + NodePerXYPlane(ZIdx);
end
for ZIdx = ZSize+1
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            if MeshMeasurements.LocalGridFineness(max([XIdx-1 1]),YIdx,ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_ZXMinusDirec = true;
            else
                AdjVolIsFiner_ZXMinusDirec = false;
            end            
            if MeshMeasurements.LocalGridFineness(XIdx,max([YIdx-1 1]),ZIdx-1)>MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx)
                AdjVolIsFiner_YZMinusDirec = true;
            else
                AdjVolIsFiner_YZMinusDirec = false;
            end
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = 1+(MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)-1)^2;
            if AdjVolIsFiner_ZXMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1)-1);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)-1);
            end
            if AdjVolIsFiner_YZMinusDirec
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1)-1);
            else
                NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    =  NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                    + (MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1)-1);
            end
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for XIdx = XSize+1
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1);
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        NodePerXYPlane(ZIdx) = NodePerXYPlane(ZIdx) + NodePerXRow(YIdx,ZIdx);
    end
    for YIdx = YSize+1
        for XIdx = 1:XSize
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1);
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        for XIdx = XSize+1
            NodePerCoarseGrid(XIdx,YIdx,ZIdx) ...
                = 1;
            NodePerXRow(YIdx,ZIdx) = NodePerXRow(YIdx,ZIdx) + NodePerCoarseGrid(XIdx,YIdx,ZIdx);
        end
        NodePerXYPlane(ZIdx) = NodePerXYPlane(ZIdx) + NodePerXRow(YIdx,ZIdx);
    end
    NodeNum = NodeNum + NodePerXYPlane(ZIdx);
end

ElemPer.VolPerCoarseGrid    = VolPerCoarseGrid;
ElemPer.VolPerXRow          = VolPerXRow;
ElemPer.VolPerXYPlane       = VolPerXYPlane;
ElemPer.YZFacePerCoarseGrid = YZFacePerCoarseGrid;
ElemPer.YZFacePerYRow       = YZFacePerYRow;
ElemPer.YZFacePerYZPlane    = YZFacePerYZPlane;
ElemPer.YZFaceNum           = YZFaceNum;
ElemPer.ZXFacePerCoarseGrid = ZXFacePerCoarseGrid;
ElemPer.ZXFacePerZRow       = ZXFacePerZRow;
ElemPer.ZXFacePerZXPlane    = ZXFacePerZXPlane;
ElemPer.ZXFaceNum           = ZXFaceNum;
ElemPer.XYFacePerCoarseGrid = XYFacePerCoarseGrid;
ElemPer.XYFacePerXRow       = XYFacePerXRow;
ElemPer.XYFacePerXYPlane    = XYFacePerXYPlane;
ElemPer.XYFaceNum           = XYFaceNum;

ElemPer.XEdgePerCoarseGrid  = XEdgePerCoarseGrid;
ElemPer.XEdgePerYRow        = XEdgePerYRow;
ElemPer.XEdgePerYZPlane     = XEdgePerYZPlane;
ElemPer.XEdgeNum            = XEdgeNum;
ElemPer.YEdgePerCoarseGrid  = YEdgePerCoarseGrid;
ElemPer.YEdgePerZRow        = YEdgePerZRow;
ElemPer.YEdgePerZXPlane     = YEdgePerZXPlane;
ElemPer.YEdgeNum            = YEdgeNum;
ElemPer.ZEdgePerCoarseGrid  = ZEdgePerCoarseGrid;
ElemPer.ZEdgePerXRow        = ZEdgePerXRow;
ElemPer.ZEdgePerXYPlane     = ZEdgePerXYPlane;
ElemPer.ZEdgeNum            = ZEdgeNum;
ElemPer.NodePerCoarseGrid   = NodePerCoarseGrid;
ElemPer.NodePerXRow         = NodePerXRow;
ElemPer.NodePerXYPlane      = NodePerXYPlane;

Num_of_Elem.SpV = VolNum;
Num_of_Elem.SpP = YZFaceNum + ZXFaceNum + XYFaceNum;
Num_of_Elem.SpS = XEdgeNum  + YEdgeNum  + ZEdgeNum;
Num_of_Elem.SpN = NodeNum;


[SpElemPositionIdx,ElemFineness] = CalcSpElemPositionIdx(MeshMeasurements,ElemPer,Num_of_Elem);

if any(ElemFineness.SpV==0) || any(ElemFineness.SpP==0) ...
        || any(ElemFineness.SpS==0)
    disp('Fineness undefined for some elements')
end
disp('RefMesh_Subgrid:Constructing DivMatrix sD')
sD = Subgrid_sD(SpElemPositionIdx,ElemFineness,ElemPer);
disp('RefMesh_Subgrid:Constructing CurlMatrix sC')
sC = Subgrid_sC(SpElemPositionIdx,ElemFineness,ElemPer);
disp('RefMesh_Subgrid:Constructing GradMatrix sG')
sG = Subgrid_sG(SpElemPositionIdx,ElemFineness,ElemPer);
MeshEndsAt = [MeshMeasurements.XCoord;MeshMeasurements.YCoord;MeshMeasurements.ZCoord];
SpDiscreteWidth_Coarse = [MeshMeasurements.dxCoarse;MeshMeasurements.dyCoarse;MeshMeasurements.dzCoarse];
SubgridsStartsFrom = [MeshMeasurements.XFineFromCoord;MeshMeasurements.YFineFromCoord;MeshMeasurements.ZFineFromCoord];
SubgridsStartsFrom(logical(SubgridsStartsFrom==0)) = -Inf;
SubgridsEndsAt = [MeshMeasurements.XFineToCoord;MeshMeasurements.YFineToCoord;MeshMeasurements.ZFineToCoord];
SubgridsEndsAt(logical(abs(SubgridsEndsAt-MeshEndsAt))<EPSILON) = Inf;
NonZeroEntityNum=0;
for SpNIdx = 1:Num_of_Elem.SpN
    if size(find(any([logical(abs(SpDiscreteWidth_Coarse.*SpElemPositionIdx.SpN(:,SpNIdx)-SubgridsStartsFrom)<EPSILON);...
        logical(abs(SpDiscreteWidth_Coarse.*SpElemPositionIdx.SpN(:,SpNIdx)-SubgridsEndsAt)<EPSILON)],2)),1)...
        >=2
    disp(['RefMesh_Subgrid: SpNIdx = ', num2str(SpNIdx),' is judged as being on the edge/corner of the Subgrid region'])
    NonZeroEntityNum = NonZeroEntityNum + 1;
    isOnFineGridCorner_NoneZeroRow(NonZeroEntityNum) = 1;
    isOnFineGridCorner_NoneZeroCol(NonZeroEntityNum) = SpNIdx;
    isOnFineGridCorner_NoneZeroVal(NonZeroEntityNum) = true;
    end
end
if NonZeroEntityNum ~= 0 &&  max(MeshMeasurements.LocalGridFineness,[],'all')>1
    isOnFineGridCorner = sparse(isOnFineGridCorner_NoneZeroRow,isOnFineGridCorner_NoneZeroCol,isOnFineGridCorner_NoneZeroVal,1,Num_of_Elem.SpN);
    if max(LocalUpdateNum,[],'all') == max(MeshMeasurements.LocalGridFineness,[],'all')
        SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion = isOnFineGridCorner;
    else
        SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion = sparse(false(1,Num_of_Elem.SpN));
    end
elseif NonZeroEntityNum ~= 0 && max(LocalUpdateNum,[],'all') >1
    isOnFineGridCorner =  sparse(false(1,Num_of_Elem.SpN));
    SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion ...
        = sparse(isOnFineGridCorner_NoneZeroRow,isOnFineGridCorner_NoneZeroCol,isOnFineGridCorner_NoneZeroVal,1,Num_of_Elem.SpN);
else
    isOnFineGridCorner =  sparse(false(1,Num_of_Elem.SpN));
    SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion = sparse(false(1,Num_of_Elem.SpN));
end

disp('RefMesh_Subgrid:Calculating NodePositions')
for SpVIdx = 1:Num_of_Elem.SpV
    NodePos.Dual(SpVIdx).Vec = ...
        [MeshMeasurements.dxCoarse;MeshMeasurements.dyCoarse;MeshMeasurements.dzCoarse] .* SpElemPositionIdx.SpV(:,SpVIdx);
end
NodePos.Prim = Subgrid_NodePos_Prim(SpElemPositionIdx,ElemPer,ElemFineness,isOnFineGridCorner,sG,sC,sD,MeshMeasurements);

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            VolIdxList = ...
                (1:ElemPer.VolPerCoarseGrid(XIdx,YIdx,ZIdx))...
                +sum(ElemPer.VolPerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
                +sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))...
                +sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
            SpElemProperties.SpV.UpdNum(VolIdxList) = LocalUpdateNum(XIdx,YIdx,ZIdx);
        end
    end
end

disp('RefMesh_Subgrid:Setting Boundary Condition: Perfect Electric Wall')
SpElemProperties.SpP.ElecWall = sparse(false(1,Num_of_Elem.SpP));
MeshEndsAt = [MeshMeasurements.XCoord;MeshMeasurements.YCoord;MeshMeasurements.ZCoord];    
for SpPIdx = [...
        (1:ElemPer.YZFacePerYZPlane(1)) ...
        (sum(ElemPer.YZFacePerYZPlane(1:XSize))+1:ElemPer.YZFaceNum) ...
        ElemPer.YZFaceNum + (1:ElemPer.ZXFacePerZXPlane(1)) ...
        ElemPer.YZFaceNum + (sum(ElemPer.ZXFacePerZXPlane(1:YSize))+1:ElemPer.ZXFaceNum) ...
        ElemPer.YZFaceNum + ElemPer.ZXFaceNum + (1:ElemPer.XYFacePerXYPlane(1)) ...
        ElemPer.YZFaceNum + ElemPer.ZXFaceNum + (sum(ElemPer.XYFacePerXYPlane(1:ZSize))+1:ElemPer.XYFaceNum)]
    NodeLogIdxs_DefFace = logical(logical(sC(SpPIdx,:))*logical(sG));
    MeanPosOfFace = mean([NodePos.Prim(NodeLogIdxs_DefFace).Vec],2);
    if size(find(any([logical(abs(MeanPosOfFace-[0;0;0])<EPSILON);...
            logical(abs(MeanPosOfFace-MeshEndsAt)<EPSILON)],2)),1)...
            >=1
        SpElemProperties.SpP.ElecWall(SpPIdx) = true;
    end
end
SpElemProperties.SpS.PEC = sparse(false(1,Num_of_Elem.SpS));
for SpPIdx = find(SpElemProperties.SpP.ElecWall)
    for IncSpSIdx = find(sC(SpPIdx,:))
        SpElemProperties.SpS.PEC(IncSpSIdx) = true;
    end
end

disp('RefMesh_Subgrid:Pre-calculating the area of the faces in the FDTD-like region.')
SpElemProperties.SpP.PrimAreaIsGiven = true(1,Num_of_Elem.SpP);
SpElemProperties.SpP.DualLengIsGiven = true(1,Num_of_Elem.SpP);
SpElemProperties.SpS.PrimLengIsGiven = true(1,Num_of_Elem.SpS);
SpElemProperties.SpS.DualAreaIsGiven = true(1,Num_of_Elem.SpS);
for SpPIdx = 1:Num_of_Elem.SpP
    IncSpVIdx = find(sD(:,SpPIdx).');
    if SpElemProperties.SpP.ElecWall(SpPIdx) == false
        if ElemFineness.SpV(IncSpVIdx(1)) ~= ElemFineness.SpV(IncSpVIdx(2))
            SpPLogIdx_tgt = logical(logical(sD(:,SpPIdx).')*logical(sD));
            SpSLogIdx_tgt = logical(sC(SpPLogIdx_tgt,:));
            SpElemProperties.SpP.PrimAreaIsGiven(SpPLogIdx_tgt) = false;
            SpElemProperties.SpP.DualLengIsGiven(SpPLogIdx_tgt) = false;
            SpElemProperties.SpS.PrimLengIsGiven(SpSLogIdx_tgt) = false;
            SpElemProperties.SpS.DualAreaIsGiven(SpSLogIdx_tgt) = false;
        end
    end
end

SpElemProperties.SpP.PrimArea = zeros(Num_of_Elem.SpP,1);
SpElemProperties.SpP.DualLeng = zeros(Num_of_Elem.SpP,1);
SpElemProperties.SpS.PrimLeng = zeros(Num_of_Elem.SpS,1);
SpElemProperties.SpS.DualArea = zeros(Num_of_Elem.SpS,1);
StandardArea = [...
    MeshMeasurements.dyCoarse*MeshMeasurements.dzCoarse;...
    MeshMeasurements.dzCoarse*MeshMeasurements.dxCoarse;...
    MeshMeasurements.dxCoarse*MeshMeasurements.dyCoarse];
StandardLength = ...
    [MeshMeasurements.dxCoarse;MeshMeasurements.dyCoarse;MeshMeasurements.dzCoarse];
for SpPIdx = find(SpElemProperties.SpP.PrimAreaIsGiven)
    PrimFaceDirec = SpElemPositionIdx.SpV*sD(:,SpPIdx);
    PrimFaceDirec = norm(PrimFaceDirec).^(-1)*PrimFaceDirec;
    SpElemProperties.SpP.PrimArea(SpPIdx) = ElemFineness.SpP(SpPIdx)^(-2)*dot(PrimFaceDirec.',StandardArea.*PrimFaceDirec);
end
for SpPIdx = find(SpElemProperties.SpP.DualLengIsGiven)
    DualEdge = StandardLength.*(SpElemPositionIdx.SpV*sD(:,SpPIdx));
    SpElemProperties.SpP.DualLeng(SpPIdx) = norm(DualEdge);
end
for SpSIdx = find(SpElemProperties.SpS.PrimLengIsGiven)
    PrimEdge = StandardLength.*(SpElemPositionIdx.SpN*sG(SpSIdx,:).');
    SpElemProperties.SpS.PrimLeng(SpSIdx) = norm(PrimEdge);
end
for SpSIdx = find(SpElemProperties.SpS.DualAreaIsGiven)
    DualFaceDirec = SpElemPositionIdx.SpN*sG(SpSIdx,:).';
    DualFaceDirec = norm(DualFaceDirec).^(-1)*DualFaceDirec;
    SpElemProperties.SpS.DualArea(SpSIdx) = ElemFineness.SpS(SpSIdx)^(-2)*dot(DualFaceDirec,StandardArea.*DualFaceDirec);
end
end
%%

function sD = Subgrid_sD(SpElemPositionIdx,ElemFineness,ElemPer)
global EPSILON
sD_Num_Nonzeros=0;
sD_Nonzeros_Row = 0;
sD_Nonzeros_Col = 0;
sD_Nonzeros_Val = 0;

for SpVIdx = 1:size(ElemFineness.SpV,2)
    XIdx = ceil(SpElemPositionIdx.SpV(1,SpVIdx));
    YIdx = ceil(SpElemPositionIdx.SpV(2,SpVIdx));
    ZIdx = ceil(SpElemPositionIdx.SpV(3,SpVIdx));
    
    IncYZSpPIdxs_Guess = [((1:ElemPer.YZFacePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.YZFacePerCoarseGrid(XIdx,1:YIdx-1,ZIdx))...
        + sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))...
        + sum(ElemPer.YZFacePerYZPlane(1:XIdx-1)))...
        ((1:ElemPer.YZFacePerCoarseGrid(XIdx+1,YIdx,ZIdx))...
        + sum(ElemPer.YZFacePerCoarseGrid(XIdx+1,1:YIdx-1,ZIdx))...
        + sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx+1))...
        + sum(ElemPer.YZFacePerYZPlane(1:XIdx  )))...
        ];
    for SpPIdx = IncYZSpPIdxs_Guess
        if ( SpElemPositionIdx.SpV(2,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - EPSILON < SpElemPositionIdx.SpP(2,SpPIdx)...
                && SpElemPositionIdx.SpV(2,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) + EPSILON > SpElemPositionIdx.SpP(2,SpPIdx))...
                &&( SpElemPositionIdx.SpV(3,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - EPSILON < SpElemPositionIdx.SpP(3,SpPIdx)...
                && SpElemPositionIdx.SpV(3,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) + EPSILON > SpElemPositionIdx.SpP(3,SpPIdx))
            if abs(SpElemPositionIdx.SpV(1,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - SpElemPositionIdx.SpP(1,SpPIdx) ) < EPSILON
                sD_Num_Nonzeros = sD_Num_Nonzeros +1;
                sD_Nonzeros_Row(sD_Num_Nonzeros) = SpVIdx;
                sD_Nonzeros_Col(sD_Num_Nonzeros) = SpPIdx;
                sD_Nonzeros_Val(sD_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpV(1,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) - SpElemPositionIdx.SpP(1,SpPIdx) ) < EPSILON
                sD_Num_Nonzeros = sD_Num_Nonzeros +1;
                sD_Nonzeros_Row(sD_Num_Nonzeros) = SpVIdx;
                sD_Nonzeros_Col(sD_Num_Nonzeros) = SpPIdx;
                sD_Nonzeros_Val(sD_Num_Nonzeros) =  1;
            end
        end
    end
    
    IncZXSpPIdxs_Guess = [...
        ((1:ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,1:ZIdx-1))...
        + sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx))...
        + sum(ElemPer.ZXFacePerZXPlane(1:YIdx-1))...
        +ElemPer.YZFaceNum)...
        ((1:ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx+1,ZIdx))...
        + sum(ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx+1,1:ZIdx-1))...
        + sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx+1))...
        + sum(ElemPer.ZXFacePerZXPlane(1:YIdx  ))...
        +ElemPer.YZFaceNum)...
        ];
    for SpPIdx = IncZXSpPIdxs_Guess
        if ( SpElemPositionIdx.SpV(3,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - EPSILON < SpElemPositionIdx.SpP(3,SpPIdx)...
                && SpElemPositionIdx.SpV(3,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) + EPSILON > SpElemPositionIdx.SpP(3,SpPIdx))...
                &&( SpElemPositionIdx.SpV(1,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - EPSILON < SpElemPositionIdx.SpP(1,SpPIdx)...
                && SpElemPositionIdx.SpV(1,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) + EPSILON > SpElemPositionIdx.SpP(1,SpPIdx))
            if abs(SpElemPositionIdx.SpV(2,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - SpElemPositionIdx.SpP(2,SpPIdx) ) < EPSILON
                sD_Num_Nonzeros = sD_Num_Nonzeros +1;
                sD_Nonzeros_Row(sD_Num_Nonzeros) = SpVIdx;
                sD_Nonzeros_Col(sD_Num_Nonzeros) = SpPIdx;
                sD_Nonzeros_Val(sD_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpV(2,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) - SpElemPositionIdx.SpP(2,SpPIdx) ) < EPSILON
                sD_Num_Nonzeros = sD_Num_Nonzeros +1;
                sD_Nonzeros_Row(sD_Num_Nonzeros) = SpVIdx;
                sD_Nonzeros_Col(sD_Num_Nonzeros) = SpPIdx;
                sD_Nonzeros_Val(sD_Num_Nonzeros) =  1;
            end
        end
    end
    
    IncXYSpPIdxs_Guess = [...
        ((1:ElemPer.XYFacePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.XYFacePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
        + sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1))...
        +ElemPer.YZFaceNum+ElemPer.ZXFaceNum)...
        ((1:ElemPer.XYFacePerCoarseGrid(XIdx,YIdx,ZIdx+1))...
        + sum(ElemPer.XYFacePerCoarseGrid(1:XIdx-1,YIdx,ZIdx+1))...
        + sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx+1))...
        + sum(ElemPer.ZXFacePerZXPlane(1:ZIdx  ))...
        +ElemPer.YZFaceNum+ElemPer.ZXFaceNum)...
        ];
    for SpPIdx = IncXYSpPIdxs_Guess
        if ( SpElemPositionIdx.SpV(1,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - EPSILON < SpElemPositionIdx.SpP(1,SpPIdx)...
                && SpElemPositionIdx.SpV(1,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) + EPSILON > SpElemPositionIdx.SpP(1,SpPIdx))...
                &&( SpElemPositionIdx.SpV(2,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - EPSILON < SpElemPositionIdx.SpP(2,SpPIdx)...
                && SpElemPositionIdx.SpV(2,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) + EPSILON > SpElemPositionIdx.SpP(2,SpPIdx))
            if abs(SpElemPositionIdx.SpV(3,SpVIdx) - 0.5/ElemFineness.SpV(SpVIdx) - SpElemPositionIdx.SpP(3,SpPIdx) ) < EPSILON
                sD_Num_Nonzeros = sD_Num_Nonzeros +1;
                sD_Nonzeros_Row(sD_Num_Nonzeros) = SpVIdx;
                sD_Nonzeros_Col(sD_Num_Nonzeros) = SpPIdx;
                sD_Nonzeros_Val(sD_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpV(3,SpVIdx) + 0.5/ElemFineness.SpV(SpVIdx) - SpElemPositionIdx.SpP(3,SpPIdx) ) < EPSILON
                sD_Num_Nonzeros = sD_Num_Nonzeros +1;
                sD_Nonzeros_Row(sD_Num_Nonzeros) = SpVIdx;
                sD_Nonzeros_Col(sD_Num_Nonzeros) = SpPIdx;
                sD_Nonzeros_Val(sD_Num_Nonzeros) =  1;
            end
        end
    end
    
end

sD = sparse(sD_Nonzeros_Row,sD_Nonzeros_Col,sD_Nonzeros_Val,size(ElemFineness.SpV,2),size(ElemFineness.SpP,2));


end

%%

function sC = Subgrid_sC(SpElemPositionIdx,ElemFineness,ElemPer)
global EPSILON
sC_Num_Nonzeros=0;
sC_Nonzeros_Row = 0;
sC_Nonzeros_Col = 0;
sC_Nonzeros_Val = 0;

for SpPIdx = 1:ElemPer.YZFaceNum
    XIdx = ceil(SpElemPositionIdx.SpP(1,SpPIdx)+EPSILON);
    YIdx = ceil(SpElemPositionIdx.SpP(2,SpPIdx)+EPSILON);
    ZIdx = ceil(SpElemPositionIdx.SpP(3,SpPIdx)+EPSILON);
    
    IncYSpSIdxs_Guess = [((1:ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,1:ZIdx-1))...
        + sum(ElemPer.YEdgePerZRow(1:XIdx-1,YIdx))...
        + sum(ElemPer.YEdgePerZXPlane(1:YIdx-1))...
        + ElemPer.XEdgeNum)...
        ((1:ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,ZIdx+1))...
        + sum(ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,1:ZIdx  ))...
        + sum(ElemPer.YEdgePerZRow(1:XIdx-1,YIdx))...
        + sum(ElemPer.YEdgePerZXPlane(1:YIdx-1))...
        + ElemPer.XEdgeNum)...
        ];
    for SpSIdx = IncYSpSIdxs_Guess
        if ( SpElemPositionIdx.SpP(2,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - EPSILON < SpElemPositionIdx.SpS(2,SpSIdx)...
                && SpElemPositionIdx.SpP(2,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) + EPSILON > SpElemPositionIdx.SpS(2,SpSIdx))...
                && abs(SpElemPositionIdx.SpP(1,SpPIdx) - SpElemPositionIdx.SpS(1,SpSIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpP(3,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(3,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) =  1;
            elseif abs(SpElemPositionIdx.SpP(3,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(3,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) = -1;
            end
        end
    end

    IncZSpSIdxs_Guess = [((1:ElemPer.ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.ZEdgePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
        + sum(ElemPer.ZEdgePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.ZEdgePerXYPlane(1:ZIdx-1))...
        + ElemPer.XEdgeNum+ElemPer.YEdgeNum)...
        ((1:ElemPer.ZEdgePerCoarseGrid(XIdx,YIdx+1,ZIdx))...
        + sum(ElemPer.ZEdgePerCoarseGrid(1:XIdx-1,YIdx+1,ZIdx))...
        + sum(ElemPer.ZEdgePerXRow(1:YIdx  ,ZIdx))...
        + sum(ElemPer.ZEdgePerXYPlane(1:ZIdx-1))...
        + ElemPer.XEdgeNum+ElemPer.YEdgeNum)...
        ];
    for SpSIdx = IncZSpSIdxs_Guess
        if ( SpElemPositionIdx.SpP(3,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - EPSILON < SpElemPositionIdx.SpS(3,SpSIdx)...
                && SpElemPositionIdx.SpP(3,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) + EPSILON > SpElemPositionIdx.SpS(3,SpSIdx))...
                && abs(SpElemPositionIdx.SpP(1,SpPIdx) - SpElemPositionIdx.SpS(1,SpSIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpP(2,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(2,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpP(2,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(2,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) =  1;
            end
        end
    end
    
end


for SpPIdx = ElemPer.YZFaceNum+1:ElemPer.YZFaceNum+ElemPer.ZXFaceNum
    XIdx = ceil(SpElemPositionIdx.SpP(1,SpPIdx)+EPSILON);
    YIdx = ceil(SpElemPositionIdx.SpP(2,SpPIdx)+EPSILON);
    ZIdx = ceil(SpElemPositionIdx.SpP(3,SpPIdx)+EPSILON);
    
    IncZSpSIdxs_Guess = [...
        ((1:ElemPer.ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.ZEdgePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
        + sum(ElemPer.ZEdgePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.ZEdgePerXYPlane(1:ZIdx-1))...
        + ElemPer.XEdgeNum+ElemPer.YEdgeNum)...
        ((1:ElemPer.ZEdgePerCoarseGrid(XIdx+1,YIdx,ZIdx))...
        + sum(ElemPer.ZEdgePerCoarseGrid(1:XIdx  ,YIdx,ZIdx))...
        + sum(ElemPer.ZEdgePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.ZEdgePerXYPlane(1:ZIdx-1))...
        + ElemPer.XEdgeNum+ElemPer.YEdgeNum)...
        ];
    for SpSIdx = IncZSpSIdxs_Guess
        if ( SpElemPositionIdx.SpP(3,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - EPSILON < SpElemPositionIdx.SpS(3,SpSIdx)...
                && SpElemPositionIdx.SpP(3,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) + EPSILON > SpElemPositionIdx.SpS(3,SpSIdx))...
                && abs(SpElemPositionIdx.SpP(2,SpPIdx) - SpElemPositionIdx.SpS(2,SpSIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpP(1,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(1,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) =  1;
            elseif abs(SpElemPositionIdx.SpP(1,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(1,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) = -1;
            end
        end
    end

    IncXSpSIdxs_Guess = [...
        ((1:ElemPer.XEdgePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.XEdgePerCoarseGrid(XIdx,1:YIdx-1,ZIdx))...
        + sum(ElemPer.XEdgePerYRow(1:ZIdx-1,XIdx))...
        + sum(ElemPer.XEdgePerYZPlane(1:XIdx-1))...
        )...
        ((1:ElemPer.XEdgePerCoarseGrid(XIdx,YIdx,ZIdx+1))...
        + sum(ElemPer.XEdgePerCoarseGrid(XIdx,1:YIdx-1,ZIdx+1))...
        + sum(ElemPer.XEdgePerYRow(1:ZIdx  ,XIdx))...
        + sum(ElemPer.XEdgePerYZPlane(1:XIdx-1))...
        )...
        ];
    for SpSIdx = IncXSpSIdxs_Guess
        if ( SpElemPositionIdx.SpP(1,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - EPSILON < SpElemPositionIdx.SpS(1,SpSIdx)...
                && SpElemPositionIdx.SpP(1,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) + EPSILON > SpElemPositionIdx.SpS(1,SpSIdx))...
                && abs(SpElemPositionIdx.SpP(2,SpPIdx) - SpElemPositionIdx.SpS(2,SpSIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpP(3,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(3,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpP(3,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(3,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) =  1;
            end
        end
    end
    
end

for SpPIdx = ElemPer.YZFaceNum+ElemPer.ZXFaceNum+1:ElemPer.YZFaceNum+ElemPer.ZXFaceNum+ElemPer.XYFaceNum
    XIdx = ceil(SpElemPositionIdx.SpP(1,SpPIdx)+EPSILON);
    YIdx = ceil(SpElemPositionIdx.SpP(2,SpPIdx)+EPSILON);
    ZIdx = ceil(SpElemPositionIdx.SpP(3,SpPIdx)+EPSILON);

    IncXSpSIdxs_Guess = [...
        ((1:ElemPer.XEdgePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.XEdgePerCoarseGrid(XIdx,1:YIdx-1,ZIdx))...
        + sum(ElemPer.XEdgePerYRow(1:ZIdx-1,XIdx))...
        + sum(ElemPer.XEdgePerYZPlane(1:XIdx-1))...
        )...
        ((1:ElemPer.XEdgePerCoarseGrid(XIdx,YIdx+1,ZIdx))...
        + sum(ElemPer.XEdgePerCoarseGrid(XIdx,1:YIdx  ,ZIdx))...
        + sum(ElemPer.XEdgePerYRow(1:ZIdx-1,XIdx))...
        + sum(ElemPer.XEdgePerYZPlane(1:XIdx-1))...
        )...
        ];
    for SpSIdx = IncXSpSIdxs_Guess
        if ( SpElemPositionIdx.SpP(1,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - EPSILON < SpElemPositionIdx.SpS(1,SpSIdx)...
                && SpElemPositionIdx.SpP(1,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) + EPSILON > SpElemPositionIdx.SpS(1,SpSIdx))...
                && abs(SpElemPositionIdx.SpP(3,SpPIdx) - SpElemPositionIdx.SpS(3,SpSIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpP(2,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(2,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) =  1;
            elseif abs(SpElemPositionIdx.SpP(2,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(2,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) = -1;
            end
        end
    end

    IncYSpSIdxs_Guess = [((1:ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,1:ZIdx-1))...
        + sum(ElemPer.YEdgePerZRow(1:XIdx-1,YIdx))...
        + sum(ElemPer.YEdgePerZXPlane(1:YIdx-1))...
        + ElemPer.XEdgeNum)...
        ((1:ElemPer.YEdgePerCoarseGrid(XIdx+1,YIdx,ZIdx))...
        + sum(ElemPer.YEdgePerCoarseGrid(XIdx+1,YIdx,1:ZIdx-1))...
        + sum(ElemPer.YEdgePerZRow(1:XIdx  ,YIdx))...
        + sum(ElemPer.YEdgePerZXPlane(1:YIdx-1))...
        + ElemPer.XEdgeNum)...
        ];
    for SpSIdx = IncYSpSIdxs_Guess
        if ( SpElemPositionIdx.SpP(2,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - EPSILON < SpElemPositionIdx.SpS(2,SpSIdx)...
                && SpElemPositionIdx.SpP(2,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) + EPSILON > SpElemPositionIdx.SpS(2,SpSIdx))...
                && abs(SpElemPositionIdx.SpP(3,SpPIdx) - SpElemPositionIdx.SpS(3,SpSIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpP(1,SpPIdx) - 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(1,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpP(1,SpPIdx) + 0.5/ElemFineness.SpP(SpPIdx) - SpElemPositionIdx.SpS(1,SpSIdx) ) < EPSILON
                sC_Num_Nonzeros = sC_Num_Nonzeros +1;
                sC_Nonzeros_Row(sC_Num_Nonzeros) = SpPIdx;
                sC_Nonzeros_Col(sC_Num_Nonzeros) = SpSIdx;
                sC_Nonzeros_Val(sC_Num_Nonzeros) =  1;
            end
        end
    end
    
end

sC = sparse(sC_Nonzeros_Row,sC_Nonzeros_Col,sC_Nonzeros_Val,size(ElemFineness.SpP,2),size(ElemFineness.SpS,2));

end
%%
function sG = Subgrid_sG(SpElemPositionIdx,ElemFineness,ElemPer)
global EPSILON
sG_Num_Nonzeros = 0;
sG_Nonzeros_Row = 0;
sG_Nonzeros_Col = 0;
sG_Nonzeros_Val = 0;

for SpSIdx = 1:ElemPer.XEdgeNum
    XIdx = ceil(SpElemPositionIdx.SpS(1,SpSIdx)+EPSILON);
    YIdx = ceil(SpElemPositionIdx.SpS(2,SpSIdx)+EPSILON);
    ZIdx = ceil(SpElemPositionIdx.SpS(3,SpSIdx)+EPSILON);
    
    IncSpNIdxs_Guess = [...
        ((1:ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.NodePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
        + sum(ElemPer.NodePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.NodePerXYPlane(1:ZIdx-1))...
        )...
        ((1:ElemPer.NodePerCoarseGrid(XIdx+1,YIdx,ZIdx))...
        + sum(ElemPer.NodePerCoarseGrid(1:XIdx  ,YIdx,ZIdx))...
        + sum(ElemPer.NodePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.NodePerXYPlane(1:ZIdx-1))...
        )...
        ];
    for SpNIdx = IncSpNIdxs_Guess
        if  abs(SpElemPositionIdx.SpS(2,SpSIdx) - SpElemPositionIdx.SpN(2,SpNIdx) ) < EPSILON ...
                && abs(SpElemPositionIdx.SpS(3,SpSIdx) - SpElemPositionIdx.SpN(3,SpNIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpS(1,SpSIdx)-0.5/ElemFineness.SpS(SpSIdx) - SpElemPositionIdx.SpN(1,SpNIdx) ) < EPSILON
                sG_Num_Nonzeros = sG_Num_Nonzeros +1;
                sG_Nonzeros_Row(sG_Num_Nonzeros) = SpSIdx;
                sG_Nonzeros_Col(sG_Num_Nonzeros) = SpNIdx;
                sG_Nonzeros_Val(sG_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpS(1,SpSIdx)+0.5/ElemFineness.SpS(SpSIdx) - SpElemPositionIdx.SpN(1,SpNIdx) ) < EPSILON
                sG_Num_Nonzeros = sG_Num_Nonzeros +1;
                sG_Nonzeros_Row(sG_Num_Nonzeros) = SpSIdx;
                sG_Nonzeros_Col(sG_Num_Nonzeros) = SpNIdx;
                sG_Nonzeros_Val(sG_Num_Nonzeros) =  1;
            end
        end
    end
end

for SpSIdx = ElemPer.XEdgeNum+1:ElemPer.XEdgeNum+ElemPer.YEdgeNum
    XIdx = ceil(SpElemPositionIdx.SpS(1,SpSIdx)+EPSILON);
    YIdx = ceil(SpElemPositionIdx.SpS(2,SpSIdx)+EPSILON);
    ZIdx = ceil(SpElemPositionIdx.SpS(3,SpSIdx)+EPSILON);
    
    IncSpNIdxs_Guess = [...
        ((1:ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.NodePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
        + sum(ElemPer.NodePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.NodePerXYPlane(1:ZIdx-1))...
        )...
        ((1:ElemPer.NodePerCoarseGrid(XIdx,YIdx+1,ZIdx))...
        + sum(ElemPer.NodePerCoarseGrid(1:XIdx-1,YIdx+1,ZIdx))...
        + sum(ElemPer.NodePerXRow(1:YIdx  ,ZIdx))...
        + sum(ElemPer.NodePerXYPlane(1:ZIdx-1))...
        )...
        ];
    for SpNIdx = IncSpNIdxs_Guess
        if  abs(SpElemPositionIdx.SpS(3,SpSIdx) - SpElemPositionIdx.SpN(3,SpNIdx) ) < EPSILON ...
                && abs(SpElemPositionIdx.SpS(1,SpSIdx) - SpElemPositionIdx.SpN(1,SpNIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpS(2,SpSIdx)-0.5/ElemFineness.SpS(SpSIdx) - SpElemPositionIdx.SpN(2,SpNIdx) ) < EPSILON
                sG_Num_Nonzeros = sG_Num_Nonzeros +1;
                sG_Nonzeros_Row(sG_Num_Nonzeros) = SpSIdx;
                sG_Nonzeros_Col(sG_Num_Nonzeros) = SpNIdx;
                sG_Nonzeros_Val(sG_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpS(2,SpSIdx)+0.5/ElemFineness.SpS(SpSIdx) - SpElemPositionIdx.SpN(2,SpNIdx) ) < EPSILON
                sG_Num_Nonzeros = sG_Num_Nonzeros +1;
                sG_Nonzeros_Row(sG_Num_Nonzeros) = SpSIdx;
                sG_Nonzeros_Col(sG_Num_Nonzeros) = SpNIdx;
                sG_Nonzeros_Val(sG_Num_Nonzeros) =  1;
            end
        end
    end
end

for SpSIdx = ElemPer.XEdgeNum+ElemPer.YEdgeNum+1:ElemPer.XEdgeNum+ElemPer.YEdgeNum+ElemPer.ZEdgeNum
    XIdx = ceil(SpElemPositionIdx.SpS(1,SpSIdx)+EPSILON);
    YIdx = ceil(SpElemPositionIdx.SpS(2,SpSIdx)+EPSILON);
    ZIdx = ceil(SpElemPositionIdx.SpS(3,SpSIdx)+EPSILON);
    
    IncSpNIdxs_Guess = [...
        ((1:ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx))...
        + sum(ElemPer.NodePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
        + sum(ElemPer.NodePerXRow(1:YIdx-1,ZIdx))...
        + sum(ElemPer.NodePerXYPlane(1:ZIdx-1))...
        )...
        ((1:ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx+1))...
        + sum(ElemPer.NodePerCoarseGrid(1:XIdx-1,YIdx,ZIdx+1))...
        + sum(ElemPer.NodePerXRow(1:YIdx-1,ZIdx+1))...
        + sum(ElemPer.NodePerXYPlane(1:ZIdx  ))...
        )...
        ];
    for SpNIdx = IncSpNIdxs_Guess
        if  abs(SpElemPositionIdx.SpS(1,SpSIdx) - SpElemPositionIdx.SpN(1,SpNIdx) ) < EPSILON ...
                && abs(SpElemPositionIdx.SpS(2,SpSIdx) - SpElemPositionIdx.SpN(2,SpNIdx) ) < EPSILON
            if abs(SpElemPositionIdx.SpS(3,SpSIdx)-0.5/ElemFineness.SpS(SpSIdx) - SpElemPositionIdx.SpN(3,SpNIdx) ) < EPSILON
                sG_Num_Nonzeros = sG_Num_Nonzeros +1;
                sG_Nonzeros_Row(sG_Num_Nonzeros) = SpSIdx;
                sG_Nonzeros_Col(sG_Num_Nonzeros) = SpNIdx;
                sG_Nonzeros_Val(sG_Num_Nonzeros) = -1;
            elseif abs(SpElemPositionIdx.SpS(3,SpSIdx)+0.5/ElemFineness.SpS(SpSIdx) - SpElemPositionIdx.SpN(3,SpNIdx) ) < EPSILON
                sG_Num_Nonzeros = sG_Num_Nonzeros +1;
                sG_Nonzeros_Row(sG_Num_Nonzeros) = SpSIdx;
                sG_Nonzeros_Col(sG_Num_Nonzeros) = SpNIdx;
                sG_Nonzeros_Val(sG_Num_Nonzeros) =  1;
            end
        end
    end
end

sG = sparse(sG_Nonzeros_Row,sG_Nonzeros_Col,sG_Nonzeros_Val,size(ElemFineness.SpS,2),size(SpElemPositionIdx.SpN,2));

end

function NodePos_Prim = Subgrid_NodePos_Prim(SpElemPositionIdx,ElemPer,ElemFineness,isOnFineGridCorner,sG,sC,sD,MeshMeasurements)
clear AdjustPrimNodePos

NodePosPrimM = zeros(3,size(sG,2));

for SpNIdx = 1:size(sG,2)
    NodePosPrimM(:,SpNIdx) = SpElemPositionIdx.SpN(:,SpNIdx);
end
ProgressPercentage = -25;
for SpVIdx = 1:size(sD,1)
    if mod(SpVIdx,round(0.25*size(sD,1))) == 1
        ProgressPercentage = ProgressPercentage + 25;
        disp(['Subgrid_NodePos_Prim : Progress - ',num2str(ProgressPercentage),' %'])
    end
    YZSpPs = 1:ElemPer.YZFaceNum;
    IncYZSpPIdx = find(sD(SpVIdx,YZSpPs))+0;
    BoundaryYZSpPIdx = IncYZSpPIdx(logical(ElemFineness.SpP(IncYZSpPIdx)/ElemFineness.SpV(SpVIdx) >1));
    if size(BoundaryYZSpPIdx,2)==0
        continue;
    end
    NodePosPrimM = NodePosPrimM+AdjustPrimNodePos(BoundaryYZSpPIdx,SpVIdx,SpElemPositionIdx,ElemFineness,isOnFineGridCorner,sG,sC,sD);
    
    ZXSpPs = ElemPer.YZFaceNum+1:ElemPer.YZFaceNum+ElemPer.ZXFaceNum;
    IncZXSpPIdx = find(sD(SpVIdx,ZXSpPs))+ElemPer.YZFaceNum;
    BoundaryZXSpPIdx = IncZXSpPIdx(logical(ElemFineness.SpP(IncZXSpPIdx)/ElemFineness.SpV(SpVIdx) >1));
    if size(BoundaryZXSpPIdx,2)==0
        continue;
    end
    NodePosPrimM = NodePosPrimM+AdjustPrimNodePos(BoundaryZXSpPIdx,SpVIdx,SpElemPositionIdx,ElemFineness,isOnFineGridCorner,sG,sC,sD);
    
    XYSpPs = ElemPer.YZFaceNum+ElemPer.ZXFaceNum+1:ElemPer.YZFaceNum+ElemPer.ZXFaceNum+ElemPer.XYFaceNum; 
    IncXYSpPIdx = find(sD(SpVIdx,XYSpPs))+ElemPer.YZFaceNum+ElemPer.ZXFaceNum;
    BoundaryXYSpPIdx = IncXYSpPIdx(logical(ElemFineness.SpP(IncXYSpPIdx)/ElemFineness.SpV(SpVIdx) >1));
    if size(BoundaryXYSpPIdx,2)==0
        continue;
    end
    NodePosPrimM = NodePosPrimM+AdjustPrimNodePos(BoundaryXYSpPIdx,SpVIdx,SpElemPositionIdx,ElemFineness,isOnFineGridCorner,sG,sC,sD);
end

for SpNIdx = 1:size(sG,2)
    NodePos_Prim(SpNIdx).Vec = [MeshMeasurements.dxCoarse;MeshMeasurements.dyCoarse;MeshMeasurements.dzCoarse].*NodePosPrimM(:,SpNIdx);
end

end
%%
function NodePosPrimDeltaM = AdjustPrimNodePos(BoundarySpPIdx,SpVIdx,SpElemPositionIdx,ElemFineness,isOnFineGridCorner,sG,sC,sD)
global EPSILON
persistent isAlreadyAdjusted 

if isempty(isAlreadyAdjusted)
    isAlreadyAdjusted = false(1,size(sG,2));
end

NodePosPrimDeltaM = sparse(3,size(sG,2));

% find incident faces with greater fineness. 
% Here we assume that all of these faces are quadrangles.
% if FRatio is even
% find the central node
% for each edges inc to the central node, compute the shift amount of the other endpoint.
% to compute the shift amount, first project the dual edge to a plane which makes
% a right angle with the prototype plane 
% (the plane containing the inc-faces when node-shifts are 0)
% and intersects the protoplane at the prototype edge.
% The shift-amount of the shifted node is determined so that the projected
% dual edge and the edge makes a right angle.
% do until all the nodes on the boundary of the faces are done
% for all edges inc to these endpoints,  compute the shift amount of the other endpoint.
% elseif FRatio is odd
% find the central face
% define all the nodes on the boundary of this face as central nodes.
% for each edges inc to the central node, compute the shift amount of the other endpoint.
% do until all the nodes on the boundary of the faces are done
% for all edges inc to these endpoints,  compute the shift amount of the other endpoint.
% end if

switch mod(ElemFineness.SpP(BoundarySpPIdx(1))/ElemFineness.SpV(SpVIdx),2)
    case 0
        BoundarySpSLogIdx = logical(sum(logical(sC(BoundarySpPIdx,:))));
        BoundarySpSIdx    = find(sum(logical(sC(BoundarySpPIdx,:))));
        BoundarySpNLogIdx = logical(sum(logical(sG(BoundarySpSLogIdx,:))));
        BoundarySpNIdx    = find(sum(logical(sG(BoundarySpSLogIdx,:))));
        %% x方向に制限されている
        CenterSpNLocalLogIdx = ...
            logical(sum((SpElemPositionIdx.SpN([2;3],BoundarySpNLogIdx)-(SpElemPositionIdx.SpV([2;3],SpVIdx))*ones(1,size(BoundarySpNIdx,2))).^2)<EPSILON);
        CenterSpNLocalIdx = ...
            find(sum((SpElemPositionIdx.SpN([2;3],BoundarySpNLogIdx)-(SpElemPositionIdx.SpV([2;3],SpVIdx))*ones(1,size(BoundarySpNIdx,2))).^2)<EPSILON);
        sC_LocalOnBoundary = sC(BoundarySpPIdx,BoundarySpSIdx);
        sG_LocalOnBoundary = sG(BoundarySpSIdx,BoundarySpNIdx);
        Next_JustDoneNodeSet_LocalIdx = CenterSpNLocalIdx;
        DoneFlag_Local = CenterSpNLocalLogIdx;
        while any(~DoneFlag_Local)
            JustDoneNodeSet_LocalIdx = Next_JustDoneNodeSet_LocalIdx;
            Next_JustDoneNodeSet_LocalIdx = [];
            for JustDoneNode_LocalIdx = JustDoneNodeSet_LocalIdx
                JustDoneNode_GlobalIdx = BoundarySpNIdx(JustDoneNode_LocalIdx);
                for IncSpS_LocalIdx = find(sG_LocalOnBoundary(:,JustDoneNode_LocalIdx).')
                    IncSpS_GlobalIdx = BoundarySpSIdx(IncSpS_LocalIdx);
                    IncSpP_LocalIdx  = find(sC_LocalOnBoundary(:,IncSpS_LocalIdx),1);
                    IncSpP_GlobalIdx = BoundarySpPIdx(IncSpP_LocalIdx);
                    PrimEdgeVec_Prototype  = SpElemPositionIdx.SpN* sG(IncSpS_GlobalIdx,:).';
                    DualEdgeVec = SpElemPositionIdx.SpV* sD(:,IncSpP_GlobalIdx);
                    SpanVec1 = norm(PrimEdgeVec_Prototype)^(-1)*PrimEdgeVec_Prototype;
                    SpanVec2 = [...
                        1*(abs(DualEdgeVec(1))>=max(abs(DualEdgeVec([2 3]))));...
                        1*(abs(DualEdgeVec(2))>=max(abs(DualEdgeVec([3 1]))));...
                        1*(abs(DualEdgeVec(3))>=max(abs(DualEdgeVec([1 2]))))...
                        ];
                    SpanVec2 = SpanVec2 - dot(SpanVec2,SpanVec1)*SpanVec1;
                    SpanVec2 = norm(SpanVec2)^(-1)*SpanVec2;
                    ProjectingPlaneVec = cross(SpanVec1,SpanVec2);
                    ProjectingPlaneVec = norm(ProjectingPlaneVec)^(-1)*ProjectingPlaneVec;
                    ProjectedDualEdge  = DualEdgeVec - dot(DualEdgeVec,ProjectingPlaneVec)*ProjectingPlaneVec;
                    PrimEdgeVec_Actual = ...
                        (norm(PrimEdgeVec_Prototype)/dot(ProjectedDualEdge,SpanVec2))...
                        *(dot(ProjectedDualEdge,SpanVec2)*SpanVec1 - dot(ProjectedDualEdge,SpanVec1)*SpanVec2);
                    PrimEdgeVec_Actual = - PrimEdgeVec_Actual*sG(IncSpS_GlobalIdx,JustDoneNode_GlobalIdx);
                    Endpoint2_LocalIdx = find(sG_LocalOnBoundary(IncSpS_LocalIdx,:));
                    Endpoint2_LocalIdx = Endpoint2_LocalIdx(logical(Endpoint2_LocalIdx~=JustDoneNode_LocalIdx));
                    Endpoint2_GlobalIdx = BoundarySpNIdx(Endpoint2_LocalIdx);
                    if DoneFlag_Local(Endpoint2_LocalIdx) == false
                        NodePosPrimDeltaM(:,Endpoint2_GlobalIdx) ...
                            = PrimEdgeVec_Actual + NodePosPrimDeltaM(:,JustDoneNode_GlobalIdx) ...
                            -(SpElemPositionIdx.SpN(:,Endpoint2_GlobalIdx)-SpElemPositionIdx.SpN(:,JustDoneNode_GlobalIdx));
                        Next_JustDoneNodeSet_LocalIdx(size(Next_JustDoneNodeSet_LocalIdx,2)+1) = Endpoint2_LocalIdx;
                        DoneFlag_Local(Endpoint2_LocalIdx) = true;
                    end
                end
            end
        end
    case 1
        
end
for SpNIdx = BoundarySpNIdx(logical(isAlreadyAdjusted(BoundarySpNLogIdx)))
    if isOnFineGridCorner(SpNIdx) == false
        NodePosPrimDeltaM(:,SpNIdx) = [0;0;0];
    end
end
isAlreadyAdjusted(BoundarySpNLogIdx) = true;

end
%%
% function sD = Subgrid_sD(MeshMeasurements,Num_of_Elem,ElemPer)
% XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
% YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
% ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
% 
% sD = sparse(Num_of_Elem.SpV,Num_of_Elem.SpP);
% 
% for ZIdx = 1:ZSize
%     for YIdx = 1:YSize
%         for XIdx = 1:XSize
%             sD = AddYZFace_XMinusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%             sD = AddYZFace_XPlusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%             sD = AddZXFace_YMinusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%             sD = AddZXFace_YPlusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%             sD = AddXYFace_ZMinusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%             sD = AddXYFace_ZPlusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%         end
%     end
% end
% end
% %%
% function sD = AddYZFace_XMinusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% if XIdx == 1
%     GridFineness_XMinus = GridFineness;
% else
%     GridFineness_XMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
% end
% FRatio_XMinus = GridFineness_XMinus/GridFineness;
% for LocalZIdx = 1:GridFineness
%     for LocalYIdx = 1:GridFineness
%         for LocalXIdx = 1
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_XMinus<=1
%                 YZSpPIdx = LocalYIdx + GridFineness*(LocalZIdx-1) ...
%                     +sum(ElemPer.YZFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx))+sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))+sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%             else
%                 YZSpPIdx = zeros(1,FRatio_XMinus^2);
%                 for LocalLocalZIdx = 1:FRatio_XMinus
%                     for LocalLocalYIdx = 1:FRatio_XMinus
%                         YZSpPIdx(LocalLocalYIdx+FRatio_XMinus*(LocalLocalZIdx-1))...
%                             = LocalLocalYIdx...
%                             +(FRatio_XMinus  )*(LocalYIdx-1)...
%                             +(FRatio_XMinus  )*GridFineness*(LocalLocalZIdx-1)...
%                             +(FRatio_XMinus^2)*GridFineness*(LocalZIdx-1) ...
%                             +sum(ElemPer.YZFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx))+sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))+sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%                     end
%                 end
%             end
%             sD(SpVIdx,YZSpPIdx) = -1;
%         end
%         for LocalXIdx = 2:GridFineness
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) ...
%                 +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_XMinus<=1
%                 YZSpPIdx = LocalYIdx + GridFineness*(LocalZIdx-1) ...
%                     +(GridFineness^2)*(LocalXIdx-1)...
%                     +sum(ElemPer.YZFacePerCoaeseGrid(Xidx,1:YIdx-1,ZIdx))+sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))+sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%             else
%                 YZSpPIdx ...
%                     = LocalYIdx + GridFineness*(LocalZIdx-1) ...
%                     +(GridFineness^2)*(LocalXIdx-1)-(GridFineness^2)+(FRatio_XMinus*GridFineness)^2 ...
%                     +sum(ElemPer.YZFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx))+sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))+sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%             end
%             sD(SpVIdx,YZSpPIdx) = -1;
%         end
%     end
% end
% end
% %%
% function sD = AddYZFace_XPlusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% XSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% 
% if XIdx == 1
%     GridFineness_XMinus = GridFineness;
% else
%     GridFineness_XMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
% end
% 
% if XIdx == XSize
%     GridFineness_XPlus = GridFineness;
% else
%     GridFineness_XPlus = MeshMeasurements.LocalGridFineness(XIdx+1,YIdx,ZIdx);
% end
% FRatio_XMinus = GridFineness_XMinus/GridFineness;
% FRatio_XPlus = GridFineness_XPlus/GridFineness;
% 
% for LocalZIdx = 1:GridFineness
%     for LocalYIdx = 1:GridFineness
%         for LocalXIdx = 1:GridFineness-1
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) + (GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             YZSpPIdx ...
%                 = LocalYIdx + GridFineness*(LocalZIdx-1) ...
%                 +(GridFineness^2)*(LocalXIdx  ) - GridFineness^2 + (FRatio_XMinus*GridFineness)^2 ...
%                 +sum(ElemPer.YZFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx))...
%                 +sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))...
%                 +sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%             sD(SpVIdx,YZSpPIdx) = 1;
%         end
%         for LocalXIdx = GridFineness
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_XPlus<=1
%                 YZSpPIdx = LocalYIdx + GridFineness*(LocalZIdx-1) ...
%                     +sum(ElemPer.YZFacePerCoaeseGrid(XIdx+1,1:YIdx-1,ZIdx))...
%                     +sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx+1))...
%                     +sum(ElemPer.YZFacePerYZPlane(1:XIdx  ));
%             else
%                 YZSpPIdx = zeros(1,FRatio_XPlus^2);
%                 for LocalLocalZIdx = 1:FRatio_XPlus
%                     for LocalLocalYIdx = 1:FRatio_XPlus
%                         YZSpPIdx(LocalLocalYIdx+FRatio_XPlus*(LocalLocalZIdx-1))...
%                             = LocalLocalYIdx ...
%                             +(FRatio_XPlus  )*(LocalYIdx-1)...
%                             +(FRatio_XPlus  )*GridFineness*(LocalLocalZIdx-1)...
%                             +(FRatio_XPlus^2)*GridFineness*(LocalZIdx-1) ...
%                             +sum(ElemPer.YZFacePerCoaeseGrid(XIdx+1,1:YIdx-1,ZIdx))...
%                             +sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx+1))...
%                             +sum(ElemPer.YZFacePerYZPlane(1:XIdx  ));
%                     end
%                 end    
%             end
%             sD(SpVIdx,YZSpPIdx) = 1;
%         end
%     end
% end
% end
% %%
% function sD = AddZXFace_YMinusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% if YIdx == 1
%     GridFineness_YMinus = GridFineness;
% else
%     GridFineness_YMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);
% end
% FRatio_YMinus = GridFineness_YMinus/GridFineness;
% for LocalXIdx = 1:GridFineness
%     for LocalZIdx = 1:GridFineness
%         for LocalYIdx = 1
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_YMinus<=1
%                 ZXSpPIdx = LocalZIdx + GridFineness*(LocalXIdx-1) ...
%                     +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,YIdx,1:ZIdx-1))+sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx))+sum(ElemPer.ZXFacePerZXPlane(1:YIdx-1));
%             else
%                 ZXSpPIdx = zeros(1,FRatio_YMinus^2);
%                 for LocalLocalXIdx = 1:FRatio_YMinus
%                     for LocalLocalZIdx = 1:FRatio_YMinus
%                         ZXSpPIdx(LocalLocalZIdx+FRatio_YMinus*(LocalLocalXIdx-1)) ...
%                             = LocalLocalZIdx...
%                             +(FRatio_YMinus  )*(LocalZIdx-1)...
%                             +(FRatio_YMinus  )*GridFineness*(LocalLocalXIdx-1)...
%                             +(FRatio_YMinus^2)*GridFineness*(LocalXIdx-1) ...
%                             +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,YIdx,1:ZIdx-1))+sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx))+sum(ElemPer.ZXFacePerZXPlane(1:YIdx-1));
%                     end
%                 end
%             end
%             sD(SpVIdx,ZXSpPIdx) = -1;
%         end
%         for LocalYIdx = 2:GridFineness
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_YMinus<=1
%                 ZXSpPIdx = LocalZIdx + GridFineness*(LocalXIdx-1) +(GridFineness^2)*(LocalYIdx-1)...
%                     +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx))+sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx))+sum(ElemPer.ZXFacePerZXPlane(1:YIdx-1));
%             else
%                 ZXSpPIdx ...
%                     = LocalZIdx + GridFineness*(LocalXIdx-1) ...
%                     +(GridFineness^2)*(LocalYIdx-1)-(GridFineness^2)+(FRatio_YMinus*GridFineness)^2 ...
%                     +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx))+sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx))+sum(ElemPer.ZXFacePerZXPlane(1:YISdx-1));
%             end
%             sD(SpVIdx,ZXSpPIdx) = -1;
%         end
%     end
% end
% end
% %%
% function sD = AddZXFace_YPlusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% YSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% 
% if YIdx == 1
%     GridFineness_YMinus = GridFineness;
% else
%     GridFineness_YMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);
% end
% 
% if YIdx == YSize
%     GridFineness_YPlus = GridFineness;
% else
%     GridFineness_YPlus = MeshMeasurements.LocalGridFineness(XIdx,YIdx+1,ZIdx);
% end
% FRatio_YMinus = GridFineness_YMinus/GridFineness;
% FRatio_YPlus = GridFineness_YPlus/GridFineness;
% 
% for LocalXIdx = 1:GridFineness
%     for LocalZIdx = 1:GridFineness
%         for LocalYIdx = 1:GridFineness-1
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) + (GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             ZXSpPIdx ...
%                 = LocalZIdx + GridFineness*(LocalXIdx-1) ...
%                 +(GridFineness^2)*(LocalYIdx  ) - GridFineness^2 + (FRatio_YMinus*GridFineness)^2 ...
%                 +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,YIdx,1:ZIdx-1))...
%                 +sum(ElemPer.ZXFacePerZRow(1:ZIdx-1,YIdx))...
%                 +sum(ElemPer.ZXFacePerZXPlane(1:YIdx-1));
%             sD(SpVIdx,ZXSpPIdx) = 1;
%         end
%         for LocalYIdx = GridFineness
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_YPlus<=1
%                 ZXSpPIdx = LocalZIdx + GridFineness*(LocalYIdx-1) ...
%                     +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,YIdx+1,1:ZIdx-1))...
%                     +sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx+1))...
%                     +sum(ElemPer.ZXFacePerZXPlane(1:YIdx  ));
%             else
%                 ZXSpPIdx = zeros(1,FRatio_YPlus^2);
%                 for LocalLocalXIdx = 1:FRatio_YPlus
%                     for LocalLocalZIdx = 1:FRatio_YPlus
%                         ZXSpPIdx(LocalLocalZIdx+FRatio_YPlus*(LocalLocalXIdx-1))...
%                             = LocalLocalZIdx ...
%                             +(FRatio_YPlus  )*(LocalZIdx-1)...
%                             +(FRatio_YPlus  )*GridFineness*(LocalLocalXIdx-1)...
%                             +(FRatio_YPlus^2)*GridFineness*(LocalXIdx-1) ...
%                             +sum(ElemPer.ZXFacePerCoaeseGrid(XIdx,YIdx+1,1:ZIdx-1))...
%                             +sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx+1))...
%                             +sum(ElemPer.ZXFacePerZXPlane(1:YIdx  ));
%                     end
%                 end    
%             end
%             sD(SpVIdx,ZXSpPIdx) = 1;
%         end
%     end
% end
% end
% 
% %%
% 
% function sD = AddXYFace_ZMinusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% if ZIdx == 1
%     GridFineness_ZMinus = GridFineness;
% else
%     GridFineness_ZMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
% end
% FRatio_ZMinus = GridFineness_ZMinus/GridFineness;
% for LocalYIdx = 1:GridFineness
%     for LocalXIdx = 1:GridFineness
%         for LocalZIdx = 1
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_ZMinus<=1
%                 XYSpPIdx = LocalXIdx + GridFineness*(LocalYIdx-1) ...
%                     +sum(ElemPer.XYFacePerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1));
%             else
%                 XYSpPIdx = zeros(1,FRatio_ZMinus^2);
%                 for LocalLocalYIdx = 1:FRatio_ZMinus
%                     for LocalLocalXIdx = 1:FRatio_ZMinus
%                         XYSpPIdx(LocalLocalXIdx+FRatio_ZMinus*(LocalLocalYIdx-1))...
%                             = LocalLocalXIdx...
%                             +(FRatio_ZMinus  )*(LocalXIdx-1)...
%                             +(FRatio_ZMinus  )*GridFineness*(LocalLocalYIdx-1)...
%                             +(FRatio_ZMinus^2)*GridFineness*(LocalYIdx-1) ...
%                             +sum(ElemPer.XYFacePerCoaeseGrid(XIdx-1,YIdx,ZIdx))+sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))+sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%                     end
%                 end
%             end
%             sD(SpVIdx,XYSpPIdx) = -1;
%         end
%         for LocalZIdx = 2:GridFineness
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) ...
%                 +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_ZMinus<=1
%                 XYSpPIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                     +sum(ElemPer.XYFacePerCoaeseGrid(1:Xidx-1,YIdx,ZIdx))+sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1));
%             else
%                 XYSpPIdx ...
%                     = LocalXIdx + GridFineness*(LocalYIdx-1) ...
%                     +(GridFineness^2)*(LocalYIdx-1)-(GridFineness^2)+(FRatio_ZMinus*GridFineness)^2 ...
%                     +sum(ElemPer.XYFacePerCoaeseGrid(1:Xidx-1,YIdx,ZIdx))+sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1));
%             end
%             sD(SpVIdx,XYSpPIdx) = -1;
%         end
%     end
% end
% end
% %%
% function sD = AddXYFace_ZPlusSide(sD,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% 
% if ZIdx == 1
%     GridFineness_ZMinus = GridFineness;
% else
%     GridFineness_ZMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
% end
% 
% if ZIdx == ZSize
%     GridFineness_ZPlus = GridFineness;
% else
%     GridFineness_ZPlus = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx+1);
% end
% FRatio_ZMinus = GridFineness_ZMinus/GridFineness;
% FRatio_ZPlus = GridFineness_ZPlus/GridFineness;
% 
% for LocalYIdx = 1:GridFineness
%     for LocalXIdx = 1:GridFineness
%         for LocalZIdx = 1:GridFineness-1
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) + (GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             XYSpPIdx ...
%                 = LocalXIdx + GridFineness*(LocalYIdx-1) ...
%                 +(GridFineness^2)*(LocalZIdx  ) - GridFineness^2 + (FRatio_ZMinus*GridFineness)^2 ...
%                 +sum(ElemPer.XYFacePerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))...
%                 +sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx  ))...
%                 +sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1));
%             sD(SpVIdx,XYSpPIdx) = 1;
%         end
%         for LocalZIdx = GridFineness
%             SpVIdx = LocalXIdx + GridFineness*(LocalYIdx-1) +(GridFineness^2)*(LocalZIdx-1)...
%                 +sum(ElemPer.VolPerCoaeseGrid(1:XIdx-1,YIdx,ZIdx))+sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))+sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
%             if FRatio_ZPlus<=1
%                 XYSpPIdx = LocalXIdx + GridFineness*(LocalYIdx-1) ...
%                     +sum(ElemPer.XYFacePerCoaeseGrid(1:XIdx-1,YIdx,ZIdx+1))...
%                     +sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx+1))...
%                     +sum(ElemPer.XYFacePerXYPlane(1:ZIdx  ));
%             else
%                 XYSpPIdx = zeros(1,FRatio_ZPlus^2);
%                 for LocalLocalYIdx = 1:FRatio_ZPlus
%                     for LocalLocalXIdx = 1:FRatio_ZPlus
%                         XYSpPIdx(LocalLocalXIdx+FRatio_ZPlus*(LocalLocalYIdx-1))...
%                             = LocalLocalXIdx ...
%                             +(FRatio_ZPlus  )*(LocalXIdx-1)...
%                             +(FRatio_ZPlus  )*GridFineness*(LocalLocalYIdx-1)...
%                             +(FRatio_ZPlus^2)*GridFineness*(LocalYIdx-1) ...
%                             +sum(ElemPer.XYFacePerCoaeseGrid(1:XIdx-1,YIdx,ZIdx+1))...
%                             +sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx+1))...
%                             +sum(ElemPer.XYFacePerXYPlane(1:ZIdx  ));
%                     end
%                 end    
%             end
%             sD(SpVIdx,XYSpPIdx) = 1;
%         end
%     end
% end
% 
% end
% 
% %%
% function sC = Subgrid_sC(MeshMeasurements,Num_of_Elem,ElemPer)
% 
% XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
% YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
% ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
% 
% sC = sparse(Num_of_Elem.SpP,Num_of_Elem.SpS);     
% 
% for ZIdx = 1:ZSize
%     for YIdx = 1:YSize
%         for XIdx = 1:XSize
%             sC = AddYEdgeIncToYZFace_YMinusSide(sC,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
%         end
%         for XIdx = XSize+1
%         end
%     end
% end
% end
% %%
% function sC = AddYEdgeIncToYZFace_YMinusSide(sC,XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
% 
% GridFineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% 
% if XIdx == 1
%     GridFineness_XMinus = GridFineness;
% else
%     GridFineness_XMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
% end
% if YIdx == 1
%     GridFineness_YMinus = GridFineness;
% else
%     GridFineness_YMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
% end
% if ZIdx == 1
%     GridFineness_ZMinus = GridFineness;
% else
%     GridFineness_ZMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
% end
% FRatio_XMinus = GridFineness_XMinus/GridFineness;
% FRatio_YMinus = GridFineness_YMinus/GridFineness;
% FRatio_ZMinus = GridFineness_ZMinus/GridFineness;
% if ZIdx == 1 || XIdx == 1
%     GridFineness_ZXMinus = GridFineness;
% else
%     GridFineness_ZXMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1);
% end
% FRatio_ZXMinus = GridFineness_ZXMinus/GridFineness;
% for LocalZIdx = 1:GridFineness
%     for LocalYIdx = 1:GridFineness
%         for LocalXIdx = 1
%             if FRatio_XMinus<=1
%                 YZSpPIdx = ...
%                     LocalYIdx + GridFineness*(LocalZIdx-1) + (GridFineness^2)*(LocalXIdx-1) ...
%                     +sum(ElemPer.YZFacePerCoaeseGrid(XIdx,1:YIdx-1,ZIdx)) ...
%                     +sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))...
%                     +sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
%             else
%                 for LocalLocalZIdx = 1:FRatio_XMinus
%                     for LocalLocalYIdx = 1:FRatio_XMinus
%                     end
%                 end
%             end
%         end
%         for LocalXIdx = 2:GridFineness           
%         end
%     end
% end