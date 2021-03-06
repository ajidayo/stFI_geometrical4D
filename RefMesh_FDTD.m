function [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties,ElemPer] ...
    = RefMesh_FDTD(MeshMeasurements,LocalUpdateNum)
% global EPSILON
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

Num_of_Elem.SpV = XSize*YSize*ZSize;
Num_of_Elem.SpP = ...
    (XSize+1)*YSize*ZSize ...
    + XSize*(YSize+1)*ZSize ...
    + XSize*YSize*(ZSize+1);
Num_of_Elem.SpS ...
    =(XSize  )*(YSize+1)*(ZSize+1) ...
    +(XSize+1)*(YSize  )*(ZSize+1) ...
    +(XSize+1)*(YSize+1)*(ZSize  );
Num_of_Elem.SpN = (XSize+1)*(YSize+1)*(ZSize+1);

sG = FDTD_sG(MeshMeasurements,Num_of_Elem);
sC = FDTD_sC(MeshMeasurements,Num_of_Elem);
sD = FDTD_sD(MeshMeasurements,Num_of_Elem);
NodePos = FDTD_NodePos(MeshMeasurements);
%%
SpElemProperties.SpV.UpdNum     = ones(Num_of_Elem.SpV,1);
VolPerXRow          = XSize;
VolPerXYPlane       = XSize*YSize;
%VolNum              = XSize*YSize*ZSize;
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = round(2*XSize/3):XSize
              SpVIdx = XIdx + (YIdx-1)*VolPerXRow + (ZIdx-1)*VolPerXYPlane;
              SpElemProperties.SpV.UpdNum(SpVIdx) = LocalUpdateNum;
        end
    end
end

%%
SpElemProperties.SpS.PEC        = FDTD_SpS_PEC(MeshMeasurements,Num_of_Elem);
SpElemProperties.SpP.ElecWall   = FDTD_SpP_ElecWall(MeshMeasurements,Num_of_Elem);

SpElemProperties.SpP.PrimAreaIsGiven = true(Num_of_Elem.SpP,1);
SpElemProperties.SpP.DualLengIsGiven = true(Num_of_Elem.SpP,1);
SpElemProperties.SpS.PrimLengIsGiven = true(Num_of_Elem.SpS,1);
SpElemProperties.SpS.DualAreaIsGiven = true(Num_of_Elem.SpS,1);
SpElemProperties.SpP.PrimArea = ones(Num_of_Elem.SpP,1);
SpElemProperties.SpP.DualLeng = ones(Num_of_Elem.SpP,1);
SpElemProperties.SpS.PrimLeng = ones(Num_of_Elem.SpS,1);
SpElemProperties.SpS.DualArea = ones(Num_of_Elem.SpS,1);

ElemPer.XEdgePerCoarseGrid  = 1*ones(XSize,YSize+1,ZSize+1);
ElemPer.XEdgePerYRow        = (YSize+1)*ones(ZSize+1,XSize);
ElemPer.XEdgePerYZPlane     = (YSize+1)*(ZSize+1)*ones(XSize);
ElemPer.XEdgeNum            = (YSize+1)*(ZSize+1)*XSize;
ElemPer.YEdgePerCoarseGrid  = 1*ones(XSize+1,YSize,ZSize+1);
ElemPer.YEdgePerZRow        = (ZSize+1)*ones(XSize+1,YSize);
ElemPer.YEdgePerZXPlane     = (ZSize+1)*(XSize+1)*ones(YSize);
ElemPer.YEdgeNum            = (ZSize+1)*(XSize+1)*YSize;
ElemPer.XEdgePerCoarseGrid  = 1*ones(XSize+1,YSize+1,ZSize);
ElemPer.ZEdgePerXRow        = (ZSize+1)*ones(XSize+1,ZSize);
ElemPer.ZEdgePerXYPlane     = (ZSize+1)*(XSize+1)*ones(ZSize);
ElemPer.ZEdgeNum            = (ZSize+1)*(XSize+1)*ZSize;

ElemPer.YZFacePerCoarseGrid = 1*ones(XSize+1,YSize,ZSize);
ElemPer.YZFacePerYRow       = (YSize)*ones(ZSize,XSize+1);
ElemPer.YZFacePerYZPlane    = (YSize)*(ZSize)*ones(XSize+1);
ElemPer.YZFaceNum           = YSize*ZSize*(XSize+1);
ElemPer.ZXFacePerCoarseGrid = 1*ones(XSize,YSize+1,ZSize);
ElemPer.ZXFacePerZRow       = (ZSize)*ones(XSize,YSize+1);
ElemPer.ZXFacePerZXPlane    = (ZSize)*(XSize)*ones(YSize+1);
ElemPer.ZXFaceNum           = ZSize*XSize*(YSize+1);
ElemPer.XYFacePerCoarseGrid = 1*ones(XSize,YSize,ZSize+1);
ElemPer.XYFacePerXRow       = (XSize)*ones(YSize,ZSize+1);
ElemPer.XYFacePerXYPlane    = (XSize)*(YSize)*ones(ZSize+1);
ElemPer.XYFaceNum           = XSize*YSize*(ZSize+1);

SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion = sparse(false(Num_of_Elem.SpN,1));
end

function sG = FDTD_sG(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = XSize+1;
YEdgePerXYPlane     = (XSize+1)*YSize;
YEdgeNum            = (XSize+1)*YSize*(ZSize+1);
ZEdgePerXRow        = XSize+1;
ZEdgePerXYPlane     = (XSize+1)*(YSize+1);
%ZEdgeNum            = (XSize+1)*(YSize+1)*ZSize;

NodePerXRow        = XSize+1;
NodePerXYPlane     = (XSize+1)*(YSize+1);
%NodeNum            = (XSize+1)*(YSize+1)*(ZSize+1);

sG = spalloc(Num_of_Elem.SpS,Num_of_Elem.SpN,2*Num_of_Elem.SpS);
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize
            SIdx    = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow  + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) = -1;
            IncNIdx = XIdx+1 + (YIdx-1)*NodePerXRow  + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize
        for XIdx = 1:XSize+1
            SIdx    = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow  + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) = -1;
            IncNIdx = XIdx   + (YIdx  )*NodePerXRow  + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize+1
            SIdx    = XIdx   + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow  + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) = -1;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow  + (ZIdx  )*NodePerXYPlane;
            sG(SIdx,IncNIdx) =  1;
        end
    end
end

end

function sC = FDTD_sC(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

YZFacePerXRow       = (XSize+1);
YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum           = (XSize+1)*YSize*ZSize;
ZXFacePerXRow       = XSize;
ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum           = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum          = XSize*YSize*(ZSize+1);

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = (XSize+1);
YEdgePerXYPlane     = (XSize+1)*YSize;
YEdgeNum            = (XSize+1)*YSize*(ZSize+1);
ZEdgePerXRow        = (XSize+1);
ZEdgePerXYPlane     = (XSize+1)*(YSize+1);
%ZEdgeNum           = (XSize+1)*(YSize+1)*ZSize;

sC = spalloc(Num_of_Elem.SpP,Num_of_Elem.SpS,4*Num_of_Elem.SpP);

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize+1
            PIdx    = XIdx   + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            IncSIdx = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane...
                +XEdgeNum;
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx  )*YEdgePerXYPlane...
                +XEdgeNum;  
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx-1)*ZEdgePerXRow  + (ZIdx-1)*ZEdgePerXYPlane...
                +XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx  )*ZEdgePerXRow  + (ZIdx-1)*ZEdgePerXYPlane...
                +XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize
            PIdx    = XIdx   + (YIdx-1)*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane ...
                +YZFaceNum;
            IncSIdx = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx  )*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx   + (YIdx-1)*ZEdgePerXRow  + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx+1 + (YIdx-1)*ZEdgePerXRow  + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) = -1;
        end
    end
end
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            PIdx    = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx-1)*XYFacePerXYPlane ...
                +YZFaceNum+ZXFaceNum;
            IncSIdx = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;  
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;
            sC(PIdx,IncSIdx) =  1;
        end
    end
end
end
function sD = FDTD_sD(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;
VolPerXRow          = XSize;
VolPerXYPlane       = XSize*YSize;
%VolNum              = XSize*YSize*ZSize;

YZFacePerXRow       = XSize+1;
YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum        = (XSize+1)*YSize*ZSize;
ZXFacePerXRow       = XSize;
ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum        = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum        = XSize*YSize*(ZSize+1);

sD = spalloc(Num_of_Elem.SpV,Num_of_Elem.SpP,6*Num_of_Elem.SpV);
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            VIdx = XIdx + (YIdx-1)*VolPerXRow + (ZIdx-1)*VolPerXYPlane;
            IncPIdx = XIdx   + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            sD(VIdx,IncPIdx) = -1;
            IncPIdx = XIdx+1 + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            sD(VIdx,IncPIdx) =  1;
            IncPIdx = XIdx   + (YIdx-1)*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane ...
                +YZFaceNum;
            sD(VIdx,IncPIdx) = -1;
            IncPIdx = XIdx   + (YIdx  )*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane ...
                +YZFaceNum;            
            sD(VIdx,IncPIdx) =  1;
            IncPIdx = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx-1)*XYFacePerXYPlane ...
                +YZFaceNum+ZXFaceNum;
            sD(VIdx,IncPIdx) = -1;
            IncPIdx = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx  )*XYFacePerXYPlane ...
                +YZFaceNum+ZXFaceNum;
            sD(VIdx,IncPIdx) =  1;
        end
    end
end
end

function NodePos = FDTD_NodePos(MeshMeasurements)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

NodePerXRow        = XSize+1;
NodePerXYPlane     = (XSize+1)*(YSize+1);
%NodeNum            = (XSize+1)*(YSize+1)*(ZSize+1);
VolPerXRow          = XSize;
VolPerXYPlane       = XSize*YSize;
%VolNum              = XSize*YSize*ZSize;



for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize+1
            SpNIdx = XIdx   + (YIdx-1)*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            XPos = MeshMeasurements.dxCoarse*(XIdx-1);
            YPos = MeshMeasurements.dyCoarse*(YIdx-1);
            ZPos = MeshMeasurements.dzCoarse*(ZIdx-1);
            NodePos.Prim(SpNIdx).Vec = [XPos;YPos;ZPos]; 
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpNIdx = XIdx   + (YIdx-1)*VolPerXRow + (ZIdx-1)*VolPerXYPlane;
            XPos = MeshMeasurements.dxCoarse*(XIdx-0.5);
            YPos = MeshMeasurements.dyCoarse*(YIdx-0.5);
            ZPos = MeshMeasurements.dzCoarse*(ZIdx-0.5);
            NodePos.Dual(SpNIdx).Vec = [XPos;YPos;ZPos];
        end
    end
end
end

function SpElemProperties_SpS_PEC = FDTD_SpS_PEC(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = XSize+1;
YEdgePerXYPlane     = (XSize+1)*YSize;
YEdgeNum            = (XSize+1)*YSize*(ZSize+1);
ZEdgePerXRow        = XSize+1;
ZEdgePerXYPlane     = (XSize+1)*(YSize+1);
%ZEdgeNum            = (XSize+1)*(YSize+1)*ZSize;

SpElemProperties_SpS_PEC = zeros(1,Num_of_Elem.SpS);
for ZIdx = 1:ZSize+1
    for YIdx = [1 YSize+1]
        for XIdx = 1:XSize
            SIdx = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            SpElemProperties_SpS_PEC(SIdx) = true;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = [1 YSize+1]
        for XIdx = 1:XSize+1
            SIdx = XIdx   + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane +XEdgeNum+YEdgeNum;
            SpElemProperties_SpS_PEC(SIdx) = true;
        end
    end
end
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize
        for XIdx = [1 XSize+1]
            SIdx = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
            SpElemProperties_SpS_PEC(SIdx) = true;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize+1
        for XIdx = [1 XSize+1]
            SIdx = XIdx   + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane +XEdgeNum+YEdgeNum;
            SpElemProperties_SpS_PEC(SIdx) = true;
        end
    end
end
for ZIdx = [1 ZSize+1]
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize
            SIdx = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            SpElemProperties_SpS_PEC(SIdx) = true;
        end
    end
end
for ZIdx = [1 ZSize+1]
    for YIdx = 1:YSize
        for XIdx = 1:XSize+1
            SIdx = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
            SpElemProperties_SpS_PEC(SIdx) = true;
        end
    end
end
end

function SpElemProperties_SpP_ElecWall = FDTD_SpP_ElecWall(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

YZFacePerXRow       = (XSize+1);
YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum           = (XSize+1)*YSize*ZSize;
ZXFacePerXRow       = XSize;
ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum           = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum          = XSize*YSize*(ZSize+1);

SpElemProperties_SpP_ElecWall = zeros(1,Num_of_Elem.SpP);
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = [1 XSize+1]
            PIdx = XIdx   + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            SpElemProperties_SpP_ElecWall(PIdx) = true;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = [1 YSize+1]
        for XIdx = 1:XSize
            PIdx = XIdx   + (YIdx-1)*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane+YZFaceNum;
            SpElemProperties_SpP_ElecWall(PIdx) = true;
        end
    end
end
for ZIdx = [1 ZSize+1]
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            PIdx = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx-1)*XYFacePerXYPlane+YZFaceNum+ZXFaceNum;
            SpElemProperties_SpP_ElecWall(PIdx) = true;
        end
    end
end

end