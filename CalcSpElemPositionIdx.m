function [SpElemPositionIdx,ElemFineness] = CalcSpElemPositionIdx(MeshMeasurements,ElemPer,Num_of_Elem)

[SpElemPositionIdx.SpV,ElemFineness.SpV] = CalcSpElemPositionIdx_SpV(MeshMeasurements,ElemPer,Num_of_Elem);
[SpElemPositionIdx.SpP,ElemFineness.SpP] = CalcSpElemPositionIdx_SpP(MeshMeasurements,ElemPer,Num_of_Elem);
[SpElemPositionIdx.SpS,ElemFineness.SpS] = CalcSpElemPositionIdx_SpS(MeshMeasurements,ElemPer,Num_of_Elem);
[SpElemPositionIdx.SpM,ElemFineness.SpN] = CalcSpElemPositionIdx_SpN(MeshMeasurements,ElemPer,Num_of_Elem);
end

%%

function [SpElemPositionIdx_SpV,ElemFineness_SpV] = CalcSpElemPositionIdx_SpV(MeshMeasurements,ElemPer,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

SpElemPositionIdx_SpV = zeros(3,Num_of_Elem.SpV);
ElemFineness_SpV = zeros(1,Num_of_Elem.SpV);

for ZIdx = 1:ZSize    
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            Offset = ...
                +sum(ElemPer.VolPerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
                +sum(ElemPer.VolPerXRow(1:YIdx-1,ZIdx))...
                +sum(ElemPer.VolPerXYPlane(1:ZIdx-1));
            SpVIdx = ...
                Offset+1:Offset+ElemPer.VolPerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpV(:,SpVIdx), ElemFineness_SpV(SpVIdx)] = LocallyCalcSpVPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end
end
%%
function [SpVPositionIdx_Local,SpVFineness_Local]=  LocallyCalcSpVPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)

SpVPositionIdx_Local = zeros(3,ElemPer.VolPerCoarseGrid(XIdx,YIdx,ZIdx));
SpVFineness_Local    = zeros(1,ElemPer.VolPerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
SpVFineness_Local(1:ElemPer.VolPerCoarseGrid(XIdx,YIdx,ZIdx)) = Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];
for LocalZIdx = 1:Fineness
    for LocalYIdx = 1:Fineness
        for LocalXIdx = 1:Fineness
            LocalSpVIdx = LocalXIdx + Fineness*(LocalYIdx-1) + Fineness^2*(LocalZIdx-1);
            SpVPositionIdx_Local(1:3,LocalSpVIdx) = PosIdxOffset ...
                +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-0.5;LocalZIdx-0.5];
            
        end
    end
end
end
%%
function [SpElemPositionIdx_SpP,ElemFineness_SpP] = CalcSpElemPositionIdx_SpP(MeshMeasurements,ElemPer,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

SpElemPositionIdx_SpP = zeros(3,Num_of_Elem.SpP);
ElemFineness_SpP      = zeros(1,Num_of_Elem.SpP);

for ZIdx = 1:ZSize    
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpPIdxOffset = ...
                +sum(ElemPer.YZFacePerCoarseGrid(XIdx,1:YIdx-1,ZIdx))...
                +sum(ElemPer.YZFacePerYRow(1:ZIdx-1,XIdx))...
                +sum(ElemPer.YZFacePerYZPlane(1:XIdx-1));
            YZSpPIdx = ...
                SpPIdxOffset+1:SpPIdxOffset+ElemPer.YZFacePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpP(:,YZSpPIdx),ElemFineness_SpP(YZSpPIdx)] ...
                = LocallyCalcYZSpPPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

for XIdx = 1:XSize    
    for ZIdx = 1:ZSize
        for YIdx = 1:YSize
            SpPIdxOffset = ...
                +sum(ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,1:ZIdx-1))...
                +sum(ElemPer.ZXFacePerZRow(1:XIdx-1,YIdx))...
                +sum(ElemPer.ZXFacePerZXPlane(1:YIdx-1))...
                +ElemPer.YZFaceNum;
            ZXSpPIdx = ...
                SpPIdxOffset+1:SpPIdxOffset+ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpP(:,ZXSpPIdx),ElemFineness_SpP(ZXSpPIdx)] ...
                = LocallyCalcZXSpPPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

for YIdx = 1:YSize    
    for XIdx = 1:XSize
        for ZIdx = 1:ZSize
            SpPIdxOffset = ...
                +sum(ElemPer.XYFacePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
                +sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx))...
                +sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1))...
                +ElemPer.YZFaceNum+ElemNum.ZXFaceNum;
            XYSpPIdx = ...
                SpPIdxOffset+1:SpPIdxOffset+ElemPer.XYFacePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpP(:,XYSpPIdx),ElemFineness_SpP(XYSpPIdx)] ...
                = LocallyCalcXYSpPPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

end

function [SpPPositionIdx_Local, SpPFineness_Local] = LocallyCalcYZSpPPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
SpPPositionIdx_Local = zeros(3,ElemPer.YZFacePerCoarseGrid(XIdx,YIdx,ZIdx));
SpPFineness_Local    = zeros(1,ElemPer.YZFacePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_XMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx);
FRatio_XMinus = Fineness_XMinus/Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];
for LocalZIdx = 1:Fineness
    for LocalYIdx = 1:Fineness
        for LocalXIdx = 1
            if FRatio_XMinus <= 1
                Local_YZSpPIdx = LocalYIdx + Fineness*(LocalZIdx-1);
                SpPPositionIdx_Local(1:3,Local_YZSpPIdx) = PosIdxOffset ...
                    +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-0.5];
                SpPFineness_Local(Local_YZSpPIdx) = Fineness;
            else
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalZIdx = 1:FRatio_XMinus
                    for LocalLocalYIdx = 1:FRatio_XMinus
                        Local_YZSpPIdx ...
                            =  LocalLocalYIdx ...
                            +  FRatio_XMinus*(LocalYIdx-1) ...
                            + (FRatio_XMinus  *Fineness)*(LocalLocalZIdx-1) ...
                            + (FRatio_XMinus^2*Fineness)*(LocalZIdx-1);
                        SpPPositionIdx_Local(1:3,Local_YZSpPIdx) = PosIdxOffset +LocalPosIdxOffset ...
                            +(FRatio_XMinus*Fineness)^(-1)*[0;LocalLocalYIdx-0.5;LocalLocalZIdx-0.5];
                        SpPFineness_Local(Local_YZSpPIdx) = FRatio_XMinus*Fineness;
                    end
                end
            end
        end
        for LocalXIdx = 2:Fineness
            Local_YZSpPIdx = LocalYIdx + Fineness*(LocalZIdx-1) + Fineness^2*(LocalXIdx-2)+(FRatio_XMinus*Fineness)^2;
            SpPPositionIdx_Local(1:3,Local_YZSpPIdx) = PosIdxOffset ...
                +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-0.5];
            SpPFineness_Local(Local_YZSpPIdx) = Fineness;
        end
    end
end
end

function [SpPPositionIdx_Local,SpPFineness_Local] = LocallyCalcZXSpPPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
SpPPositionIdx_Local = zeros(3,ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx));
SpPFineness_Local    = zeros(1,ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_YMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx);
FRatio_YMinus = Fineness_YMinus/Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];
for LocalXIdx = 1:Fineness
    for LocalZIdx = 1:Fineness
        for LocalYIdx = 1
            if FRatio_YMinus <= 1
                Local_ZXSpPIdx = LocalZIdx + Fineness*(LocalXIdx-1);
                SpPPositionIdx_Local(1:3,Local_ZXSpPIdx) = PosIdxOffset ...
                    +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-0.5];
                SpPFineness_Local(Local_ZXSpPIdx) = Fineness;
            else
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalXIdx = 1:FRatio_YMinus
                    for LocalLocalZIdx = 1:FRatio_YMinus
                        Local_ZXSpPIdx ...
                            =  LocalLocalZIdx ...
                            +  FRatio_YMinus*(LocalZIdx-1) ...
                            + (FRatio_YMinus  *Fineness)*(LocalLocalXIdx-1) ...
                            + (FRatio_YMinus^2*Fineness)*(LocalXIdx-1);
                        SpPPositionIdx_Local(1:3,Local_ZXSpPIdx) = PosIdxOffset + LocalPosIdxOffset ...
                            +(FRatio_YMinus*Fineness)^(-1)*[LocalLocalXIdx-0.5;0;LocalLocalZIdx-0.5];
                        SpPFineness_Local(Local_ZXSpPIdx) = FRatio_YMinus*Fineness;
                    end
                end
            end
        end
        for LocalYIdx = 2:Fineness
            Local_ZXSpPIdx = LocalZIdx + Fineness*(LocalXIdx-1) + Fineness^2*(LocalYIdx-2)+(FRatio_YMinus*Fineness)^2;
            SpPPositionIdx_Local(1:3,Local_ZXSpPIdx) = PosIdxOffset ...
                +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-0.5];
            SpPFineness_Local(Local_ZXSpPIdx) = Fineness;
        end
    end
end
end

function [SpPPositionIdx_Local,SpPFineness_Local] = LocallyCalcXYSpPPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
SpPPositionIdx_Local = zeros(3,ElemPer.XYFacePerCoarseGrid(XIdx,YIdx,ZIdx));
SpPFineness_Local    = zeros(1,ElemPer.ZXFacePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_ZMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx-1);
FRatio_ZMinus = Fineness_ZMinus/Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];
for LocalYIdx = 1:Fineness
    for LocalXIdx = 1:Fineness
        for LocalZIdx = 1
            if FRatio_ZMinus <= 1
                Local_XYSpPIdx = LocalXIdx + Fineness*(LocalYIdx-1);
                SpPPositionIdx_Local(1:3,Local_XYSpPIdx) = PosIdxOffset ...
                    +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-0.5;LocalZIdx-1];
                SpPFineness_Local(Local_XYSpPIdx) = Fineness;
            else
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalYIdx = 1:FRatio_ZMinus
                    for LocalLocalXIdx = 1:FRatio_ZMinus
                        Local_XYSpPIdx ...
                            =  LocalLocalXIdx ...
                            +  FRatio_ZMinus*(LocalXIdx-1) ...
                            + (FRatio_ZMinus  *Fineness)*(LocalLocalYIdx-1) ...
                            + (FRatio_ZMinus^2*Fineness)*(LocalYIdx-1);
                        SpPPositionIdx_Local(1:3,Local_XYSpPIdx) = PosIdxOffset + LocalPosIdxOffset ...
                            +(FRatio_ZMinus*Fineness)^(-1)*[LocalLocalXIdx-0.5;LocalLocalYIdx-0.5;0];
                        SpPFineness_Local(Local_XYSpPIdx) = FRatio_ZMinus*Fineness;
                    end
                end
            end
        end
        for LocalZIdx = 2:Fineness
            Local_XYSpPIdx = LocalXIdx + Fineness*(LocalYIdx-1) + Fineness^2*(LocalZIdx-2)+(FRatio_ZMinus*Fineness)^2;
            SpPPositionIdx_Local(1:3,Local_XYSpPIdx) = PosIdxOffset ...
                +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-0.5;LocalZIdx-1];
            SpPFineness_Local(Local_XYSpPIdx) = Fineness;
        end
    end
end
end

%%

function [SpElemPositionIdx_SpS,ElemFineness_SpS] = CalcSpElemPositionIdx_SpS(MeshMeasurements,ElemPer,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

SpElemPositionIdx_SpS = zeros(3,Num_of_Elem.SpS);
ElemFineness_SpS      = zeros(1,Num_of_Elem.SpS);

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpSIdxOffset = ...
                +sum(ElemPer.XEdgePerCoarseGrid(XIdx,1:YIdx-1,ZIdx))...
                +sum(ElemPer.XEdgePerYRow(1:ZIdx-1,XIdx))...
                +sum(ElemPer.XEdgePerYZPlane(1:XIdx-1));
            XSpSIdx = ...
                SpSIdxOffset+1:SpSIdxOffset+ElemPer.XEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpS(:,XSpSIdx),ElemFineness_SpS(XSpSIdx)] ...
                = LocallyCalcXSpSPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpSIdxOffset = ...
                +sum(ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,1:ZIdx-1))...
                +sum(ElemPer.YEdgePerZRow(1:XIdx-1,YIdx))...
                +sum(ElemPer.YEdgePerZXPlane(1:YIdx-1));
            YSpSIdx = ...
                SpSIdxOffset+1:SpSIdxOffset+ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpS(:,YSpSIdx),ElemFineness_SpS(YSpSIdx)] ...
                = LocallyCalcYSpSPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpSIdxOffset = ...
                +sum(ElemPer.ZEdgePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
                +sum(ElemPer.ZEdgePerXRow(1:YIdx-1,ZIdx))...
                +sum(ElemPer.ZEdgePerXYPlane(1:ZIdx-1));
            ZSpSIdx = ...
                SpSIdxOffset+1:SpSIdxOffset+ElemPer.ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpS(:,ZSpSIdx),ElemFineness_SpS(ZSpSIdx)] ...
                = LocallyCalcZSpSPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

end

function [SpSPositionIdxM_Local,SpSFineness_Local] = LocallyCalcXSpSPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
SpSPositionIdxM_Local = zeros(3,ElemNum.XEdgePerCoarseGrid(XIdx,YIdx,ZIdx));
SpSFineness_Local     = zeros(1,ElemNum.XEdgePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_YMinus  = MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx  );
Fineness_ZMinus  = MeshMeasurements.LocalGridFineness(XIdx,YIdx  ,ZIdx-1);
Fineness_YZMinus = MeshMeasurements.LocalGridFineness(XIdx,YIdx-1,ZIdx-1);
FRatio_YMinus  = Fineness_YMinus/Fineness;
FRatio_ZMinus  = Fineness_ZMinus/Fineness;
FRatio_YZMinus = Fineness_YZMinus/Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];
if FRatio_ZMinus<=1 && FRatio_YMinus<=1
    if FRatio_YZMinus<=1       
        for LocalZIdx = 1:Fineness
            for LocalYIdx = 1:Fineness
                for LocalXIdx = 1:Fineness
                    Local_XSpSIdx = LocalYIdx +Fineness*(LocalZIdx-1) +Fineness^2*(LocalXIdx-1)...
                        +sum(ElemPer.XEdgePerCoarseGrid(XIdx,1:YIdx-1,ZIdx))...
                        +sum(ElemPer.XEdgePerYRow(1:ZIdx-1,XIdx))...
                        +sum(ElemPer.XEdgePerYZPlane(1:XIdx-1));
                    SpSPositionIdxM_Local(:,Local_XSpSIdx) = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-1];
                    SpSFineness_Local(Local_XSpSIdx) = Fineness;
                end
            end
        end
    else % if FRatio_YZMinus> 1 
        for LocalXIdx = 1:Fineness
            for LocalZIdx = 1
                for LocalYIdx = 1
                    Local_XSpSIdx = 1:FRatio_YZMinus ...
                        + (Fineness^2-1+FRatio_YZMinus)*(LocalXIdx-1);
                    SpSPositionIdxM_Local(2:3,Local_XSpSIdx) ...
                        = PosIdxOffset(2:3);
                    SpSPositionIdxM_Local(1,Local_XSpSIdx) ...
                        = (FRatio_YZMinus*Fineness)^(-1)*([1:FRatio_YZMinus]-0.5) + (Fineness)^(-1)*(LocalXIdx-1);   
                    SpSFineness_Local(Local_XSpSIdx) = FRatio_YZMinus*Fineness;
                end
                for LocalYIdx = 2:Fineness
                    Local_XSpSIdx ...
                        = LocalYIdx-1+FRatio_YZMinus...
                        +(Fineness^2-1+FRatio_YZMinus)*(LocalXIdx-1);
                    SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                        = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-1];
                    SpSFineness_Local(Local_XSpSIdx) = Fineness;
                end
            end
            for LocalZIdx = 2:Fineness
                for LocalYIdx = 1:Fineness
                    Local_XSpSIdx ...
                        = LocalYIdx ...
                        + Fineness*(LocalZIdx-2)+FRatio_YZMinus+Fineness-1 ...
                        +(Fineness^2-1+FRatio_YZMinus)*(LocalXIdx-1);
                    SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                        = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-1];
                    SpSFineness_Local(Local_XSpSIdx) = Fineness;
                end
            end
        end
    end
elseif FRatio_ZMinus> 1 && FRatio_YMinus<=1
    for LocalXIdx = 1:Fineness
        for LocalYIdx = 1:Fineness
            for LocalZIdx = 1
                LocalXSpSIdx_Offset = (FRatio_ZMinus^2)*(LocalYIdx-1)+0*(LocalZIdx-1)...
                    +(Fineness*(Fineness-1)+FRatio_ZMinus^2*Fineness)*(LocalXIdx);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalXIdx = 1:FRatio_ZMinus
                    for LocalLocalYIdx = 1:FRatio_ZMinus
                       Local_XSpSIdx = LocalXSpSIdx_Offset...
                           +LocalLocalYIdx+FRatio_ZMinus*(LocalLocalXIdx-1);
                       SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                           = PosIdxOffset + LocalPosIdxOffset...
                           +(FRatio_ZMinus*Fineness)^(-1)*[LocalLocalXIdx-0.5;LocalLocalYIdx-1;0];
                       SpSFineness_Local(Local_XSpSIdx) = FRatio_ZMinus*Fineness;
                    end
                end
            end
            for LocalZIdx = 2:Fineness
                Local_XSpSIdx = LocalYIdx + Fineness*(LocalZIdx-2)+(FRatio_ZMinus^2*Fineness)...
                    + (Fineness*(Fineness-1)+FRatio_ZMinus^2*Fineness)*(LocalXIdx-1);
                SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                    = PosIdxOffset...
                    +Fineness^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-1];
                SpSFineness_Local(Local_XSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_ZMinus<=1 && FRatio_YMinus> 1
    for LocalXIdx = 1:Fineness
        for LocalZIdx = 1:Fineness
            for LocalYIdx = 1
                LocalXSpSIdx_Offset = (Fineness-1+FRatio_YMinus)*(LocalZIdx-1)+(Fineness*(Fineness-1)+FRatio_YMinus*Fineness)*(LocalXIdx);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalXIdx = 1:FRatio_YMinus
                    for LocalLocalZIdx = 1:FRatio_YMinus
                        Local_XSpSIdx = LocalXSpSIdx_Offset...
                            + LocalLocalZIdx +  FRatio_YMinus*(LocalLocalXIdx-1);
                        SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio_YMinus*Fineness)^(-1)*[LocalLocalXIdx-0.5;0;LocalLocalZIdx-1];
                        SpSFineness_Local(Local_XSpSIdx) = FRatio_YMinus*Fineness;
                    end
                end
            end
            for LocalYIdx = 2:Fineness
                Local_XSpSIdx  = LocalYIdx-1+FRatio_YMinus^2 ...
                    +(Fineness-1+FRatio_YMinus)*(LocalZIdx-1)...
                    +(Fineness*(Fineness-1)+FRatio_YMinus*Fineness)*(LocalXIdx);
                SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                    = PosIdxOffset ...
                    + Fineness^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-1];
                SpSFineness_Local(Local_XSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_ZMinus> 1 && FRatio_YMinus> 1 && FRatio_ZMinus==FRatio_ZMinus
    FRatio = FRatio_YMinus;
    for LocalXIdx = 1:Fineness
        for LocalZIdx = 1
            for LocalYIdx = 1
                LocalXSpSIdx_Offset = (2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalXIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalXIdx = 1:FRatio   
                    for LocalLocalYIdx = 1:FRatio
                        Local_XSpSIdx = LocalXSpSIdx_Offset...
                            + LocalLocalYIdx + (FRatio+FRatio-1)*(LocalLocalXIdx-1);
                        SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-0.5;LocalLocalYIdx-1;0];
                        SpSFineness_Local(Local_XSpSIdx) = FRatio*Fineness;
                    end
                    for LocalLocalZIdx = 2:FRatio
                        Local_XSpSIdx = LocalXSpSIdx_Offset...
                            + LocalLocalZIdx-1 + FRatio + (FRatio+FRatio-1)*(LocalLocalXIdx-1);
                        SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-0.5;0;LocalLocalZIdx-1];
                        SpSFineness_Local(Local_XSpSIdx) = Fineness;
                    end
                end
            end
            for LocalYIdx = 2:Fineness
                LocalXSpSIdx_Offset = (Fineness^2)*(LocalYIdx-2)+(2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalXIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalXIdx = 1:FRatio
                    for LocalLocalYIdx = 1:FRatio
                        Local_XSpSIdx = LocalXSpSIdx_Offset...
                            + LocalLocalYIdx + FRatio*(LocalLocalXIdx-1);
                        SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-0.5;LocalLocalYIdx-1;0];
                        SpSFineness_Local(Local_XSpSIdx) = FRatio*Fineness;
                    end
                end
            end
        end
        for LocalZIdx = 2:Fineness
            for LocalYIdx = 1
                LocalXSpSIdx_Offset = (FRatio^2+Fineness-1)*(LocalZIdx-2)+(FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalXIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalXIdx = 1:FRatio
                    for LocalLocalZIdx = 1:FRatio
                        Local_XSpSIdx = LocalXSpSIdx_Offset...
                            + LocalLocalZIdx + FRatio*(LocalLocalXIdx-1);
                        SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-0.5;0;LocalLocalZIdx-1];
                        SpSFineness_Local(Local_XSpSIdx) = FRatio*Fineness;
                    end
                end
            end
            for LocalYIdx = 2:Fineness
                Local_XSpSIdx = LocalYIdx-1+FRatio^2 ...
                    +(FRatio^2+Fineness-1)*(LocalZIdx-2) +(FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalXIdx-1);
                SpSPositionIdxM_Local(:,Local_XSpSIdx) ...
                            = PosIdxOffset + (Fineness)^(-1)*[LocalXIdx-0.5;LocalYIdx-1;LocalZIdx-1];
                        SpSFineness_Local(Local_XSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_ZMinus< 1 && FRatio_YMinus< 1 && FRatio_ZMinus~=FRatio_ZMinus
    disp('LocallyCalcXSpSPositionIdx:Unexpected fineness distribution')
    disp('(FRatio_YMinus< 1 && FRatio_ZMinus< 1 && FRatio_YMinus~=FRatio_YMinus)')
    disp(['at [',num2str(XIdx),num2str(YIdx),num2str(ZIdx),']'])
end

end

function [SpSPositionIdxM_Local,SpSFineness_Local] = LocallyCalcYSpSPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
SpSPositionIdxM_Local = zeros(3,ElemNum.YEdgePerCoarseGrid(XIdx,YIdx,ZIdx));
SpSFineness_Local     = zeros(1,ElemNum.YEdgePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_XMinus  = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx  );
Fineness_ZMinus  = MeshMeasurements.LocalGridFineness(XIdx  ,YIdx,ZIdx-1);
Fineness_ZXMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx,ZIdx-1);
FRatio_XMinus  = Fineness_XMinus/Fineness;
FRatio_ZMinus  = Fineness_ZMinus/Fineness;
FRatio_ZXMinus = Fineness_ZXMinus/Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];

if FRatio_XMinus<=1 && FRatio_ZMinus<=1
    if FRatio_ZXMinus<=1       
        for LocalZIdx = 1:Fineness
            for LocalYIdx = 1:Fineness
                for LocalXIdx = 1:Fineness
                    Local_YSpSIdx = LocalZIdx +Fineness*(LocalXIdx-1) +Fineness^2*(LocalYIdx-1)...
                        +sum(ElemPer.YEdgePerCoarseGrid(XIdx,YIdx,1:ZIdx-1))...
                        +sum(ElemPer.YEdgePerZRow(1:XIdx-1,YIdx))...
                        +sum(ElemPer.YEdgePerZXPlane(1:YIdx-1));
                    SpSPositionIdxM_Local(:,Local_YSpSIdx) = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-1];
                    SpSFineness_Local(Local_YSpSIdx) = Fineness;
                end
            end
        end
    else % if FRatio_ZXMinus> 1 
        for LocalYIdx = 1:Fineness
            for LocalXIdx = 1
                for LocalZIdx = 1
                    Local_YSpSIdx = 1:FRatio_ZXMinus ...
                        + (Fineness^2-1+FRatio_ZXMinus)*(LocalYIdx-1);
                    SpSPositionIdxM_Local([1 3],Local_YSpSIdx) ...
                        = PosIdxOffset([1 3]);
                    SpSPositionIdxM_Local(2,Local_YSpSIdx) ...
                        = (FRatio_ZXMinus*Fineness)^(-1)*([1:FRatio_ZXMinus]-0.5) + (Fineness)^(-1)*(LocalYIdx-1);   
                    SpSFineness_Local(Local_YSpSIdx) = FRatio_ZXMinus*Fineness;
                end
                for LocalZIdx = 2:Fineness
                    Local_YSpSIdx ...
                        = LocalZIdx-1+FRatio_ZXMinus...
                        +(Fineness^2-1+FRatio_ZXMinus)*(LocalYIdx-1);
                    SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                        = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-1];
                    SpSFineness_Local(Local_YSpSIdx) = Fineness;
                end
            end
            for LocalXIdx = 2:Fineness
                for LocalZIdx = 1:Fineness
                    Local_YSpSIdx ...
                        = LocalZIdx ...
                        + Fineness*(LocalXIdx-2)+FRatio_ZXMinus+Fineness-1 ...
                        +(Fineness^2-1+FRatio_ZXMinus)*(LocalYIdx-1);
                    SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                        = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-1];
                    SpSFineness_Local(Local_YSpSIdx) = Fineness;
                end
            end
        end
    end
elseif FRatio_XMinus> 1 && FRatio_ZMinus<=1
    for LocalYIdx = 1:Fineness
        for LocalZIdx = 1:Fineness
            for LocalXIdx = 1
                LocalYSpSIdx_Offset = (FRatio_XMinus^2)*(LocalZIdx-1)+0*(LocalXIdx-1)...
                    +(Fineness*(Fineness-1)+FRatio_XMinus^2*Fineness)*(LocalYIdx);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalYIdx = 1:FRatio_XMinus
                    for LocalLocalZIdx = 1:FRatio_XMinus
                       Local_YSpSIdx = LocalYSpSIdx_Offset...
                           +LocalLocalZIdx+FRatio_XMinus*(LocalLocalYIdx-1);
                       SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                           = PosIdxOffset + LocalPosIdxOffset...
                           +(FRatio_XMinus*Fineness)^(-1)*[0;LocalLocalYIdx-0.5;LocalLocalZIdx-1];
                       SpSFineness_Local(Local_YSpSIdx) = FRatio_XMinus*Fineness;
                    end
                end
            end
            for LocalXIdx = 2:Fineness
                Local_YSpSIdx = LocalZIdx + Fineness*(LocalXIdx-2)+(FRatio_XMinus^2*Fineness)...
                    + (Fineness*(Fineness-1)+FRatio_XMinus^2*Fineness)*(LocalYIdx-1);
                SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                    = PosIdxOffset...
                    +Fineness^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-1];
                SpSFineness_Local(Local_YSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_XMinus<=1 && FRatio_ZMinus> 1
    for LocalYIdx = 1:Fineness
        for LocalXIdx = 1:Fineness
            for LocalZIdx = 1
                LocalYSpSIdx_Offset = (Fineness-1+FRatio_ZMinus)*(LocalXIdx-1)+(Fineness*(Fineness-1)+FRatio_ZMinus*Fineness)*(LocalYIdx);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalYIdx = 1:FRatio_ZMinus
                    for LocalLocalXIdx = 1:FRatio_ZMinus
                        Local_YSpSIdx = LocalYSpSIdx_Offset...
                            + LocalLocalXIdx +  FRatio_ZMinus*(LocalLocalYIdx-1);
                        SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio_ZMinus*Fineness)^(-1)*[LocalLocalXIdx-1;LocalLocalYIdx-0.5;0];
                        SpSFineness_Local(Local_YSpSIdx) = FRatio_ZMinus*Fineness;
                    end
                end
            end
            for LocalZIdx = 2:Fineness
                Local_YSpSIdx  = LocalZIdx-1+FRatio_ZMinus^2 ...
                    +(Fineness-1+FRatio_ZMinus)*(LocalXIdx-1)...
                    +(Fineness*(Fineness-1)+FRatio_ZMinus*Fineness)*(LocalYIdx);
                SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                    = PosIdxOffset ...
                    + Fineness^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-1];
                SpSFineness_Local(Local_YSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_XMinus> 1 && FRatio_ZMinus> 1 && FRatio_XMinus==FRatio_XMinus
    FRatio = FRatio_ZMinus;
    for LocalYIdx = 1:Fineness
        for LocalXIdx = 1
            for LocalZIdx = 1
                LocalYSpSIdx_Offset = (2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalYIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalYIdx = 1:FRatio   
                    for LocalLocalXIdx = 1:FRatio
                        Local_YSpSIdx = LocalYSpSIdx_Offset...
                            + LocalLocalXIdx + (FRatio+FRatio-1)*(LocalLocalYIdx-1);
                        SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[0;LocalLocalYIdx-0.5;LocalLocalXIdx-1];
                        SpSFineness_Local(Local_YSpSIdx) = FRatio*Fineness;
                    end
                    for LocalLocalXIdx = 2:FRatio
                        Local_YSpSIdx = LocalYSpSIdx_Offset...
                            + LocalLocalXIdx-1 + FRatio + (FRatio+FRatio-1)*(LocalLocalYIdx-1);
                        SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-1;LocalLocalYIdx-0.5;0];
                        SpSFineness_Local(Local_YSpSIdx) = Fineness;
                    end
                end
            end
            for LocalZIdx = 2:Fineness
                LocalYSpSIdx_Offset = (Fineness^2)*(LocalZIdx-2)+(2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalYIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalYIdx = 1:FRatio
                    for LocalLocalZIdx = 1:FRatio
                        Local_YSpSIdx = LocalYSpSIdx_Offset...
                            + LocalLocalZIdx + FRatio*(LocalLocalYIdx-1);
                        SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[0;LocalLocalYIdx-0.5;LocalLocalZIdx-1];
                        SpSFineness_Local(Local_YSpSIdx) = FRatio*Fineness;
                    end
                end
            end
        end
        for LocalXIdx = 2:Fineness
            for LocalZIdx = 1
                LocalYSpSIdx_Offset = (FRatio^2+Fineness-1)*(LocalXIdx-2)+(FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalYIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalYIdx = 1:FRatio
                    for LocalLocalXIdx = 1:FRatio
                        Local_YSpSIdx = LocalYSpSIdx_Offset...
                            + LocalLocalXIdx + FRatio*(LocalLocalYIdx-1);
                        SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-1;LocalLocalYIdx-0.5;0];
                        SpSFineness_Local(Local_YSpSIdx) = FRatio*Fineness;
                    end
                end
            end
            for LocalZIdx = 2:Fineness
                Local_YSpSIdx = LocalZIdx-1+FRatio^2 ...
                    +(FRatio^2+Fineness-1)*(LocalXIdx-2) +(FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalYIdx-1);
                SpSPositionIdxM_Local(:,Local_YSpSIdx) ...
                            = PosIdxOffset + (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-0.5;LocalZIdx-1];
                        SpSFineness_Local(Local_YSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_XMinus< 1 && FRatio_ZMinus< 1 && FRatio_XMinus~=FRatio_XMinus
    disp('LocallyCalcYSpSPositionIdx:Unexpected fineness distribution')
    disp('(FRatio_ZMinus< 1 && FRatio_XMinus< 1 && FRatio_ZMinus~=FRatio_ZMinus)')
    disp(['at [',num2str(XIdx),num2str(YIdx),num2str(ZIdx),']'])
end

end



function [SpSPositionIdxM_Local,SpSFineness_Local] = LocallyCalcZSpSPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)
SpSPositionIdxM_Local = zeros(3,ElemNum.ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx));
SpSFineness_Local     = zeros(1,ElemNum.ZEdgePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_XMinus  = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx  ,ZIdx);
Fineness_YMinus  = MeshMeasurements.LocalGridFineness(XIdx  ,YIdx-1,ZIdx);
Fineness_XYMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx);
FRatio_XMinus  = Fineness_XMinus/Fineness;
FRatio_YMinus  = Fineness_YMinus/Fineness;
FRatio_XYMinus = Fineness_XYMinus/Fineness;

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];

if FRatio_YMinus<=1 && FRatio_XMinus<=1
    if FRatio_XYMinus<=1       
        for LocalZIdx = 1:Fineness
            for LocalYIdx = 1:Fineness
                for LocalXIdx = 1:Fineness
                    Local_ZSpSIdx = LocalXIdx +Fineness*(LocalYIdx-1) +Fineness^2*(LocalZIdx-1)...
                        +sum(ElemPer.ZEdgePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
                        +sum(ElemPer.ZEdgePerXRow(1:YIdx-1,ZIdx))...
                        +sum(ElemPer.ZEdgePerXYPlane(1:ZIdx-1));
                    SpSPositionIdxM_Local(:,Local_ZSpSIdx) = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-0.5];
                    SpSFineness_Local(Local_ZSpSIdx) = Fineness;
                end
            end
        end
    else % if FRatio_XYMinus> 1 
        for LocalZIdx = 1:Fineness
            for LocalYIdx = 1
                for LocalXIdx = 1
                    Local_ZSpSIdx = 1:FRatio_XYMinus ...
                        + (Fineness^2-1+FRatio_XYMinus)*(LocalZIdx-1);
                    SpSPositionIdxM_Local([1 2],Local_ZSpSIdx) ...
                        = PosIdxOffset([1 2]);
                    SpSPositionIdxM_Local(3,Local_ZSpSIdx) ...
                        = (FRatio_XYMinus*Fineness)^(-1)*([1:FRatio_XYMinus]-0.5) + (Fineness)^(-1)*(LocalYIdx-1);   
                    SpSFineness_Local(Local_ZSpSIdx) = FRatio_XYMinus*Fineness;
                end
                for LocalXIdx = 2:Fineness
                    Local_ZSpSIdx ...
                        = LocalXIdx-1+FRatio_XYMinus...
                        +(Fineness^2-1+FRatio_XYMinus)*(LocalZIdx-1);
                    SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                        = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-0.5];
                    SpSFineness_Local(Local_ZSpSIdx) = Fineness;
                end
            end
            for LocalYIdx = 2:Fineness
                for LocalXIdx = 1:Fineness
                    Local_ZSpSIdx ...
                        = LocalXIdx ...
                        + Fineness*(LocalYIdx-2)+FRatio_XYMinus+Fineness-1 ...
                        +(Fineness^2-1+FRatio_XYMinus)*(LocalZIdx-1);
                    SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                        = PosIdxOffset...
                        +(Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-0.5];
                    SpSFineness_Local(Local_ZSpSIdx) = Fineness;
                end
            end
        end
    end
elseif FRatio_YMinus> 1 && FRatio_XMinus<=1
    for LocalXIdx = 1:Fineness
        for LocalYIdx = 1:Fineness
            for LocalZIdx = 1
                LocalZSpSIdx_Offset = (FRatio_YMinus^2)*(LocalYIdx-1)+0*(LocalZIdx-1)...
                    +(Fineness*(Fineness-1)+FRatio_YMinus^2*Fineness)*(LocalXIdx);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalZIdx = 1:FRatio_YMinus
                    for LocalLocalXIdx = 1:FRatio_YMinus
                       Local_ZSpSIdx = LocalZSpSIdx_Offset...
                           +LocalLocalXIdx+FRatio_YMinus*(LocalLocalZIdx-1);
                       SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                           = PosIdxOffset + LocalPosIdxOffset...
                           +(FRatio_YMinus*Fineness)^(-1)*[LocalLocalXIdx-1;0;LocalLocalZIdx-0.5];
                       SpSFineness_Local(Local_ZSpSIdx) = FRatio_YMinus*Fineness;
                    end
                end
            end
            for LocalZIdx = 2:Fineness
                Local_ZSpSIdx = LocalXIdx + Fineness*(LocalYIdx-2)+(FRatio_YMinus^2*Fineness)...
                    + (Fineness*(Fineness-1)+FRatio_YMinus^2*Fineness)*(LocalZIdx-1);
                SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                    = PosIdxOffset...
                    +Fineness^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-0.5];
                SpSFineness_Local(Local_ZSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_YMinus<=1 && FRatio_XMinus> 1
    for LocalZIdx = 1:Fineness
        for LocalYIdx = 1:Fineness
            for LocalXIdx = 1
                LocalZSpSIdx_Offset = (Fineness-1+FRatio_XMinus)*(LocalYIdx-1)+(Fineness*(Fineness-1)+FRatio_XMinus*Fineness)*(LocalZIdx);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalZIdx = 1:FRatio_XMinus
                    for LocalLocalYIdx = 1:FRatio_XMinus
                        Local_ZSpSIdx = LocalZSpSIdx_Offset...
                            + LocalLocalYIdx +  FRatio_XMinus*(LocalLocalZIdx-1);
                        SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio_XMinus*Fineness)^(-1)*[0;LocalLocalYIdx-1;LocalLocalZIdx-0.5];
                        SpSFineness_Local(Local_ZSpSIdx) = FRatio_XMinus*Fineness;
                    end
                end
            end
            for LocalXIdx = 2:Fineness
                Local_ZSpSIdx  = LocalXIdx-1+FRatio_XMinus^2 ...
                    +(Fineness-1+FRatio_XMinus)*(LocalYIdx-1)...
                    +(Fineness*(Fineness-1)+FRatio_XMinus*Fineness)*(LocalZIdx);
                SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                    = PosIdxOffset ...
                    + Fineness^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-0.5];
                SpSFineness_Local(Local_ZSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_YMinus> 1 && FRatio_XMinus> 1 && FRatio_YMinus==FRatio_YMinus
    FRatio = FRatio_XMinus;
    for LocalZIdx = 1:Fineness
        for LocalYIdx = 1
            for LocalXIdx = 1
                LocalZSpSIdx_Offset = (2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalZIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalZIdx = 1:FRatio   
                    for LocalLocalYIdx = 1:FRatio
                        Local_ZSpSIdx = LocalZSpSIdx_Offset...
                            + LocalLocalYIdx + (FRatio+FRatio-1)*(LocalLocalZIdx-1);
                        SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[0;LocalLocalYIdx-1;LocalLocalZIdx-0.5];
                        SpSFineness_Local(Local_ZSpSIdx) = FRatio*Fineness;
                    end
                    for LocalLocalYIdx = 2:FRatio
                        Local_ZSpSIdx = LocalZSpSIdx_Offset...
                            +0 +LocalLocalYIdx-1+FRatio +(FRatio+FRatio-1)*(LocalLocalZIdx-1);
                        SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[0;LocalLocalYIdx-1;LocalLocalZIdx-0.5];
                        SpSFineness_Local(Local_ZSpSIdx) = Fineness;
                    end
                end
            end
            for LocalXIdx = 2:Fineness
                LocalZSpSIdx_Offset = (Fineness^2)*(LocalXIdx-2)+(2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalZIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalZIdx = 1:FRatio
                    for LocalLocalXIdx = 1:FRatio
                        Local_ZSpSIdx = LocalZSpSIdx_Offset...
                            + LocalLocalXIdx + FRatio*(LocalLocalZIdx-1);
                        SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[LocalLocalXIdx-1;LocalLocalZIdx-0.5];
                        SpSFineness_Local(Local_ZSpSIdx) = FRatio*Fineness;
                    end
                end
            end
        end
        for LocalYIdx = 2:Fineness
            for LocalXIdx = 1
                LocalZSpSIdx_Offset = (FRatio^2+Fineness-1)*(LocalYIdx-2)+(FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalZIdx-1);
                LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
                for LocalLocalZIdx = 1:FRatio
                    for LocalLocalYIdx = 1:FRatio
                        Local_ZSpSIdx = LocalZSpSIdx_Offset...
                            + LocalLocalYIdx + FRatio*(LocalLocalZIdx-1);
                        SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                            = PosIdxOffset + LocalPosIdxOffset...
                            + (FRatio*Fineness)^(-1)*[0;LocalLocalYIdx-1;LocalLocalZIdx-0.5];
                        SpSFineness_Local(Local_ZSpSIdx) = FRatio*Fineness;
                    end
                end
            end
            for LocalXIdx = 2:Fineness
                Local_ZSpSIdx = LocalXIdx-1+FRatio^2 ...
                    +(FRatio^2+Fineness-1)*(LocalYIdx-2) +(FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)...
                    +(2*FRatio^2*(Fineness-1)+2*(FRatio-0.5)*FRatio)*(LocalZIdx-1);
                SpSPositionIdxM_Local(:,Local_ZSpSIdx) ...
                            = PosIdxOffset + (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-0.5];
                        SpSFineness_Local(Local_ZSpSIdx) = Fineness;
            end
        end
    end
elseif FRatio_YMinus< 1 && FRatio_XMinus< 1 && FRatio_YMinus~=FRatio_YMinus
    disp('LocallyCalcZSpSPositionIdx:Unexpected fineness distribution')
    disp('(FRatio_XMinus< 1 && FRatio_YMinus< 1 && FRatio_XMinus~=FRatio_XMinus)')
    disp(['at [',num2str(XIdx),num2str(YIdx),num2str(ZIdx),']'])
end

end

%% 
function [SpElemPositionIdx_SpN,ElemFineness_SpN] = CalcSpElemPositionIdx_SpN(MeshMeasurements,ElemPer,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

SpElemPositionIdx_SpN = zeros(3,Num_of_Elem.SpN);
ElemFineness_SpN      = zeros(1,Num_of_Elem.SpN);

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpNIdxOffset = ...
                +sum(ElemPer.NodePerCoarseGrid(1:XIdx-1,YIdx,ZIdx))...
                +sum(ElemPer.NodePerXRow(1:YIdx-1,ZIdx))...
                +sum(ElemPer.NodePerXYPlane(1:ZIdx-1));
            SpNIdx = ...
                SpNIdxOffset+1:SpNIdxOffset+ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx);
            [SpElemPositionIdx_SpN(:,SpNIdx),ElemFineness_SpN(SpNIdx)] ...
                = LocallyCalcSpNPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer);
        end
    end
end

end

function  [SpNPositionIdx_Local,SpNFineness_Local] = LocallyCalcSpNPositionIdx(XIdx,YIdx,ZIdx,MeshMeasurements,ElemPer)

SpNPositionIdx_Local = zeros(3,ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx));
SpNFineness_Local    = zeros(1,ElemPer.NodePerCoarseGrid(XIdx,YIdx,ZIdx));

Fineness = MeshMeasurements.LocalGridFineness(XIdx,YIdx,ZIdx);
Fineness_XMinus  = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx  ,ZIdx  );
Fineness_YMinus  = MeshMeasurements.LocalGridFineness(XIdx  ,YIdx-1,ZIdx  );
Fineness_ZMinus  = MeshMeasurements.LocalGridFineness(XIdx  ,YIdx  ,ZIdx-1);
Fineness_YZMinus = MeshMeasurements.LocalGridFineness(XIdx  ,YIdx-1,ZIdx-1);
Fineness_ZXMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx  ,ZIdx-1);
Fineness_XYMinus = MeshMeasurements.LocalGridFineness(XIdx-1,YIdx-1,ZIdx  );

FRatio_XMinus  = Fineness_XMinus/Fineness;
FRatio_YMinus  = Fineness_YMinus/Fineness;
FRatio_ZMinus  = Fineness_ZMinus/Fineness;
FRatio_YZMinus = Fineness_YZMinus/Fineness;
FRatio_ZXMinus = Fineness_ZXMinus/Fineness;
FRatio_XYMinus = Fineness_XYMinus/Fineness;

FRatioOnLocalPosYeqZeq0 = max([FRatio_YMinus FRatio_ZMinus FRatio_YZMinus]);
FRatioOnLocalPosZeqXeq0 = max([FRatio_ZMinus FRatio_XMinus FRatio_ZXMinus]);
FRatioOnLocalPosXeqYeq0 = max([FRatio_XMinus FRatio_YMinus FRatio_XYMinus]);

FRatio_OnEdges = [FRatioOnLocalPosYeqZeq0 FRatioOnLocalPosZeqXeq0 FRatioOnLocalPosXeqYeq0];

if any(FRatio_OnEdges(logical(FRatio_OnEdges~=1))~=FRatio_OnEdges(find(FRatio_OnEdges~=1,1)))
    disp('LocallyCalcSpNPositionIdx:Unexpected Fineness pattern')
end

if abs(FRatioOnLocalPosXeqYeq0/FRatio_YMinus-FRatioOnLocalPosXeqYeq0)<EPSILON
    disp('LocallyCalcSpNPositionIdx:Unexpected Fineness pattern ')
end

PosIdxOffset(1:3,1) = [XIdx-1;YIdx-1;ZIdx-1];
Local_SpNNum = 0; 
for LocalZIdx = 1
    for LocalYIdx = 1
        for LocalXIdx = 1
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            for LocalLocalZIdx = 1
                if FRatio_ZMinus<=1
                    Local_SpNIdx = Local_SpNNum + (1:FRatioOnLocalPosYeqZeq0);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*FRatioOnLocalPosYeqZeq0)^(-1)*((1:FRatioOnLocalPosYeqZeq0)-1);
                    SpNPositionIdx_Local([2 3],Local_SpNIdx)...
                        = PosIdxOffset([2 3]) + LocalPosIdxOffset([2 3]);
                    Local_SpNNum = Local_SpNNum + FRatioOnLocalPosYeqZeq0;
                    
                    Local_SpNIdx = Local_SpNNum + (2:FRatioOnLocalPosZeqXeq0);
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                        + (Fineness*FRatioOnLocalPosZeqXeq0)^(-1)*((2:FRatioOnLocalPosZeqXeq0)-1);
                    SpNPositionIdx_Local([1 3],Local_SpNIdx)...
                        = PosIdxOffset([1 3]) + LocalPosIdxOffset([1 3]);
                    Local_SpNNum = Local_SpNNum + FRatioOnLocalPosZeqXeq0-1;
                else
                    Local_SpNIdx = Local_SpNNum + (1:FRatio_ZMinus^2);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*FRatio_ZMinus)^(-1)*(mod((1:FRatio_ZMinus^2)-1,FRatio_ZMinus));
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                        + (Fineness*FRatio_ZMinus)^(-1)*(floor(((1:FRatio_ZMinus^2)-1)/FRatio_ZMinus));
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3);
                    Local_SpNNum = Local_SpNNum + FRatio_ZMinus^2;
                end
            end
            for LocalLocalZIdx = 2:FRatioOnLocalPosXeqYeq0
                if FRatioOnLocalPosXeqYeq0 ~= FRatio_YMinus
                    Local_SpNIdx = Local_SpNNum+1;
                    SpNPositionIdx_Local(3,Local_SpNIdx) ...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                        + (Fineness*FRatioOnLocalPosXeqYeq0)^(-1)*(LocalLocalZIdx-1);
                    SpNPositionIdx_Local([1 2],Local_SpNIdx)...
                        = PosIdxOffset([1 2]) + LocalPosIdxOffset([1 2]);
                    Local_SpNNum = Local_SpNNum + 1;
                else
                    Local_SpNIdx = Local_SpNNum + (1:FRatio_YMinus);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*Ratio_YMinus)^(-1)*((1:Ratio_YMinus)-1);
                    SpNPositionIdx_Local(2,Local_SpNIdx)...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2);
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                        + (Fineness*Ratio_YMinus)^(-1)*(LocalLocalZIdx-1);
                    Local_SpNNum = Local_SpNNum + Ratio_YMinus;
                end
                if FRatioOnLocalPosXeqYeq0 ~= FRatio_XMinus
                    % do nothing
                else
                    Local_SpNIdx = Local_SpNNum + (2:FRatio_YMinus);
                    SpNPositionIdx_Local(1,Local_SpNIdx)...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1);
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset + LocalPosIdxOffset ...
                        + (Fineness*Ratio_YMinus)^(-1)*((2:Ratio_YMinus)-1);
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                        + (Fineness*Ratio_YMinus)^(-1)*(LocalLocalZIdx-1);
                    Local_SpNNum = Local_SpNNum + Ratio_YMinus-1;
                end
            end
        end
        for LocalXIdx = 2:Fineness
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            for LocalLocalZIdx = 1
                if FRatioOnLocalPosYeqZeq0 ~= FRatio_ZMinus
                    Local_SpNIdx = Local_SpNNum + (1:FRatioOnLocalPosYeqZeq0);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*FRatioOnLocalPosYeqZeq0)^(-1)*((1:FRatioOnLocalPosYeqZeq0)-1);
                    SpNPositionIdx_Local([2 3],Local_SpNIdx)...
                        = PosIdxOffset([2 3]) + LocalPosIdxOffset([2 3]);
                    Local_SpNNum = Local_SpNNum + FRatioOnLocalPosYeqZeq0;
                else
                    Local_SpNIdx = Local_SpNNum + (1:FRatio_ZMinus^2);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*FRatio_ZMinus)^(-1)*(mod((1:FRatio_ZMinus^2)-1,FRatio_ZMinus));
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                        + (Fineness*FRatio_ZMinus)^(-1)*(floor(((1:FRatio_ZMinus^2)-1)/FRatio_ZMinus));
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3);
                    Local_SpNNum = Local_SpNNum + FRatio_ZMinus^2;
                end
            end
            for LocalLocalZIdx = 2:FRatio_YMinus
                Local_SpNIdx = Local_SpNNum + (1:FRatio_YMinus);
                SpNPositionIdx_Local(1,Local_SpNIdx) ...
                    = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                    + (Fineness*FRatio_YMinus)^(-1)*((1:FRatio_YMinus)-1);
                SpNPositionIdx_Local(2,Local_SpNIdx)...
                    = PosIdxOffset(2) + LocalPosIdxOffset(2);
                SpNPositionIdx_Local(3,Local_SpNIdx) ...
                    = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                    + (Fineness*FRatio_YMinus)^(-1)*(LocalLocalZIdx-1);
                Local_SpNNum = Local_SpNNum + FRatio_ZMinus^2;
            end
        end
    end
    for LocalYIdx = 2:Fineness
        for LocalXIdx = 1
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            for LocalLocalZIdx = 1
                if FRatioOnLocalPosZeqXeq0 ~= FRatio_ZMinus
                    Local_SpNIdx = Local_SpNNum + (1:FRatioOnLocalPosZeqXeq0);
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                        + (Fineness*FRatioOnLocalPosZeqXeq0)^(-1)*((1:FRatioOnLocalPosZeqXeq0)-1);
                    SpNPositionIdx_Local([1 3],Local_SpNIdx)...
                        = PosIdxOffset([1 3]) + LocalPosIdxOffset([1 3]);
                    Local_SpNNum = Local_SpNNum + FRatioOnLocalPosZeqXeq0;
                else
                    Local_SpNIdx = Local_SpNNum + (1:FRatio_ZMinus^2);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*FRatio_ZMinus)^(-1)*(mod((1:FRatio_ZMinus^2)-1,FRatio_ZMinus));
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                        + (Fineness*FRatio_ZMinus)^(-1)*(floor(((1:FRatio_ZMinus^2)-1)/FRatio_ZMinus));
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3);
                    Local_SpNNum = Local_SpNNum + FRatio_ZMinus^2;
                end
            end
            for LocalLocalZIdx = 2:FRatio_XMinus
                Local_SpNIdx = Local_SpNNum + (1:FRatio_XMinus);
                SpNPositionIdx_Local(1,Local_SpNIdx)...
                    = PosIdxOffset(1) + LocalPosIdxOffset(1);
                SpNPositionIdx_Local(2,Local_SpNIdx) ...
                    = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                    + (Fineness*FRatio_XMinus)^(-1)*((1:FRatio_XMinus)-1);
                SpNPositionIdx_Local(3,Local_SpNIdx) ...
                    = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                    + (Fineness*FRatio_XMinus)^(-1)*(LocalLocalZIdx-1);
                Local_SpNNum = Local_SpNNum + FRatio_ZMinus^2;
            end
        end
        for LocalXIdx = 2:Fineness
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            if FRatio_ZMinus<=1
                Local_SpNIdx = Local_SpNNum + 1;
                SpNPositionIdx_Local(:,Local_SpNIdx) ...
                    = PosIdxOffset + LocalPosIdxOffset;
                Local_SpNNum = Local_SpNNum + 1;
            else
                Local_SpNIdx = Local_SpNNum + (1:FRatio_ZMinus^2);
                SpNPositionIdx_Local(1,Local_SpNIdx) ...
                    = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                    + (Fineness*FRatio_ZMinus)^(-1)*(mod((1:FRatio_ZMinus^2)-1,FRatio_ZMinus));
                SpNPositionIdx_Local(2,Local_SpNIdx) ...
                    = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                    + (Fineness*FRatio_ZMinus)^(-1)*(floor(((1:FRatio_ZMinus^2)-1)/FRatio_ZMinus));
                SpNPositionIdx_Local(3,Local_SpNIdx)...
                    = PosIdxOffset(3) + LocalPosIdxOffset(3);
                Local_SpNNum = Local_SpNNum + FRatio_ZMinus^2;
            end
        end
    end
end
for LocalZIdx = 2:Fineness
    for LocalYIdx = 1
        for LocalXIdx = 1
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            for LocalLocalZIdx = 1:FRatioOnLocalPosXeqYeq0
                if FRatioOnLocalPosXeqYeq0 ~= FRatio_YMinus
                    Local_SpNIdx = Local_SpNNum+1;
                    SpNPositionIdx_Local(3,Local_SpNIdx) ...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                        + (Fineness*FRatioOnLocalPosXeqYeq0)^(-1)*(LocalLocalZIdx-1);
                    SpNPositionIdx_Local([1 2],Local_SpNIdx)...
                        = PosIdxOffset([1 2]) + LocalPosIdxOffset([1 2]);
                    Local_SpNNum = Local_SpNNum + 1;
                else
                    Local_SpNIdx = Local_SpNNum + (1:FRatio_YMinus);
                    SpNPositionIdx_Local(1,Local_SpNIdx) ...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                        + (Fineness*Ratio_YMinus)^(-1)*((1:Ratio_YMinus)-1);
                    SpNPositionIdx_Local(2,Local_SpNIdx)...
                        = PosIdxOffset(2) + LocalPosIdxOffset(2);
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                        + (Fineness*Ratio_YMinus)^(-1)*(LocalLocalZIdx-1);
                    Local_SpNNum = Local_SpNNum + Ratio_YMinus;
                end
                if FRatioOnLocalPosXeqYeq0 ~= FRatio_XMinus
                    % do nothing
                else
                    Local_SpNIdx = Local_SpNNum + (2:FRatio_YMinus);
                    SpNPositionIdx_Local(1,Local_SpNIdx)...
                        = PosIdxOffset(1) + LocalPosIdxOffset(1);
                    SpNPositionIdx_Local(2,Local_SpNIdx) ...
                        = PosIdxOffset + LocalPosIdxOffset ...
                        + (Fineness*Ratio_YMinus)^(-1)*((2:Ratio_YMinus)-1);
                    SpNPositionIdx_Local(3,Local_SpNIdx)...
                        = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                        + (Fineness*Ratio_YMinus)^(-1)*(LocalLocalZIdx-1);
                    Local_SpNNum = Local_SpNNum + Ratio_YMinus-1;
                end                    
            end  
        end
        for LocalXIdx = 2:Fineness
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            if FRatio_YMinus<=1
                Local_SpNIdx = Local_SpNNum + 1;
                SpNPositionIdx_Local(:,Local_SpNIdx) ...
                    = PosIdxOffset + LocalPosIdxOffset;
                Local_SpNNum = Local_SpNNum + 1;
            else
                Local_SpNIdx = Local_SpNNum + (1:FRatio_YMinus^2);
                SpNPositionIdx_Local(1,Local_SpNIdx) ...
                    = PosIdxOffset(1) + LocalPosIdxOffset(1) ...
                    + (Fineness*FRatio_YMinus)^(-1)*(mod((1:FRatio_YMinus^2)-1,FRatio_YMinus));
                SpNPositionIdx_Local(2,Local_SpNIdx)...
                    = PosIdxOffset(2) + LocalPosIdxOffset(2);
                SpNPositionIdx_Local(3,Local_SpNIdx) ...
                    = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                    + (Fineness*FRatio_YMinus)^(-1)*(floor(((1:FRatio_YMinus^2)-1)/FRatio_YMinus));
                Local_SpNNum = Local_SpNNum + FRatio_YMinus^2;
            end
        end
    end
    for LocalYIdx = 2:Fineness
        for LocalXIdx = 1
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            if FRatio_XMinus <= 1
                Local_SpNIdx = Local_SpNNum + 1;
                SpNPositionIdx_Local(:,Local_SpNIdx) ...
                    = PosIdxOffset + LocalPosIdxOffset;
                Local_SpNNum = Local_SpNNum + 1;
            else
                Local_SpNIdx = Local_SpNNum + (1:FRatio_XMinus^2);
                SpNPositionIdx_Local(1,Local_SpNIdx)...
                    = PosIdxOffset(1) + LocalPosIdxOffset(1);
                SpNPositionIdx_Local(2,Local_SpNIdx) ...
                    = PosIdxOffset(2) + LocalPosIdxOffset(2) ...
                    + (Fineness*FRatio_XMinus)^(-1)*(mod((1:FRatio_XMinus^2)-1,FRatio_XMinus));
                 SpNPositionIdx_Local(3,Local_SpNIdx) ...
                    = PosIdxOffset(3) + LocalPosIdxOffset(3) ...
                    + (Fineness*FRatio_XMinus)^(-1)*(floor(((1:FRatio_XMinus^2)-1)/FRatio_XMinus));
                Local_SpNNum = Local_SpNNum + FRatio_YMinus^2;
            end
        end
        for LocalXIdx = 2:Fineness
            LocalPosIdxOffset(1:3,1) = (Fineness)^(-1)*[LocalXIdx-1;LocalYIdx-1;LocalZIdx-1];
            Local_SpNIdx = Local_SpNNum + 1;
            SpNPositionIdx_Local(:,Local_SpNIdx) ...
                = PosIdxOffset + LocalPosIdxOffset;
            Local_SpNNum = Local_SpNNum + 1;
        end
    end
end
end