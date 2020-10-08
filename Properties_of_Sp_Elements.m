function [SpElemProperties,STElemProperties,Num_of_Elem,PrimFacePos] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos)
%% UpdNum
UpdNum_SpV = SpElemProperties.SpV.UpdNum;
UpdNum_SpP = zeros(Num_of_Elem.SpP,1);
UpdNum_SpS = zeros(Num_of_Elem.SpS,1);
UpdNum_SpN = zeros(Num_of_Elem.SpN,1);
for SpPIdx=1:Num_of_Elem.SpP
    UpdNum_SpP(SpPIdx) = max(UpdNum_SpV(logical(sD(:,SpPIdx))));
end
for SpSIdx=1:Num_of_Elem.SpS
    UpdNum_SpS(SpSIdx) = max(UpdNum_SpP(logical(sC(:,SpSIdx))));
end
for SpNIdx=1:Num_of_Elem.SpN
    UpdNum_SpN(SpNIdx) = max(UpdNum_SpS(logical(sG(:,SpNIdx))));
end
SpElemProperties.SpP.UpdNum=UpdNum_SpP;
SpElemProperties.SpS.UpdNum=UpdNum_SpS;
SpElemProperties.SpN.UpdNum=UpdNum_SpN;
%% FirstSTNIdx, RefSpN
CurrentSTNIdx=1;
for SpNIdx=1:Num_of_Elem.SpN
    SpElemProperties.SpN.FirstSTNIdx(SpNIdx) = CurrentSTNIdx;
    for STNIdx = SpElemProperties.SpN.FirstSTNIdx(SpNIdx):...
            SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+UpdNum_SpN(SpNIdx)
        STElemProperties.STN.RefSpN(STNIdx) = SpNIdx;
    end
    CurrentSTNIdx = CurrentSTNIdx +UpdNum_SpN(SpNIdx)+1;
end
Num_of_Elem.STN = CurrentSTNIdx -1;
%% FirstSpSIdx
CurrentSTSIdx = 1;
for SpSIdx = 1:Num_of_Elem.SpS
    SpElemProperties.SpS.FirstSTSIdx(SpSIdx) = CurrentSTSIdx;
    CurrentSTSIdx = CurrentSTSIdx + UpdNum_SpS(SpSIdx) +1;
end
Num_of_Elem.SplikeSTS = CurrentSTSIdx -1;
for SpNIdx=1:Num_of_Elem.SpN
    SpElemProperties.SpN.FirstSTSIdx(SpNIdx) = CurrentSTSIdx;
    CurrentSTSIdx = CurrentSTSIdx + UpdNum_SpN(SpNIdx) +1;    
end
Num_of_Elem.STS = CurrentSTSIdx -1;
Num_of_Elem.TilikeSTS = Num_of_Elem.STS-Num_of_Elem.SplikeSTS;

%% FirstSTPIdx
CurrentSTPIdx=1;
for SpPIdx=1:Num_of_Elem.SpP
    SpElemProperties.SpP.FirstSTPIdx(SpPIdx) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpP(SpPIdx) +1;
end
Num_of_Elem.SplikeSTP = CurrentSTPIdx -1;
for SpSIdx=1:Num_of_Elem.SpS
    SpElemProperties.SpS.FirstSTPIdx(SpSIdx) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpS(SpSIdx) +1;
end
Num_of_Elem.STP = CurrentSTPIdx -1;
Num_of_Elem.TilikeSTP=Num_of_Elem.STP-Num_of_Elem.SplikeSTP;

%% FirstSTOmegaIdx
CurrentSTOmegaIdx =1;
for SpVIdx =1:Num_of_Elem.SpV
    SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx) = CurrentSTOmegaIdx;
    CurrentSTOmegaIdx = CurrentSTOmegaIdx + SpElemProperties.SpV.UpdNum(SpVIdx) +1;
end
Num_of_Elem.SplikeSTOmega = CurrentSTOmegaIdx-1;

for SpSIdx =1:Num_of_Elem.SpP
    SpElemProperties.SpP.FirstSTOmegaIdx(SpSIdx) = CurrentSTOmegaIdx;
    CurrentSTOmegaIdx = CurrentSTOmegaIdx + SpElemProperties.SpP.UpdNum(SpSIdx) +1;
end
Num_of_Elem.STOmega = CurrentSTOmegaIdx-1;
Num_of_Elem.TilikeSTOmega = Num_of_Elem.STOmega-Num_of_Elem.SplikeSTOmega;
%% STV
CurrentSTVIdx=0;
for SpVIdx=1:Num_of_Elem.SpV
    UpdNum = UpdNum_SpV(SpVIdx);
    for TimeIdx = -0.5/UpdNum
        CurrentSTVIdx = CurrentSTVIdx+1;
        SpElemProperties.SpV.FirstSTVIdx(SpVIdx)    = CurrentSTVIdx;
        STElemProperties.STV.TimeIdx(CurrentSTVIdx) = TimeIdx;
        STElemProperties.STV.RefSpV(CurrentSTVIdx)  = SpVIdx;
    end
    for TimeIdx = 0.5/UpdNum:1/UpdNum:(UpdNum-0.5)/UpdNum
        CurrentSTVIdx = CurrentSTVIdx+1;
        STElemProperties.STV.TimeIdx(CurrentSTVIdx) = TimeIdx;
        STElemProperties.STV.RefSpV(CurrentSTVIdx)  = SpVIdx;
    end
end
Num_of_Elem.STV = CurrentSTVIdx;

STElemProperties.STP.TaskIdx = zeros(1,Num_of_Elem.STP);
%% tasktype
% Usage Guide: SpElemPropreties.SpP.Belong_to_ST_FI(SpPIdx)
% First check if SpPs are on dt-Interfaces. Interface SpPs belongs to ST_FI-region.
% SpPs adjacent to interface-SpPs belongs to ST_FI-region.
SpElemProperties.SpP.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpP));
SpElemProperties.SpS.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpS));
SpElemProperties.SpP.UpdNumBoundary  = logical(sparse(1,Num_of_Elem.SpP));
SpElemProperties.SpS.UpdNumBoundary  = logical(sparse(1,Num_of_Elem.SpS));
% AdjM_SpP_via_SpV = sD.'*sD;
for SpPIdx = 1:Num_of_Elem.SpP
    Sum = 0;
    IncSpV_List = find(sD(:,SpPIdx)).';
    for IncSpV = IncSpV_List
        Sum = Sum + SpElemProperties.SpV.UpdNum(IncSpV);
    end
    if Sum == size(IncSpV_List,2)*SpElemProperties.SpV.UpdNum(IncSpV_List(1))
        SpElemProperties.SpP.Belong_to_ST_FI(SpPIdx) = false;
    else
        SpElemProperties.SpP.Belong_to_ST_FI(SpPIdx) = true;
        SpElemProperties.SpP.UpdNumBoundary(SpPIdx)  = true;
        %        SpElemProperties.SpP.Belong_to_ST_FI(logical(AdjM_SpP_via_SpV(SpPIdx,:))) = true;
        for IncSpS = find(sC(SpPIdx,:))
            SpElemProperties.SpS.UpdNumBoundary(IncSpS)  = true;
            SpElemProperties.SpS.Belong_to_ST_FI(IncSpS) = true;
        end
    end
end
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI==true)
    for IncSpP = find(sC(:,SpSIdx)).'
        SpElemProperties.SpP.Belong_to_ST_FI(IncSpP) = true;
    end
    for AdjSpS = find(logical(sG(SpSIdx,:))*logical(sG).')
        SpElemProperties.SpS.Belong_to_ST_FI(AdjSpS) = true;
    end
end
% SpElemProperties.SpN.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpN));
% for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
%     for IncSpN = find(sG(SpSIdx,:))
%         SpElemProperties.SpN.Belong_to_ST_FI(IncSpN) = true;
%     end
% end
% SpElemProperties.SpV.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpV));
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpV = find(sD(:,SpPIdx).')
        SpElemProperties.SpV.Belong_to_ST_FI(IncSpV) = true;
    end
end
%dummy 
SpElemProperties.SpS.UpdNumCorner = logical(sparse(1,Num_of_Elem.SpS));

%% PML, Dummy
SpElemProperties.SpS.PML        = false(1,Num_of_Elem.SpS);
SpElemProperties.SpP.PML        = false(1,Num_of_Elem.SpP);
VolPerXRow          = XSize;
VolPerXYPlane       = XSize*YSize;
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            x = (XIdx-0.5)*MeshMeasurements.dx;
            y = (YIdx-0.5)*MeshMeasurements.dy;
            z = (ZIdx-0.5)*MeshMeasurements.dz;
            if  abs(func_sigma_123(1,x,y,z))>EPSILON || abs(func_sigma_123(2,x,y,z))>EPSILON || abs(func_sigma_123(3,x,y,z))>EPSILON
                SpVIdx = XIdx + (YIdx-1)*VolPerXRow + (ZIdx-1)*VolPerXYPlane;
                for IncSpPIdx = find(sD(SpVIdx,:))
                    SpElemProperties.SpP.PML(IncSpPIdx)=true;
                    for IncSpSIdx = find(sC(IncSpPIdx,:))
                        SpElemProperties.SpS.PML(IncSpSIdx)=true;
                    end
                end
            end
        end
    end
end

SpElemProperties.SpP.FirstPMLImagDualSTP=sparse(1,Num_of_Elem.SpP);
SpElemProperties.SpS.FirstPMLImagDualSTP=sparse(1,Num_of_Elem.SpS);

Num_of_PMLTilikeImagDualSTP=0;
for SpPIdx = find(SpElemProperties.SpP.PML)
    Num_of_PMLTilikeImagDualSTP=Num_of_PMLTilikeImagDualSTP+1+SpElemProperties.SpP.UpdNum(SpPIdx);
end
Num_of_Elem.TilikePMLImagDualSTP = Num_of_PMLTilikeImagDualSTP;
Num_of_PMLSplikeImagDualSTP=0;
for SpSIdx = find(SpElemProperties.SpS.PML)
    Num_of_PMLSplikeImagDualSTP=Num_of_PMLSplikeImagDualSTP+1+SpElemProperties.SpS.UpdNum(SpSIdx);
end
Num_of_Elem.SplikePMLImagDualSTP = Num_of_PMLSplikeImagDualSTP;
Num_of_Elem.PMLImagDualSTP= Num_of_Elem.TilikePMLImagDualSTP+Num_of_Elem.SplikePMLImagDualSTP;


PMLImagDualSTPIdx = 1;
for SpPIdx = find(SpElemProperties.SpP.PML==true)
    SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx) = PMLImagDualSTPIdx;
    PMLImagDualSTPIdx = PMLImagDualSTPIdx + SpElemProperties.SpP.UpdNum(SpPIdx)+1;
end
for SpSIdx = find(SpElemProperties.SpS.PML==true)
    SpElemProperties.SpS.FirstPMLImagDualSTP(SpSIdx) = PMLImagDualSTPIdx;
    PMLImagDualSTPIdx = PMLImagDualSTPIdx + SpElemProperties.SpS.UpdNum(SpSIdx)+1;
end


%% position of Primal Faces
PrimFacePos(Num_of_Elem.SpP).Vec = [0;0;0];
for SpPIdx = 1:Num_of_Elem.SpP
    Nodes = find(( logical(sG).'*logical(sC(SpPIdx,:)).' ).');
    PosVec_SpP = [0;0;0];
    for SpNIdx = Nodes
        PosVec_SpP = PosVec_SpP+NodePos.Prim(SpNIdx).Vec;
    end
    PrimFacePos(SpPIdx).Vec = PosVec_SpP/size(Nodes,2);
end
end