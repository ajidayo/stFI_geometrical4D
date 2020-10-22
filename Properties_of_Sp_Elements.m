function [SpElemProperties,STElemProperties,Num_of_Elem,SigmaSpV] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos)
% global EPSILON
disp('Properties_of_Sp_Elements: Determing Updating Numbers per timestep for each Sp-Element')
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

disp('Properties_of_Sp_Elements: Defining reference relations between Sp/ST Elements, ')
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
disp('Properties_of_Sp_Elements: Determing Dof-calculating tasks by which each Field-DoF is updated.')
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

SpElemProperties.SpV.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpV));
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpV = find(sD(:,SpPIdx).')
        SpElemProperties.SpV.Belong_to_ST_FI(IncSpV) = true;
        SpElemProperties.SpV.IncToSTFISpP(IncSpV) = true;
    end
end

%% PML
NodePosDual_M = zeros(3,Num_of_Elem.SpV) ;
for SpVIdx = 1:Num_of_Elem.SpV
    NodePosDual_M(:,SpVIdx) = NodePos.Dual(SpVIdx).Vec;
end

SigmaSpV = zeros(3,Num_of_Elem.SpV); 
SpElemProperties_SpS_PML        = false(1,Num_of_Elem.SpS);
SpElemProperties_SpP_PML        = false(1,Num_of_Elem.SpP);
for SpVIdx = 1:Num_of_Elem.SpV
    if mod(SpVIdx,round(0.2*Num_of_Elem.SpV)) == 0
        disp(['Properties_of_Sp_Elements: Computing SigmaSpV. Current progress - ', num2str(100*SpVIdx/Num_of_Elem.SpV), "%"])
    end
    SigmaSpV(:,SpVIdx) = func_sigma_123(NodePosDual_M(:,SpVIdx));
    if norm(SigmaSpV(:,SpVIdx))>0
        for IncSpPIdx = find(sD(SpVIdx,:))
            SpElemProperties_SpP_PML(IncSpPIdx)=true;
            for IncSpSIdx = find(sC(IncSpPIdx,:))
                SpElemProperties_SpS_PML(IncSpSIdx)=true;
            end
        end
    end
end
SpElemProperties.SpS.PML = SpElemProperties_SpS_PML;
SpElemProperties.SpP.PML = SpElemProperties_SpP_PML;

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

Num_of_Elem.PMLSpP = size(find(SpElemProperties.SpP.PML==true),2);
Num_of_Elem.PMLSpS = size(find(SpElemProperties.SpS.PML==true),2);


%% Check Exclusivity of TaskTypes
if any(SpElemProperties.SpP.Belong_to_ST_FI.* SpElemProperties.SpP.PML)
    warning('More than one Sp-faces have both STFI and PML property (prohibited)')
end
if any(SpElemProperties.SpS.Belong_to_ST_FI.*SpElemProperties.SpS.PML)
    warning('More than one Sp-edges have both STFI and PML property (prohibited)')
end

%% position of Primal Faces
% disp('Properties_of_Sp_Elements: Roughly calculating the position of each primal face.')
% NodePos_Prim_Matrix = zeros(3,Num_of_Elem.SpN) ;
% for SpNIdx = 1:Num_of_Elem.SpN
%     NodePos_Prim_Matrix(:,SpNIdx) = NodePos.Prim(SpNIdx).Vec;
% end
% PrimFacePos(Num_of_Elem.SpP).Vec = [0;0;0];
% for SpPIdx = 1:Num_of_Elem.SpP
%     if mod(SpPIdx,round(0.1*Num_of_Elem.SpP)) == 1
%         disp(['Processing SpPIdx = ', num2str(SpPIdx),'/',num2str(Num_of_Elem.SpP)])
%     end
%     Nodes_LogIdx = logical(logical(sC(SpPIdx,:))*logical(sG));
%     PrimFacePos(SpPIdx).Vec = mean(NodePos_Prim_Matrix(:,Nodes_LogIdx),2);
% end
end