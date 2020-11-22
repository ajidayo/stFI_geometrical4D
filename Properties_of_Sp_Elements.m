function [SpElemProperties,STElemProperties,Num_of_Elem,SigmaSpV] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePosDual_M,FuncHandle_sigma123)
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
SpElemProperties_SpP_Belong_to_ST_FI = false(1,Num_of_Elem.SpP);
SpElemProperties_SpS_Belong_to_ST_FI = false(1,Num_of_Elem.SpS);
SpElemProperties_SpP_UpdNumBoundary  = false(1,Num_of_Elem.SpP);
SpElemProperties_SpS_UpdNumBoundary  = false(1,Num_of_Elem.SpS);
SpElemProperties_SpV_Belong_to_ST_FI = false(1,Num_of_Elem.SpV);
SpElemProperties_SpV_IncToSTFISpP    = false(1,Num_of_Elem.SpV);
for SpPIdx = 1:Num_of_Elem.SpP
    if mod(SpPIdx,round(0.2*Num_of_Elem.SpP)) == 1
        disp(['Properties_of_Sp_Elements: Detecting STFI SpPs, current progress - ', num2str(100*SpPIdx/Num_of_Elem.SpP), '%'])
    end
    if size(unique(SpElemProperties.SpV.UpdNum(logical(sD(:,SpPIdx).'))),2)~=1
        SpElemProperties_SpP_Belong_to_ST_FI(SpPIdx) = true;
        SpElemProperties_SpP_UpdNumBoundary(SpPIdx)  = true;
        IncSpSLogIdx = logical(sC(SpPIdx,:));
        SpElemProperties_SpS_UpdNumBoundary(IncSpSLogIdx)  = true;
        SpElemProperties_SpS_Belong_to_ST_FI(IncSpSLogIdx) = true;
    end
end
SpSsAdjToSpS_viaSpNs = logical(logical(sG)*logical(sG).');
for SpSIdx = find(SpElemProperties_SpS_Belong_to_ST_FI==true)
    SpElemProperties_SpP_Belong_to_ST_FI(logical(sC(:,SpSIdx).')) = true;
    SpElemProperties_SpS_Belong_to_ST_FI(SpSsAdjToSpS_viaSpNs(SpSIdx,:)) = true;
end
for SpPIdx = find(SpElemProperties_SpP_Belong_to_ST_FI)
    IncSpVLogIdx = logical(sD(:,SpPIdx).');
    SpElemProperties_SpV_Belong_to_ST_FI(IncSpVLogIdx) = true;
    SpElemProperties_SpV_IncToSTFISpP(IncSpVLogIdx) = true;
end
SpElemProperties.SpP.Belong_to_ST_FI = sparse(SpElemProperties_SpP_Belong_to_ST_FI);
SpElemProperties.SpS.Belong_to_ST_FI = sparse(SpElemProperties_SpS_Belong_to_ST_FI);
SpElemProperties.SpP.UpdNumBoundary  = sparse(SpElemProperties_SpP_UpdNumBoundary);
SpElemProperties.SpS.UpdNumBoundary  = sparse(SpElemProperties_SpS_UpdNumBoundary);
SpElemProperties.SpV.Belong_to_ST_FI = sparse(SpElemProperties_SpV_Belong_to_ST_FI); 
SpElemProperties.SpV.IncToSTFISpP    = sparse(SpElemProperties_SpV_IncToSTFISpP);

%% PML
SigmaSpV = zeros(3,Num_of_Elem.SpV); 
SpElemProperties_SpS_PML        = false(1,Num_of_Elem.SpS);
SpElemProperties_SpP_PML        = false(1,Num_of_Elem.SpP);
IncIncSpS = logical(logical(sD)*logical(sC));
for SpVIdx = 1:Num_of_Elem.SpV
    if mod(SpVIdx,round(0.2*Num_of_Elem.SpV)) == 1
        disp(['Properties_of_Sp_Elements: Computing SigmaSpV, current progress - ', num2str(100*SpVIdx/Num_of_Elem.SpV), '%'])
    end
    SigmaSpV(:,SpVIdx) = FuncHandle_sigma123(NodePosDual_M(:,SpVIdx));
    if norm(SigmaSpV(:,SpVIdx))>0
        SpElemProperties_SpP_PML(logical(sD(SpVIdx,:)))=true;
        SpElemProperties_SpS_PML(IncIncSpS(SpVIdx,:))=true;
    end
end
SpElemProperties.SpS.PML = SpElemProperties_SpS_PML;
SpElemProperties.SpP.PML = SpElemProperties_SpP_PML;

Num_of_Elem.PMLSpP = size(find(SpElemProperties.SpP.PML==true),2);
Num_of_Elem.PMLSpS = size(find(SpElemProperties.SpS.PML==true),2);
Num_of_Elem.TilikePMLImagDualSTP = sum((SpElemProperties.SpP.UpdNum(logical(SpElemProperties.SpP.PML))+1),1);
Num_of_Elem.SplikePMLImagDualSTP = sum((SpElemProperties.SpS.UpdNum(logical(SpElemProperties.SpS.PML))+1),1);
Num_of_Elem.PMLImagDualSTP       = Num_of_Elem.TilikePMLImagDualSTP+Num_of_Elem.SplikePMLImagDualSTP;

SpElemProperties.SpP.FirstPMLImagDualSTP=sparse(1,Num_of_Elem.SpP);
SpElemProperties.SpS.FirstPMLImagDualSTP=sparse(1,Num_of_Elem.SpS);
PMLImagDualSTPIdx = 1;
for SpPIdx = find(SpElemProperties.SpP.PML==true)
    SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx) = PMLImagDualSTPIdx;
    PMLImagDualSTPIdx = PMLImagDualSTPIdx + SpElemProperties.SpP.UpdNum(SpPIdx)+1;
end
for SpSIdx = find(SpElemProperties.SpS.PML==true)
    SpElemProperties.SpS.FirstPMLImagDualSTP(SpSIdx) = PMLImagDualSTPIdx;
    PMLImagDualSTPIdx = PMLImagDualSTPIdx + SpElemProperties.SpS.UpdNum(SpSIdx)+1;
end

%% find TF-SF Boundary
disp('Properties_of_Sp_Elements: Detecting TF-SF Boundaries')
SpElemProperties.SpP.isInSFRegion     = true(1,Num_of_Elem.SpP);
SpElemProperties.SpS.isInSFRegion     = true(1,Num_of_Elem.SpS);
for SpVIdx = find(SpElemProperties.SpV.isInSFRegion==false)
    SpElemProperties.SpP.isInSFRegion(logical(sD(SpVIdx,:))) = false;
end
for SpPIdx = find(SpElemProperties.SpP.isInSFRegion==false)
    SpElemProperties.SpS.isInSFRegion(logical(sC(SpPIdx,:))) = false;
end
SpElemProperties.SpP.isOnTFSFBoundary = false(1,Num_of_Elem.SpP);
SpElemProperties.SpS.isOnTFSFBoundary = false(1,Num_of_Elem.SpS);
for SpPIdx = 1:Num_of_Elem.SpP
    if size(unique(SpElemProperties.SpV.isInSFRegion(logical(sD(:,SpPIdx).'))),2) ~=1
        SpElemProperties.SpP.isOnTFSFBoundary(SpPIdx) = true;
        SpElemProperties.SpS.isOnTFSFBoundary(logical(sC(SpPIdx,:))) = true;
    end
end


%% Check Exclusivity of TaskTypes
if any(SpElemProperties.SpP.Belong_to_ST_FI.* SpElemProperties.SpP.PML)
    warning('More than one Sp-faces have both STFI and PML property (prohibited)')
end
if any(SpElemProperties.SpS.Belong_to_ST_FI.*SpElemProperties.SpS.PML)
    warning('More than one Sp-edges have both STFI and PML property (prohibited)')
end
if any(SpElemProperties.SpP.isOnTFSFBoundary.*SpElemProperties.SpP.PML)
    warning('More than one Sp-faces have both isOnTFSFBoundary and PML property (prohibited)')
end
if any(SpElemProperties.SpS.isOnTFSFBoundary.*SpElemProperties.SpS.PML)
    warning('More than one Sp-edges have both isOnTFSFBoundary and PML property (prohibited)')
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