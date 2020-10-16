function [kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties,Num_of_Elem)
disp('ComputeKappa_4D_ST: Computing the value of kappa for each ST-face pair in duality relations.')
global SpDIM
kappa = zeros(Num_of_Elem.STP,1);

FaceArea.Prim = zeros(Num_of_Elem.SpP,1);
FaceArea.Dual = zeros(Num_of_Elem.SpS,1);
DeltaSTNPos   = sparse(SpDIM, Num_of_Elem.STN);
[kappa,FaceArea] = ComputeKappa_for_STFI_SpPs_and_SpSs(kappa,FaceArea,DeltaSTNPos,cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties);
%disp('call ComputeKappa_for_NonSTFI_SpSs')
[kappa,FaceArea.Dual] = ComputeKappa_for_NonSTFI_SpSs(kappa,FaceArea.Dual,cdt,sG,sC,sD,NodePos,SpElemProperties);
%disp('call ComputeKappa_for_NonSTFI_SpPs')
[kappa,FaceArea.Prim] = ComputeKappa_for_NonSTFI_SpPs(kappa,FaceArea.Prim,cdt,sG,sC,sD,NodePos,SpElemProperties);
end
%%
function [kappa,FaceArea] = ComputeKappa_for_STFI_SpPs_and_SpSs(kappa,FaceArea,DeltaSTNPos,cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties)
disp('ComputeKappa_for_STFI_SpPs_and_SpSs: Computing Kappa For STFI edges')
for SpSIdx = find(SpElemProperties.SpS.UpdNumBoundary)
    if any(SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion(logical(sG(SpSIdx,:))))
        continue;
    elseif SpElemProperties.SpS.PEC(SpSIdx) == true
        STPIdx_RefferingSpS = ...
            SpElemProperties.SpS.FirstSTPIdx(SpSIdx):SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
        kappa(STPIdx_RefferingSpS) = -1;
        FaceArea.Dual(SpSIdx) = 1;
        continue;
    end
    PrimEdgeTgtLeng = SpEdgeLength(SpSIdx,sG,NodePos.Prim);
    DualFaceTgtArea = SpFaceArea(SpSIdx,sD.',sC.',NodePos.Dual);
    Area_PrimSTP = PrimEdgeTgtLeng*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
    Area_DualSTP = DualFaceTgtArea;
    STPIdxList_RefferingSpS = ...
        SpElemProperties.SpS.FirstSTPIdx(SpSIdx):SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    kappa(STPIdxList_RefferingSpS) = -1*Area_DualSTP/Area_PrimSTP;
    FaceArea.Dual(SpSIdx) = Area_DualSTP;
    DeltaSTNPos = ...
        Compute_DeltaSTNPos_for_SingleBoundaryEdge(SpSIdx,DeltaSTNPos,NodePos,PrimEdgeTgtLeng,DualFaceTgtArea,sG,sC,D2,D3,SpElemProperties,STElemProperties,cdt);
end

for SpPIdx = find(SpElemProperties.SpP.UpdNumBoundary) % compute Delta for SpSs on outer boundaries
    if ~any(SpElemProperties.SpS.PEC(logical(sC(SpPIdx,:))))
        continue;
    end
    SpNIdx_DefPrimFace = find(logical(sC(SpPIdx,:))*logical(sG));
    for Time = 1:SpElemProperties.SpP.UpdNum(SpPIdx)-1
        DeltaSTN_Fetch = [0;0;0];
        for SpNIdx = SpNIdx_DefPrimFace
            STNIdx = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+Time;
            if norm(DeltaSTNPos(:,STNIdx))~=0
                DeltaSTN_Fetch = DeltaSTNPos(:,STNIdx);
                break;
            end
        end
        if norm(DeltaSTN_Fetch) == 0
            warning('ComputeKappa_4D_ST:norm(DeltaSTN_Fetch) == 0')
           % pause;
        end
        for SpNIdx = SpNIdx_DefPrimFace
            STNIdx = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+Time;
            DeltaSTNPos(:,STNIdx) = DeltaSTN_Fetch;
        end
    end
end

for SpNIdx = find(SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion)
    UpdNum_tgt = SpElemProperties.SpN.UpdNum(SpNIdx);
    for Time = 1:UpdNum_tgt-1
        STNIdx_tgt = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+Time;
        DeltaSTNPos(:,STNIdx_tgt) = [0;0;0];
        AdjSpNIdx_List = find(logical(sG(:,SpNIdx).')*logical(sG));
        AdjSpNIdx_List = AdjSpNIdx_List(logical(SpElemProperties.SpN.UpdNum(AdjSpNIdx_List)==UpdNum_tgt));
        AdjSpN_NotOnCorner_List =  AdjSpNIdx_List(~logical(SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion(AdjSpNIdx_List)));
        if size(AdjSpN_NotOnCorner_List,2) > 0
            for AdjSpN_NotOnCorner_Idx = AdjSpN_NotOnCorner_List
                STNIdx_Fetch = SpElemProperties.SpN.FirstSTNIdx(AdjSpN_NotOnCorner_Idx)+Time;
                DeltaSTNPos(:,STNIdx_tgt) = DeltaSTNPos(:,STNIdx_tgt) + DeltaSTNPos(:,STNIdx_Fetch);
            end
        else
            AdjAdjSpN_List = find(logical(sG(:,SpNIdx).')*logical(sG)*logical(sG.')*logical(sG));
            for AdjAdjSpN_NotOnCorner_Idx = unique(AdjAdjSpN_List(~logical(SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion(AdjAdjSpN_List))))
                STNIdx_Fetch = SpElemProperties.SpN.FirstSTNIdx(AdjAdjSpN_NotOnCorner_Idx)+Time;
                DeltaSTNPos(:,STNIdx_tgt) = DeltaSTNPos(:,STNIdx_tgt) + DeltaSTNPos(:,STNIdx_Fetch);
            end
        end
    end
end

for SpSIdx = find(SpElemProperties.SpS.UpdNumBoundary)
    if ~any(SpElemProperties.SpN.isOn_EdgesOf_HomoUpdNumRegion(logical(sG(SpSIdx,:)))) ...
            || SpElemProperties.SpS.PEC(SpSIdx) == true
        continue;
    end  
    STPIdxList_RefferingSpS = ...
        SpElemProperties.SpS.FirstSTPIdx(SpSIdx):SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    [kappa(STPIdxList_RefferingSpS), FaceArea.Dual(SpSIdx)] = ...
        ComputeKappa_for_SingleSTFISpS(SpSIdx,cdt,sG,sC,sD,DeltaSTNPos,NodePos,SpElemProperties);
end

STFISpSList = find(SpElemProperties.SpS.Belong_to_ST_FI);
for SpSIdx = STFISpSList(~SpElemProperties.SpS.UpdNumBoundary(logical(SpElemProperties.SpS.Belong_to_ST_FI)))
    if SpElemProperties.SpS.PEC(SpSIdx) == true
        STPIdx_RefferingSpS = ...
            SpElemProperties.SpS.FirstSTPIdx(SpSIdx):SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
        kappa(STPIdx_RefferingSpS) = -1;
        FaceArea.Dual(SpSIdx) = 1;
        continue;
    end
    STPIdxList_RefferingSpS = ...
        SpElemProperties.SpS.FirstSTPIdx(SpSIdx):SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    [kappa(STPIdxList_RefferingSpS), FaceArea.Dual(SpSIdx)] = ...
        ComputeKappa_for_SingleSTFISpS(SpSIdx,cdt,sG,sC,sD,DeltaSTNPos,NodePos,SpElemProperties);
end

disp('ComputeKappa_for_STFI_SpPs_and_SpSs: Computing Kappa For STFI faces')
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    if SpElemProperties.SpP.ElecWall(SpPIdx) == true
        STPIdx_RefferingSpP = ...
            SpElemProperties.SpP.FirstSTPIdx(SpPIdx):SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);
        kappa(STPIdx_RefferingSpP) = 1;
        FaceArea.Prim(SpPIdx) = 1;
        continue;
    end
    STPIdxList_RefferingSpP = ...
        SpElemProperties.SpP.FirstSTPIdx(SpPIdx):SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);   
    [kappa(STPIdxList_RefferingSpP),FaceArea.Prim(SpPIdx)] = ...
        ComputeKappa_for_SingleSTFISpP(SpPIdx,cdt,sD,D0,D1,D2,D3,DeltaSTNPos,NodePos,SpElemProperties,STElemProperties);
 end
end
%%
function DeltaSTNPos = Compute_DeltaSTNPos_for_SingleBoundaryEdge(SpSIdx,DeltaSTNPos,NodePos,PrimEdgeTgtLeng,DualFaceTgtArea,sG,sC,D2,D3,SpElemProperties,STElemProperties,cdt)
global EPSILON
kappa       = -1*DualFaceTgtArea/(PrimEdgeTgtLeng*cdt/SpElemProperties.SpS.UpdNum(SpSIdx));
LocalBasis1 = findDeltaDirection(SpSIdx,sG,sC,NodePos.Prim,SpElemProperties);
if norm(LocalBasis1)<EPSILON
    disp(SpSIdx)
end
SpS_Direc   = [NodePos.Prim(logical(sG(SpSIdx,:))).Vec]*sG(SpSIdx,logical(sG(SpSIdx,:))).';
SpS_Direc   = norm(SpS_Direc).^(-1)*SpS_Direc;
LocalBasis2 = SpS_Direc - dot(SpS_Direc,LocalBasis1)*LocalBasis1;
LocalBasis2 = LocalBasis2/norm(LocalBasis2);
LocalBasis3 = cross(LocalBasis1,LocalBasis2);

SpEndPoints = find(sG(SpSIdx,:));
STNIdx_Past = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)-1;
STNIdx_Futr = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)  ;
for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+1:...
        SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)-1
    STNIdx_Past = STNIdx_Past+1;
    STNIdx_Futr = STNIdx_Futr+1;
    
    DualSTNIdx_DefDualFace = find(logical(D2(:,STPIdx)).'*logical(D3).');
    DualSTNIdx_DefDualFace = SortNodesAlongFaceOrientation(DualSTNIdx_DefDualFace,STPIdx,D3.' ,D2.' );
    PosVec_Projected_to_03Plane = zeros(2,size(DualSTNIdx_DefDualFace,2));
    
    Row0 = 1;
    Row3 = 2;
    ColIdx = 0;
    for STNTildeIdx = DualSTNIdx_DefDualFace 
        ColIdx = ColIdx +1;
        PosVec_Projected_to_03Plane(Row0,ColIdx) = cdt *STElemProperties.STV.TimeIdx(STNTildeIdx);
        SpNTildeIdxTemp                          =      STElemProperties.STV.RefSpV(STNTildeIdx);
        PosVec_Projected_to_03Plane(Row3,ColIdx) = dot(NodePos.Dual(SpNTildeIdxTemp).Vec, LocalBasis3);
    end
    Area03_STPtilde       = SignedArea2D(PosVec_Projected_to_03Plane);
    Delta_of_DeltaNodePos = -(Area03_STPtilde/(PrimEdgeTgtLeng*dot(SpS_Direc, LocalBasis2)))/kappa;
    DeltaSTNPos(:,STNIdx_Futr) = DeltaSTNPos(:,STNIdx_Past) ...
        + Delta_of_DeltaNodePos*LocalBasis1;
end
%STNIdx_Last  = STNIdx_Futr;
%STNIdx_First = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints);
%DeltaSTNPos(STNIdx_First) = DeltaSTNPos(STNIdx_Last);
end
%%
function DeltaDirection = findDeltaDirection(SpSTgt,sG,sC,NodePosPrim,SpElemProperties)
global EPSILON
for FaceIdx = find(sC(:,SpSTgt)).'
    GuessedDirection = [0;0;0];
    for EdgeIdx = find(sC(FaceIdx,:))
        if SpElemProperties.SpS.UpdNumBoundary(EdgeIdx) == true
            continue;
        elseif  logical(sG(SpSTgt,:))*(logical(sG(EdgeIdx,:)).') ~= 0
            if norm(GuessedDirection) == 0 
                GuessedDirection = [NodePosPrim(logical(sG(EdgeIdx,:))).Vec]*sG(EdgeIdx,logical(sG(EdgeIdx,:))).';
                GuessedDirection = GuessedDirection/norm(GuessedDirection);
            elseif norm(cross(GuessedDirection, ...
                    [NodePosPrim(logical(sG(EdgeIdx,:))).Vec]*sG(EdgeIdx,logical(sG(EdgeIdx,:))).'))<EPSILON
                DeltaDirection = GuessedDirection;
                
                return;
            end
        end
    end
end
end

%%
function SortedNodesIdx = SortNodesAlongFaceOrientation(NodesDefiningFace_Idxs,FaceTar,GradM,CurlM)
InputNodesNum = size(NodesDefiningFace_Idxs,2);
SortedNodesIdx = zeros(1,InputNodesNum);
LeftNodes = NodesDefiningFace_Idxs;
SortedNodesIdx(1) = NodesDefiningFace_Idxs(1);
LeftNodes(1) = [];
ProgressIdx=1;
for Counter = 2:InputNodesNum
    for EdgeIdx = find(CurlM(FaceTar,:))
        if GradM(EdgeIdx,SortedNodesIdx(ProgressIdx)) ~= 0 && ...
                GradM(EdgeIdx,SortedNodesIdx(ProgressIdx))*CurlM(FaceTar,EdgeIdx)==-1
            GuideEdgeIdx = EdgeIdx;
            break;
        end
    end
    for LeftNodeScan = 1:size(LeftNodes,2)
        GuessedNodeIdx = LeftNodes(LeftNodeScan);
        if GradM(GuideEdgeIdx,GuessedNodeIdx)~=0
            ProgressIdx = ProgressIdx+1;
            SortedNodesIdx(ProgressIdx) = GuessedNodeIdx;
            LeftNodes(LeftNodeScan) = [];
            break;
        end
    end
end
end
%%
function AreaReturn = SignedArea2D(NodePos2D)
Pos0 = NodePos2D(:,1);
Vec1 = zeros('like',Pos0);
Vec2 = NodePos2D(:,2) - NodePos2D(:,1);
NodePos2D(:,1:2) = [];
Num_of_TriangleDivision = size(NodePos2D,2);
AreaReturn = 0;
for TriDivIdx = 1:Num_of_TriangleDivision
    Vec1 = Vec1 + Vec2;
    Vec2 = NodePos2D(:,1) - Pos0;
    NodePos2D(:,1) = [];
    TriDivArea_Signed = 0.5 * Vec1.'*[0 1;-1 0]*Vec2;
    AreaReturn = AreaReturn + TriDivArea_Signed;
end
end
%%
function AreaReturn_Vec = VectorArea3D(NodePos3D)
Pos0 = NodePos3D(:,1);
Vec1 = zeros('like',Pos0);
Vec2 = NodePos3D(:,2) - NodePos3D(:,1);
NodePos3D(:,1:2) = [];
Num_of_TriangleDivision = size(NodePos3D,2);
AreaReturn_Vec = [0;0;0];
for TriDivIdx = 1:Num_of_TriangleDivision
    Vec1 = Vec2;
    Vec2 = NodePos3D(:,1) - Pos0;
    NodePos3D(:,1) = [];
    TriDivArea_Signed = 0.5 * cross(Vec1,Vec2);
    AreaReturn_Vec = AreaReturn_Vec + TriDivArea_Signed;
end
end

%%
function [kappa_partial, FaceArea_Dual] = ComputeKappa_for_SingleSTFISpS(SpS_tgt,cdt,sG,sC,sD,DeltaSTNPos,NodePos,SpElemProperties)

TerminalPrimSpNs    = find(sG(SpS_tgt,:));
UpdNum_SpSTgt       = SpElemProperties.SpS.UpdNum(SpS_tgt);
UpdNum_TermSpNs     = SpElemProperties.SpN.UpdNum(TerminalPrimSpNs).';
UpdNumRatio = UpdNum_TermSpNs/UpdNum_SpSTgt;
TermSTNIdx_SpLikeSTS_Offset  = SpElemProperties.SpN.FirstSTNIdx(TerminalPrimSpNs)-UpdNumRatio;

kappa_partial = zeros(1+UpdNum_SpSTgt,1);

for TimeSec = 1:UpdNum_SpSTgt
    STPTgt_LocalIdx = 1 + TimeSec;
    TermSTNIdx_SpLikeSTS_Past  = TermSTNIdx_SpLikeSTS_Offset + UpdNumRatio*TimeSec;

    SpLikePrimSTSVec_Past = [NodePos.Prim(logical(sG(SpS_tgt,:))).Vec]*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
    SpLikePrimSTSVec_Past = SpLikePrimSTSVec_Past + DeltaSTNPos(:,TermSTNIdx_SpLikeSTS_Past            )*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
    SpLikePrimSTSVec_Futr = [NodePos.Prim(logical(sG(SpS_tgt,:))).Vec]*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
    SpLikePrimSTSVec_Futr = SpLikePrimSTSVec_Futr + DeltaSTNPos(:,TermSTNIdx_SpLikeSTS_Past+UpdNumRatio)*sG(SpS_tgt,logical(sG(SpS_tgt,:))).'; 

    LocalBasis1 = 0.5*(SpLikePrimSTSVec_Past+SpLikePrimSTSVec_Futr);
    LocalBasis1 = norm(LocalBasis1)^(-1)*LocalBasis1;
    
    %     BulgeArea_EP1Side = ...
    %         CalcBulgeArea(STNIdx_forEPs_Past(1),SpEndPoints(1),SpS_tgt,DeltaSTNPos,NodePos.Prim,sG,SpElemProperties,cdt);
    %     BulgeArea_EP2Side = ...
    %         CalcBulgeArea(STNIdx_forEPs_Past(2),SpEndPoints(2),SpS_tgt,DeltaSTNPos,NodePos.Prim,sG,SpElemProperties,cdt);
    BulgeArea_EP1Side = ...
        CalcBulgeArea(TermSTNIdx_SpLikeSTS_Past(1),TerminalPrimSpNs(1),SpS_tgt,DeltaSTNPos,LocalBasis1,sG,SpElemProperties,cdt);
    BulgeArea_EP2Side = ...
        CalcBulgeArea(TermSTNIdx_SpLikeSTS_Past(2),TerminalPrimSpNs(2),SpS_tgt,DeltaSTNPos,LocalBasis1,sG,SpElemProperties,cdt);

    Area_PrimSTP  = ...
        0.5*(dot(SpLikePrimSTSVec_Past,LocalBasis1)+dot(SpLikePrimSTSVec_Futr,LocalBasis1)) *cdt/UpdNum_SpSTgt + (BulgeArea_EP1Side+BulgeArea_EP2Side);
    %     PrimEdgeTgtLeng = SpEdgeLength(SpS_tgt,sG,NodePos.Prim);
    %     Area_PrimSTP  = ...
    %         PrimEdgeTgtLeng *cdt/UpdNum_SpSTgt + (BulgeArea_EP1Side+BulgeArea_EP2Side);
    
    % compute dualface-area as vector and project on perp-plane of the average prim edge direction
    NodeList_DefDualSpP = find(logical(sC(:,SpS_tgt).')*logical(sD.'));
    NodeList_DefDualSpP = SortNodesAlongFaceOrientation(NodeList_DefDualSpP,SpS_tgt,sD.',sC.');
    NodePos_DefDualSpP_SpPart = zeros(3,size(NodeList_DefDualSpP,2));
    ColIdx = 0;
    for SpNIdx = NodeList_DefDualSpP
        ColIdx = ColIdx +1;
        NodePos_DefDualSpP_SpPart(:,ColIdx) = NodePos.Dual(SpNIdx).Vec;
        % project dual face to plane perp to LocalBasis1
        NodePos_DefDualSpP_SpPart(:,ColIdx) = NodePos_DefDualSpP_SpPart(:,ColIdx)-dot(NodePos_DefDualSpP_SpPart(:,ColIdx),LocalBasis1)*LocalBasis1;
    end
    FaceArea_Dual  = norm(VectorArea3D(NodePos_DefDualSpP_SpPart));
    Area_DualSTP  = FaceArea_Dual;
    
    kappa_partial(STPTgt_LocalIdx) = -1*Area_DualSTP/Area_PrimSTP;  
    if kappa_partial(STPTgt_LocalIdx)>0
        disp('kappa>0 for Edgelike-STP')
    end
end
kappa_partial(1) = kappa_partial(1+UpdNum_SpSTgt);
end

%%
function BulgeArea = CalcBulgeArea(STNTgt_Oldest,SpN_tgt,SpS_tgt,DeltaSTNPos,LocalBasis1,sG,SpElemProperties,cdt)

% LocalBasis1 = ...
%     [NodePosPrim(logical(sG(SpSTgt,:))).Vec]*sG(SpSTgt,logical(sG(SpSTgt,:))).';
% LocalBasis1 = LocalBasis1/norm(LocalBasis1);

BulgeArea = 0;
for LocalTimeSec = 1:SpElemProperties.SpN.UpdNum(SpN_tgt)/SpElemProperties.SpS.UpdNum(SpS_tgt)
    EdgeLengPast = ...
        sG(SpS_tgt,SpN_tgt) * dot(DeltaSTNPos(:,STNTgt_Oldest+LocalTimeSec-1), LocalBasis1);
    EdgeLengFutr = ...
        sG(SpS_tgt,SpN_tgt) * dot(DeltaSTNPos(:,STNTgt_Oldest+LocalTimeSec  ), LocalBasis1);
    BulgeArea = BulgeArea ...
        + (cdt/SpElemProperties.SpN.UpdNum(SpN_tgt)) * (EdgeLengPast + EdgeLengFutr)/2;
end
end

%%

% function [kappa_partial, FaceArea_Dual] = ComputeKappa_for_SingleBoundarySpS_OnCorners(SpS_tgt,cdt,sG,sC,sD,DeltaSTNPos,NodePos,SpElemProperties)
% 
% TerminalPrimSpNs    = find(sG(SpS_tgt,:));
% UpdNum_SpSTgt       = SpElemProperties.SpS.UpdNum(SpS_tgt);
% UpdNum_TermSpNs     = SpElemProperties.SpN.UpdNum(TerminalPrimSpNs).';
% UpdNumRatio = UpdNum_TermSpNs/UpdNum_SpSTgt;
% TermSTNIdx_SpLikeSTS_Offset  = SpElemProperties.SpN.FirstSTNIdx(TerminalPrimSpNs)-UpdNumRatio;
% 
% kappa_partial = zeros(1+UpdNum_SpSTgt,1);
% 
% for TimeSec = 1:UpdNum_SpSTgt
%     STPTgt_LocalIdx = 1 + TimeSec;
%     TermSTNIdx_SpLikeSTS_Past  = TermSTNIdx_SpLikeSTS_Offset + UpdNumRatio*TimeSec;
%     % compute sp-part of the average direction of the two SpLike prim ST-edge
%     SpLikePrimSTSVec_Past = [NodePos.Prim(logical(sG(SpS_tgt,:))).Vec]*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
%     SpLikePrimSTSVec_Past = SpLikePrimSTSVec_Past + DeltaSTNPos(:,TermSTNIdx_SpLikeSTS_Past            )*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
%     SpLikePrimSTSVec_Futr = [NodePos.Prim(logical(sG(SpS_tgt,:))).Vec]*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
%     SpLikePrimSTSVec_Futr = SpLikePrimSTSVec_Futr + DeltaSTNPos(:,TermSTNIdx_SpLikeSTS_Past+UpdNumRatio)*sG(SpS_tgt,logical(sG(SpS_tgt,:))).';
%     
%     LocalBasis1 = 0.5*(SpLikePrimSTSVec_Past+SpLikePrimSTSVec_Futr);
%     LocalBasis1 = norm(LocalBasis1)^(-1)*LocalBasis1;
%     
%     % project sp-part of the prim ST-edges on the average direction
%     % vector
%     Area_PrimSTP  = ...
%         0.5*(dot(SpLikePrimSTSVec_Past,LocalBasis1)+dot(SpLikePrimSTSVec_Futr,LocalBasis1))*cdt/UpdNum_SpSTgt;
%     NodeList_DefDualSpP = find(logical(sC(:,SpS_tgt).')*logical(sD.'));
%     NodeList_DefDualSpP = SortNodesAlongFaceOrientation(NodeList_DefDualSpP,SpS_tgt,sD.',sC.');
%     NodePos_DefDualSpP_SpPart = zeros(3,size(NodeList_DefDualSpP,2));
%     ColIdx = 0;
%     for SpNIdx = NodeList_DefDualSpP
%         ColIdx = ColIdx +1;
%         NodePos_DefDualSpP_SpPart(:,ColIdx) = NodePos.Dual(SpNIdx).Vec;
%         % project sp-part of the dual ST-face on the perp-plane of the
%         % average direction vector of the prim SpS
%         NodePos_DefDualSpP_SpPart(:,ColIdx) = NodePos_DefDualSpP_SpPart(:,ColIdx)-dot(NodePos_DefDualSpP_SpPart(:,ColIdx),LocalBasis1)*LocalBasis1;
%     end
%     FaceArea_Dual  = norm(VectorArea3D(NodePos_DefDualSpP_SpPart));
%     Area_DualSTP  = FaceArea_Dual;
%     
%     kappa_partial(STPTgt_LocalIdx) = -1*Area_DualSTP/Area_PrimSTP;
%     if kappa_partial(STPTgt_LocalIdx)>0
%         disp('kappa>0 for Edge-associated-STP')
%     end
% end
% kappa_partial(1) = kappa_partial(1+UpdNum_SpSTgt);
% end

%%
function [kappa_Local,FaceArea_Prim] = ComputeKappa_for_SingleSTFISpP(SpPIdx,cdt,sD,D0,D1,D2,D3,DeltaSTNPos,NodePos,SpElemProperties,STElemProperties)
kappa_Local = zeros(1+SpElemProperties.SpP.UpdNum(SpPIdx),1);
for Time = 0:SpElemProperties.SpP.UpdNum(SpPIdx)-1
    STPIdx_Local  = 1+Time;
    STPIdx_Global = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+Time;
    LocalBasis3 = [NodePos.Dual(logical(sD(:,SpPIdx).')).Vec]*sD(logical(sD(:,SpPIdx)),SpPIdx);
    LocalBasis3 = LocalBasis3/norm(LocalBasis3);
    STNIdx_DefPrimFace = find(logical(D1(STPIdx_Global,:))*logical(D0));
    STNIdx_DefPrimFace = SortNodesAlongFaceOrientation(STNIdx_DefPrimFace,STPIdx_Global,D0   ,D1   );
    PosVec_SpatialPart = zeros(3,size(STNIdx_DefPrimFace,2));
    ColIdx = 0;
    for STNIdx = STNIdx_DefPrimFace
        ColIdx = ColIdx +1;
        SpNIdx = STElemProperties.STN.RefSpN(STNIdx);
        PosVec_SpatialPart(:,ColIdx) = NodePos.Prim(SpNIdx).Vec + DeltaSTNPos(:,STNIdx);
        % project prim face to perp-plane of the dual Sp-edge (perp to LocalBasis3)
        PosVec_SpatialPart(:,ColIdx) = PosVec_SpatialPart(:,ColIdx)-dot(PosVec_SpatialPart(:,ColIdx),LocalBasis3)*LocalBasis3;
    end
    Area_PrimSTP  = norm(VectorArea3D(PosVec_SpatialPart));
    if Time == 0
        FaceArea_Prim = Area_PrimSTP;
    end
    STNTilIdx_DefDualFace = find(logical(D2(:,STPIdx_Global)).'*logical(D3).');
    STNTilIdx_DefDualFace = SortNodesAlongFaceOrientation(STNTilIdx_DefDualFace,STPIdx_Global,D3.' ,D2.');
    PosVec_Projected_to_03Plane  = zeros(2,size(STNTilIdx_DefDualFace,2));
    Row0 = 1;
    Row3 = 2;
    ColIdx = 0;
    for STNTildeIdx = STNTilIdx_DefDualFace
        ColIdx = ColIdx +1;
        PosVec_Projected_to_03Plane(Row0,ColIdx) = cdt *STElemProperties.STV.TimeIdx(STNTildeIdx);
        SpNTildeIdx = STElemProperties.STV.RefSpV(STNTildeIdx);
        PosVec_Projected_to_03Plane(Row3,ColIdx) = dot(NodePos.Dual(SpNTildeIdx).Vec, LocalBasis3);
    end
    Area_DualSTP = abs(SignedArea2D(PosVec_Projected_to_03Plane));
    kappa_Local(STPIdx_Local) = Area_DualSTP/Area_PrimSTP;
end
kappa_Local(1+SpElemProperties.SpP.UpdNum(SpPIdx)) = kappa_Local(1);
end
%%
function [kappa,FaceAreaPrim] = ComputeKappa_for_NonSTFI_SpPs(kappa,FaceAreaPrim,cdt,sG,sC,sD,NodePos,SpElemProperties)
disp('ComputeKappa_for_NonSTFI_SpPs: Computing kappa for SpFI faces')
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI==false)
    switch SpElemProperties.SpP.ElecWall(SpPIdx)
        case true
            for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                    SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
                kappa(STPIdx) = 1;
            end
            FaceAreaPrim(SpPIdx) = 1;
        case false
            if SpElemProperties.SpP.PrimAreaIsGiven(SpPIdx) == true && SpElemProperties.SpP.DualLengIsGiven(SpPIdx) == true
                Area_PrimSTP = SpElemProperties.SpP.PrimArea(SpPIdx);
                Area_DualSTP = SpElemProperties.SpP.DualLeng(SpPIdx)*cdt/SpElemProperties.SpP.UpdNum(SpPIdx);
            else
                Area_PrimSTP = SpFaceArea(SpPIdx,sG,sC,NodePos.Prim);
                Area_DualSTP = SpEdgeLength(SpPIdx,sD.',NodePos.Dual)*cdt/SpElemProperties.SpP.UpdNum(SpPIdx);
            end
            for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                    SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
                kappa(STPIdx) = Area_DualSTP/Area_PrimSTP;
                if kappa(STPIdx) == 0
                    disp(['STP=',num2str(STPIdx),'SpPIdx=',num2str(SpPIdx)])
                end
            end
            FaceAreaPrim(SpPIdx) = Area_PrimSTP;
    end
end
end
function [kappa,FaceAreaDual] = ComputeKappa_for_NonSTFI_SpSs(kappa,FaceAreaDual,cdt,sG,sC,sD,NodePos,SpElemProperties)
disp('ComputeKappa_for_NonSTFI_SpSs: Computing kappa for SpFI edges')
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI==false)
    switch SpElemProperties.SpS.PEC(SpSIdx)
        case true
            for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                    SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
                kappa(STPIdx) = -1;
            end
            FaceAreaDual(SpSIdx) = 1;
        case false
            if SpElemProperties.SpS.PrimLengIsGiven(SpSIdx) == true && SpElemProperties.SpS.DualAreaIsGiven(SpSIdx) == true  
                Area_PrimSTP = SpElemProperties.SpS.PrimLeng(SpSIdx)*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
                Area_DualSTP = SpElemProperties.SpS.DualArea(SpSIdx);               
            else
                Area_PrimSTP = SpEdgeLength(SpSIdx,sG,NodePos.Prim)*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
                Area_DualSTP = SpFaceArea(SpSIdx,sD.',sC.',NodePos.Dual);
            end
            for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                    SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
                kappa(STPIdx) = -1*Area_DualSTP/Area_PrimSTP;
                if kappa(STPIdx) == 0
                    disp(['STP=',num2str(STPIdx),'SpSIdx=',num2str(SpSIdx)])
                end
            end
            FaceAreaDual(SpSIdx) = Area_DualSTP;
    end
end
end

%%
function Area_Face = SpFaceArea(SpPIdx,sG,sC,NodePos)
Nodes = find(logical(sC(SpPIdx,:))*logical(sG));
Area_Face = ThreeD_SpP_Area(Nodes,sG,NodePos);
end
function Length_Edge = SpEdgeLength(SpSIdx,sG,NodePos)
Endpoints = find(sG(SpSIdx,:));
Length_Edge = norm(NodePos(Endpoints(2)).Vec - NodePos(Endpoints(1)).Vec);
end
function Area_3DFace = ThreeD_SpP_Area(Nodes,sG,NodePos)
StartNode = Nodes(1);
Nodes(1) = [];
LastNode = StartNode;
IsAdj_to_StartNode = logical(sG).'*logical(sG(:,LastNode));
for LeftNodeScan = 1:size(Nodes,2)
    NodeIdx = Nodes(LeftNodeScan);
    if IsAdj_to_StartNode(NodeIdx)==true
        NextNode = NodeIdx;
        Nodes(LeftNodeScan) = [];
        break;
    else
        NextNode = 0;
    end
end
clearvars IsAdj_to_StartNode
Vec2 = NodePos(NextNode).Vec - NodePos(LastNode).Vec;
Vec1 = zeros('like',Vec2);
Num_of_TriangleDivision = size(Nodes,2);
FaceArea_Signed = 0;
for TriDivIdx = 1:Num_of_TriangleDivision
    IsAdj_to_NextNode = logical(sG).'*logical(sG(:,NextNode));
    for LeftNodeScan = 1:size(Nodes,2)
        NodeIdx = Nodes(LeftNodeScan);
        if IsAdj_to_NextNode(NodeIdx)==true && NodeIdx ~= LastNode
            LastNode = NextNode;
            NextNode = NodeIdx;
            Vec1 = Vec1+Vec2;
            Vec2 = NodePos(NextNode).Vec - NodePos(LastNode).Vec;
            Nodes(LeftNodeScan) = [];
            break;
        end
    end
    TriDivArea = 0.5*norm(cross(Vec1,Vec2));
    FaceArea_Signed = FaceArea_Signed + TriDivArea;
end
Area_3DFace = abs(FaceArea_Signed);
end
