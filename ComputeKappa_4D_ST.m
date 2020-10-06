function [kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties,Num_of_Elem)
global SpDIM
kappa = zeros(Num_of_Elem.STP,1);

FaceArea.Prim = zeros(Num_of_Elem.SpP,1);
FaceArea.Dual = zeros(Num_of_Elem.SpS,1);
DeltaSTNPos   = sparse(SpDIM, Num_of_Elem.STN);
[kappa,FaceArea] = ComputeKappa_for_STFI_SpPs_and_SpSs(kappa,FaceArea,DeltaSTNPos,cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties);
disp('call ComputeKappa_for_SpFI_SpSs')
[kappa,FaceArea.Dual]   = ComputeKappa_for_SpFI_SpSs(kappa,FaceArea.Prim,cdt,sG,sC,sD,NodePos,SpElemProperties);
disp('call ComputeKappa_for_SpFI_SpPs')
[kappa,FaceArea.Prim]   = ComputeKappa_for_SpFI_SpPs(kappa,FaceArea.Dual,cdt,sG,sC,sD,NodePos,SpElemProperties);
end
%%
function [kappa,FaceArea] = ComputeKappa_for_STFI_SpPs_and_SpSs(kappa,FaceArea,DeltaSTNPos,cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,STElemProperties)

for SpSIdx = find(SpElemProperties.SpS.UpdNumBoundary)
    if SpElemProperties.SpS.PEC(SpSIdx) == true
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            kappa(STPIdx) = -1;
        end
        FaceArea.Dual(SpSIdx) = 1;
        continue;
    elseif SpElemProperties.SpS.UpdNumCorner == true
        continue;
    else
        PrimEdgeTgtLeng = SpEdgeLength(SpSIdx,sG,NodePos.Prim);
        DualFaceTgtArea = SpFaceArea(SpSIdx,sD.',sC.',NodePos.Dual);
        Area_PrimSTP = PrimEdgeTgtLeng*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
        Area_DualSTP = DualFaceTgtArea;
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            kappa(STPIdx) = -1*Area_DualSTP/Area_PrimSTP;
        end
        FaceArea.Dual(SpSIdx) = Area_DualSTP;
        DeltaSTNPos = ...
            Compute_DeltaSTNPos_for_SingleBoundaryEdge(SpSIdx,DeltaSTNPos,NodePos,PrimEdgeTgtLeng,DualFaceTgtArea,sG,sC,D2,D3,SpElemProperties,STElemProperties,cdt);
    end 
end

for SpPIdx = find(SpElemProperties.SpP.UpdNumBoundary)
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
            disp('norm(DeltaSTN_Fetch) == 0')
            pause;
        end
        for SpNIdx = SpNIdx_DefPrimFace
            STNIdx = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+Time;
            DeltaSTNPos(:,STNIdx) = DeltaSTN_Fetch;
        end
    end
end

for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI==true)
    if SpElemProperties.SpS.UpdNumBoundary(SpSIdx)==true
        continue;
    elseif SpElemProperties.SpS.PEC(SpSIdx) == true
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            kappa(STPIdx) = -1;
        end
        FaceArea.Dual(SpSIdx) = 1;
        continue;
    end
    [kappa, FaceArea.Dual(SpSIdx)] = ...
        ComputeKappa_for_SingleNonBoundarySpS(kappa,SpSIdx,cdt,sG,sC,sD,DeltaSTNPos,NodePos,SpElemProperties);
end

for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    if SpElemProperties.SpP.ElecWall(SpPIdx) == true
        for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
            kappa(STPIdx) = 1;
        end
        FaceArea.Prim(SpPIdx) = 1;
        continue;
    end
    
    for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
            SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)-1
        STNIdx_DefPrimFace = find(logical(D1(STPIdx,:))*logical(D0));
        STNIdx_DefPrimFace = SortNodesAlongFaceOrientation(STNIdx_DefPrimFace,STPIdx,D0   ,D1   );
        PosVec_SpatialPart = zeros(3,size(STNIdx_DefPrimFace,2));
        ColIdx = 0;
        for STNIdx = STNIdx_DefPrimFace
            ColIdx = ColIdx +1;
            SpNIdx = STElemProperties.STN.RefSpN(STNIdx);
            PosVec_SpatialPart(:,ColIdx) = NodePos.Prim(SpNIdx).Vec + DeltaSTNPos(:,STNIdx);
        end
%         if STPIdx == 16
%             disp('hoge')
%         end
        Area_PrimSTP  = norm(VectorArea3D(PosVec_SpatialPart));
        if STPIdx == SpElemProperties.SpP.FirstSTPIdx(SpPIdx)
            FaceArea.Prim(SpPIdx) = Area_PrimSTP;
        end
%        FaceAreaPrimSTP(STPIdx) = Area_PrimSTP;
        LocalBasis3 = [NodePos.Dual(logical(sD(:,SpPIdx).')).Vec]*sD(logical(sD(:,SpPIdx)),SpPIdx);
        LocalBasis3 = LocalBasis3/norm(LocalBasis3);
        STNTilIdx_DefDualFace = find(logical(D2(:,STPIdx)).'*logical(D3).');
        %         SignFlipper_STP(1:Num_of_Elem.SplikeSTP,1) ...
        %             =    ones(Num_of_Elem.SplikeSTP,1);
        %         SignFlipper_STP(Num_of_Elem.SplikeSTP+1:Num_of_Elem.STP,1) ...
        %             = -1*ones(Num_of_Elem.TilikeSTP,1);
        %         SignFlipper_STOmega(1:Num_of_Elem.SplikeSTOmega,1) ...
        %             =    ones(Num_of_Elem.SplikeSTOmega,1);
        %         SignFlipper_STOmega(Num_of_Elem.SplikeSTOmega+1:Num_of_Elem.STOmega,1) ...
        %             = -1*ones(Num_of_Elem.TilikeSTOmega,1);
        %         STNTildeIdx_DefiningDualFace = SortSTNsAlongFaceOrientation(STNTildeIdx_DefiningDualFace,STPIdx,...
        %             (D3...
        %             *spdiags(SignFlipper_STOmega,0,Num_of_Elem.STOmega,Num_of_Elem.STOmega)...
        %             ).' ,...
        %             (spdiags(SignFlipper_STOmega,0,Num_of_Elem.STOmega,Num_of_Elem.STOmega)...
        %             *D2...
        %             *spdiags(SignFlipper_STP,0,Num_of_Elem.STP,Num_of_Elem.STP)...
        %             ).' ...
        %             );
        STNTilIdx_DefDualFace = SortNodesAlongFaceOrientation(STNTilIdx_DefDualFace,STPIdx,D3.' ,D2.');
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
        
        kappa(STPIdx) = Area_DualSTP/Area_PrimSTP;
    end
    STP_First = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
    STP_Last  = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);
    kappa(STP_Last) = kappa(STP_First);
end
end
%%
function DeltaSTNPos = Compute_DeltaSTNPos_for_SingleBoundaryEdge(SpSIdx,DeltaSTNPos,NodePos,PrimEdgeTgtLeng,DualFaceTgtArea,sG,sC,D2,D3,SpElemProperties,STElemProperties,cdt)
kappa       = -1*DualFaceTgtArea/(PrimEdgeTgtLeng*cdt/SpElemProperties.SpS.UpdNum(SpSIdx));
LocalBasis1 = findDeltaDirection(SpSIdx,sG,sC,NodePos.Prim,SpElemProperties);
% if LocalBasis1(1)==0
%     disp(SpSIdx)
% end
SpS_Direc   = [NodePos.Prim(logical(sG(SpSIdx,:))).Vec]*sG(SpSIdx,logical(sG(SpSIdx,:))).';
LocalBasis2 = SpS_Direc - dot(SpS_Direc,LocalBasis1)*LocalBasis1;
LocalBasis2 = LocalBasis2/norm(LocalBasis2);
LocalBasis3 = cross(LocalBasis1,LocalBasis2);

SpEndPoints = find(sG(SpSIdx,:));
STNIdx_Past = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)-1;
STNIdx_Futr = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)  ;
for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+1:...
        SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
    STNIdx_Past = STNIdx_Past+1;
    STNIdx_Futr = STNIdx_Futr+1;
    
    STNTilIdx_DefDualFace = find(logical(D2(:,STPIdx)).'*logical(D3).');
    STNTilIdx_DefDualFace = SortNodesAlongFaceOrientation(STNTilIdx_DefDualFace,STPIdx,D3.' ,D2.' );
    PosVec_Projected_to_03Plane = zeros(2,size(STNTilIdx_DefDualFace,2));
    
    Row0 = 1;
    Row3 = 2;
    ColIdx = 0;
    for STNTildeIdx = STNTilIdx_DefDualFace 
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
STNIdx_Last  = STNIdx_Futr;
STNIdx_First = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints);
DeltaSTNPos(STNIdx_First) = DeltaSTNPos(STNIdx_Last);
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
function SortedNodesIdx = SortNodesAlongFaceOrientation(NodesDefiningFace_Idxs,STPTar,GradM,CurlM)
InputNodesNum = size(NodesDefiningFace_Idxs,2);
SortedNodesIdx = zeros(1,InputNodesNum);
LeftNodes = NodesDefiningFace_Idxs;
SortedNodesIdx(1) = NodesDefiningFace_Idxs(1);
LeftNodes(1) = [];
ProgressIdx=1;
%STPTar
for Counter = 2:InputNodesNum
    %Counter
    %ProgressIdx 
    for EdgeIdx = find(CurlM(STPTar,:))
        if GradM(EdgeIdx,SortedNodesIdx(ProgressIdx)) ~= 0 && ...
                GradM(EdgeIdx,SortedNodesIdx(ProgressIdx))*CurlM(STPTar,EdgeIdx)==-1
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
function [kappa, FaceArea_Dual] = ComputeKappa_for_SingleNonBoundarySpS(kappa,SpSTgt,cdt,sG,sC,sD,DeltaSTNPos,NodePos,SpElemProperties)
PrimEdgeTgtLeng         = SpEdgeLength(SpSTgt,sG,NodePos.Prim);
DualFaceTgtArea         = SpFaceArea(SpSTgt,sD.',sC.',NodePos.Dual);
FaceArea_Dual   = DualFaceTgtArea;

SpEndPoints         = find(sG(SpSTgt,:));
UpdNum_SpSTgt       = SpElemProperties.SpS.UpdNum(SpSTgt);
UpdNumEP            = SpElemProperties.SpN.UpdNum(SpEndPoints);

STNIdx_forEPs_Past  = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)-UpdNumEP/UpdNum_SpSTgt;
STPTgt              = SpElemProperties.SpS.FirstSTPIdx(SpSTgt);
for TimeSecTgt=1:UpdNum_SpSTgt
    STPTgt = STPTgt + 1;
    STNIdx_forEPs_Past  = STNIdx_forEPs_Past + UpdNumEP/UpdNum_SpSTgt;
    
    BulgeArea_EP1Side = ...
        CalcBulgeArea(STNIdx_forEPs_Past(1),SpEndPoints(1),SpSTgt,DeltaSTNPos,NodePos.Prim,sG,SpElemProperties,cdt);
    BulgeArea_EP2Side = ...
        CalcBulgeArea(STNIdx_forEPs_Past(2),SpEndPoints(2),SpSTgt,DeltaSTNPos,NodePos.Prim,sG,SpElemProperties,cdt);
   
    Area_PrimSTP  = ...
        PrimEdgeTgtLeng *cdt/UpdNum_SpSTgt + (BulgeArea_EP1Side+BulgeArea_EP2Side);
    Area_DualSTP  = DualFaceTgtArea;
    
    kappa(STPTgt) = -1*Area_DualSTP/Area_PrimSTP;  
    if kappa(STPTgt)>0
        disp('kappa(STPTgt)>0')
    end
end
kappa(STPTgt-UpdNum_SpSTgt) = kappa(STPTgt);
end

%%
function BulgeArea = ...
    CalcBulgeArea(STNTgt_Oldest,SpNTgt,SpSTgt,DeltaSTNPos,NodePosPrim,sG,SpElemProperties,cdt)
Basis1 = ...
    [NodePosPrim(logical(sG(SpSTgt,:))).Vec]*sG(SpSTgt,logical(sG(SpSTgt,:))).';
Basis1 = Basis1/norm(Basis1);

BulgeArea = 0;
for LocalTimeSec = 1:SpElemProperties.SpN.UpdNum(SpNTgt)/SpElemProperties.SpS.UpdNum(SpSTgt)
    EdgeLengPast = ...
        sG(SpSTgt,SpNTgt) * dot(DeltaSTNPos(:,STNTgt_Oldest+LocalTimeSec-1), Basis1);
    EdgeLengFutr = ...
        sG(SpSTgt,SpNTgt) * dot(DeltaSTNPos(:,STNTgt_Oldest+LocalTimeSec  ), Basis1);
    BulgeArea = BulgeArea ...
        + (cdt/SpElemProperties.SpN.UpdNum(SpNTgt)) * (EdgeLengPast + EdgeLengFutr)/2;
end
end

%%

%%
function [kappa,FaceAreaPrim] = ComputeKappa_for_SpFI_SpPs(kappa,FaceAreaPrim,cdt,sG,sC,sD,NodePos,SpElemProperties)
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
function [kappa,FaceAreaDual] = ComputeKappa_for_SpFI_SpSs(kappa,FaceAreaDual,cdt,sG,sC,sD,NodePos,SpElemProperties)
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