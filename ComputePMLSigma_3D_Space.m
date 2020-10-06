function [Sigma] = ComputePMLSigma_3D_Space(cdt,NodePos,sG,sD,SpElemProperties,Num_of_Elem)

Sigma.SpP.OneTwoThree = sparse(Num_of_Elem.SpP,1);
Sigma.SpP.TwoThreeOne = sparse(Num_of_Elem.SpP,1);
Sigma.SpP.ThreeOneTwo = sparse(Num_of_Elem.SpP,1);
Sigma.SpS.OneTwoThree = sparse(Num_of_Elem.SpS,1);
Sigma.SpS.TwoThreeOne = sparse(Num_of_Elem.SpS,1);
Sigma.SpS.ThreeOneTwo = sparse(Num_of_Elem.SpS,1);

for SpPIdx = find(SpElemProperties.SpP.PML)
    if SpElemProperties.SpP.ElecWall(SpPIdx)
        continue
    end
    DualNode_Sta = find(sD(:,SpPIdx).'==-1,1);
    DualNode_Tgt = find(sD(:,SpPIdx).'== 1,1);
    ApproxPoint = 0.5* (NodePos.Dual(DualNode_Tgt).Vec + NodePos.Dual(DualNode_Sta).Vec);
    LineVec     =       NodePos.Dual(DualNode_Tgt).Vec - NodePos.Dual(DualNode_Sta).Vec;
    Sigma.SpP.OneTwoThree(SpPIdx) = ...
        (LineVec(1)*func_sigma_123(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(1)...
        +LineVec(2)*func_sigma_123(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(2)...
        +LineVec(3)*func_sigma_123(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(3)...
        )/norm(LineVec);
    Sigma.SpP.TwoThreeOne(SpPIdx) = ...
        (LineVec(1)*func_sigma_231(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(1)...
        +LineVec(2)*func_sigma_231(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(2)...
        +LineVec(3)*func_sigma_231(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(3)...
        )/norm(LineVec);
    Sigma.SpP.ThreeOneTwo(SpPIdx) = ...
        (LineVec(1)*func_sigma_312(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(1)...
        +LineVec(2)*func_sigma_312(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(2)...
        +LineVec(3)*func_sigma_312(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(3)...
        )/norm(LineVec);
end
for SpSIdx = find(SpElemProperties.SpS.PML)
    if SpElemProperties.SpS.PEC(SpSIdx)
        continue
    end
    PrimNode_Sta = find(sG(SpSIdx,:)==-1,1);
    PrimNode_Tgt = find(sG(SpSIdx,:)== 1,1);
    ApproxPoint = 0.5* (NodePos.Prim(PrimNode_Tgt).Vec + NodePos.Prim(PrimNode_Sta).Vec);
    LineVec     =       NodePos.Prim(PrimNode_Tgt).Vec - NodePos.Prim(PrimNode_Sta).Vec;
    Sigma.SpS.OneTwoThree(SpSIdx) = ...
        (LineVec(1)*func_sigma_123(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(1)...
        +LineVec(2)*func_sigma_123(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(2)...
        +LineVec(3)*func_sigma_123(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(3)...
        )/norm(LineVec);
    Sigma.SpS.TwoThreeOne(SpSIdx) = ...
        (LineVec(1)*func_sigma_231(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(1)...
        +LineVec(2)*func_sigma_231(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(2)...
        +LineVec(3)*func_sigma_231(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(3)...
        )/norm(LineVec);
    Sigma.SpS.ThreeOneTwo(SpSIdx) = ...
        (LineVec(1)*func_sigma_312(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(1)...
        +LineVec(2)*func_sigma_312(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(2)...
        +LineVec(3)*func_sigma_312(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineVec(3)...
        )/norm(LineVec);
end

Sigma.SpP.OneTwoThree = cdt * Sigma.SpP.OneTwoThree;
Sigma.SpP.TwoThreeOne = cdt * Sigma.SpP.TwoThreeOne;
Sigma.SpP.ThreeOneTwo = cdt * Sigma.SpP.ThreeOneTwo;
Sigma.SpS.OneTwoThree = cdt * Sigma.SpS.OneTwoThree;
Sigma.SpS.TwoThreeOne = cdt * Sigma.SpS.TwoThreeOne;
Sigma.SpS.ThreeOneTwo = cdt * Sigma.SpS.ThreeOneTwo;

% for SpPIdx = find(SpElemProperties.SpP.PML)
%     DualNode_Sta = find(sD(:,SpPIdx).'==-1);
%     DualNode_Tgt = find(sD(:,SpPIdx).'== 1);
%     Sigma.SpP.OneTwoThree(SpPIdx) = VectorStraightLineIntegral(@func_sigma_123,NodePos.Dual(DualNode_Sta).Vec,NodePos.Dual(DualNode_Tgt).Vec);
%     Sigma.SpP.OneTwoThree(SpPIdx) = Sigma.SpP.OneTwoThree(SpPIdx)/norm(NodePos.Dual(logical(sD(:,SpPIdx).')).Vec*sD(SpPIdx,logical(sD(:,SpPIdx))));
%     Sigma.SpP.TwoThreeOne(SpPIdx) = VectorStraightLineIntegral(@func_sigma_231,NodePos.Dual(DualNode_Sta).Vec,NodePos.Dual(DualNode_Tgt).Vec);
%     Sigma.SpP.TwoThreeOne(SpPIdx) = Sigma.SpP.TwoThreeOne(SpPIdx)/norm(NodePos.Dual(logical(sD(:,SpPIdx).')).Vec*sD(SpPIdx,logical(sD(:,SpPIdx))));
%     Sigma.SpP.ThreeOneTwo(SpPIdx) = VectorStraightLineIntegral(@func_sigma_312,NodePos.Dual(DualNode_Sta).Vec,NodePos.Dual(DualNode_Tgt).Vec);
%     Sigma.SpP.ThreeOneTwo(SpPIdx) = Sigma.SpP.ThreeOneTwo(SpPIdx)/norm(NodePos.Dual(logical(sD(:,SpPIdx).')).Vec*sD(SpPIdx,logical(sD(:,SpPIdx))));
% end
% 
% for SpSIdx = find(SpElemProperties.SpS.PML)
%     NodesDefiningDualFace_Idxs = find(logical(sC(:,SpSIdx).')*logical(sD.'));
%     NodesDefiningDualFace_Idxs = SortNodesAlongFaceOrientation(NodesDefiningDualFace_Idxs,SpSPIdx,sD.',sC.');
%     BoundaryDualEdges = logical(sC(:,SpSIdx).');
%     LocalGradM=sD(NodesDefiningDualFace_Idxs,BoundaryDualEdges).';
%     LocalNodePosM= [NodePos.Dual(NodesDefiningDualFace_Idxs).Vec];
%     
%     Sigma.SpS.OneTwoThree(SpSIdx) = VectorRectangleSurfaceIntegral(@func_sigma_123,LocalGradM,LocalNodePosM);
%     Sigma.SpS.OneTwoThree(SpSIdx) = Sigma.SpS.OneTwoThree(SpPIdx)/FaceArea.Dual(SpSIdx);
%     Sigma.SpS.TwoThreeOne(SpSIdx) = VectorRectangleSurfaceIntegral(@func_sigma_231,LocalGradM,LocalNodePosM);
%     Sigma.SpS.TwoThreeOne(SpSIdx) = Sigma.SpS.TwoThreeOne(SpPIdx)/FaceArea.Dual(SpSIdx);   
%     Sigma.SpS.ThreeOneTwo(SpSIdx) = VectorRectangleSurfaceIntegral(@func_sigma_312,LocalGradM,LocalNodePosM);
%     Sigma.SpS.ThreeOneTwo(SpSIdx) = Sigma.SpS.ThreeOneTwo(SpPIdx)/FaceArea.Dual(SpSIdx);    
% end

end
function IntegratedValue = VectorStraightLineIntegral(OneFormFuncHandle,StaNodePos,TgtNodePos)
LineVec = TgtNodePos-StaNodePos;

FuncAlongLine = @(x) ...
    (abs(LineVec(1))*OneFormFuncHandle(1,LineVec(1)*x+StaNodePos(1),LineVec(2)*x+StaNodePos(2),LineVec(3)*x+StaNodePos(3))...
    +abs(LineVec(2))*OneFormFuncHandle(2,LineVec(1)*x+StaNodePos(1),LineVec(2)*x+StaNodePos(2),LineVec(3)*x+StaNodePos(3))...
    +abs(LineVec(3))*OneFormFuncHandle(3,LineVec(1)*x+StaNodePos(1),LineVec(2)*x+StaNodePos(2),LineVec(3)*x+StaNodePos(3))...
    )/norm(LineVec);

IntegratedValue = integral(FuncAlongLine,0,1);

end

function IntegratedValue = VectorRectangleSurfaceIntegral(TwoFormFuncHandle,LocalGradM,LocalNodePosM)
global EPSILON

if size(LocalNodePosM,2)~=4
    disp('VectorRectangleSurfaceIntegral:Non-quadrangle integration domains are currently not available.')
    return;
end
% N.B. Nodes must be sorted along the orientation of the integration area.
LocalOrigin = LocalNodePosM(:,1);
for EdgeIdx = find(LocalGradM(:,1).')
    if LocalGradM(EdgeIdx,2)~=0
        LocalBasisMatrix(:,1) = -LocalGradM(EdgeIdx,1)*(LocalGradM(EdgeIdx,:)*LocalNodePosM.').';
    elseif LocalGradM(EdgeIdx,4)~=0
        LocalBasisMatrix(:,2) = -LocalGradM(EdgeIdx,1)*(LocalGradM(EdgeIdx,:)*LocalNodePosM.').';
    end
end
for LocalAxis = 1
    LocalBasisMatrix(:,LocalAxis) = LocalBasisMatrix(:,LocalAxis)/norm(LocalBasisMatrix(:,LocalAxis));
end
for LocalAxis = 2
    LocalBasisMatrix(:,LocalAxis) = LocalBasisMatrix(:,LocalAxis)-dot(LocalBasisMatrix(:,LocalAxis-1),LocalBasisMatrix(:,LocalAxis))*LocalBasisMatrix(:,LocalAxis-1);
    LocalBasisMatrix(:,LocalAxis) = LocalBasisMatrix(:,LocalAxis)/norm(LocalBasisMatrix(:,LocalAxis));
end

x1min=0;
x2min=0;
LocalAxis1 = 1;
x1max = dot(LocalNodePosM(:,2)-LocalOrigin,LocalBasisMatrix(:,LocalAxis1));
LocalAxis2 = 2;
x2max = dot(LocalNodePosM(:,4)-LocalOrigin,LocalBasisMatrix(:,LocalAxis2));
if abs(dot(LocalNodePosM(:,4)-LocalOrigin,LocalBasisMatrix(:,LocalAxis1))-x1max)<EPSILON ||...
        abs(dot(LocalNodePosM(:,4)-LocalOrigi,LocalBasisMatrix(:,LocalAxis2))-x2max)<EPSILON
    disp('VectorRectangleSurfaceIntegral:Non-rectangle integration domains are currently not available.')
    return;
end

NormLineVec = cross(LocalBasisMatrix(:,1),LocalBasisMatrix(:,2));
FuncOnSurface = @(x1,x2) ...
    (1/norm(NormLineVec))*...
    (NormLineVec(1)*TwoFormFuncHandle(1,...
    x1*LocalBasisMatrix(1,LocalAxis1)+x2*LocalBasisMatrix(1,LocalAxis2)+LocalOrigin(1),...
    x1*LocalBasisMatrix(1,LocalAxis1)+x2*LocalBasisMatrix(1,LocalAxis2)+LocalOrigin(1),...
    x1*LocalBasisMatrix(1,LocalAxis1)+x2*LocalBasisMatrix(1,LocalAxis2)+LocalOrigin(1))...
    +NormLineVec(3)*TwoFormFuncHandle(2,...
    x1*LocalBasisMatrix(2,LocalAxis1)+x2*LocalBasisMatrix(2,LocalAxis12)+LocalOrigin(2),...
    x1*LocalBasisMatrix(2,LocalAxis1)+x2*LocalBasisMatrix(2,LocalAxis2)+LocalOrigin(2),...
    x1*LocalBasisMatrix(2,LocalAxis1)+x2*LocalBasisMatrix(2,LocalAxis2)+LocalOrigin(2))...
    +NormLineVec(3)*TwoFormFuncHandle(3,...
    x1*LocalBasisMatrix(3,LocalAxis1)+x2*LocalBasisMatrix(3,LocalAxis2)+LocalOrigin(3),...
    x1*LocalBasisMatrix(3,LocalAxis1)+x2*LocalBasisMatrix(3,LocalAxis2)+LocalOrigin(3),...
    x1*LocalBasisMatrix(3,LocalAxis1)+x2*LocalBasisMatrix(3,LocalAxis2)+LocalOrigin(3))...
);

IntegratedValue = integral2(FuncOnSurface,x1min,x1max,x2min,x2max);
end

function value = func_sigma_231(dim,x,y,z)

switch dim 
    case 1
        value =func_sigma_123(2,x,y,z);
    case 2
        value =func_sigma_123(3,x,y,z);
    case 3
        value =func_sigma_123(1,x,y,z);
end
end

function value = func_sigma_312(dim,x,y,z)

switch dim 
    case 1
        value =func_sigma_123(3,x,y,z);
    case 2
        value =func_sigma_123(1,x,y,z);
    case 3
        value =func_sigma_123(2,x,y,z);
end
end