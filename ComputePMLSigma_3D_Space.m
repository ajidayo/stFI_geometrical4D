function [Sigma] = ComputePMLSigma_3D_Space(SigmaSpV,cdt,NodePosPrim_M,NodePosDual_M,sG,sC,sD,SpElemProperties)
disp('ComputePMLSigma_3D_Space: Evaluating the value of sigma for each element in the PML region.')
% compute RefSigma.{SpP, SpS}.{OneTwoThree, TwoThreeOne, ThreeOneTwo}
% For primal reference-faces(dual   reference-edges), weight sigma by relative lengths.
% For dual   reference-faces(primal reference-edges), weight sigma by relative areas. 
SpPNum = size(SpElemProperties.SpP.PML,2);
SpSNum = size(SpElemProperties.SpS.PML,2);
Sigma_SpP_OneTwoThree = zeros(SpPNum,1);
Sigma_SpP_TwoThreeOne = zeros(SpPNum,1);
Sigma_SpP_ThreeOneTwo = zeros(SpPNum,1);
Sigma_SpS_OneTwoThree = zeros(SpSNum,1);
Sigma_SpS_TwoThreeOne = zeros(SpSNum,1);
Sigma_SpS_ThreeOneTwo = zeros(SpSNum,1);
SigmaCalcMethod = "New"
% SigmaCalcMethod = "Old"
switch SigmaCalcMethod
    case "Old"
        for SpPIdx = find(SpElemProperties.SpP.PML)
            if SpElemProperties.SpP.ElecWall(SpPIdx)
                continue
            end
            DualNode_Sta = find(sD(:,SpPIdx).'==-1,1);
            DualNode_Tgt = find(sD(:,SpPIdx).'== 1,1);
            ApproxPoint = 0.5* (NodePosDual_M(:,DualNode_Tgt) + NodePosDual_M(:,DualNode_Sta));
            LineDirec     =       NodePosDual_M(:,DualNode_Tgt) - NodePosDual_M(:,DualNode_Sta);
            LineDirec     = norm(LineDirec)^(-1)*LineDirec;
            Sigma_SpP_OneTwoThree(SpPIdx) = ...
                (LineDirec(1)*func_sigma_123(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(1)...
                +LineDirec(2)*func_sigma_123(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(2)...
                +LineDirec(3)*func_sigma_123(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(3)...
                );%/norm(LineDirec);
            Sigma_SpP_TwoThreeOne(SpPIdx) = ...
                (LineDirec(1)*func_sigma_231(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(1)...
                +LineDirec(2)*func_sigma_231(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(2)...
                +LineDirec(3)*func_sigma_231(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(3)...
                );%/norm(LineDirec);
            Sigma_SpP_ThreeOneTwo(SpPIdx) = ...
                (LineDirec(1)*func_sigma_312(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(1)...
                +LineDirec(2)*func_sigma_312(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(2)...
                +LineDirec(3)*func_sigma_312(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(3)...
                );%/norm(LineDirec);
            clearvars PrimNode_Sta PrimNode_Tgt
        end
        for SpSIdx = find(SpElemProperties.SpS.PML)
            if SpElemProperties.SpS.PEC(SpSIdx)
                continue
            end
            PrimNode_Sta = find(sG(SpSIdx,:)==-1,1);
            PrimNode_Tgt = find(sG(SpSIdx,:)== 1,1);
            ApproxPoint = 0.5* (NodePosPrim_M(:,PrimNode_Tgt) + NodePosPrim_M(:,PrimNode_Sta).Vec);
            LineDirec     =       NodePosPrim_M(:,PrimNode_Tgt) - NodePosPrim_M(:,PrimNode_Sta);
            LineDirec     = norm(LineDirec)^(-1)*LineDirec;
            Sigma_SpS_OneTwoThree(SpSIdx) = ...
                (LineDirec(1)*func_sigma_123(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(1)...
                +LineDirec(2)*func_sigma_123(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(2)...
                +LineDirec(3)*func_sigma_123(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(3)...
                );
            Sigma_SpS_TwoThreeOne(SpSIdx) = ...
                (LineDirec(1)*func_sigma_231(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(1)...
                +LineDirec(2)*func_sigma_231(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(2)...
                +LineDirec(3)*func_sigma_231(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(3)...
                );
            Sigma_SpS_ThreeOneTwo(SpSIdx) = ...
                (LineDirec(1)*func_sigma_312(1,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(1)...
                +LineDirec(2)*func_sigma_312(2,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(2)...
                +LineDirec(3)*func_sigma_312(3,ApproxPoint(1),ApproxPoint(2),ApproxPoint(3))*LineDirec(3)...
                );
            clearvars PrimNode_Sta PrimNode_Tgt
        end
    case "New"
        for SpPIdx = find(SpElemProperties.SpP.PML)
            if SpElemProperties.SpP.ElecWall(SpPIdx)
                continue
            end
            DualEdgeDirec = NodePosDual_M*sD(:,SpPIdx);
            DualEdgeDirec = norm(DualEdgeDirec).^(-1)*DualEdgeDirec;
            SigmaTemp123 = 0.5*(SigmaSpV*logical(sD(:,SpPIdx)));
            SigmaTemp231 = [SigmaTemp123(2);SigmaTemp123(3);SigmaTemp123(1)];
            SigmaTemp312 = [SigmaTemp123(3);SigmaTemp123(1);SigmaTemp123(2)];
            Sigma_SpP_OneTwoThree(SpPIdx) = dot(DualEdgeDirec.*SigmaTemp123,DualEdgeDirec);
            Sigma_SpP_TwoThreeOne(SpPIdx) = dot(DualEdgeDirec.*SigmaTemp231,DualEdgeDirec);
            Sigma_SpP_ThreeOneTwo(SpPIdx) = dot(DualEdgeDirec.*SigmaTemp312,DualEdgeDirec);
        end
        sCLogical = logical(sC);
        sDLogical = logical(sD);
        for SpSIdx = find(SpElemProperties.SpS.PML)
            if SpElemProperties.SpS.PEC(SpSIdx)
                continue
            end
            EdgeDirec = NodePosPrim_M*sG(SpSIdx,:).';
            EdgeDirec = norm(EdgeDirec).^(-1)*EdgeDirec;
            IncIncSpVLogIdx = logical(sDLogical*sCLogical(:,SpSIdx)).';
            SigmaTemp123 = size(find(IncIncSpVLogIdx),2).^(-1)*(SigmaSpV*IncIncSpVLogIdx.');
            SigmaTemp231 = [SigmaTemp123(2);SigmaTemp123(3);SigmaTemp123(1)];
            SigmaTemp312 = [SigmaTemp123(3);SigmaTemp123(1);SigmaTemp123(2)];
            Sigma_SpS_OneTwoThree(SpSIdx) = dot(EdgeDirec.*SigmaTemp123,EdgeDirec);
            Sigma_SpS_TwoThreeOne(SpSIdx) = dot(EdgeDirec.*SigmaTemp231,EdgeDirec);
            Sigma_SpS_ThreeOneTwo(SpSIdx) = dot(EdgeDirec.*SigmaTemp312,EdgeDirec);
        end
    otherwise 
        warning('Sigma not calculated')
end

Sigma.SpP.OneTwoThree = cdt * Sigma_SpP_OneTwoThree;
Sigma.SpP.TwoThreeOne = cdt * Sigma_SpP_TwoThreeOne;
Sigma.SpP.ThreeOneTwo = cdt * Sigma_SpP_ThreeOneTwo;
Sigma.SpS.OneTwoThree = cdt * Sigma_SpS_OneTwoThree;
Sigma.SpS.TwoThreeOne = cdt * Sigma_SpS_TwoThreeOne;
Sigma.SpS.ThreeOneTwo = cdt * Sigma_SpS_ThreeOneTwo;


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
% function IntegratedValue = VectorStraightLineIntegral(OneFormFuncHandle,StaNodePos,TgtNodePos)
% LineVec = TgtNodePos-StaNodePos;
% 
% FuncAlongLine = @(x) ...
%     (abs(LineVec(1))*OneFormFuncHandle(1,LineVec(1)*x+StaNodePos(1),LineVec(2)*x+StaNodePos(2),LineVec(3)*x+StaNodePos(3))...
%     +abs(LineVec(2))*OneFormFuncHandle(2,LineVec(1)*x+StaNodePos(1),LineVec(2)*x+StaNodePos(2),LineVec(3)*x+StaNodePos(3))...
%     +abs(LineVec(3))*OneFormFuncHandle(3,LineVec(1)*x+StaNodePos(1),LineVec(2)*x+StaNodePos(2),LineVec(3)*x+StaNodePos(3))...
%     )/norm(LineVec);
% 
% IntegratedValue = integral(FuncAlongLine,0,1);
% 
% end
% 
% function IntegratedValue = VectorRectangleSurfaceIntegral(TwoFormFuncHandle,LocalGradM,LocalNodePosM)
% global EPSILON
% 
% if size(LocalNodePosM,2)~=4
%     disp('VectorRectangleSurfaceIntegral:Non-quadrangle integration domains are currently not available.')
%     return;
% end
% % N.B. Nodes must be sorted along the orientation of the integration area.
% LocalOrigin = LocalNodePosM(:,1);
% for EdgeIdx = find(LocalGradM(:,1).')
%     if LocalGradM(EdgeIdx,2)~=0
%         LocalBasisMatrix(:,1) = -LocalGradM(EdgeIdx,1)*(LocalGradM(EdgeIdx,:)*LocalNodePosM.').';
%     elseif LocalGradM(EdgeIdx,4)~=0
%         LocalBasisMatrix(:,2) = -LocalGradM(EdgeIdx,1)*(LocalGradM(EdgeIdx,:)*LocalNodePosM.').';
%     end
% end
% for LocalAxis = 1
%     LocalBasisMatrix(:,LocalAxis) = LocalBasisMatrix(:,LocalAxis)/norm(LocalBasisMatrix(:,LocalAxis));
% end
% for LocalAxis = 2
%     LocalBasisMatrix(:,LocalAxis) = LocalBasisMatrix(:,LocalAxis)-dot(LocalBasisMatrix(:,LocalAxis-1),LocalBasisMatrix(:,LocalAxis))*LocalBasisMatrix(:,LocalAxis-1);
%     LocalBasisMatrix(:,LocalAxis) = LocalBasisMatrix(:,LocalAxis)/norm(LocalBasisMatrix(:,LocalAxis));
% end
% 
% x1min=0;
% x2min=0;
% LocalAxis1 = 1;
% x1max = dot(LocalNodePosM(:,2)-LocalOrigin,LocalBasisMatrix(:,LocalAxis1));
% LocalAxis2 = 2;
% x2max = dot(LocalNodePosM(:,4)-LocalOrigin,LocalBasisMatrix(:,LocalAxis2));
% if abs(dot(LocalNodePosM(:,4)-LocalOrigin,LocalBasisMatrix(:,LocalAxis1))-x1max)<EPSILON ||...
%         abs(dot(LocalNodePosM(:,4)-LocalOrigi,LocalBasisMatrix(:,LocalAxis2))-x2max)<EPSILON
%     disp('VectorRectangleSurfaceIntegral:Non-rectangle integration domains are currently not available.')
%     return;
% end
% 
% NormLineVec = cross(LocalBasisMatrix(:,1),LocalBasisMatrix(:,2));
% FuncOnSurface = @(x1,x2) ...
%     (1/norm(NormLineVec))*...
%     (NormLineVec(1)*TwoFormFuncHandle(1,...
%     x1*LocalBasisMatrix(1,LocalAxis1)+x2*LocalBasisMatrix(1,LocalAxis2)+LocalOrigin(1),...
%     x1*LocalBasisMatrix(1,LocalAxis1)+x2*LocalBasisMatrix(1,LocalAxis2)+LocalOrigin(1),...
%     x1*LocalBasisMatrix(1,LocalAxis1)+x2*LocalBasisMatrix(1,LocalAxis2)+LocalOrigin(1))...
%     +NormLineVec(3)*TwoFormFuncHandle(2,...
%     x1*LocalBasisMatrix(2,LocalAxis1)+x2*LocalBasisMatrix(2,LocalAxis12)+LocalOrigin(2),...
%     x1*LocalBasisMatrix(2,LocalAxis1)+x2*LocalBasisMatrix(2,LocalAxis2)+LocalOrigin(2),...
%     x1*LocalBasisMatrix(2,LocalAxis1)+x2*LocalBasisMatrix(2,LocalAxis2)+LocalOrigin(2))...
%     +NormLineVec(3)*TwoFormFuncHandle(3,...
%     x1*LocalBasisMatrix(3,LocalAxis1)+x2*LocalBasisMatrix(3,LocalAxis2)+LocalOrigin(3),...
%     x1*LocalBasisMatrix(3,LocalAxis1)+x2*LocalBasisMatrix(3,LocalAxis2)+LocalOrigin(3),...
%     x1*LocalBasisMatrix(3,LocalAxis1)+x2*LocalBasisMatrix(3,LocalAxis2)+LocalOrigin(3))...
% );
% 
% IntegratedValue = integral2(FuncOnSurface,x1min,x1max,x2min,x2max);
% end

function value = func_sigma_231(PosVec)
value = func_sigma_123(PosVec);
value = [value(2);value(3);value(1)];
end

function value = func_sigma_312(PosVec)
value =func_sigma_123(PosVec);
value = [value(3);value(1);value(2)];
end