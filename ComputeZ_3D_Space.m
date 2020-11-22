function ZDiscret = ComputeZ_3D_Space(FuncHandle_Impedance,NodePosDual_M,sC,sD,SpElemProperties,Num_of_Elem)
disp('ComputeZ_3D_Space: Calculating Z for each STP')
% global EPSILON
sCLogical = logical(sC);
sDLogical = logical(sD);

ZDiscret = ones(Num_of_Elem.STP,1);

Z_SpV = zeros(1,size(sD,1));
for SpVIdx = 1:size(sD,1)
    Z_SpV(SpVIdx) = FuncHandle_Impedance(NodePosDual_M(:,SpVIdx));
end
SpVIdx_NotVacuum = logical(Z_SpV~=1);
SpPIdx_NotVacuum = find(any(sDLogical(SpVIdx_NotVacuum,:),1));
SpSIdx_NotVacuum = find(any(sCLogical(SpPIdx_NotVacuum,:),1));

for SpPIdx = SpPIdx_NotVacuum
    STPIdx_RefElemIsSpPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)...
        :SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);
    IncSpVLogIdx = logical(sD(:,SpPIdx).');
    ZDiscretInv = size(find(IncSpVLogIdx),2).^(-1)*(Z_SpV.^(-1)*IncSpVLogIdx.');
    ZDiscret(STPIdx_RefElemIsSpPIdx) = ZDiscretInv^(-1);
%    PosVec = NodePosDual_M(:,IncSpVLogIdx);
%    ZDiscret(STPIdx_RefElemIsSpPIdx) ...
%         = (...
%         0.5*FuncHandle_Impedance(PosVec).^(-1)+0.5*FuncHandle_Impedance(PosVec).^(-1)...
%         ).^(-1);
end
IncIncSpVLogIdxMatrix = logical(sCLogical.'*sDLogical.');
for SpSIdx = SpSIdx_NotVacuum
   % ZDiscretInv = 0;
    IncIncSpVLogIdx = IncIncSpVLogIdxMatrix(SpSIdx,:);
    STPIdx_RefElemIsSpSIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)...
        :SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    ZDiscretInv = size(find(IncIncSpVLogIdx),2).^(-1)*(Z_SpV.^(-1)*IncIncSpVLogIdx.');
    ZDiscret(STPIdx_RefElemIsSpSIdx) = ZDiscretInv^(-1);
end
    
end