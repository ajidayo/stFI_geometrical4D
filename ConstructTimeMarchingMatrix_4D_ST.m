function TMM = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Sigma,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem)
disp('ConstructTimeMarchingMatrix_4D_ST: Composes the explicit time-marching matrix from each DoF-calculating tasks')
TMM_Redundant = speye(Num_of_Elem.STP+Num_of_Elem.STSource+Num_of_Elem.PMLImagDualSTP);
for nthTask = 1:size(TaskOrder,2)
    if mod(nthTask-1,10^(-1+round(log10(size(TaskOrder,2)))) ) == 0
        disp(['Processing ', num2str(nthTask),'/',...
            num2str(size(TaskOrder,2)),' th task'])
    end
    TaskIdx = TaskOrder(nthTask);
    switch Task(TaskIdx).Type
        case "InitialCondition"
        case "BoundaryCondition"
        case "STFIPrim_Block"
            disp('TaskType:STFIPrim_Block')
            TMM_Redundant = ...
                Task(TaskIdx).TMMOnestepGenerFunc(Task(TaskIdx),D2,size(TMM_Redundant,1))*TMM_Redundant;
        case "STFIDual_Block"
            disp('TaskType:STFIDual_Block')
            TMM_Redundant = ...
                Task(TaskIdx).TMMOnestepGenerFunc(Task(TaskIdx),D1,Kappa_over_Z,size(TMM_Redundant,1))*TMM_Redundant;
        case "PML"
            disp('TaskType:PML')
             TMM_Redundant = SingleTask_PML(Task(TaskIdx),sC,Kappa_over_Z,Sigma,SpElemProperties,Num_of_Elem,TMM_Redundant);
        case "ST_FI_Prim"
            TMM_Redundant = SingleTask_STFI_Prim(Task(TaskIdx),D2,TMM_Redundant);
        case "ST_FI_Dual"
            TMM_Redundant = SingleTask_STFI_Dual(Task(TaskIdx),D1,Kappa_over_Z,TMM_Redundant);
        case "Sp_FI"
            disp('TaskType:SpFI')
            TMM_Redundant = SingleTask_SpFI(Task(TaskIdx),sC,Kappa_over_Z,Source,SpElemProperties,Num_of_Elem,TMM_Redundant);
        otherwise
            disp("No such Tasktype defined.")
    end
end

TMM = RemoveRedundantDoF(TMM_Redundant,SpElemProperties,Num_of_Elem);

end

%%
function TMM_Redundant = SingleTask_STFI_Prim(STFI_PrimTask,D2,TMM_Redundant)
%disp('SingleTask_STFI_Prim CALLED')
TMM_Onestep=speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
STPNum = size(D2,2);
TMM_Onestep(STFI_PrimTask.STPtgt,1:STPNum) = ...
    -(1.0/D2(STFI_PrimTask.STHtgt,STFI_PrimTask.STPtgt))*D2(STFI_PrimTask.STHtgt,:);
TMM_Onestep(STFI_PrimTask.STPtgt,STFI_PrimTask.STPtgt)=0;
TMM_Redundant=TMM_Onestep*TMM_Redundant;
end

%%
function TMM_Redundant= SingleTask_STFI_Dual(STFI_DualTask,D1,Kappa_over_Z,TMM_Redundant)
%disp('SingleTask_STFI_Dual CALLED')
TMM_Onestep=speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
STPNum = size(D1,1);
TMM_Onestep(STFI_DualTask.STPtgt,1:STPNum) = ...
    -(1.0/(D1(STFI_DualTask.STPtgt,STFI_DualTask.STStgt)*Kappa_over_Z(STFI_DualTask.STPtgt)))...
    *(D1(:,STFI_DualTask.STStgt).'*spdiags(Kappa_over_Z,0,STPNum,STPNum));
TMM_Onestep(STFI_DualTask.STPtgt,STFI_DualTask.STPtgt)=0;
TMM_Redundant=TMM_Onestep*TMM_Redundant;
end

%%
function TMM_Redundant = ...
    SingleTask_SpFI(Task,sC,Kappa_over_Z,Source,SpElemProperties,Num_of_Elem,TMM_Redundant)
%disp('SingleTask_SpFI CALLED')
%TMM_Onestep      = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));

%Sp_FI_Region_tgt = SpFI_TaskInfo.Sp_FI_Region_tgt;
TimeSection_tgt  = Task.TimeSection_tgt;

ConservedSTP_Faraday = true(Num_of_Elem.STP,1);
ConservedSTP_Ampere  = true(Num_of_Elem.STP,1);
Map_STPFutr_SpP         = logical(sparse(Num_of_Elem.STP,Num_of_Elem.SpP));
Map_SpP_STPPast         = logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));
Map_STPFutr_SpS         = logical(sparse(Num_of_Elem.STP,Num_of_Elem.SpS));
Map_SpS_STPPast         = logical(sparse(Num_of_Elem.SpS,Num_of_Elem.STP));
Map_SpPinctoSpS_STPPast = logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));
Map_SpSinctoSpP_STPFutr = logical(sparse(Num_of_Elem.SpS,Num_of_Elem.STP));

for SpPIdx  = Task.SpPtar
    STPFutr = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSection_tgt;
    STPPast = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSection_tgt-1;
    Map_STPFutr_SpP(STPFutr,SpPIdx)    = true;
    Map_SpP_STPPast(SpPIdx,STPPast)    = true;
    ConservedSTP_Faraday(STPFutr)      = false;
    for SpSIdx  = find(sC(SpPIdx,:))
        STPFutr = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt;
        Map_SpSinctoSpP_STPFutr(SpSIdx,STPFutr)    = true;
    end
end
for SpSIdx  = Task.SpStar
    STPFutr = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt;
    STPPast = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt-1;
    Map_STPFutr_SpS(STPFutr,SpSIdx)    = true;
    Map_SpS_STPPast(SpSIdx,STPPast)    = true;
    ConservedSTP_Ampere(STPFutr)      = false;
    for SpPIdx  = find(sC(:,SpSIdx)).'
        STPPast = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSection_tgt-1;
        Map_SpPinctoSpS_STPPast(SpPIdx,STPPast)    = true;
    end
end

z       = spdiags(Kappa_over_Z.^(-1),0,Num_of_Elem.STP,Num_of_Elem.STP);
zinv    = spdiags(Kappa_over_Z      ,0,Num_of_Elem.STP,Num_of_Elem.STP);

STSourceNum=0;
for SpSourceIdx = 1:size(Source,2) 
    STSourceNum=STSourceNum+Source(SpSourceIdx).UpdNum;
end

SourcePattern = sparse(Num_of_Elem.SpS,STSourceNum);
ST_SourceIdx = 0;
for SpSourceIdx = 1:size(Source,2)
    for CurrentTimeSec = 1:Source(SpSourceIdx).UpdNum 
        ST_SourceIdx = ST_SourceIdx+1;
        if ismember(Source(SpSourceIdx).DualFace_tgt,Task.SpStar) && CurrentTimeSec == TimeSection_tgt
            SourcePattern(Source(SpSourceIdx).DualFace_tgt,ST_SourceIdx)=1;
        end
    end
end
% the minus for the dual grid is caused by the sign [z] has.
TMM_Onestep = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(1:Num_of_Elem.STP,1:Num_of_Elem.STP) = ...
    z*Map_STPFutr_SpS*(Map_SpS_STPPast - sC.'*Map_SpPinctoSpS_STPPast)*zinv...
    +spdiags(ConservedSTP_Ampere,0,Num_of_Elem.STP,Num_of_Elem.STP);
TMM_Onestep(1:Num_of_Elem.STP,Num_of_Elem.STP+1:Num_of_Elem.STP+STSourceNum)...
    = z*Map_STPFutr_SpS*SourcePattern;
TMM_Redundant   = TMM_Onestep * TMM_Redundant;

TMM_Onestep = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(1:Num_of_Elem.STP,1:Num_of_Elem.STP) = ...
    Map_STPFutr_SpP*(Map_SpP_STPPast - sC  *Map_SpSinctoSpP_STPFutr)...
    +spdiags(ConservedSTP_Faraday,0,Num_of_Elem.STP,Num_of_Elem.STP);
TMM_Redundant   = TMM_Onestep * TMM_Redundant;

end

%%
function TMM_Redundant = ...
    SingleTask_PML(PML_TaskInfo,sC,Kappa_over_Z,Sigma,SpElemProperties,Num_of_Elem,TMM_Redundant)

TimeSection_tgt   = PML_TaskInfo.TimeSection_tgt;

ConservedPrimDoF_Faraday        = true(Num_of_Elem.STP,1);
ConservedPrimDoF_ConstiG2F      = true(Num_of_Elem.STP,1);
ConservedPMLDualDoF_ConstiF2G   = true(Num_of_Elem.PMLImagDualSTP,1);
ConservedPMLDualDoF_Ampere      = true(Num_of_Elem.PMLImagDualSTP,1);

% Num_of_DoF = Num_of_Elem.STP + Num_of_Elem.STSource + Num_of_Elem.PMLImagDualSTP;

Map_DoFFutr_SpP             = logical(sparse(Num_of_Elem.STP            ,Num_of_Elem.SpP));
Map_DoFPast_SpP             = logical(sparse(Num_of_Elem.STP            ,Num_of_Elem.SpP));
Map_DoFFutr_SpS             = logical(sparse(Num_of_Elem.STP            ,Num_of_Elem.SpS));
Map_DoFPast_SpS             = logical(sparse(Num_of_Elem.STP            ,Num_of_Elem.SpS));
Map_PMLDualDoFPast_SpP      = logical(sparse(Num_of_Elem.PMLImagDualSTP ,Num_of_Elem.SpP));
Map_PMLDualDoFPast_SpS      = logical(sparse(Num_of_Elem.PMLImagDualSTP ,Num_of_Elem.SpS));
Map_PMLDualDoFFutr_SpP   = logical(sparse(Num_of_Elem.PMLImagDualSTP ,Num_of_Elem.SpP));
Map_PMLDualDoFFutr_SpS   = logical(sparse(Num_of_Elem.PMLImagDualSTP ,Num_of_Elem.SpS));
Map_DoFFutr_IncSpS          = logical(sparse(Num_of_Elem.STP            ,Num_of_Elem.SpS));
Map_PMLDualDoFPast_IncSpP   = logical(sparse(Num_of_Elem.PMLImagDualSTP ,Num_of_Elem.SpP));
Map_DoFPast_IncSpP          = logical(sparse(Num_of_Elem.STP            ,Num_of_Elem.SpP));

for SpPIdx  = PML_TaskInfo.SpPtar
    DoFFutr         = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)          +TimeSection_tgt  ;
    DoFPast         = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)          +TimeSection_tgt-1;
    PMLDualDoFFutr  = SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx)  +TimeSection_tgt  ;
    PMLDualDoFPast  = SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx)  +TimeSection_tgt-1;
    
    Map_DoFFutr_SpP(DoFFutr,SpPIdx)                     = true;
    Map_DoFPast_SpP(DoFPast,SpPIdx)                     = true;
    Map_PMLDualDoFFutr_SpP(PMLDualDoFFutr,SpPIdx)       = true;
    Map_PMLDualDoFPast_SpP(PMLDualDoFPast,SpPIdx)       = true;
    ConservedPrimDoF_Faraday(DoFFutr)                   = false;
    ConservedPMLDualDoF_ConstiF2G(PMLDualDoFFutr)       = false;    
    for IncSpSIdx  = find(sC(SpPIdx,:))
        DoFFutr = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx)+TimeSection_tgt;
        Map_DoFFutr_IncSpS(DoFFutr,IncSpSIdx)    = true;
    end
end

for SpSIdx  = PML_TaskInfo.SpStar
    DoFFutr         = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)          +TimeSection_tgt  ;
    DoFPast         = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)          +TimeSection_tgt-1;
    PMLDualDoFFutr  = SpElemProperties.SpS.FirstPMLImagDualSTP(SpSIdx)  +TimeSection_tgt  ;
    PMLDualDoFPast  = SpElemProperties.SpS.FirstPMLImagDualSTP(SpSIdx)  +TimeSection_tgt-1;
    
    Map_DoFFutr_SpS(DoFFutr,SpSIdx)                     = true;
    Map_DoFPast_SpS(DoFPast,SpSIdx)                     = true;
    Map_PMLDualDoFFutr_SpS(PMLDualDoFFutr,SpSIdx)       = true;
    Map_PMLDualDoFPast_SpS(PMLDualDoFPast,SpSIdx)       = true;
    ConservedPrimDoF_ConstiG2F(DoFFutr)                 = false;
    ConservedPMLDualDoF_Ampere(PMLDualDoFFutr)          = false;
    for IncSpPIdx  = find(sC(:,SpSIdx).')
        switch SpElemProperties.SpP.PML(IncSpPIdx)
            case true
                PMLDualDoFPast  = SpElemProperties.SpP.FirstPMLImagDualSTP(IncSpPIdx)   +TimeSection_tgt-1;
                Map_PMLDualDoFPast_IncSpP(PMLDualDoFPast,IncSpPIdx)     = true;
            case false
                PrimSTPPast = SpElemProperties.SpP.FirstSTPIdx(IncSpPIdx) +TimeSection_tgt-1;
                Map_DoFPast_IncSpP(PrimSTPPast,IncSpPIdx) = true;
        end
    end
end

PMLDualDoF_Idxs = ...
    Num_of_Elem.STP+Num_of_Elem.STSource+1:...
    Num_of_Elem.STP+Num_of_Elem.STSource+Num_of_Elem.PMLImagDualSTP;
PrimDoF_Idxs = 1:Num_of_Elem.STP;

PMLON = true

Z0 = 1;

% Ampere
% disp('Ampere')
PMLCoeff_Ampere_gSpPast = (1-0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.ThreeOneTwo)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.ThreeOneTwo).^(-1);
PMLCoeff_Ampere_gSpPast = spdiags(PMLCoeff_Ampere_gSpPast,0,Num_of_Elem.SpS,Num_of_Elem.SpS);
PMLCoeff_Ampere_gTi     = (1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.ThreeOneTwo).^(-1);
PMLCoeff_Ampere_gTi     = spdiags(PMLCoeff_Ampere_gTi,0,Num_of_Elem.SpP,Num_of_Elem.SpP);
if ~PMLON
    PMLCoeff_Ampere_gSpPast = speye(Num_of_Elem.SpS);
    PMLCoeff_Ampere_gTi     = speye(Num_of_Elem.SpP);
end

TMM_Onestep = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(PMLDualDoF_Idxs,PMLDualDoF_Idxs) = ...
    Map_PMLDualDoFFutr_SpS* ...
    (PMLCoeff_Ampere_gSpPast*Map_PMLDualDoFPast_SpS.' -sC.'*PMLCoeff_Ampere_gTi*Map_PMLDualDoFPast_IncSpP.');
TMM_Onestep(PMLDualDoF_Idxs,PrimDoF_Idxs) = ...
    Map_PMLDualDoFFutr_SpS* ...
    (-sC.'*PMLCoeff_Ampere_gTi*Map_DoFPast_IncSpP.'*spdiags(Kappa_over_Z,0,Num_of_Elem.STP,Num_of_Elem.STP));
TMM_Onestep(PMLDualDoF_Idxs, PMLDualDoF_Idxs) ...
    = TMM_Onestep(PMLDualDoF_Idxs, PMLDualDoF_Idxs) +spdiags(ConservedPMLDualDoF_Ampere,0,Num_of_Elem.PMLImagDualSTP,Num_of_Elem.PMLImagDualSTP);
TMM_Redundant = TMM_Onestep*TMM_Redundant;

% Constitutive
% disp('Constitutive')

PMLCoeff_ConstiG2F_fTiPast = (1-0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.TwoThreeOne)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.TwoThreeOne).^(-1);
PMLCoeff_ConstiG2F_fTiPast = spdiags(PMLCoeff_ConstiG2F_fTiPast,0,Num_of_Elem.SpS,Num_of_Elem.SpS);
PMLCoeff_ConstiG2F_gSpFutr = (1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.OneTwoThree)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.TwoThreeOne).^(-1);
PMLCoeff_ConstiG2F_gSpFutr = spdiags(PMLCoeff_ConstiG2F_gSpFutr,0,Num_of_Elem.SpS,Num_of_Elem.SpS);
PMLCoeff_ConstiG2F_gSpPast= (1-0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.OneTwoThree)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.TwoThreeOne).^(-1);
PMLCoeff_ConstiG2F_gSpPast = spdiags(PMLCoeff_ConstiG2F_gSpPast,0,Num_of_Elem.SpS,Num_of_Elem.SpS);
if ~PMLON
    PMLCoeff_ConstiG2F_fTiPast = speye(Num_of_Elem.SpS);
    PMLCoeff_ConstiG2F_gSpFutr = speye(Num_of_Elem.SpS);
    PMLCoeff_ConstiG2F_gSpPast = speye(Num_of_Elem.SpS);
end

ZOverKappa_PrimDoFFutr = sparse(Num_of_Elem.STP,Num_of_Elem.STP);
for SpSIdx = find(SpElemProperties.SpS.PML)
    for PrimDoFIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt
        ZOverKappa_PrimDoFFutr(PrimDoFIdx,PrimDoFIdx) = Kappa_over_Z(PrimDoFIdx).^(-1);
    end
end
TMM_Onestep = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(PrimDoF_Idxs,PrimDoF_Idxs) = Map_DoFFutr_SpS*PMLCoeff_ConstiG2F_fTiPast*Map_DoFPast_SpS.';
TMM_Onestep(PrimDoF_Idxs,PMLDualDoF_Idxs)...
    = ZOverKappa_PrimDoFFutr*Map_DoFFutr_SpS*...
    (PMLCoeff_ConstiG2F_gSpFutr*Map_PMLDualDoFFutr_SpS.' -PMLCoeff_ConstiG2F_gSpPast*Map_PMLDualDoFPast_SpS.');
TMM_Onestep(PrimDoF_Idxs,PrimDoF_Idxs) = ...
    TMM_Onestep(PrimDoF_Idxs,PrimDoF_Idxs)+spdiags(ConservedPrimDoF_ConstiG2F,0,Num_of_Elem.STP,Num_of_Elem.STP);
TMM_Redundant = TMM_Onestep*TMM_Redundant;

% Faraday
% disp('Faraday')
PMLCoeff_Faraday_fSpPast = (1-0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.ThreeOneTwo)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.ThreeOneTwo).^(-1);
PMLCoeff_Faraday_fSpPast = spdiags(PMLCoeff_Faraday_fSpPast,0,Num_of_Elem.SpP,Num_of_Elem.SpP);
PMLCoeff_Faraday_fTi     = (1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpS.ThreeOneTwo).^(-1);
PMLCoeff_Faraday_fTi = spdiags(PMLCoeff_Faraday_fTi,0,Num_of_Elem.SpS,Num_of_Elem.SpS);

if ~PMLON
    PMLCoeff_Faraday_fSpPast    = speye(Num_of_Elem.SpP);
    PMLCoeff_Faraday_fTi        = speye(Num_of_Elem.SpS);
end

TMM_Onestep = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(PrimDoF_Idxs,PrimDoF_Idxs) = ...
    Map_DoFFutr_SpP*(PMLCoeff_Faraday_fSpPast*Map_DoFPast_SpP.' -sC*PMLCoeff_Faraday_fTi*Map_DoFFutr_IncSpS.');
TMM_Onestep(PrimDoF_Idxs,PrimDoF_Idxs) = ...
    TMM_Onestep(PrimDoF_Idxs,PrimDoF_Idxs) +spdiags(ConservedPrimDoF_Faraday,0,Num_of_Elem.STP,Num_of_Elem.STP);
TMM_Redundant = TMM_Onestep*TMM_Redundant;

% Constitutive
% disp('Constitutive')
PMLCoeff_ConstiF2G_gTiPast = (1-0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.TwoThreeOne)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.TwoThreeOne).^(-1);
PMLCoeff_ConstiF2G_gTiPast = spdiags(PMLCoeff_ConstiF2G_gTiPast,0,Num_of_Elem.SpP,Num_of_Elem.SpP);
PMLCoeff_ConstiF2G_fSpFutr = (1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.OneTwoThree)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.TwoThreeOne).^(-1);
PMLCoeff_ConstiF2G_fSpFutr = spdiags(PMLCoeff_ConstiF2G_fSpFutr,0,Num_of_Elem.SpP,Num_of_Elem.SpP);
PMLCoeff_ConstiF2G_fSpPast = (1-0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.OneTwoThree)...
    .*(1+0.5*Z0/PML_TaskInfo.UpdNum*Sigma.SpP.TwoThreeOne).^(-1);
PMLCoeff_ConstiF2G_fSpPast = spdiags(PMLCoeff_ConstiF2G_fSpPast,0,Num_of_Elem.SpP,Num_of_Elem.SpP);

if ~PMLON
    PMLCoeff_ConstiF2G_gTiPast = speye(Num_of_Elem.SpP);
    PMLCoeff_ConstiF2G_fSpFutr = speye(Num_of_Elem.SpP);
    PMLCoeff_ConstiF2G_fSpPast = speye(Num_of_Elem.SpP);
end

KappaOverZ_PMLDualDoFFutr = sparse(Num_of_Elem.PMLImagDualSTP,Num_of_Elem.PMLImagDualSTP);
for SpPIdx = find(SpElemProperties.SpP.PML)
    for DualDoFIdx = SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx)   +TimeSection_tgt
        STPIdx     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)           +TimeSection_tgt;
        KappaOverZ_PMLDualDoFFutr(DualDoFIdx,DualDoFIdx) = Kappa_over_Z(STPIdx);
    end
end
TMM_Onestep = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(PMLDualDoF_Idxs,PMLDualDoF_Idxs) = Map_PMLDualDoFFutr_SpP*PMLCoeff_ConstiF2G_gTiPast*Map_PMLDualDoFPast_SpP.';
TMM_Onestep(PMLDualDoF_Idxs,PrimDoF_Idxs)    = ...
    KappaOverZ_PMLDualDoFFutr*Map_PMLDualDoFFutr_SpP...
    *(PMLCoeff_ConstiF2G_fSpFutr*Map_DoFFutr_SpP.' -PMLCoeff_ConstiF2G_fSpPast*Map_DoFPast_SpP.');
TMM_Onestep(PMLDualDoF_Idxs,PMLDualDoF_Idxs) = ...
    TMM_Onestep(PMLDualDoF_Idxs,PMLDualDoF_Idxs) +spdiags(ConservedPMLDualDoF_ConstiF2G,0,Num_of_Elem.PMLImagDualSTP,Num_of_Elem.PMLImagDualSTP);
TMM_Redundant = TMM_Onestep*TMM_Redundant;


end

%%
function TMM = RemoveRedundantDoF(TMM_Redundant,SpElemProperties,Num_of_Elem)

Num_of_PMLDualDoF_Last = size(find(SpElemProperties.SpP.PML==true),2)+size(find(SpElemProperties.SpS.PML==true),2);
Num_of_PMLDualDoF_Init = Num_of_PMLDualDoF_Last;

Store_Last=logical(sparse([],[],[],...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_PMLDualDoF_Last,...
    Num_of_Elem.STP+Num_of_Elem.STSource+Num_of_Elem.PMLImagDualSTP,...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_PMLDualDoF_Last));
for SpPIdx=1:Num_of_Elem.SpP
    STP_tgt=SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);
    Store_Last(SpPIdx,STP_tgt)=true;
end
for SpSIdx=1:Num_of_Elem.SpS
    STP_tgt=SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    Store_Last(Num_of_Elem.SpP+SpSIdx,STP_tgt)=true;
end
DualDoFLast_Idx = 0;
for SpPIdx = find(SpElemProperties.SpP.PML)
    DualDoFLast_Idx = DualDoFLast_Idx+1;
    PMLDualDoF_tgt = SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);
    Store_Last(Num_of_Elem.SpP+Num_of_Elem.SpS+DualDoFLast_Idx,...
        Num_of_Elem.STP+Num_of_Elem.STSource+PMLDualDoF_tgt)=true;
end
for SpSIdx = find(SpElemProperties.SpS.PML)
    DualDoFLast_Idx = DualDoFLast_Idx+1;
    PMLDualDoF_tgt = SpElemProperties.SpS.FirstPMLImagDualSTP(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    Store_Last(Num_of_Elem.SpP+Num_of_Elem.SpS+DualDoFLast_Idx,...
        Num_of_Elem.STP+Num_of_Elem.STSource+PMLDualDoF_tgt)=true;
end

Store_Init=logical(sparse([],[],[],...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource+Num_of_PMLDualDoF_Init,...
    Num_of_Elem.STP+Num_of_Elem.STSource+Num_of_Elem.PMLImagDualSTP,...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource+Num_of_PMLDualDoF_Init));
for SpPIdx=1:Num_of_Elem.SpP
    STP_tgt=SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
    Store_Init(SpPIdx,STP_tgt)=true; 
end
for SpSIdx=1:Num_of_Elem.SpS
    STP_tgt=SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
    Store_Init(Num_of_Elem.SpP+SpSIdx,STP_tgt)=true;
end
for ST_SourceIdx =1:Num_of_Elem.STSource
    Store_Init(Num_of_Elem.SpP+Num_of_Elem.SpS+ST_SourceIdx,Num_of_Elem.STP+ST_SourceIdx)=true;
end
DualDoFInit_Idx = 0;
for SpPIdx =find(SpElemProperties.SpP.PML)
    DualDoFInit_Idx = DualDoFInit_Idx+1;
    PMLDualDoF_tgt =  SpElemProperties.SpP.FirstPMLImagDualSTP(SpPIdx);
    Store_Init(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource+DualDoFInit_Idx,...
        Num_of_Elem.STP+Num_of_Elem.STSource+PMLDualDoF_tgt)=true;
end
for SpSIdx = find(SpElemProperties.SpS.PML)
    DualDoFInit_Idx = DualDoFInit_Idx+1;
    PMLDualDoF_tgt = SpElemProperties.SpS.FirstPMLImagDualSTP(SpSIdx);
    Store_Init(Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource+DualDoFInit_Idx,...
        Num_of_Elem.STP+Num_of_Elem.STSource+PMLDualDoF_tgt)=true;
end

TMM = Store_Last*TMM_Redundant*Store_Init.';
end