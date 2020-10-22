function [Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateTask_STFI_4DST_Block(SpElemProperties,STElemProperties,Task,TaskDepGraph,sC,D2,D1)

if any(SpElemProperties.SpP.Belong_to_ST_FI) == false ...
        && any(SpElemProperties.SpS.Belong_to_ST_FI) == false
    return;
end

TaskFieldNames = fieldnames(Task(size(Task,2)));
switch isempty(Task(size(Task,2)).(TaskFieldNames{1}))
    case 1
        GlobalTaskOffset = 0;
    otherwise
        GlobalTaskOffset = size(Task,2);
end

StaTask = [];
TgtTask = [];
%%

[Task_Local,STElemProperties.STP.TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIPrim_Block(SpElemProperties,GlobalTaskOffset,sC,D2,STElemProperties.STP.TaskIdx); 
Task = [Task Task_Local];
StaTask = [StaTask StaTask_Local];
TgtTask = [TgtTask TgtTask_Local];

%%

GlobalTaskOffset = size(Task,2);

[Task_Local,STElemProperties.STP.TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIDual_Block(SpElemProperties,GlobalTaskOffset,sC,D1,STElemProperties.STP.TaskIdx);
Task = [Task Task_Local];
StaTask = [StaTask StaTask_Local];
TgtTask = [TgtTask TgtTask_Local];

if isempty(StaTask) == false
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask


end


function [Task_Local,STElemProperties_STP_TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIPrim_Block(SpElemProperties,GlobalTaskOffset,sC,D2,STElemProperties_STP_TaskIdx) 
%disp('GenerateTask_STFIPrim_Block Called')
LocalTaskOffset = 0;

StaTask_Local = [];
TgtTask_Local = [];
TaskDepEdgeNum_Local = 0;

Num_of_STP = size(D2,2);

sCPattern_Modified = logical(sC);
sCPattern_Modified(~logical(SpElemProperties.SpP.Belong_to_ST_FI),:) = false;
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    for IncSpPIdx = find(sC(:,SpSIdx).')
        if SpElemProperties.SpP.UpdNum(IncSpPIdx) ~= SpElemProperties.SpS.UpdNum(SpSIdx)
            sCPattern_Modified(IncSpPIdx,SpSIdx) = false;
        end
    end
end
Graph_SpPviaSpS_SomeGraphEdgeExcluded = graph(sCPattern_Modified*sCPattern_Modified.','omitselfloops');
[SubgraphBin_SpP,SpPNumInSubgraph] = conncomp(Graph_SpPviaSpS_SomeGraphEdgeExcluded);
for SubgraphIdx = 1:size(SpPNumInSubgraph,2)
    SpPIdxListInSubgraph = find(SubgraphBin_SpP==SubgraphIdx);
    RepresentiveSpP = SpPIdxListInSubgraph(1);
    if SpElemProperties.SpP.Belong_to_ST_FI(RepresentiveSpP) == false
        continue;
    elseif SpElemProperties.SpP.ElecWall(RepresentiveSpP) == true
        continue;
    end
    %disp(['SpPs in SubgraphIdx = ',num2str(SubgraphIdx),' is STFI Subgraph'])
    STPtgt_OffSet = SpElemProperties.SpP.FirstSTPIdx(SpPIdxListInSubgraph);
    STHtgt_OffSet = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdxListInSubgraph);
    for TimeSection = 1:SpElemProperties.SpP.UpdNum(RepresentiveSpP)
        GlobalTaskIdx = GlobalTaskOffset + TimeSection;
        LocalTaskIdx  = LocalTaskOffset  + TimeSection;
        Task_Local(LocalTaskIdx) = GenerEmptyTask;
        Task_Local(LocalTaskIdx).Type = "STFIPrim_Block";
        Task_Local(LocalTaskIdx).TMMOnestepGenerFunc = @TMMOnestep_STFIPrim_Block;
        Task_Local(LocalTaskIdx).STPtgt = STPtgt_OffSet + TimeSection;
        Task_Local(LocalTaskIdx).STHtgt = STHtgt_OffSet + TimeSection;
        STPtgt_LogIdx = false(1,Num_of_STP);
        STPtgt_LogIdx(Task_Local(LocalTaskIdx).STPtgt) = true;
        Task_Local(LocalTaskIdx).DepSTP ...
            = find(logical(sum(logical(D2(Task_Local(LocalTaskIdx).STHtgt,:)),1)).*~STPtgt_LogIdx);
        STElemProperties_STP_TaskIdx(Task_Local(LocalTaskIdx).STPtgt) ...
            = GlobalTaskIdx;
        if TimeSection > 1
            TaskDepEdgeNum_Local = TaskDepEdgeNum_Local+1;
            StaTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx -1;
            TgtTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx;
        end
    end
    GlobalTaskOffset = GlobalTaskOffset + SpElemProperties.SpP.UpdNum(RepresentiveSpP);
    LocalTaskOffset  = LocalTaskOffset  + SpElemProperties.SpP.UpdNum(RepresentiveSpP);
end
end

function [Task_Local,STElemProperties_STP_TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIDual_Block(SpElemProperties,GlobalTaskOffset,sC,D1,STElemProperties_STP_TaskIdx)
%disp('GenerateTask_STFIDual_Block Called')
LocalTaskOffset = 0;

StaTask_Local = [];
TgtTask_Local = [];
TaskDepEdgeNum_Local = 0;

Num_of_STP = size(D1,1);

sCPattern_Modified = logical(sC);
sCPattern_Modified(:,~logical(SpElemProperties.SpS.Belong_to_ST_FI)) = false;
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpSIdx = find(sC(SpPIdx,:))
        if SpElemProperties.SpS.UpdNum(IncSpSIdx) ~= SpElemProperties.SpP.UpdNum(SpPIdx)
            sCPattern_Modified(SpPIdx,IncSpSIdx) = false;
        end
    end
end
Graph_SpSviaSpP_SomeGraphEdgeExcluded = graph(sCPattern_Modified.'*sCPattern_Modified,'omitselfloops');
[SubgraphBin_SpS,SpSNumInSubgraph] = conncomp(Graph_SpSviaSpP_SomeGraphEdgeExcluded);
for SubgraphIdx = 1:size(SpSNumInSubgraph,2)
    SpSIdxListInSubgraph = find(SubgraphBin_SpS==SubgraphIdx);
    RepresentiveSpS = SpSIdxListInSubgraph(1);
    if SpElemProperties.SpS.Belong_to_ST_FI(RepresentiveSpS) == false
        continue;
    elseif SpElemProperties.SpS.PEC(RepresentiveSpS) == true
        continue;
    end
   % disp(['SpSs in SubgraphIdx = ',num2str(SubgraphIdx),' is STFI Subgraph'])
    STPtgt_OffSet = SpElemProperties.SpS.FirstSTPIdx(SpSIdxListInSubgraph);
    STStgt_OffSet = -1+SpElemProperties.SpS.FirstSTSIdx(SpSIdxListInSubgraph);
    for TimeSection = 1:SpElemProperties.SpS.UpdNum(RepresentiveSpS)
        GlobalTaskIdx = GlobalTaskOffset + TimeSection;
        LocalTaskIdx  = LocalTaskOffset + TimeSection;
        Task_Local(LocalTaskIdx) = GenerEmptyTask;
        Task_Local(LocalTaskIdx).Type = "STFIDual_Block";
        Task_Local(LocalTaskIdx).TMMOnestepGenerFunc = @TMMOnestep_STFIDual_Block;
        Task_Local(LocalTaskIdx).STPtgt = STPtgt_OffSet + TimeSection;
        Task_Local(LocalTaskIdx).STStgt = STStgt_OffSet + TimeSection;
        STPtgt_LogIdx = false(1,Num_of_STP);
        STPtgt_LogIdx(Task_Local(LocalTaskIdx).STPtgt) = true;
        Task_Local(LocalTaskIdx).DepSTP ...
            = find(logical(sum(logical(D1(:,Task_Local(LocalTaskIdx).STStgt).'),1)).*~STPtgt_LogIdx);
        STElemProperties_STP_TaskIdx(Task_Local(LocalTaskIdx).STPtgt) = GlobalTaskIdx;
        if TimeSection > 1
            TaskDepEdgeNum_Local = TaskDepEdgeNum_Local+1;
            StaTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx -1;
            TgtTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx;
        end
    end
    GlobalTaskOffset = GlobalTaskOffset + SpElemProperties.SpS.UpdNum(RepresentiveSpS);
    LocalTaskOffset  = LocalTaskOffset  + SpElemProperties.SpS.UpdNum(RepresentiveSpS);
end
end

function TMM_Onestep = TMMOnestep_STFIPrim_Block(Task,D2,SizeOfTMMRedundant)
STPNum = size(D2,2);
DoFIdxAssignedToSTP = 1:STPNum;
STPtgt_2DLogIdx = false(SizeOfTMMRedundant,1);
STPtgt_2DLogIdx(Task.STPtgt) = true;
STPtgt_2DLogIdx = spdiags(STPtgt_2DLogIdx,0,SizeOfTMMRedundant,SizeOfTMMRedundant);
TMM_Onestep = speye(SizeOfTMMRedundant,SizeOfTMMRedundant);

% assuming D2(Task.STHtgt(i),Task.STPtgt(i)) =
% 1 for all i = 1:size(Task.STPtgt,2)
TMM_Onestep(Task.STPtgt,DoFIdxAssignedToSTP) ...
    = -D2(Task.STHtgt,:);
TMM_Onestep(STPtgt_2DLogIdx) = 0;
end

function TMM_Onestep = TMMOnestep_STFIDual_Block(Task,D1,Kappa_over_Z,SizeOfTMMRedundant)
STPtgtNum = size(Task.STPtgt,2);
STPNum = size(D1,1);
DoFIdxAssignedToSTP = 1:STPNum;
STPtgt_2DLogIdx = false(SizeOfTMMRedundant,1);
STPtgt_2DLogIdx(Task.STPtgt) = true;
STPtgt_2DLogIdx = spdiags(STPtgt_2DLogIdx,0,SizeOfTMMRedundant,SizeOfTMMRedundant);

TMM_Onestep = speye(SizeOfTMMRedundant,SizeOfTMMRedundant);
% assuming D1(Task.STPtgt(i),Task.STStgt(i)) =
% 1 for all i = 1:size(Task.STPtgt,2)
TMM_Onestep(Task.STPtgt,DoFIdxAssignedToSTP) ...
    = -spdiags(Kappa_over_Z(Task.STPtgt).^(-1),0,STPtgtNum,STPtgtNum)*D1(:,Task.STStgt).'*spdiags(Kappa_over_Z,0,STPNum,STPNum);
TMM_Onestep(STPtgt_2DLogIdx) = 0;
end