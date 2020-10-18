function [Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateTask_STFI_4DST_Block(SpElemProperties,STElemProperties,Task,TaskDepGraph,sC)

switch size(fieldnames(Task(size(Task,2))),1)
    case 0
        GlobalTaskOffset = 0;
    otherwise
        GlobalTaskOffset = size(Task,2);
end

StaTask = [];
TgtTask = [];
%%
[Task_Local,STElemProperties.STP.TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIPrim_Block(SpElemProperties,GlobalTaskOffset,sC,STElemProperties.STP.TaskIdx); 
Task = [Task Task_Local];
StaTask = [StaTask StaTask_Local];
TgtTask = [TgtTask TgtTask_Local];

%%
GlobalTaskOffset = size(Task,2);
[Task_Local,STElemProperties.STP.TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIDual_Block(SpElemProperties,GlobalTaskOffset,sC,STElemProperties_STP_TaskIdx);
Task = [Task Task_Local];
StaTask = [StaTask StaTask_Local];
TgtTask = [TgtTask TgtTask_Local];

if isempty(StaTask) == false
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask


end


function [Task_STFIPrim_Block_Local,STElemProperties_STP_TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIPrim_Block(SpElemProperties,GlobalTaskOffset,sC,STElemProperties_STP_TaskIdx) 
StaTask_Local = [];
TgtTask_Local = [];
TaskDepEdgeNum_Local = 0;

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
for SubgraphIdx = size(SpPNumInSubgraph,2)
    SpPIdxListInSubgraph = find(SubgraphBin_SpP==SubgraphIdx);
    RepresentiveSpP = SpPIdxListInSubgraph(1);
    if SpElemProperties.SpP.Belong_to_ST_FI(RepresentiveSpP) == false
        continue;
    elseif SpElemProperties.SpP.ElecWall(RepresentiveSpP) == true
        continue;
    end
    STPtgt_OffSet = SpElemProperties.SpP.FirstSTPIdx(SpPIdxListInSubgraph);
    STHtgt_OffSet = -1+SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdxListInSubgraph);
    LocalTaskOffset = 0;
    for TimeSection = 1:SpElemProperties.SpP.UpdNum(RepresentiveSpP)
        GlobalTaskIdx = GlobalTaskOffset + TimeSection;
        LocalTaskIdx  = LocalTaskOffset  + TimeSection;
        Task_STFIPrim_Block_Local(LocalTaskIdx).STPtgt = STPtgt_OffSet + TimeSection;
        Task_STFIPrim_Block_Local(LocalTaskIdx).STHtgt = STHtgt_OffSet + TimeSection;
        Task_STFIPrim_Block_Local(LocalTaskIdx).Type = "STFIPrim_Block";
        STElemProperties_STP_TaskIdx(Task_STFIPrim_Block_Local(LocalTaskIdx).STPtgt) ...
            = GlobalTaskIdx;
        if TimeSection >= 1
            TaskDepEdgeNum_Local = TaskDepEdgeNum_Local+1;
            StaTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx -1;
            TgtTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx;
        end
    end
    GlobalTaskOffset = GlobalTaskOffset + SpElemProperties.SpP.UpdNum(RepresentiveSpP);
end
end

function [Task_Local,STElemProperties_STP_TaskIdx,StaTask_Local,TgtTask_Local] ...
    = GenerateTask_STFIDual_Block(SpElemProperties,GlobalTaskOffset,sC,STElemProperties_STP_TaskIdx) 
StaTask_Local = [];
TgtTask_Local = [];
TaskDepEdgeNum_Local = 0;

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
for SubgraphIdx = size(SpSNumInSubgraph,2)
    SpSIdxListInSubgraph = find(SubgraphBin_SpS==SubgraphIdx);
    RepresentiveSpS = SpSIdxListInSubgraph(1);
    if SpElemProperties.SpS.Belong_to_ST_FI(RepresentiveSpS) == false
        continue;
    elseif SpElemProperties.SpS.PEC(RepresentiveSpS) == true
        continue;
    end
    STPtgt_OffSet = SpElemProperties.SpS.FirstSTPIdx(SpSIdxListInSubgraph);
    STHtildetgt_OffSet = -1+SpElemProperties.SpS.FirstSTSIdx(SpSIdxListInSubgraph);
    LocalTaskOffset = 0;
    for TimeSection = 1:SpElemProperties.SpS.UpdNum(RepresentiveSpS)
        GlobalTaskIdx = GlobalTaskOffset + Timesection;
        LocalTaskIdx  = LocalTaskOffset + TimeSection;
        Task_Local(LocalTaskIdx).STPtgt = STPtgt_OffSet + TimeSection;
        Task_Local(LocalTaskIdx).STHtildetgt = STHtildetgt_OffSet + TimeSection;
        Task_Local(LocalTaskIdx).Type = "STFIDual_Block";
        STElemProperties_STP_TaskIdx(Task_Local(LocalTaskIdx).STPtgt) ...
            = GlobalTaskIdx;
        if TimeSection >= 1
            TaskDepEdgeNum_Local = TaskDepEdgeNum_Local+1;
            StaTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx -1;
            TgtTask_Local(TaskDepEdgeNum_Local) = GlobalTaskIdx;
        end
    end
    GlobalTaskOffset = GlobalTaskOffset + SpElemProperties.SpS.UpdNum(RepresentiveSpS);
end

end