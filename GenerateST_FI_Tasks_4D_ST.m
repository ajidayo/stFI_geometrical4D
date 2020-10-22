function [Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateST_FI_Tasks_4D_ST(SpElemProperties,STElemProperties,Task,TaskDepGraph)

switch size(fieldnames(Task(size(Task,2))),1)
    case 0
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end

% sDPattern_Modified = logical(sD);
% sDPattern_Modified(:,~logical(SpElemProperties.SpP.Belong_to_ST_FI)) = false;
% for SpVIdx = find(SpElemProperties.SpV.IncToSTFISpP)
%     for IncSpPIdx = find(sD(SpVIdx,:))
%         if SpElemProperties.SpV.UpdNum(SpVIdx) ~= SpElemProperties.SpP.UpdNum(IncSpPIdx)
%             sD(SpVIdx,IncSpPIdx) = false;
%         end
%     end
% end
% Graph_SpPviaSpV_SomeGraphEdgeExcluded = graph(sDPattern_Modified.'*sDPattern_Modified,'omitselfloops');
% [SubgraphBin_SpP,SpPNumInSubgraph] = conncomp(Graph_SpPviaSpV_SomeGraphEdgeExcluded);
% for SubgraphIdx = size(SpPNumInSubgraph,2)
%     SpPIdxListInSubgraph = find(SubgraphBin_SpP==SubgraphIdx);
%     
% end

for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    if SpElemProperties.SpS.PEC(SpSIdx)==true
        continue;
    end
    STP_Tgt             = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
    STOmegaTilde_Tgt    = SpElemProperties.SpS.FirstSTSIdx(SpSIdx)-1;
    for CurrentTimeSec = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
        GlobalTaskIdx = GlobalTaskIdx +1;
        Task(GlobalTaskIdx) = GenerEmptyTask;
        STP_Tgt = STP_Tgt +1;
        STOmegaTilde_Tgt = STOmegaTilde_Tgt+1;
        Task(GlobalTaskIdx).Type = "ST_FI_Dual";
        if CurrentTimeSec == 1
            Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)=GlobalTaskIdx;
        end
        Task(GlobalTaskIdx).STPtgt            = STP_Tgt;
        Task(GlobalTaskIdx).STStgt   = STOmegaTilde_Tgt;
        STElemProperties.STP.TaskIdx(STP_Tgt)  = GlobalTaskIdx;
    end
end
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    if SpElemProperties.SpP.ElecWall(SpPIdx)==true
        continue;
    end
    STP_Tgt       = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
    STOmega_Tgt   = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
    for CurrentTimeSec = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
        GlobalTaskIdx = GlobalTaskIdx +1;
        Task(GlobalTaskIdx) = GenerEmptyTask;
        STP_Tgt = STP_Tgt+1;
        STOmega_Tgt = STOmega_Tgt+1;
        Task(GlobalTaskIdx).Type = "ST_FI_Prim";
        if CurrentTimeSec == 1
            Map_SpElem_to_FirstGlobTask.SpPIdx(SpPIdx)=GlobalTaskIdx;
        end
        Task(GlobalTaskIdx).STPtgt     = STP_Tgt;
        Task(GlobalTaskIdx).STHtgt = STOmega_Tgt;
        STElemProperties.STP.TaskIdx(STP_Tgt)  = GlobalTaskIdx;
    end
end

EdgeIdx =0;
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    if SpElemProperties.SpS.PEC(SpSIdx)==true
        continue;
    end
    for CurrentTimeSec = 2:SpElemProperties.SpS.UpdNum(SpSIdx)
        EdgeIdx = EdgeIdx+1;
        StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)-1+CurrentTimeSec-1;
        TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)-1+CurrentTimeSec  ;
    end
end
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    if SpElemProperties.SpP.ElecWall(SpPIdx)==true
        continue;
    end
    for CurrentTimeSec = 2:SpElemProperties.SpP.UpdNum(SpPIdx)
        EdgeIdx = EdgeIdx+1;
        StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpPIdx(SpPIdx)-1+CurrentTimeSec-1;
        TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpPIdx(SpPIdx)-1+CurrentTimeSec  ;
    end
end

if exist('StaTask','var') == 1
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask


end