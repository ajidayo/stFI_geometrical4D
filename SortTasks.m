function TaskOrder = SortTasks(Task,TaskDepGraph,STElemProperties,D1,D2)

%% add edges to task dependency graph: Dependencies between SpFI tasks and STFI tasks
EdgeIdx =0;

for TgtTaskIdx = 1:size(Task,2)
    switch Task(TgtTaskIdx).Type
        case "InitialCondition"
        case "BoundaryCondition"
        case "ST_FI_Prim"
            STOmega_Tgt = Task(TgtTaskIdx).STOmega_Tgt;
            STP_Tgt = Task(TgtTaskIdx).STP_Tgt;
            for IncSTPIdx = find(D2(STOmega_Tgt,:))
                if IncSTPIdx == STP_Tgt
                    continue;
                end
                StaTaskIdx = STElemProperties.STP.TaskIdx(IncSTPIdx);
                EdgeIdx = EdgeIdx +1;
                StaTask(EdgeIdx) = StaTaskIdx;
                TgtTask(EdgeIdx) = TgtTaskIdx;
            end
        case "ST_FI_Dual"
            STOmegaTilde_Tgt = Task(TgtTaskIdx).STOmegaTilde_Tgt;
            STP_Tgt = Task(TgtTaskIdx).STP_Tgt;
            for IncSTPIdx = find(D1(:,STOmegaTilde_Tgt)).'
                if IncSTPIdx == STP_Tgt
                    continue;
                end
                StaTaskIdx = STElemProperties.STP.TaskIdx(IncSTPIdx);
                EdgeIdx = EdgeIdx +1;
                StaTask(EdgeIdx) = StaTaskIdx;
                TgtTask(EdgeIdx) = TgtTaskIdx;
            end
        case "Sp_FI"
            for DepSTPIdx = Task(TgtTaskIdx).DepSTPIdx
                StaTaskIdx = STElemProperties.STP.TaskIdx(DepSTPIdx);
                EdgeIdx = EdgeIdx +1;
                StaTask(EdgeIdx) = StaTaskIdx;
                TgtTask(EdgeIdx) = TgtTaskIdx;
            end
        case "PML"
            for DepSTPIdx = Task(TgtTaskIdx).DepSTPIdx
                StaTaskIdx = STElemProperties.STP.TaskIdx(DepSTPIdx);
                EdgeIdx = EdgeIdx +1;
                StaTask(EdgeIdx) = StaTaskIdx;
                TgtTask(EdgeIdx) = TgtTaskIdx;
            end
        otherwise
            disp(["Sorttasks:No such tasktype defined.(",Task(TgtTaskIdx).Type,")"])
    end
end

% for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
%     for TimeSec = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
%         STPTgt = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSec;
%         TgtGlobTaskIdx = STElemProperties.STP.Task_GlobIdx(STPTgt);
%         STOmegaTilde_Tgt = Task(TgtGlobTaskIdx).STOmegaTilde_Tgt;
%         for IncSTPIdx = find(D1(:,STOmegaTilde_Tgt)).'
%             if IncSTPIdx == STPTgt
%                 continue;
%             end
%             StaGlobTaskIdx = STElemProperties.STP.Task_GlobIdx(IncSTPIdx);
%             EdgeIdx = EdgeIdx +1;
%             StaTask(EdgeIdx) = StaGlobTaskIdx;
%             TgtTask(EdgeIdx) = TgtGlobTaskIdx;
%         end
%     end
% end
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%      for TimeSec = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
%         STOmegaTar      = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx) ...
%             +TimeSec;
%         for IncSTPIdx = find(D2(STOmegaTar,:))
%             
%         end
%      end
% end

% for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
%     for IncSpPIdx = find(sC(:,SpSIdx)).'
%         UpdRatio = SpElemProperties.SpS.UpdNum(SpSIdx)/SpElemProperties.SpP.UpdNum(IncSpPIdx);
%         for CurrentTimeSection = 2:SpElemProperties.SpS.UpdNum(SpSIdx)
%             EdgeIdx = EdgeIdx +1;
%             StaTask(EdgeIdx) = ...
%                 Map_SpElem_to_FirstGlobTask.SpPIdx(IncSpPIdx)-1 ...
%                 +CurrentTimeSection-1 +UpdRatio*(CurrentTimeSection-1);
%             TgtTask(EdgeIdx) = ...
%                 Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx   )-1 ...
%                 +CurrentTimeSection  ;
%         end
%     end
% end
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%     for IncSpSIdx = find(sC(SpPIdx,:))
%         UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
%         for CurrentTimeSection = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
%             EdgeIdx = EdgeIdx +1;
%             StaTask(EdgeIdx) = ...
%                 Map_SpElem_to_FirstGlobTask.SpSIdx(IncSpSIdx)-1 ...
%                 +CurrentTimeSection-1 +UpdRatio*(CurrentTimeSection-1);
%             TgtTask(EdgeIdx) = ...
%                 Map_SpElem_to_FirstGlobTask.SpPIdx(   SpPIdx)-1 ...
%                 +CurrentTimeSection  ;
%         end
%     end
% end

% for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
%     SpFI_RegionIdx = max(SpElemProperties.SpP.SpFI_TaskIdx(logical(sC(:,SpSIdx)).'));
%     if SpFI_RegionIdx > 0
%         for CurrentTimeSec = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
%             EdgeIdx = EdgeIdx+1;
%             StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)-1+CurrentTimeSec;
%             TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_RegionIdx)-1+CurrentTimeSec  ;
%         end
%         for CurrentTimeSec = 2:SpElemProperties.SpS.UpdNum(SpSIdx)
%             EdgeIdx = EdgeIdx+1;
%             StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_RegionIdx)-1+CurrentTimeSec-1;
%             TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)-1+CurrentTimeSec;
%         end
%     end
% end
if exist('StaTask','var') == 1
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask
if numedges(TaskDepGraph) == 0
    TaskOrder = 1;
else
    TaskOrder = toposort(TaskDepGraph);
end
end