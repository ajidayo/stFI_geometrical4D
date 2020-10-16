function TaskOrder = SortTasks(Task,TaskDepGraph,STElemProperties,D1,D2)
disp('SortTasks: Sorting DoF-calculating tasks in a viable order')
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
                if StaTask(EdgeIdx) == 0
                    disp('StaTask == 0')
                end
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
                if StaTask(EdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        case "Sp_FI"
            for DepSTPIdx = Task(TgtTaskIdx).DepSTPIdx
                StaTaskIdx = STElemProperties.STP.TaskIdx(DepSTPIdx);
                EdgeIdx = EdgeIdx +1;
                StaTask(EdgeIdx) = StaTaskIdx;
                TgtTask(EdgeIdx) = TgtTaskIdx;
                if StaTask(EdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        case "PML"
            for DepSTPIdx = Task(TgtTaskIdx).DepSTPIdx
                StaTaskIdx = STElemProperties.STP.TaskIdx(DepSTPIdx);
                EdgeIdx = EdgeIdx +1;
                StaTask(EdgeIdx) = StaTaskIdx;
                TgtTask(EdgeIdx) = TgtTaskIdx;
                if StaTask(EdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        otherwise
            disp(["Sorttasks:No such tasktype defined.(",Task(TgtTaskIdx).Type,")"])
    end
end
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