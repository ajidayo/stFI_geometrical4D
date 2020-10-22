function TaskOrder = SortTasks(Task,TaskDepGraph,STElemProperties,D1,D2)
disp('SortTasks: Sorting DoF-calculating tasks in a viable order')
%% add edges to task dependency graph: Dependencies between SpFI tasks and STFI tasks
TaskDepEdgeIdx =0;
TaskDepEdgeNum =0;

for TaskIdx = 1:size(Task,2)
    switch Task(TaskIdx).Type
        case "InitialCondition"
        case "BoundaryCondition"
        case "STFIPrim_Block"
            StaTaskList = unique(STElemProperties.STP.TaskIdx(Task(TaskIdx).DepSTP));
            TaskDepEdgeIdx = TaskDepEdgeNum + (1:size(StaTaskList,2));
            TaskDepEdgeNum = TaskDepEdgeNum + size(StaTaskList,2);
            StaTask(TaskDepEdgeIdx) =  StaTaskList;
            TgtTask(TaskDepEdgeIdx) = TaskIdx;
            if StaTask(TaskDepEdgeIdx) == 0
                disp('StaTask == 0')
            end
        case "STFIDual_Block"
            StaTaskList = unique(STElemProperties.STP.TaskIdx(Task(TaskIdx).DepSTP));
            TaskDepEdgeIdx = TaskDepEdgeNum + (1:size(StaTaskList,2));
            TaskDepEdgeNum = TaskDepEdgeNum + size(StaTaskList,2);
            StaTask(TaskDepEdgeIdx) =  StaTaskList;
            TgtTask(TaskDepEdgeIdx) = TaskIdx;
            if StaTask(TaskDepEdgeIdx) == 0
                disp('StaTask == 0')
            end
        case "ST_FI_Prim"
            STOmega_Tgt = Task(TaskIdx).STHtgt;
            STP_Tgt = Task(TaskIdx).STPtgt;
            for IncSTPIdx = find(D2(STOmega_Tgt,:))
                if IncSTPIdx == STP_Tgt
                    continue;
                end
                StaTaskIdx = STElemProperties.STP.TaskIdx(IncSTPIdx);
                TaskDepEdgeIdx = TaskDepEdgeIdx +1;
                TaskDepEdgeNum = TaskDepEdgeNum +1;
                StaTask(TaskDepEdgeIdx) = StaTaskIdx;
                TgtTask(TaskDepEdgeIdx) = TaskIdx;
                if StaTask(TaskDepEdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        case "ST_FI_Dual"
            STOmegaTilde_Tgt = Task(TaskIdx).STStgt;
            STP_Tgt = Task(TaskIdx).STPtgt;
            for IncSTPIdx = find(D1(:,STOmegaTilde_Tgt)).'
                if IncSTPIdx == STP_Tgt
                    continue;
                end
                StaTaskIdx = STElemProperties.STP.TaskIdx(IncSTPIdx);
                TaskDepEdgeIdx = TaskDepEdgeIdx +1;
                TaskDepEdgeNum = TaskDepEdgeNum +1;
                StaTask(TaskDepEdgeIdx) = StaTaskIdx;
                TgtTask(TaskDepEdgeIdx) = TaskIdx;
                if StaTask(TaskDepEdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        case "Sp_FI"
            for DepSTPIdx = Task(TaskIdx).DepSTP
                StaTaskIdx = STElemProperties.STP.TaskIdx(DepSTPIdx);
                TaskDepEdgeIdx = TaskDepEdgeIdx +1;
                TaskDepEdgeNum = TaskDepEdgeNum +1;
                StaTask(TaskDepEdgeIdx) = StaTaskIdx;
                TgtTask(TaskDepEdgeIdx) = TaskIdx;
                if StaTask(TaskDepEdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        case "PML"
            for DepSTPIdx = Task(TaskIdx).DepSTP
                StaTaskIdx = STElemProperties.STP.TaskIdx(DepSTPIdx);
                TaskDepEdgeIdx = TaskDepEdgeIdx +1;
                TaskDepEdgeNum = TaskDepEdgeNum +1;
                StaTask(TaskDepEdgeIdx) = StaTaskIdx;
                TgtTask(TaskDepEdgeIdx) = TaskIdx;
                if StaTask(TaskDepEdgeIdx) == 0
                    disp('StaTask == 0')
                end
            end
        otherwise
            disp(["Sorttasks:No such tasktype defined.(",Task(TaskIdx).Type,")"])
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