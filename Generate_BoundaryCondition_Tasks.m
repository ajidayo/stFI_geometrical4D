function [Task,STElemProperties] = Generate_BoundaryCondition_Tasks(Task,SpElemProperties,STElemProperties)
TaskFieldNames = fieldnames(Task(size(Task,2)));
switch isempty(Task(size(Task,2)).(TaskFieldNames{1}))
    case 1
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end
GlobalTaskIdx = GlobalTaskIdx+1;
Task(GlobalTaskIdx) = GenerEmptyTask;
Task(GlobalTaskIdx).Type = "BoundaryCondition";

for PECSpSIdx = find(SpElemProperties.SpS.PEC)
    for TimeSec = 0:SpElemProperties.SpS.UpdNum(PECSpSIdx)
        STPTgt = SpElemProperties.SpS.FirstSTPIdx(PECSpSIdx)+TimeSec;
        STElemProperties.STP.TaskIdx(STPTgt) = GlobalTaskIdx;
    end
end
for EWLSpPIdx = find(SpElemProperties.SpP.ElecWall)
    for TimeSec = 0:SpElemProperties.SpP.UpdNum(EWLSpPIdx)
        STPTgt = SpElemProperties.SpP.FirstSTPIdx(EWLSpPIdx)+TimeSec;
        STElemProperties.STP.TaskIdx(STPTgt) = GlobalTaskIdx;
    end
end
end