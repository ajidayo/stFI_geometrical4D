function [Task,STElemProperties] = Generate_InitialCondition_Tasks(Task,Num_of_Elem,SpElemProperties,STElemProperties)
TaskFieldNames = fieldnames(Task(size(Task,2)));
switch isempty(Task(size(Task,2)).(TaskFieldNames{1}))
    case 1
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end
GlobalTaskIdx = GlobalTaskIdx+1;
Task(GlobalTaskIdx) = GenerEmptyTask;
Task(GlobalTaskIdx).Type = "InitialCondition";

for SpSIdx = 1:Num_of_Elem.SpS
    STPTgt = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
    STElemProperties.STP.TaskIdx(STPTgt) = GlobalTaskIdx;
end
for SpPIdx = 1:Num_of_Elem.SpP
    STPTgt = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
    STElemProperties.STP.TaskIdx(STPTgt) = GlobalTaskIdx;
end
end