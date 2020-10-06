function [Task,STElemProperties] = Generate_InitialCondition_Tasks(Task,Num_of_Elem,SpElemProperties,STElemProperties)
switch size(fieldnames(Task(size(Task,2))),1)
    case 0
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end
GlobalTaskIdx = GlobalTaskIdx+1;
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