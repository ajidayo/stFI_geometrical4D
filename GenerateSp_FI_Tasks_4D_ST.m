function [Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph)

TaskFieldNames = fieldnames(Task(size(Task,2)));
switch isempty(Task(size(Task,2)).(TaskFieldNames{1}))
    case 1
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end

sDPattern_PartiallyOmitted=logical(sD);


for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    sDPattern_PartiallyOmitted(:,SpPIdx)=false;
end
for SpPIdx = find(SpElemProperties.SpP.PML)
    sDPattern_PartiallyOmitted(:,SpPIdx)=false;
end
NoSpFIFlag = true;
for SpPIdx = 1:Num_of_Elem.SpP
    if ~SpElemProperties.SpP.Belong_to_ST_FI(SpPIdx) && ~SpElemProperties.SpP.PML(SpPIdx) 
        NoSpFIFlag = false;
    end
end

if NoSpFIFlag == true
    return
end

Adj_SpP_PartiallyOmitted = graph(sDPattern_PartiallyOmitted.' * sDPattern_PartiallyOmitted,'omitselfloops');


[SG_bin_SpP,SG_ElemNum_SpP] = conncomp(Adj_SpP_PartiallyOmitted);
clearvars Adj_SpP_PartiallyOmitted

[SpFI_TaskInfo,SpElemProperties] ...
     = Eliminate_STFIandPMLSpPs(SG_bin_SpP,SG_ElemNum_SpP,SpElemProperties,Num_of_Elem);
clearvars SG_bin_SpP SG_ElemNum_SpP
for SpFI_TaskIdx =1:size(SpFI_TaskInfo,2)
    SpFI_TaskInfo(SpFI_TaskIdx).SpStar = 0;
    SpFI_TaskInfo(SpFI_TaskIdx).SpStarNum = 0;
end
SpElemProperties.SpS.SpFI_TaskIdx=zeros(Num_of_Elem.SpS,1);
for SpS_tgt=1:Num_of_Elem.SpS
    if SpElemProperties.SpS.Belong_to_ST_FI(SpS_tgt) || SpElemProperties.SpS.PEC(SpS_tgt) ...
            || SpElemProperties.SpS.PML(SpS_tgt)
        continue;
    end
    IncSpSIdx = find(sC(:,SpS_tgt)).';
    for Counter = 1:size(IncSpSIdx,2)
        SpP_fetch = IncSpSIdx(Counter);
        SpFI_TaskIdx = SpElemProperties.SpP.SpFI_TaskIdx(SpP_fetch);
        if  SpFI_TaskIdx>0
            break
        end
    end
    if SpFI_TaskIdx == 0
        disp(['SpFI_TaskIdx == 0, SpS_tgt = ',num2str(SpS_tgt)])
    end
    
    Num_of_IncludedSpS = SpFI_TaskInfo(SpFI_TaskIdx).SpStarNum +1;
    SpFI_TaskInfo(SpFI_TaskIdx).SpStarNum = Num_of_IncludedSpS;
    SpFI_TaskInfo(SpFI_TaskIdx).SpStar(Num_of_IncludedSpS) = SpS_tgt;
    SpElemProperties.SpS.SpFI_TaskIdx(SpS_tgt)=SpFI_TaskIdx;
end

for SpFI_TaskIdx = 1:size(SpFI_TaskInfo,2)
    DepSpSNum = 0;
    SpSCheckedFlag = false(1,Num_of_Elem.SpS);
    for SpPIdx = SpFI_TaskInfo(SpFI_TaskIdx).SpPtar
        for IncSpSIdx = find(sC(SpPIdx,:))
            if SpSCheckedFlag(IncSpSIdx) == true
                continue;
            else
                SpSCheckedFlag(IncSpSIdx) = true;
                if SpElemProperties.SpS.SpFI_TaskIdx(IncSpSIdx)~=SpFI_TaskIdx
                    DepSpSNum = DepSpSNum+1;
                    SpFI_TaskInfo(SpFI_TaskIdx).DepSpSIdx(DepSpSNum) ...
                        = IncSpSIdx;
                end
            end
        end
    end
end
for SpFI_TaskIdx = 1:size(SpFI_TaskInfo,2)
    DepSpPNum = 0;
    SpPCheckedFlag = false(1,Num_of_Elem.SpP);
    for SpSIdx = SpFI_TaskInfo(SpFI_TaskIdx).SpStar
        for IncSpPIdx = find(sC(:,SpSIdx)).'
            if SpPCheckedFlag(IncSpPIdx) == true
                continue;
            else
                SpPCheckedFlag(IncSpPIdx) = true;
                if SpElemProperties.SpP.SpFI_TaskIdx(IncSpPIdx)~=SpFI_TaskIdx
                    DepSpPNum = DepSpPNum+1;
                    SpFI_TaskInfo(SpFI_TaskIdx).DepSpPIdx(DepSpPNum) ...
                        = IncSpPIdx;
                end
            end
        end
    end
end


switch size(fieldnames(Task(size(Task,2))),1)
    case 0
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end
for SpFI_TaskIdx=1:size(SpFI_TaskInfo,2)
    SpFI_TaskInfo(SpFI_TaskIdx).UpdNum = ...
        SpElemProperties.SpP.UpdNum(SpFI_TaskInfo(SpFI_TaskIdx).SpPtar(1));
    for CurrentTimeSec = 1:SpFI_TaskInfo(SpFI_TaskIdx).UpdNum
        GlobalTaskIdx = GlobalTaskIdx+1;
        Task(GlobalTaskIdx) = GenerEmptyTask;
        if CurrentTimeSec == 1
            Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_TaskIdx) = GlobalTaskIdx;
            for SpPIdx = SpFI_TaskInfo(SpFI_TaskIdx).SpPtar
                Map_SpElem_to_FirstGlobTask.SpPIdx(SpPIdx)=GlobalTaskIdx;
            end
            for SpSIdx = SpFI_TaskInfo(SpFI_TaskIdx).SpStar
                Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)=GlobalTaskIdx;
            end
        end
        names = fieldnames(SpFI_TaskInfo(SpFI_TaskIdx));
        for i = 1:size(names,1)
            if names{i}~="DepSpPIdx" && names{i}~="DepSpSIdx"
                Task(GlobalTaskIdx).(names{i}) = SpFI_TaskInfo(SpFI_TaskIdx).(names{i});
            end
        end
        Task(GlobalTaskIdx).Type = "Sp_FI";
        Task(GlobalTaskIdx).TimeSection_tgt = CurrentTimeSec;
        DepSTPNum = 0;
        if isfield(SpFI_TaskInfo(SpFI_TaskIdx),'DepSpPIdx')
            for DepSpPIdx = SpFI_TaskInfo(SpFI_TaskIdx).DepSpPIdx
                DepSTPIdx = SpElemProperties.SpP.FirstSTPIdx(DepSpPIdx)-1 ...
                    +CurrentTimeSec;
                DepSTPNum = DepSTPNum +1;
                Task(GlobalTaskIdx).DepSTP(DepSTPNum) = DepSTPIdx;
            end
        end
        if isfield(SpFI_TaskInfo(SpFI_TaskIdx),'DepSpSIdx')
            for DepSpSIdx = SpFI_TaskInfo(SpFI_TaskIdx).DepSpSIdx
                DepSTPIdx = SpElemProperties.SpS.FirstSTPIdx(DepSpSIdx) ...
                    +CurrentTimeSec;
                DepSTPNum = DepSTPNum +1;
                Task(GlobalTaskIdx).DepSTP(DepSTPNum) = DepSTPIdx;
            end
        end
        for SpPIdx = Task(GlobalTaskIdx).SpPtar
            STP_Tgt = SpElemProperties.SpP.FirstSTPIdx(SpPIdx) ...
                +CurrentTimeSec;
            STElemProperties.STP.TaskIdx(STP_Tgt) = GlobalTaskIdx;
        end
        for SpSIdx = Task(GlobalTaskIdx).SpStar
            STP_Tgt = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)...
                +CurrentTimeSec;
            STElemProperties.STP.TaskIdx(STP_Tgt) = GlobalTaskIdx;
        end
    end
end

EdgeIdx =0;
for SpFI_TaskIdx = 1:size(SpFI_TaskInfo,2)
    for CurrentTimeSec = 2:SpFI_TaskInfo(SpFI_TaskIdx).UpdNum
        EdgeIdx = EdgeIdx+1;
        StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_TaskIdx)-1+CurrentTimeSec-1;
        TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_TaskIdx)-1+CurrentTimeSec;
    end
end
if exist('StaTask','var') == 1
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask

end

%% 
function [SpFI_TaskInfo,SpElemProperties] = Eliminate_STFIandPMLSpPs(SG_bin_SpP,SG_ElemNum_SpP,SpElemProperties,Num_of_Elem)
SpElemProperties.SpP.SpFI_TaskIdx=zeros(Num_of_Elem.SpP,1);
SpFI_TaskIdx=0;
for SGIdx=1:size(SG_ElemNum_SpP,2)
    LogiIdx=find(SG_bin_SpP==SGIdx);
    SpP_test=LogiIdx(1);
    if SpElemProperties.SpP.Belong_to_ST_FI(SpP_test) || SpElemProperties.SpP.PML(SpP_test)
    else
        SpFI_TaskIdx=SpFI_TaskIdx+1;
       % SpFI_TaskInfo(SpFI_TaskIdx).Sp_FI_Region_tgt = SpFI_TaskIdx;
        SpFI_TaskInfo(SpFI_TaskIdx).SpPtar = LogiIdx;
        SpFI_TaskInfo(SpFI_TaskIdx).SpPtarNum = SG_ElemNum_SpP(SGIdx);
        SpElemProperties.SpP.SpFI_TaskIdx(LogiIdx)=SpFI_TaskIdx;
    end
end
end