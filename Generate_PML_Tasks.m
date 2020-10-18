function [Task,TaskDepGraph,SpElemProperties,STElemProperties] ...
    = Generate_PML_Tasks(sC,sD,SpElemProperties,STElemProperties,Num_of_Elem,Task,TaskDepGraph)

if size(find(SpElemProperties.SpP.PML),2)==0
    return
end

sDPattern_PartiallyOmitted=logical(sD);
for SpPIdx = find(SpElemProperties.SpP.PML==false)
    sDPattern_PartiallyOmitted(:,SpPIdx)=false;
end

Adj_SpP_PartiallyOmitted = graph(sDPattern_PartiallyOmitted.' * sDPattern_PartiallyOmitted,'omitselfloops');

[SG_bin_SpP,SG_ElemNum_SpP] = conncomp(Adj_SpP_PartiallyOmitted);
clearvars Adj_SpP_PartiallyOmitted

[PML_TaskInfo,SpElemProperties] = Eliminate_NonPML_SpPs(SG_bin_SpP,SG_ElemNum_SpP,SpElemProperties,Num_of_Elem);
clearvars SG_bin_SpP SG_ElemNum_SpP
for PMLRegionIdx =1:size(PML_TaskInfo,2)
    PML_TaskInfo(PMLRegionIdx).ElemIdx.SpS = 0;
    PML_TaskInfo(PMLRegionIdx).ElemNum.SpS = 0;
end
SpElemProperties.SpS.PML_TaskIdx=zeros(Num_of_Elem.SpS,1);
for SpS_tgt = find(SpElemProperties.SpS.PML)
    if SpElemProperties.SpS.PEC(SpS_tgt) == true
        continue;
    end
    for SpP_fetch = find(sC(:,SpS_tgt).')
        PMLRegionIdx = SpElemProperties.SpP.PML_TaskIdx(SpP_fetch);
        if  PMLRegionIdx>0
            break
        end
    end
    if PMLRegionIdx == 0
        disp(['Cannot fetch PML_TaskIdx for SpS_tgt = ',num2str(SpS_tgt)])
    end
    
    Num_of_IncludedSpS = PML_TaskInfo(PMLRegionIdx).ElemNum.SpS +1;
    PML_TaskInfo(PMLRegionIdx).ElemNum.SpS = Num_of_IncludedSpS;
    PML_TaskInfo(PMLRegionIdx).ElemIdx.SpS(Num_of_IncludedSpS) = SpS_tgt;
    SpElemProperties.SpS.PML_TaskIdx(SpS_tgt) = PMLRegionIdx;
end

for PMLRegionIdx = 1:size(PML_TaskInfo,2)
    DepSpSNum = 0;
    SpSCheckedFlag = false(1,Num_of_Elem.SpS);
    for SpPIdx = PML_TaskInfo(PMLRegionIdx).ElemIdx.SpP
        for IncSpSIdx = find(sC(SpPIdx,:))
            if SpSCheckedFlag(IncSpSIdx) == true
                continue;
            else
                SpSCheckedFlag(IncSpSIdx) = true;
                if SpElemProperties.SpS.PML_TaskIdx(IncSpSIdx) ~= PMLRegionIdx
                    DepSpSNum = DepSpSNum+1;
                    PML_TaskInfo(PMLRegionIdx).DepSpSIdx(DepSpSNum) = IncSpSIdx;
                end
            end
        end
    end
end
for PMLRegionIdx = 1:size(PML_TaskInfo,2)
    DepSpPNum = 0;
    SpPCheckedFlag = false(1,Num_of_Elem.SpP);
    for SpSIdx = PML_TaskInfo(PMLRegionIdx).ElemIdx.SpS
        for IncSpPIdx = find(sC(:,SpSIdx).')
            if SpPCheckedFlag(IncSpPIdx) == true
                continue;
            else
                SpPCheckedFlag(IncSpPIdx) = true;
                if SpElemProperties.SpP.PML_TaskIdx(IncSpPIdx) ~= PMLRegionIdx
                    DepSpPNum = DepSpPNum+1;
                    PML_TaskInfo(PMLRegionIdx).DepSpPIdx(DepSpPNum) = IncSpPIdx;
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
for PMLRegionIdx=1:size(PML_TaskInfo,2)
    PML_TaskInfo(PMLRegionIdx).UpdNum = ...
        SpElemProperties.SpP.UpdNum(PML_TaskInfo(PMLRegionIdx).ElemIdx.SpP(1));
    for CurrentTimeSec = 1:PML_TaskInfo(PMLRegionIdx).UpdNum
        GlobalTaskIdx = GlobalTaskIdx+1;
        if CurrentTimeSec == 1
            Map_SpElem_to_FirstGlobTask.PML_RegionIdx(PMLRegionIdx) = GlobalTaskIdx;
            for SpPIdx = PML_TaskInfo(PMLRegionIdx).ElemIdx.SpP
                Map_SpElem_to_FirstGlobTask.SpPIdx(SpPIdx)=GlobalTaskIdx;
            end
            for SpSIdx = PML_TaskInfo(PMLRegionIdx).ElemIdx.SpS
                Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)=GlobalTaskIdx;
            end
        end
        names = fieldnames(PML_TaskInfo(PMLRegionIdx));
        for i = 1:size(names,1)
            if names{i}~="DepSpPIdx" && names{i}~="DepSpSIdx"
                Task(GlobalTaskIdx).(names{i}) = PML_TaskInfo(PMLRegionIdx).(names{i});
            end
        end
        Task(GlobalTaskIdx).Type = "PML";
        Task(GlobalTaskIdx).TimeSection_tgt = CurrentTimeSec;
        DepSTPNum = 0;
        if isfield(PML_TaskInfo(PMLRegionIdx),'DepSpPIdx')
            for DepSpPIdx = PML_TaskInfo(PMLRegionIdx).DepSpPIdx
                DepSTPIdx = SpElemProperties.SpP.FirstSTPIdx(DepSpPIdx) -1+CurrentTimeSec;
                DepSTPNum = DepSTPNum +1;
                Task(GlobalTaskIdx).DepSTPIdx(DepSTPNum) = DepSTPIdx;
            end
        end
        if isfield(PML_TaskInfo(PMLRegionIdx),'DepSpSIdx')
            for DepSpSIdx = PML_TaskInfo(PMLRegionIdx).DepSpSIdx
                DepSTPIdx = SpElemProperties.SpS.FirstSTPIdx(DepSpSIdx) +CurrentTimeSec;
                DepSTPNum = DepSTPNum +1;
                Task(GlobalTaskIdx).DepSTPIdx(DepSTPNum) = DepSTPIdx;
            end
        end
        for SpPIdx = Task(GlobalTaskIdx).ElemIdx.SpP
            STP_Tgt = SpElemProperties.SpP.FirstSTPIdx(SpPIdx) +CurrentTimeSec;
            STElemProperties.STP.TaskIdx(STP_Tgt) = GlobalTaskIdx;
        end
        for SpSIdx = Task(GlobalTaskIdx).ElemIdx.SpS
            STP_Tgt = SpElemProperties.SpS.FirstSTPIdx(SpSIdx) +CurrentTimeSec;
            STElemProperties.STP.TaskIdx(STP_Tgt) = GlobalTaskIdx;
        end
    end
end

EdgeIdx =0;
for PMLRegionIdx = 1:size(PML_TaskInfo,2)
    for CurrentTimeSec = 2:PML_TaskInfo(PMLRegionIdx).UpdNum
        EdgeIdx = EdgeIdx+1;
        StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.PML_RegionIdx(PMLRegionIdx)-1+CurrentTimeSec-1;
        TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.PML_RegionIdx(PMLRegionIdx)-1+CurrentTimeSec;
    end
end
if exist('StaTask','var') == 1
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask

end

%% 
function [PML_TaskInfo,SpElemProperties] = Eliminate_NonPML_SpPs(SG_bin_SpP,SG_ElemNum_SpP,SpElemProperties,Num_of_Elem)
SpElemProperties.SpP.PML_TaskIdx=zeros(Num_of_Elem.SpP,1);
PML_TaskIdx=0;
for SGIdx=1:size(SG_ElemNum_SpP,2)
    IncludedSpPIdx=find(SG_bin_SpP==SGIdx);
    if SpElemProperties.SpP.PML(IncludedSpPIdx(1))==false 
        continue;
    end
    PML_TaskIdx=PML_TaskIdx+1;
    PML_TaskInfo(PML_TaskIdx).PML_Region_tgt = PML_TaskIdx;
    PML_TaskInfo(PML_TaskIdx).ElemIdx.SpP=IncludedSpPIdx;
    PML_TaskInfo(PML_TaskIdx).ElemNum.SpP=SG_ElemNum_SpP(SGIdx);
    SpElemProperties.SpP.PML_TaskIdx(IncludedSpPIdx)=PML_TaskIdx;
end
end