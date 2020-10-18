function [TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Num_of_Elem,SpElemProperties)

Num_of_PMLDualDoF_Init = size(find(SpElemProperties.SpP.PML==true),2)+size(find(SpElemProperties.SpS.PML==true),2);


FieldDoF_Prim_Init_Idxs = 1:Num_of_Elem.SpP+Num_of_Elem.SpS;
FieldDoF_PMLDual_Init_Idxs    = ...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource+1:...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource+Num_of_PMLDualDoF_Init;

DoF_STSource_Idx = Num_of_Elem.SpP+Num_of_Elem.SpS+1:Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_Elem.STSource;


FieldDoF_Prim_Last_Idxs = 1:Num_of_Elem.SpP+Num_of_Elem.SpS;
FieldDoF_PMLDual_Last_Idxs    = ...
    Num_of_Elem.SpP+Num_of_Elem.SpS+1:...
    Num_of_Elem.SpP+Num_of_Elem.SpS+Num_of_PMLDualDoF_Init;

TMM_Fields  =[...
    TMM(FieldDoF_Prim_Last_Idxs      ,FieldDoF_Prim_Init_Idxs) TMM(FieldDoF_Prim_Last_Idxs      ,FieldDoF_PMLDual_Init_Idxs);...
    TMM(FieldDoF_PMLDual_Last_Idxs   ,FieldDoF_Prim_Init_Idxs) TMM(FieldDoF_PMLDual_Last_Idxs   ,FieldDoF_PMLDual_Init_Idxs)];

TMM_Sources = [...
    TMM(FieldDoF_Prim_Last_Idxs      ,DoF_STSource_Idx);...
    TMM(FieldDoF_PMLDual_Last_Idxs   ,DoF_STSource_Idx)];

end