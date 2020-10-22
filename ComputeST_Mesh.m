function [D0,D1,D2,D3] = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem)
disp('ComputeST_Mesh: Constructing ST-Mesh')
Num_of_Elem.STV     = sum(SpElemProperties.SpV.UpdNum) + Num_of_Elem.SpV;
Num_of_Elem.STOmega = sum(SpElemProperties.SpV.UpdNum) + Num_of_Elem.SpV + sum(SpElemProperties.SpP.UpdNum) + Num_of_Elem.SpP;
Num_of_Elem.STS     = sum(SpElemProperties.SpS.UpdNum) + Num_of_Elem.SpS + sum(SpElemProperties.SpN.UpdNum) + Num_of_Elem.SpN;
Num_of_Elem.STN     = sum(SpElemProperties.SpN.UpdNum) + Num_of_Elem.SpN;
D0 = sparse(Num_of_Elem.STS     ,Num_of_Elem.STN    );
D1 = sparse(Num_of_Elem.STP     ,Num_of_Elem.STS    );
D2 = sparse(Num_of_Elem.STOmega ,Num_of_Elem.STP    );
D3 = sparse(Num_of_Elem.STV     ,Num_of_Elem.STOmega);

if size(find(SpElemProperties.SpS.Belong_to_ST_FI),2)==0
    return
end
D1_Nonzero_Num = 0;
D1_Nonzero_Row = [];
D1_Nonzero_Col = [];
D1_Nonzero_Val = [];
D2_Nonzero_Num = 0;
D2_Nonzero_Row = [];
D2_Nonzero_Col = [];
D2_Nonzero_Val = [];
D3_Nonzero_Num = 0;
D3_Nonzero_Row = [];
D3_Nonzero_Col = [];
D3_Nonzero_Val = [];
                    
% %% D0
% for SpVIdx = find(SpElemProperties.SpV.Belong_to_ST_FI) 
%     for SpPIdx = find(sD(SpVIdx,:))
%         for IncSpSIdx = find(sC(SpPIdx,:))
%             for IncSpNIdx = find(sG(IncSpSIdx,:))
%                 UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(IncSpSIdx);
%                 for CurrentTime = 0:SpElemProperties.SpS.UpdNum(IncSpSIdx)
%                     STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)   +CurrentTime;
%                     STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+UpdRatio*CurrentTime;
%                     D0(STS_tar,STN_tar) = sG(IncSpSIdx,IncSpNIdx);
%                 end
%                 for TimeSect = 0
%                     STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx);
%                     STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx);
%                     D0(STS_tar,STN_tar) =  1;
%                 end
%                 for TimeSect = 1:SpElemProperties.SpN.UpdNum(IncSpNIdx)
%                     STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx)+TimeSect;
%                     STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+TimeSect-1;
%                     D0(STS_tar,STN_tar) = -1;
%                     STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+TimeSect  ;
%                     D0(STS_tar,STN_tar) =  1;
%                 end
%             end
%         end
%     end
% end

STFISpVLogIdxs = logical(SpElemProperties.SpV.Belong_to_ST_FI);
STFISpPIdxs = logical(sum(logical(sD(STFISpVLogIdxs,:)),1));
STFISpSIdxs = find(sum(logical(sC(STFISpPIdxs,:)),1));
STFISpNIdxs = find(sum(logical(sG(STFISpSIdxs,:)),1));
D0_Nonzero_Num = 0;
D0_Nonzero_Row = [];
D0_Nonzero_Col = [];
D0_Nonzero_Val = [];

for IncSpSIdx = STFISpSIdxs
    for IncSpNIdx = find(sG(IncSpSIdx,:))
        UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(IncSpSIdx);
        for CurrentTime = 0:SpElemProperties.SpS.UpdNum(IncSpSIdx)
            STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)   +CurrentTime;
            STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+UpdRatio*CurrentTime;
            D0(STS_tar,STN_tar) = sG(IncSpSIdx,IncSpNIdx);
            D0_Nonzero_Num = D0_Nonzero_Num + 1;
            D0_Nonzero_Row(D0_Nonzero_Num) = STS_tar;
            D0_Nonzero_Col(D0_Nonzero_Num) = STN_tar;
            D0_Nonzero_Val(D0_Nonzero_Num) = sG(IncSpSIdx,IncSpNIdx);
        end
    end
end
for SpNIdx = STFISpNIdxs
    for TimeSect = 0
        STS_tar = SpElemProperties.SpN.FirstSTSIdx(SpNIdx);
        STN_tar = SpElemProperties.SpN.FirstSTNIdx(SpNIdx);
        D0(STS_tar,STN_tar) =  1;
        D0_Nonzero_Num = D0_Nonzero_Num + 1;
        D0_Nonzero_Row(D0_Nonzero_Num) = STS_tar;
        D0_Nonzero_Col(D0_Nonzero_Num) = STN_tar;
        D0_Nonzero_Val(D0_Nonzero_Num) = 1;
    end
    for TimeSect = 1:SpElemProperties.SpN.UpdNum(SpNIdx)
        STS_tar = SpElemProperties.SpN.FirstSTSIdx(SpNIdx)+TimeSect;
        STN_tar = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+TimeSect-1;
        D0(STS_tar,STN_tar) = -1;
        D0_Nonzero_Num = D0_Nonzero_Num + 1;
        D0_Nonzero_Row(D0_Nonzero_Num) = STS_tar;
        D0_Nonzero_Col(D0_Nonzero_Num) = STN_tar;
        D0_Nonzero_Val(D0_Nonzero_Num) = -1;
        STN_tar = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+TimeSect  ;
        D0(STS_tar,STN_tar) =  1;
        D0_Nonzero_Num = D0_Nonzero_Num + 1;
        D0_Nonzero_Row(D0_Nonzero_Num) = STS_tar;
        D0_Nonzero_Col(D0_Nonzero_Num) = STN_tar;
        D0_Nonzero_Val(D0_Nonzero_Num) = 1;
    end
end
%% D1
% for SpVIdx = find(SpElemProperties.SpV.Belong_to_ST_FI)
%     for SpPIdx = find(sD(SpVIdx,:))
%         for IncSpSIdx = find(sC(SpPIdx,:))
%             UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
%             for CurrentTime = 0:SpElemProperties.SpP.UpdNum(SpPIdx)
%                 STP_tar = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+CurrentTime;
%                 STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)+UpdRatio*CurrentTime;
%                 D1(STP_tar,STS_tar) = sC(SpPIdx,IncSpSIdx);
%             end
%             for TimeSect = 0
%                 STP_tar = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx);
%                 STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx);
%                 D1(STP_tar,STS_tar) = -1;
%             end
%             for TimeSect = 1:SpElemProperties.SpS.UpdNum(IncSpSIdx)
%                 STP_tar = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx)+TimeSect;
%                 STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)+TimeSect-1;
%                 D1(STP_tar,STS_tar) =  1;
%                 STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)+TimeSect  ;
%                 D1(STP_tar,STS_tar) = -1;
%             end
%             for IncSpNIdx = find(sG(IncSpSIdx,:))
%                 UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(IncSpSIdx);
%                 for TimeSect = 0
%                     STP_tar = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx);
%                     STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx);
%                     D1(STP_tar,STS_tar) =  sG(IncSpSIdx,IncSpNIdx);
%                 end
%                 for TimeSect = 1:SpElemProperties.SpS.UpdNum(IncSpSIdx)
%                     for LocalTimeSec_for_SpN = 1:UpdRatio
%                         STP_tar = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx)+TimeSect;
%                         STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx)+LocalTimeSec_for_SpN+UpdRatio*(TimeSect-1);
%                         D1(STP_tar,STS_tar) =  sG(IncSpSIdx,IncSpNIdx);
%                     end
%                 end
%             end
%         end
%     end
% end


STFISpVLogIdxs = logical(SpElemProperties.SpV.Belong_to_ST_FI);
STFISpPIdxs = find(sum(logical(sD(STFISpVLogIdxs,:)),1));
STFISpSIdxs = find(sum(logical(sC(STFISpPIdxs,:)),1));
%STFISpNIdxs = find(sum(logical(sG(STFISpSIdxs,:)),1));
D1_Nonzero_Num = 0;
D1_Nonzero_Row = [];
D1_Nonzero_Col = [];
D1_Nonzero_Val = [];
for SpPIdx = STFISpPIdxs
    for IncSpSIdx = find(sC(SpPIdx,:))
        UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
        for CurrentTime = 0:SpElemProperties.SpP.UpdNum(SpPIdx)
            STP_tar = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+CurrentTime;
            STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)+UpdRatio*CurrentTime;
            D1_Nonzero_Num = D1_Nonzero_Num+1;
            D1_Nonzero_Row(D1_Nonzero_Num) = STP_tar;
            D1_Nonzero_Col(D1_Nonzero_Num) = STS_tar;
            D1_Nonzero_Val(D1_Nonzero_Num) = sC(SpPIdx,IncSpSIdx);
        end
    end
end
for SpSIdx = STFISpSIdxs
    for TimeSect = 0
        STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
        STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx);
        D1_Nonzero_Num = D1_Nonzero_Num+1;
        D1_Nonzero_Row(D1_Nonzero_Num) = STP_tar;
        D1_Nonzero_Col(D1_Nonzero_Num) = STS_tar;
        D1_Nonzero_Val(D1_Nonzero_Num) = -1;
    end
    for TimeSect = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
        STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSect;
        STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx)+TimeSect-1;
        D1_Nonzero_Num = D1_Nonzero_Num+1;
        D1_Nonzero_Row(D1_Nonzero_Num) = STP_tar;
        D1_Nonzero_Col(D1_Nonzero_Num) = STS_tar;
        D1_Nonzero_Val(D1_Nonzero_Num) = 1;
        STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx)+TimeSect  ;
        D1_Nonzero_Num = D1_Nonzero_Num+1;
        D1_Nonzero_Row(D1_Nonzero_Num) = STP_tar;
        D1_Nonzero_Col(D1_Nonzero_Num) = STS_tar;
        D1_Nonzero_Val(D1_Nonzero_Num) = -1;
    end
    for IncSpNIdx = find(sG(SpSIdx,:))
        UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(SpSIdx);
        for TimeSect = 0
            STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
            STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx);
            D1_Nonzero_Num = D1_Nonzero_Num+1;
            D1_Nonzero_Row(D1_Nonzero_Num) = STP_tar;
            D1_Nonzero_Col(D1_Nonzero_Num) = STS_tar;
            D1_Nonzero_Val(D1_Nonzero_Num) = sG(SpSIdx,IncSpNIdx);
        end
        for TimeSect = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
            for LocalTimeSec_for_SpN = 1:UpdRatio
                STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSect;
                STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx)+LocalTimeSec_for_SpN+UpdRatio*(TimeSect-1);
                D1_Nonzero_Num = D1_Nonzero_Num+1;
                D1_Nonzero_Row(D1_Nonzero_Num) = STP_tar;
                D1_Nonzero_Col(D1_Nonzero_Num) = STS_tar;
                D1_Nonzero_Val(D1_Nonzero_Num) = sG(SpSIdx,IncSpNIdx);
            end
        end
    end
end


%% D2
% for SpVIdx = find(SpElemProperties.SpV.Belong_to_ST_FI)
%     for SpPIdx = find(sD(SpVIdx,:))
%         UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(SpVIdx);
%         for CurrentTime = 0:SpElemProperties.SpV.UpdNum(SpVIdx)
%             STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx)+CurrentTime;
%             STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+UpdRatio*CurrentTime;
%             D2(STH_tar,STP_tar) = sD(SpVIdx,SpPIdx);
%         end
%         for TimeSect = 0
%             STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
%             STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
%             D2(STH_tar,STP_tar) =  1;
%         end
%         for TimeSect = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
%             STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+TimeSect;
%             STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSect-1;
%             D2(STH_tar,STP_tar) = -1;
%             STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSect;
%             D2(STH_tar,STP_tar) =  1;
%         end
%         for IncSpSIdx = find(sC(SpPIdx,:))
%             UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
%             for TimeSect = 0
%                 STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
%                 STP_tar     = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx);
%                 D2(STH_tar,STP_tar) =  sC(SpPIdx,IncSpSIdx);
%             end
%             for TimeSect = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
%                 for LocalTimeSectForSpS = 1:UpdRatio
%                     STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+TimeSect;
%                     STP_tar     = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx)+LocalTimeSectForSpS+UpdRatio*(TimeSect-1);
%                     D2(STH_tar,STP_tar) = sC(SpPIdx,IncSpSIdx);
%                 end
%             end
%         end
%     end
% end

STFISpVIdxs = find(SpElemProperties.SpV.Belong_to_ST_FI);
STFISpPIdxs = find(sum(logical(sD(STFISpVIdxs,:)),1));
%STFISpSIdxs = find(sum(logical(sC(STFISpPIdxs,:)),1));
%STFISpNIdxs = find(sum(logical(sG(STFISpSIdxs,:)),1));
D2_Nonzero_Num = 0;
D2_Nonzero_Row = [];
D2_Nonzero_Col = [];
D2_Nonzero_Val = [];

for SpVIdx = STFISpVIdxs
    for IncSpPIdx = find(sD(SpVIdx,:))
        UpdRatio = SpElemProperties.SpP.UpdNum(IncSpPIdx)/SpElemProperties.SpV.UpdNum(SpVIdx);
        for CurrentTime = 0:SpElemProperties.SpV.UpdNum(SpVIdx)
            STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx)+CurrentTime;
            STP_tar     = SpElemProperties.SpP.FirstSTPIdx(IncSpPIdx)+UpdRatio*CurrentTime;
            D2_Nonzero_Num = D2_Nonzero_Num +1;
            D2_Nonzero_Row(D2_Nonzero_Num) = STH_tar;
            D2_Nonzero_Col(D2_Nonzero_Num) = STP_tar;
            D2_Nonzero_Val(D2_Nonzero_Num) = sD(SpVIdx,IncSpPIdx);
        end
    end
end
for SpPIdx = STFISpPIdxs
    for TimeSect = 0
        STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
        STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
        D2_Nonzero_Num = D2_Nonzero_Num +1;
        D2_Nonzero_Row(D2_Nonzero_Num) = STH_tar;
        D2_Nonzero_Col(D2_Nonzero_Num) = STP_tar;
        D2_Nonzero_Val(D2_Nonzero_Num) = 1;
    end
    for TimeSect = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
        STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+TimeSect;
        STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSect-1;
        D2_Nonzero_Num = D2_Nonzero_Num +1;
        D2_Nonzero_Row(D2_Nonzero_Num) = STH_tar;
        D2_Nonzero_Col(D2_Nonzero_Num) = STP_tar;
        D2_Nonzero_Val(D2_Nonzero_Num) = -1;
        STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSect;
        D2_Nonzero_Num = D2_Nonzero_Num +1;
        D2_Nonzero_Row(D2_Nonzero_Num) = STH_tar;
        D2_Nonzero_Col(D2_Nonzero_Num) = STP_tar;
        D2_Nonzero_Val(D2_Nonzero_Num) = 1;
    end
    for IncSpSIdx = find(sC(SpPIdx,:))
        UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
        for TimeSect = 0
            STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
            STP_tar     = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx);
            D2_Nonzero_Num = D2_Nonzero_Num +1;
            D2_Nonzero_Row(D2_Nonzero_Num) = STH_tar;
            D2_Nonzero_Col(D2_Nonzero_Num) = STP_tar;
            D2_Nonzero_Val(D2_Nonzero_Num) = sC(SpPIdx,IncSpSIdx);
        end
        for TimeSect = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
            for LocalTimeSectForSpS = 1:UpdRatio
                STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+TimeSect;
                STP_tar     = SpElemProperties.SpS.FirstSTPIdx(IncSpSIdx)+LocalTimeSectForSpS+UpdRatio*(TimeSect-1);
                D2_Nonzero_Num = D2_Nonzero_Num +1;
                D2_Nonzero_Row(D2_Nonzero_Num) = STH_tar;
                D2_Nonzero_Col(D2_Nonzero_Num) = STP_tar;
                D2_Nonzero_Val(D2_Nonzero_Num) = sC(SpPIdx,IncSpSIdx);
            end
        end
    end
end

%% D3

% for SpVIdx = find(SpElemProperties.SpV.Belong_to_ST_FI)
%     for SpPIdx = find(sD(SpVIdx,:))
%         UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(SpVIdx);
%         for TimeSect = 0
%              STV_tar     = SpElemProperties.SpV.FirstSTVIdx(SpVIdx);
%              STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx);
%              D3(STV_tar,STH_tar) = -1;
%              STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
%              D3(STV_tar,STH_tar) =  sD(SpVIdx,SpPIdx);
%         end
%         for TimeSect = 1:SpElemProperties.SpV.UpdNum(SpVIdx)
%             STV_tar     = SpElemProperties.SpV.FirstSTVIdx(SpVIdx)     +TimeSect;
%             STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx) +TimeSect-1;
%             D3(STV_tar,STH_tar) =  1;
%             STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx) +TimeSect  ;
%             D3(STV_tar,STH_tar) = -1;
%             for LocalTimeSec_for_SpP = 1:UpdRatio
%                 STH_tar = ...
%                     SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+LocalTimeSec_for_SpP+UpdRatio*(TimeSect-1);
%                 D3(STV_tar,STH_tar) = sD(SpVIdx,SpPIdx);
%             end
%         end
%     end
% end

STFISpVIdxs = find(SpElemProperties.SpV.Belong_to_ST_FI);
%STFISpPIdxs = find(sum(logical(sD(STFISpVIdxs,:)),1));
%STFISpSIdxs = find(sum(logical(sC(STFISpPIdxs,:)),1));
%STFISpNIdxs = find(sum(logical(sG(STFISpSIdxs,:)),1));
D3_Nonzero_Num = 0;
D3_Nonzero_Row = [];
D3_Nonzero_Col = [];
D3_Nonzero_Val = [];

for SpVIdx = STFISpVIdxs
    for TimeSect = 0
        STV_tar     = SpElemProperties.SpV.FirstSTVIdx(SpVIdx);
        STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx);
        D3_Nonzero_Num = D3_Nonzero_Num+1;
        D3_Nonzero_Row(D3_Nonzero_Num) = STV_tar;
        D3_Nonzero_Col(D3_Nonzero_Num) = STH_tar;
        D3_Nonzero_Val(D3_Nonzero_Num) = -1;
    end
    for TimeSect = 1:SpElemProperties.SpV.UpdNum(SpVIdx)
        STV_tar     = SpElemProperties.SpV.FirstSTVIdx(SpVIdx) +TimeSect;
        STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx) +TimeSect-1;
        D3_Nonzero_Num = D3_Nonzero_Num+1;
        D3_Nonzero_Row(D3_Nonzero_Num) = STV_tar;
        D3_Nonzero_Col(D3_Nonzero_Num) = STH_tar;
        D3_Nonzero_Val(D3_Nonzero_Num) = 1;
        STH_tar = SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx) +TimeSect  ;
        D3_Nonzero_Num = D3_Nonzero_Num+1;
        D3_Nonzero_Row(D3_Nonzero_Num) = STV_tar;
        D3_Nonzero_Col(D3_Nonzero_Num) = STH_tar;
        D3_Nonzero_Val(D3_Nonzero_Num) = -1;
    end
    for SpPIdx = find(sD(SpVIdx,:))
        UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(SpVIdx);
        for TimeSect = 0
            STV_tar     = SpElemProperties.SpV.FirstSTVIdx(SpVIdx);
            STH_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
            D3_Nonzero_Num = D3_Nonzero_Num+1;
            D3_Nonzero_Row(D3_Nonzero_Num) = STV_tar;
            D3_Nonzero_Col(D3_Nonzero_Num) = STH_tar;
            D3_Nonzero_Val(D3_Nonzero_Num) = sD(SpVIdx,SpPIdx);
        end
        for TimeSect = 1:SpElemProperties.SpV.UpdNum(SpVIdx)
            STV_tar     = SpElemProperties.SpV.FirstSTVIdx(SpVIdx) +TimeSect;
            for LocalTimeSec_for_SpP = 1:UpdRatio
                STH_tar = ...
                    SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+LocalTimeSec_for_SpP+UpdRatio*(TimeSect-1);
                D3_Nonzero_Num = D3_Nonzero_Num+1;
                D3_Nonzero_Row(D3_Nonzero_Num) = STV_tar;
                D3_Nonzero_Col(D3_Nonzero_Num) = STH_tar;
                D3_Nonzero_Val(D3_Nonzero_Num) = sD(SpVIdx,SpPIdx);
            end
        end
    end
end



D0 = sparse(D0_Nonzero_Row,D0_Nonzero_Col,D0_Nonzero_Val,Num_of_Elem.STS,Num_of_Elem.STN);
D1 = sparse(D1_Nonzero_Row,D1_Nonzero_Col,D1_Nonzero_Val,Num_of_Elem.STP,Num_of_Elem.STS);
D2 = sparse(D2_Nonzero_Row,D2_Nonzero_Col,D2_Nonzero_Val,Num_of_Elem.STOmega,Num_of_Elem.STP);
D3 = sparse(D3_Nonzero_Row,D3_Nonzero_Col,D3_Nonzero_Val,Num_of_Elem.STV,Num_of_Elem.STOmega);

%%
% %% D0
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%     for IncSpSIdx = find(sC(SpPIdx,:))
%         for IncSpNIdx = find(sG(IncSpSIdx,:))
%             UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(IncSpSIdx);
%             for CurrentTime = 0:SpElemProperties.SpS.UpdNum(IncSpSIdx)
%                 STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)   +CurrentTime;
%                 STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+UpdRatio*CurrentTime;
%                 D0(STS_tar,STN_tar) = sG(IncSpSIdx,IncSpNIdx);
%             end
%             for TimeSect = 0
%                 STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx);
%                 STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx);
%                 D0(STS_tar,STN_tar) =  1;
%             end
%             for TimeSect = 1:SpElemProperties.SpN.UpdNum(IncSpNIdx)
%                 STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx)+TimeSect;
%                 STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+TimeSect-1;
%                 D0(STS_tar,STN_tar) = -1;
%                 STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+TimeSect  ;
%                 D0(STS_tar,STN_tar) =  1;
%             end
%         end
%     end
% end
% 
% %% D1
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%     for IncSpSIdx = find(sC(SpPIdx,:))
%         UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
%         for CurrentTime = 0:SpElemProperties.SpP.UpdNum(SpPIdx)
%             STP_tar = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+CurrentTime;
%             STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpSIdx)+UpdRatio*CurrentTime;
%             D1(STP_tar,STS_tar) = sC(SpPIdx,IncSpSIdx);
%         end
%     end
% end
% for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
%     for TimeSect = 0
%         STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
%         STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx);
%         D1(STP_tar,STS_tar) = -1;
%     end
%     for TimeSect = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
%         STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSect;
%         STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx)+TimeSect-1;   
%         D1(STP_tar,STS_tar) =  1;
%         STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx)+TimeSect  ;   
%         D1(STP_tar,STS_tar) = -1;
%     end
%     for IncSpNIdx = find(sG(SpSIdx,:))
%         UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(SpSIdx);
%         for TimeSect = 0
%             STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
%             STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx);
%             D1(STP_tar,STS_tar) = sG(SpSIdx,IncSpNIdx);
%         end
%         for TimeSect = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
%             for LocalTimeSec_for_SpN = 1:UpdRatio
%                 STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSect;
%                 STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx)+LocalTimeSec_for_SpN+UpdRatio*(TimeSect-1);
%                 D1(STP_tar,STS_tar) = sG(SpSIdx,IncSpNIdx);
%             end
%         end
%     end
% end
% 
% %% D2
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%     for IncSpVIdx = find(sD(:,SpPIdx)).'
%         UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(IncSpVIdx);
%         for CurrentTime = 0:SpElemProperties.SpV.UpdNum(IncSpVIdx)
%             STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx)+CurrentTime;
%             STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+UpdRatio*CurrentTime;
%             D2(STOmega_tar,STP_tar) = sD(IncSpVIdx,SpPIdx);
%         end
%     end
% end
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%     for TimeSect = 0
%         STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
%         STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
%         D2(STOmega_tar,STP_tar) = -1;
%     end
%     
%     for TimeSect = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
%         STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+TimeSect;
%         STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSect-1;
%         D2(STOmega_tar,STP_tar) = -1;
%         STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSect;
%         D2(STOmega_tar,STP_tar) =  1;
%     end
% end
% for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
%     for IncSpPIdx = find(sC(:,SpSIdx)).'
%         UpdRatio = SpElemProperties.SpS.UpdNum(SpSIdx)/SpElemProperties.SpP.UpdNum(IncSpPIdx);
%         for TimeSect = 0
%             STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(IncSpPIdx);
%             STP_tar     = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
%             D2(STOmega_tar,STP_tar) =  sC(IncSpPIdx,SpSIdx);
%         end
%         for TimeSect = 1:SpElemProperties.SpP.UpdNum(SpPIdx)
%             for LocalTimeSectForSpS = 1:UpdRatio
%                 STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(IncSpPIdx)+TimeSect;
%                 STP_tar     = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+LocalTimeSectForSpS+UpdRatio*(TimeSect-1);
%                 D2(STOmega_tar,STP_tar) = sC(IncSpPIdx,SpSIdx);
%             end
%         end
%     end
% end
% 
% %% D3
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI==true)
%     for IncSpVIdx = find(sD(:,SpPIdx)).'
%         UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(IncSpVIdx);
%         for TimeSect = 0
%              STV_tar     = SpElemProperties.SpV.FirstSTVIdx(IncSpVIdx);
%              STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx);
%              D3(STV_tar,STOmega_tar) =  1;
%              STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
%              D3(STV_tar,STOmega_tar) = sD(IncSpVIdx,SpPIdx);
%         end
%         for TimeSect = 1:SpElemProperties.SpV.UpdNum(IncSpVIdx)
%             STV_tar     = SpElemProperties.SpV.FirstSTVIdx(IncSpVIdx)     +TimeSect;
%             STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx) +TimeSect-1;
%             D3(STV_tar,STOmega_tar) = -1;
%             STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx) +TimeSect  ;
%             D3(STV_tar,STOmega_tar) =  1;
%             for LocalTimeSec_for_SpP = 1:UpdRatio
%                 STOmega_tar = ...
%                     SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+LocalTimeSec_for_SpP+UpdRatio*(TimeSect-1);
%                 D3(STV_tar,STOmega_tar) = sD(IncSpVIdx,SpPIdx);
%             end
%         end
%     end
% end
%%

end