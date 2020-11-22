function TFSFSource = Generate_TFSFSource(IncWaveEFuncHandle,IncWaveHFuncHandle,IncWaveEDirec,IncWaveHDirec,WaveNumVec,SpElemProperties,sG,sC,sD,NodePosPrim_M,NodePosDual_M,cdt)

TFSFSourceIdx = 0;
FirstST_SourceIdx = 1;
WaveNumVecDirec = norm(WaveNumVec)^(-1)*WaveNumVec;
for SpPIdx = find(SpElemProperties.SpP.isInSFRegion)
    for IncSpSIdx = find(sC(SpPIdx,:))
        if SpElemProperties.SpS.isOnTFSFBoundary(IncSpSIdx)
            TFSFSourceIdx = TFSFSourceIdx + 1;
            IncSpSPosVec = 0.5*NodePosPrim_M*abs(sG(IncSpSIdx,:).');
            IncSpSDirec  = NodePosPrim_M*sG(IncSpSIdx,:).';
            SourceDoFFuncHandle = @(w) -sC(SpPIdx,IncSpSIdx)*(cdt/SpElemProperties.SpS.UpdNum(IncSpSIdx))*dot(IncWaveEDirec,IncSpSDirec)*IncWaveEFuncHandle(w - dot(WaveNumVecDirec, IncSpSPosVec) );
            TFSFSource(TFSFSourceIdx).Elec_Magn      = 'Magn';
            TFSFSource(TFSFSourceIdx).UpdNum         = SpElemProperties.SpP.UpdNum(SpPIdx);
            TFSFSource(TFSFSourceIdx).PrimFace_tgt   = SpPIdx;
            TFSFSource(TFSFSourceIdx).SourceDoFFuncHandle = SourceDoFFuncHandle;
            TFSFSource(TFSFSourceIdx).FirstST_SourceIdx = FirstST_SourceIdx;
            FirstST_SourceIdx                           = FirstST_SourceIdx + TFSFSource(TFSFSourceIdx).UpdNum;
        end
    end
end
for SpSIdx = find(SpElemProperties.SpS.isOnTFSFBoundary)
    for IncSpPIdx = find(sC(:,SpSIdx).')
        if SpElemProperties.SpP.isInSFRegion(IncSpPIdx)
            TFSFSourceIdx = TFSFSourceIdx + 1;
            IncSpPPosVec = 0.5*NodePosDual_M*abs(sD(:,IncSpPIdx));
            IncSpPDirec  = NodePosDual_M*sD(:,IncSpPIdx);
            SourceDoFFuncHandle = @(w) sC(IncSpPIdx,SpSIdx)*(cdt/SpElemProperties.SpP.UpdNum(IncSpPIdx))*dot(IncWaveHDirec,IncSpPDirec)*IncWaveHFuncHandle(w - dot(WaveNumVecDirec, IncSpPPosVec) );
            TFSFSource(TFSFSourceIdx).Elec_Magn      = 'Elec';
            TFSFSource(TFSFSourceIdx).UpdNum         = SpElemProperties.SpS.UpdNum(SpSIdx);
            TFSFSource(TFSFSourceIdx).DualFace_tgt   = IncSpPIdx;
            TFSFSource(TFSFSourceIdx).SourceDoFFuncHandle = SourceDoFFuncHandle;
            TFSFSource(TFSFSourceIdx).FirstST_SourceIdx = FirstST_SourceIdx;
            FirstST_SourceIdx = FirstST_SourceIdx + TFSFSource(TFSFSourceIdx).UpdNum;
        end
    end
end



end