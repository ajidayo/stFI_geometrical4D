function Source = GenerateSource(MeshMeasurements,SpElemProperties,ElemPer,sC)
global SourcePeriod
SourcePeriod = 20;
WaveformPreset = 1;
switch WaveformPreset
    case 1
        disp('GenerateSource: Waveform - Sinusoidal')
        WaveformFunctionHandle = @sinewave;
    case 2
        disp('GenerateSource: Waveform - Constantly Zero')
        WaveformFunctionHandle = @zerowave;
end
Lightspeed = 1
wavelength = Lightspeed/(1/SourcePeriod);
disp(['Wavelength/meshsize = ', num2str(wavelength/max([MeshMeasurements.dxCoarse MeshMeasurements.dyCoarse MeshMeasurements.dzCoarse]) ) ])

XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

X_SourceCenter = round(1*XSize/2);
Y_SourceCenter = round(1*YSize/2);
Z_SourceCenter = round(1*ZSize/2);
SourceIdx = 0;
for ZIdx = Z_SourceCenter+1
    for YIdx = [Y_SourceCenter  Y_SourceCenter+1]
        for XIdx = [X_SourceCenter  X_SourceCenter+1]
            SpP_tgt =  1 ...
                + sum(ElemPer.XYFacePerCoarseGrid(1:XIdx-1,YIdx,ZIdx)) ...
                + sum(ElemPer.XYFacePerXRow(1:YIdx-1,ZIdx)) ...
                + sum(ElemPer.XYFacePerXYPlane(1:ZIdx-1)) ...
                + ElemPer.YZFaceNum + ElemPer.ZXFaceNum;
            for IncSpSIdx = find(sC(SpP_tgt,:))
                SourceIdx = SourceIdx+1;
                SourceSpSIdx(SourceIdx) = IncSpSIdx;
                LoopCoilSpPIdx(SourceIdx) = SpP_tgt;
            end
        end
    end
end

Source = struct;
FirstST_SourceIdx = 1;
for SourceIdx = 1:size(SourceSpSIdx,2)
    SpSIdx = SourceSpSIdx(SourceIdx);
    Source(SourceIdx).UpdNum                    = SpElemProperties.SpS.UpdNum(SpSIdx);
    Source(SourceIdx).DualFace_tgt              = SpSIdx;
    Source(SourceIdx).WaveformFunctionHandle    = WaveformFunctionHandle;
    Source(SourceIdx).WaveformSign              = sC(LoopCoilSpPIdx(SourceIdx),SpSIdx);
    Source(SourceIdx).Area_TargetDualFace       = 1;
    Source(SourceIdx).FirstST_SourceIdx         = FirstST_SourceIdx;
    FirstST_SourceIdx                           = FirstST_SourceIdx + Source(SourceIdx).UpdNum;
end

% 
% XIdx = X_SourceCenter;
% YIdx = Y_SourceCenter;
% ZIdx = Z_SourceCenter;
% SourceSIdx(1) = XIdx   + (YIdx-1)*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(2) = XIdx   + (YIdx  )*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(3) = XIdx   + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*ElemPer.YEdgePerXYPlane +ElemPer.XEdgeNum;
% SourceSIdx(4) = XIdx+1 + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*ElemPer.YEdgePerXYPlane +ElemPer.XEdgeNum;
% SpP_tgt(1:4) =  XIdx   + (YIdx-1)*ElemPer.XYFacePerXRow  + (ZIdx-1)*ElemPer.XYFacePerXYPlane + ElemPer.YZFaceNum +ElemPer.ZXFaceNum;
% 
% XIdx = X_SourceCenter+1;
% YIdx = Y_SourceCenter;
% ZIdx = Z_SourceCenter;
% SourceSIdx(5) = XIdx   + (YIdx-1)*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(6) = XIdx   + (YIdx  )*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(7) = XIdx   + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*ElemPer.YEdgePerXYPlane +ElemPer.XEdgeNum;
% SourceSIdx(8) = XIdx+1 + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*ElemPer.YEdgePerXYPlane +ElemPer.XEdgeNum;
% SpP_tgt(5:8) =  XIdx   + (YIdx-1)*ElemPer.XYFacePerXRow  + (ZIdx-1)*ElemPer.XYFacePerXYPlane + ElemPer.YZFaceNum +ElemPer.ZXFaceNum;
% 
% XIdx = X_SourceCenter;
% YIdx = Y_SourceCenter+1;
% ZIdx = Z_SourceCenter;
% SourceSIdx(9) = XIdx   + (YIdx-1)*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(10) = XIdx   + (YIdx  )*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(11) = XIdx   + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +ElemPer.XEdgeNum;
% SourceSIdx(12) = XIdx+1 + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +ElemPer.XEdgeNum;
% SpP_tgt(9:12) =  XIdx   + (YIdx-1)*ElemPer.XYFacePerXRow  + (ZIdx-1)*ElemPer.XYFacePerXYPlane + ElemPer.YZFaceNum +ElemPer.ZXFaceNum;
% 
% XIdx = X_SourceCenter+1;
% YIdx = Y_SourceCenter+1;
% ZIdx = Z_SourceCenter;
% SourceSIdx(13) = XIdx   + (YIdx-1)*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(14) = XIdx   + (YIdx  )*ElemPer.XEdgePerXRow  + (ZIdx-1)*ElemPer.XEdgePerXYPlane;
% SourceSIdx(15) = XIdx   + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*ElemPer.YEdgePerXYPlane +ElemPer.XEdgeNum;
% SourceSIdx(16) = XIdx+1 + (YIdx-1)*ElemPer.YEdgePerXRow  + (ZIdx-1)*ElemPer.YEdgePerXYPlane +ElemPer.XEdgeNum;
% SpP_tgt(13:16) =  XIdx   + (YIdx-1)*ElemPer.XYFacePerXRow  + (ZIdx-1)*ElemPer.XYFacePerXYPlane + ElemPer.YZFaceNum +ElemPer.ZXFaceNum;
% 
% for SourceIdx = 1:size(SourceSIdx,2)
%     SpSIdx = SourceSIdx(SourceIdx);
%     Source(SourceIdx).UpdNum                    = SpElemProperties.SpS.UpdNum(SpSIdx);
%     Source(SourceIdx).DualFace_tgt              = SpSIdx;
%     Source(SourceIdx).WaveformFunctionHandle    = @sinewave;
%     %Source(SourceIdx).WaveformFunctionHandle    = @zerowave;
%     Source(SourceIdx).WaveformSign              = sC(SpP_tgt(SourceIdx),SpSIdx);
%     Source(SourceIdx).Area_TargetDualFace       = 1;
%     Source(SourceIdx).FirstST_SourceIdx         = FirstST_SourceIdx;
%     FirstST_SourceIdx                           = FirstST_SourceIdx + Source(SourceIdx).UpdNum;
% end

end

function val = sinewave(x)
global SourcePeriod
    val = sin(x/SourcePeriod);
end

function val = zerowave(x)
    val = 0*x;
end