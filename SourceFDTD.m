function Source = SourceFDTD(MeshMeasurements,SpElemProperties,sC)
global SourcePeriod
SourcePeriod = 20;

disp("Waveform: Sinusoidal")
Lightspeed = 1;
wavelength = Lightspeed/(1/SourcePeriod);
disp(['Wavelength/meshsize = ', num2str(wavelength/max([MeshMeasurements.dx MeshMeasurements.dy MeshMeasurements.dz]) ) ])

XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = (XSize+1);
YEdgePerXYPlane     = (XSize+1)*YSize;
%YEdgeNum           = (XSize+1)*YSize*(ZSize+1);
%ZEdgePerXRow       = (XSize+1);
%ZEdgePerXYPlae     = (XSize+1)*(YSize+1);
%ZEdgeNum           = (XSize+1)*(YSize+1)*ZSize;

% YZFacePerXRow       = (XSize+1);
% YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum           = (XSize+1)*YSize*ZSize;
% ZXFacePerXRow       = XSize;
% ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum           = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum          = XSize*YSize*(ZSize+1);

X_SourceCenter = round(1*XSize/2);
Y_SourceCenter = round(1*YSize/2);
Z_SourceCenter = round(1*ZSize/2);

XIdx = X_SourceCenter;
YIdx = Y_SourceCenter;
ZIdx = Z_SourceCenter;
SourceSIdx(1) = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(2) = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(3) = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SourceSIdx(4) = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SpPIdx(1:4) =  XIdx   + (YIdx-1)*XYFacePerXRow  + (ZIdx-1)*XYFacePerXYPlane + YZFaceNum +ZXFaceNum;
Source = struct;
FirstST_SourceIdx = 1;

XIdx = X_SourceCenter+1;
YIdx = Y_SourceCenter;
ZIdx = Z_SourceCenter;
SourceSIdx(5) = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(6) = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(7) = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SourceSIdx(8) = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SpPIdx(5:8) =  XIdx   + (YIdx-1)*XYFacePerXRow  + (ZIdx-1)*XYFacePerXYPlane + YZFaceNum +ZXFaceNum;

XIdx = X_SourceCenter;
YIdx = Y_SourceCenter+1;
ZIdx = Z_SourceCenter;
SourceSIdx(9) = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(10) = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(11) = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SourceSIdx(12) = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SpPIdx(9:12) =  XIdx   + (YIdx-1)*XYFacePerXRow  + (ZIdx-1)*XYFacePerXYPlane + YZFaceNum +ZXFaceNum;

XIdx = X_SourceCenter+1;
YIdx = Y_SourceCenter+1;
ZIdx = Z_SourceCenter;
SourceSIdx(13) = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(14) = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(15) = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SourceSIdx(16) = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SpPIdx(13:16) =  XIdx   + (YIdx-1)*XYFacePerXRow  + (ZIdx-1)*XYFacePerXYPlane + YZFaceNum +ZXFaceNum;

for SourceIdx = 1:size(SourceSIdx,2)
    SpSIdx = SourceSIdx(SourceIdx);
    Source(SourceIdx).UpdNum                    = SpElemProperties.SpS.UpdNum(SpSIdx);
    Source(SourceIdx).DualFace_tgt              = SpSIdx;
    Source(SourceIdx).WaveformFunctionHandle    = @sinewave;
    %Source(SourceIdx).WaveformFunctionHandle    = @zerowave;
    Source(SourceIdx).WaveformSign              = sC(SpPIdx(SourceIdx),SpSIdx);
    Source(SourceIdx).Area_TargetDualFace       = 1;
    Source(SourceIdx).FirstST_SourceIdx         = FirstST_SourceIdx;
    FirstST_SourceIdx                           = FirstST_SourceIdx + Source(SourceIdx).UpdNum;
end


end

function val = sinewave(x)
global SourcePeriod
    val = sin(x/SourcePeriod);
end