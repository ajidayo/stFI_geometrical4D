function PlotMagneticFluxDensity2D_atYEquals(FieldDoFs,YConst,FaceArea,PrimFacePos,MeshMeasurements)
global SpDIM EPSILON

dx = MeshMeasurements.dxCoarse;
dy = MeshMeasurements.dyCoarse;
dz = MeshMeasurements.dzCoarse;
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

% Disp_dx = dx;
% Disp_dy = dy;
% Disp_dz = dz;

B_SquareoidMesh = zeros(XSize+1,YSize+1,ZSize+1,SpDIM);
Area_Squareoid  = ones(XSize+1,YSize+1,ZSize+1,SpDIM);
for SpPIdx = 1:size(FaceArea.Prim,1)
    xIdx = PrimFacePos(SpPIdx).Vec(1)/dx;
    yIdx = PrimFacePos(SpPIdx).Vec(2)/dy;
    zIdx = PrimFacePos(SpPIdx).Vec(3)/dz;
    if abs(xIdx - round(xIdx)) < EPSILON
        dimIdx = 1;
        xIdx = round(xIdx)+1;
        yIdx =  ceil(yIdx+EPSILON);
        zIdx =  ceil(zIdx+EPSILON);
    elseif abs(yIdx - round(yIdx)) < EPSILON
        dimIdx = 2;
        xIdx =  ceil(xIdx+EPSILON);
        yIdx = round(yIdx)+1;
        zIdx =  ceil(zIdx+EPSILON);
    elseif abs(zIdx - round(zIdx)) < EPSILON
        dimIdx = 3;
        xIdx =  ceil(xIdx+EPSILON);
        yIdx =  ceil(yIdx+EPSILON);
        zIdx = round(zIdx)+1;
    else
        disp("hoge")
        continue;
    end
    B_SquareoidMesh(xIdx,yIdx,zIdx,dimIdx) ...
        = B_SquareoidMesh(xIdx,yIdx,zIdx,dimIdx) + FieldDoFs(SpPIdx);
    Area_Squareoid(xIdx,yIdx,zIdx,dimIdx)  ...
        =  Area_Squareoid(xIdx,yIdx,zIdx,dimIdx) +       FaceArea.Prim(SpPIdx);
end

B_SquareoidMesh=B_SquareoidMesh./Area_Squareoid;


%%

YConst = round(YConst);
PlottingVectorIdx = 0;
for zIdx = 1:ZSize
    yIdx = YConst;
    for xIdx = 1:XSize
        PlottingVectorIdx = PlottingVectorIdx+1;
        
        VecPosX(PlottingVectorIdx) = (xIdx-0.5)*dx;
        VecPosZ(PlottingVectorIdx) = (zIdx-0.5)*dz;
        B_Vol_X(PlottingVectorIdx) = ...
            (1/6)*...
            (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,1)...
            +B_SquareoidMesh(xIdx+1,yIdx  ,zIdx  ,1));
        B_Vol_Z(PlottingVectorIdx) = ...
            (1/6)*...
            (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,3)...
            +B_SquareoidMesh(xIdx  ,yIdx  ,zIdx+1,3));
    end
end

%%
figure('name',['Bx*i+Bz*k at Y =',num2str(YConst-0.5) 'Calculated by New Method'])
xa = gca;
scale = 20;
quiver(VecPosX,VecPosZ,B_Vol_X,B_Vol_Z,scale)
xlabel('x','FontSize',30)
ylabel('y','FontSize',30)
zlabel('z','FontSize',30)
xlim([0 XSize*dx])
ylim([0 YSize*dy])
zlim([0 ZSize*dz])
xticks([0 dx*XSize/4 dx*XSize/2 dx*XSize*3/4 dx*XSize])
yticks([0 dy*YSize/4 dy*YSize/2 dy*YSize*3/4 dy*YSize])
zticks([0 dz*ZSize/4 dz*ZSize/2 dz*ZSize*3/4 dz*ZSize])
xa.FontSize = 20;
xticks([0 dx*XSize/4 dx*XSize/2 dx*XSize*3/4 dx*XSize])
yticks([0 dy*YSize/4 dy*YSize/2 dy*YSize*3/4 dy*YSize])
zticks([0 dz*ZSize/4 dz*ZSize/2 dz*ZSize*3/4 dz*ZSize])
pbaspect([1 1 1])

end