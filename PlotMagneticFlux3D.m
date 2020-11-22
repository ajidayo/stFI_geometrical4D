function PlotMagneticFlux3D(FieldDoFs,FaceArea,SpElemPositionIdx_SpP,MeshMeasurements)

global SpDIM EPSILON

dx = MeshMeasurements.dxCoarse;
dy = MeshMeasurements.dyCoarse;
dz = MeshMeasurements.dzCoarse;
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dzCoarse;

B_SquareoidMesh = zeros(XSize+1,YSize+1,ZSize+1,SpDIM);
Area_Squareoid  = ones(XSize+1,YSize+1,ZSize+1,SpDIM);
for SpPIdx = 1:size(FaceArea.Prim,1)
    xIdx = SpElemPositionIdx_SpP(1,SpPIdx)/dx;
    yIdx = SpElemPositionIdx_SpP(2,SpPIdx)/dy;
    zIdx = SpElemPositionIdx_SpP(3,SpPIdx)/dz;
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
        continue;
    end
    B_SquareoidMesh(xIdx,yIdx,zIdx,dimIdx) ...
        = B_SquareoidMesh(xIdx,yIdx,zIdx,dimIdx) + FieldDoFs(SpPIdx);
    Area_Squareoid(xIdx,yIdx,zIdx,dimIdx)  ...
        =  Area_Squareoid(xIdx,yIdx,zIdx,dimIdx) +       FaceArea.Prim(SpPIdx);
end

B_SquareoidMesh=B_SquareoidMesh./Area_Squareoid;


B_Vol_X = zeros(YSize,XSize,ZSize);
B_Vol_Y = zeros(YSize,XSize,ZSize);
B_Vol_Z = zeros(YSize,XSize,ZSize);
[PosX,PosY,PosZ] ...
    = meshgrid((0.5:1:MeshMeasurements.YCoord-0.5)*dy,(0.5:1:MeshMeasurements.XCoord-0.5)*dx,(0.5:1:MeshMeasurements.ZCoord-0.5)*dz);
for zIdx = 1:ZSize 
    for yIdx = 1:YSize
        for xIdx = 1:XSize
            B_Vol_X(yIdx,xIdx,zIdx) = ...
                (1/6)*...
                (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,1)...
                +B_SquareoidMesh(xIdx+1,yIdx  ,zIdx  ,1));
            B_Vol_Y(yIdx,xIdx,zIdx) = ...
                (1/6)*...
                (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,2)...
                +B_SquareoidMesh(xIdx  ,yIdx+1,zIdx  ,2));
            B_Vol_Z(yIdx,xIdx,zIdx) = ...
                (1/6)*...
                (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,3)...
                +B_SquareoidMesh(xIdx  ,yIdx  ,zIdx+1,3));
            normB(yIdx,xIdx,zIdx) = norm([B_Vol_X(yIdx,xIdx,zIdx);B_Vol_Y(yIdx,xIdx,zIdx);B_Vol_Z(yIdx,xIdx,zIdx)]);
            B_Vol_X(yIdx,xIdx,zIdx) = normB(yIdx,xIdx,zIdx).^(-1)*B_Vol_X(yIdx,xIdx,zIdx);
            B_Vol_Y(yIdx,xIdx,zIdx) = normB(yIdx,xIdx,zIdx).^(-1)*B_Vol_Y(yIdx,xIdx,zIdx);
            B_Vol_Z(yIdx,xIdx,zIdx) = normB(yIdx,xIdx,zIdx).^(-1)*B_Vol_Z(yIdx,xIdx,zIdx);
            NormCheck(yIdx,xIdx,zIdx) =  norm([B_Vol_X(yIdx,xIdx,zIdx);B_Vol_Y(yIdx,xIdx,zIdx);B_Vol_Z(yIdx,xIdx,zIdx)]);
        end
    end
end
%%

figure('name',['B, Calculated by New Method, snapshot at Y = ',num2str(0.5*MeshMeasurements.YCoord)])
xa = gca;
%[startx,starty,startz] = meshgrid(0.5:2:MeshMeasurements.XCoord-0.5,0.5:2:MeshMeasurements.YCoord-0.5,0.5:2:MeshMeasurements.ZCoord-0.5);
[startx,starty,startz] = meshgrid(0.5:1:MeshMeasurements.XCoord-0.5,0.5*MeshMeasurements.YCoord,0.5:1:MeshMeasurements.ZCoord-0.5);
%streamline(PosX,PosY,PosZ,B_Vol_X,B_Vol_Y,B_Vol_Z,startx,starty,startz)
Scale = 0.5;
hcone = coneplot(B_Vol_X,B_Vol_Y,B_Vol_Z,startx,starty,startz,Scale,normB);
%hcone.FaceColor = 'red';
hcone.EdgeColor = 'none';
%streamslice(B_Vol_X,B_Vol_Y,B_Vol_Z,[],10,[])
%streamline(stream3(PosX,PosY,PosZ,B_Vol_X,B_Vol_Y,B_Vol_Z,startx,starty,startz))
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
pbaspect([XSize/XSize YSize/XSize ZSize/XSize])
%color = colorbar('southoutside');
colorB = colorbar;
colorB.Label.String = 'B';
view([0 -1 0]) 


figure('name',['B, Calculated by New Method, snapshot at Y = ',num2str(0.5*MeshMeasurements.YCoord)])
xa = gca;
%[startx,starty,startz] = meshgrid(0.5:2:MeshMeasurements.XCoord-0.5,0.5:2:MeshMeasurements.YCoord-0.5,0.5:2:MeshMeasurements.ZCoord-0.5);
[startx,starty,startz] = meshgrid(5.5:1:MeshMeasurements.XCoord-5.5,0.5*MeshMeasurements.YCoord,5.5:1:MeshMeasurements.ZCoord-5.5);
%streamline(PosX,PosY,PosZ,B_Vol_X,B_Vol_Y,B_Vol_Z,startx,starty,startz)
Scale = 0.5;
hcone = coneplot(B_Vol_X,B_Vol_Y,B_Vol_Z,startx,starty,startz,Scale,normB);
%hcone.FaceColor = 'red';
hcone.EdgeColor = 'none';
%streamslice(B_Vol_X,B_Vol_Y,B_Vol_Z,[],10,[])
%streamline(stream3(PosX,PosY,PosZ,B_Vol_X,B_Vol_Y,B_Vol_Z,startx,starty,startz))
xlabel('x','FontSize',50)
ylabel('y','FontSize',50)
zlabel('z','FontSize',50)
xlim([5 (XSize-5)*dx])
ylim([5 (YSize-5)*dy])
zlim([5 (ZSize-5)*dz])
xticks([5 10 15 20 25 30 35])
yticks([5 10 15 20 25])
zticks([5 10 15 20 25])
xa.FontSize = 40;
pbaspect([1 (YSize-10)/(XSize-10) (ZSize-10)/(XSize-10)])
%color = colorbar('southoutside');
colorB = colorbar;
colorB.Label.String = 'Magnetic Flux Density';
view([0 -1 0]) 

end