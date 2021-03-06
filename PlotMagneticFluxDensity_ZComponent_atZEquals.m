function PlotMagneticFluxDensity_ZComponent_atZEquals(ZConst,FieldDoFs,FaceArea,SpElemPositionIdx_SpP,MeshMeasurements)
ZConst = round(ZConst);
XSize = MeshMeasurements.XCoord/MeshMeasurements.dxCoarse;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dyCoarse;
%ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

B_SquareoidMesh = zeros(XSize,YSize);
Area_Squareoid  = zeros(XSize,YSize);
% ZPos_of_SpPs = [PrimFacePos(:).Vec];
ZPos_of_SpPs = SpElemPositionIdx_SpP(3,:);

for SpPIdx = find(ZPos_of_SpPs==ZConst)
    i = ceil(SpElemPositionIdx_SpP(1,SpPIdx));
    j = ceil(SpElemPositionIdx_SpP(2,SpPIdx));
    B_SquareoidMesh(i,j) = B_SquareoidMesh(i,j) + FieldDoFs(SpPIdx);
    Area_Squareoid(i,j)  =  Area_Squareoid(i,j) + FaceArea.Prim(SpPIdx);
end

B_SquareoidMesh=B_SquareoidMesh./Area_Squareoid;

%%

figure('name',['Bz, at Z=', num2str(ZConst),' Calculated by New Method'])
xa = gca;
mesh(B_SquareoidMesh.')
xlabel('x','FontSize',30)
ylabel('y','FontSize',30)
xlim([0+0.5 XSize+0.5])
ylim([0+0.5 YSize+0.5])
xticks([0+0.5 XSize/4+0.5 XSize/2+0.5 XSize*3/4+0.5 XSize/+0.5])
yticks([0+0.5 YSize/4+0.5 YSize/2+0.5 YSize*3/4+0.5 YSize/+0.5])
xa.FontSize = 20;
xticklabels([0 XSize/4 XSize/2 XSize*3/4 XSize])
yticklabels([0 YSize/4 YSize/2 YSize*3/4 YSize])
zlabel(['Bz, at Z=', num2str(ZConst),' Calculated by New Method'],'FontSize',15)
pbaspect([XSize/XSize YSize/XSize 1])
%color = colorbar('southoutside');
color = colorbar;
color.Label.String = 'B_{z}';


end