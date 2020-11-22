function Impedance = Func_SphereImpedance(PosVec,CenterPosVec,Radius)
if norm(PosVec-CenterPosVec)<Radius
    Impedance = 0.01;
else
    Impedance = 1;
end

end