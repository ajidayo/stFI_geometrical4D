function Impedance = Func_Impedance(PosVec)
Radius = 0;
CenterPosVec = [25;15;15];
if norm(PosVec-CenterPosVec)<Radius
    Impedance = 0.01;
else
    Impedance = 1;
end

end