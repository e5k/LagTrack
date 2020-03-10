function Cd = getCdHighMach(mach, machMeth)

dragC = load('drag.mat');
drag = dragC.drag;

if strcmp(machMeth, 'high cube')
    fld = 'highcube';
elseif strcmp(machMeth, 'low cube')
    fld = 'lowcube';
elseif strcmp(machMeth, 'sphere')
    fld = 'sphere';
elseif strcmp(machMeth, 'cylinder')
    fld = 'cylinder';  
end

if mach >= drag.(fld)(end,1)
    Cd = drag.(fld)(end,2);
else
    Cd = interp1(drag.(fld)(:,1), drag.(fld)(:,2), mach, 'Spline');
end
