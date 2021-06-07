function[macfracdiff] = maxfracdiff(in1, in2)

L1 = length(in1);
L2 = length(in2);

if L1 == L2
    macfracdiff = max((in2-in1)/in1);
elseif L2 < L1
    x1 = linspace(0,1,L1);
    x2 = linspace(0,1,L2);    
    interpolated = interp1(x2,in2,x1);
    macfracdiff = max((interpolated-in1)/in1);
elseif L2 > L1
    x1 = linspace(0,1,L1);
    x2 = linspace(0,1,L2);
    interpolated = interp1(x1,in1,x2);
    macfracdiff = max((in2-interpolated)/interpolated);
end
