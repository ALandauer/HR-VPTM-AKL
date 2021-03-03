% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function tube2mla = LFM_computeTube2MLA(lensPitch, mla2sensor, deltaOT, objRad, ftl)

% tube2mla satisfies 3 conditions: 
%1. matching effective f number with the
%2. effective tubeLensRadius
%3. the above when focused on the MLA

% solve quadratic equation 
a = ftl*lensPitch/(2*mla2sensor);
% b = delta_ot*objRad - ftl*objRad - mla2sensor*ftl;
b = deltaOT*objRad - ftl*objRad;
c = -ftl*deltaOT*objRad;

delta_equation = sqrt(b^2 - 4*a*c);

tube2mla = [(-b+delta_equation)/(2*a), (-b-delta_equation)/(2*a)];
tube2mla = tube2mla(tube2mla>0);