%Garrison Hommer, 27JAN2014
%Computes dyad of 2 vectors 
%Rudnicki, Fundamentals of Continuum Mechanics, eqn. 3.23 pg. 24

function tensout = dyad11(vecin1, vecin2)

tensout=[0 0 0;0 0 0;0 0 0];
for i=1:3
    for j=1:3
        tensout(i,j) = tensout(i,j) + vecin1(i) * vecin2(j);
    end
end
