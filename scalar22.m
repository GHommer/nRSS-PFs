%Garrison Hommer, 23JAN2014
%Computes scalar product of 2 2nd order tensors
%Rudnicki, Fundamentals of Continuum Mechanics, eqn. 3.21 pg. 28

function scalout = scalar22(tens2in1, tens2in2)

scalout=0;
for i=1:3
    for j=1:3
        scalout = scalout + tens2in1(i,j)*tens2in2(i,j);
    end
end
