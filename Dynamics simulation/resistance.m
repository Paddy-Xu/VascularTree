% ok, just to calculate those parameters
% questions: 1. what does this "data" mean, data probably in size N*2
%            2. what does 12 mean
%            3. ok, maybe write equation done to see the equations
% 

for i=1:12
    d=data(i,1);
    n=(1+(6*exp(-0.085*d)+2.2-2.44*exp(-0.06*(d^0.645))) ...
    *((d/(d-1.1))^2))*((d/(d-1.1))^2);
    C(i)=(pi*d*d*d*d)/(128*n*data(i,2));
    R(i)=(128*n*data(i,2))/(pi*d*d*d*d);
    Cn(i)=(3.141592 * d * d * d * d)/(8 * n * data(i,2));
end

plot(data(:,1),C,data(:,1),Cn)