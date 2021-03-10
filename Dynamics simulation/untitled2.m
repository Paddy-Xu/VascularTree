test_fun(1)
x = randn(1000,1);  
nbins = 50;
hist(x,nbins)
load hospital
x = hospital.Weight;
pd = fitdist(x,'Normal')
