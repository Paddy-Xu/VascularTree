% f: number of points inside each bin
% x: center values
% pd: a struct of Exp
% trapz: integrals: areas under curve
% QUESTION: what does f/trapz(x,f) mean

[f,x]= hist(aaFull2',20);
pd=fitdist(aaFull2','Exponential');
res=1/pd.mu*exp(-(0:1:500)/pd.mu);

%mu of clraity*0.8= 42.325850369432885

figure
subplot(1,2,1)
hold on
bar(x,f/trapz(x,f));
plot((0:1:500),res,'r');
hold off

subplot(1,2,2)
image(clarity)

%%