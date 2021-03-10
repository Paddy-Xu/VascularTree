figure
subplot(1,2,1)
hold on
bar(x,f/trapz(x,f));
plot((0:1:500),res,'r');
hold off

subplot(1,2,2)
image(clarity)
