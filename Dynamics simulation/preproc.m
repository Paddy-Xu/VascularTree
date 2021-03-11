% nodesMatrix: at least N*10
% vesselsMatrix: at vesselsN *4
% nephronsMatrix: at least N*1

% psiMatrix = zeros(nephronsN,nephronsN)
function psiMatrix=preproc(nodesMatrix,vesselsMatrix,nephronsMatrix)
R_resistivity=.1;
R_leak=0.1;
R_neph=1;
current=1;


nephronsN=4;
vesselsN=7;
nodesN=3;

electricalResistance = zeros(1,vesselsN);
for j=1:1:vesselsN
    length = vesselsMatrix(j,3);
    inDiameter = vesselsMatrix(j,4);
    outDiameter=inDiameter*1.1;
    electricalResistance(j)=(R_resistivity*length)/(pi/4*(outDiameter^2 - inDiameter^2));
end


psiMatrix=zeros(nephronsN,nephronsN);

yZero=zeros(nephronsN+nodesN,1);

tSpan = [0 500];

for number=1:1:nephronsN
    
sol=ode23s(@model, tSpan, yZero);
t = linspace(tSpan(1), tSpan(2), 1500);
y = deval(sol,t);


for k=1:1:nephronsN
    psiMatrix(number,k)=y(k,tSpan(2));
    
end
psiScale=1/sum(psiMatrix(number,:));
psiMatrix(number,:)=psiMatrix(number,:).*psiScale;

end


    function yPrime=model(t,y)
        %yPrime=zeros(nephronsN+nephronsN,1);
        
        for k=1:1:nephronsN
            cur=0;
            cnidx = find ( (nodesMatrix(:,9) == nephronsMatrix(k,1) | nodesMatrix(:,10) == nephronsMatrix(k,1) ), 1, 'first');
            cvidx = find ( (vesselsMatrix(:,1) == nephronsMatrix(k,1)), 1, 'first');
            
            if k==number
                cur=current;
            end
            % question: what is the size of yPrime? 
            % ok, probably nephronsN + nodesN, so;
            % yPrime = zeros(nephronsN + nodesN,1)
            yPrime(k,1)=cur+(y(nephronsN+cnidx)-y(k))/electricalResistance(cvidx) - y(k)/R_neph;
            
        end
        previousEq=nephronsN;
        for k=1:1:nodesN
            cn1idx = find (nodesMatrix(:,1) == nodesMatrix(k,3), 1, 'first');
            cn2idx = find (nodesMatrix(:,1) == nodesMatrix(k,4), 1, 'first');
            cn3idx = find (nodesMatrix(:,1) == nodesMatrix(k,5), 1, 'first');
            
            cv1idx = find (vesselsMatrix(:,1) == nodesMatrix(k,6), 1, 'first');
            cv2idx = find (vesselsMatrix(:,1) == nodesMatrix(k,7),1, 'first');
            cv3idx = find (vesselsMatrix(:,1) == nodesMatrix(k,8),1, 'first');
            
            
            if (nodesMatrix (k,3)==-1)
                cpsi1=y(previousEq+k);
            else
                cpsi1=y(previousEq+cn1idx);
            end
            if (nodesMatrix (k,4)==-1)
                cneph1idx = find (nephronsMatrix(:,1) == nodesMatrix(k,9), 1, 'first');
                cpsi2=y(cneph1idx);
            else
                cpsi2=y(previousEq+cn2idx);
            end
            if (nodesMatrix (k,5)==-1)
                cneph2idx = find (nephronsMatrix(:,1) == nodesMatrix(k,10), 1, 'first');
                cpsi3=y(cneph2idx);
            else
                cpsi3=y(previousEq+cn3idx);
            end
                                   
            yPrime (previousEq+k,1) =(cpsi1-y(previousEq+k))/electricalResistance(cv1idx)+(cpsi2-y(previousEq+k))/electricalResistance(cv2idx)+(cpsi3-y(previousEq+k))/electricalResistance(cv3idx) - y(previousEq+k)/R_leak;
            
        end
        
        
    end

end