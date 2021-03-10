%Created by Postnov D.D. 07.04.2014
%Vessel Tree model

function y=vesselTree(nodesMatrix,vesselsMatrix)
mu=0.004; %blood viscosity
pIn=7000; %artery border condition
pOut=2000; %venous border condition
el=3000;% Elasticity

vesselsN=148;
nodesN=98;
vesselEq=0;
nodeEq=1;


tspan = [0 100];
for i=1:1:(nodesN*nodeEq+vesselsN*vesselEq)
   yzero(1,i)=1; 
end


sol=ode23s(@diff_equ, tspan, yzero);


t = linspace(0, 100, 150000);

y = deval(sol,t);  
format long

out(1)=pIn;
i=2;
ii=1;
while ii~=-1
    out(i)=  y(ii,900);
    if (nodesMatrix(ii,2)==0)
        ii= find (nodesMatrix(:,1) == nodesMatrix(ii,4), 1, 'first');
    else
        ii= find (nodesMatrix(:,1) == nodesMatrix(ii,3), 1, 'first');  
    end
    i=i+1;
end

plot (out)


%-------------------------------
function yprime = diff_equ(t,y)
 for nodeIdx=1:1:nodesN  
     
     nodeType=nodesMatrix(nodeIdx,2);
     
     
     cv1=nodesMatrix(nodeIdx,6); %connected vessel1
     cv1idx=find(vesselsMatrix(:,1)==cv1,1,'first');
     
     
     length1=vesselsMatrix(cv1idx,3);
     radius1=vesselsMatrix(cv1idx,4);
     
     cv2=nodesMatrix(nodeIdx,7); %connected vessel1
     cv2idx=find(vesselsMatrix(:,1)==cv2,1,'first');
     
     length2=vesselsMatrix(cv2idx,3);
     radius2=vesselsMatrix(cv2idx,4);
     
     cv3=nodesMatrix(nodeIdx,8); %connected vessel1
     cv3idx=find(vesselsMatrix(:,1)==cv3,1,'first');
     
     length3=vesselsMatrix(cv3idx,3);
     radius3=vesselsMatrix(cv3idx,4);

     
     cond1=3.141592*(radius1.^4)/(8*mu*length1);
     cond2=3.141592*(radius2.^4)/(8*mu*length2);
     cond3=3.141592*(radius3.^4)/(8*mu*length3);
     
     cn1=nodesMatrix(nodeIdx,3);
     cn1idx=find(nodesMatrix(:,1)==cn1,1,'first');
      
     cn2=nodesMatrix(nodeIdx,4);
     cn2idx=find(nodesMatrix(:,1)==cn2,1,'first');
     
     cn3=nodesMatrix(nodeIdx,5);
     cn3idx=find(nodesMatrix(:,1)==cn3,1,'first');
     

     
    if nodeType==0  %arterie side
        if (cn1==-1)
        yprime(nodeIdx*nodeEq+0,1)=((pIn-y(nodeIdx*nodeEq+0))*cond1-(y(nodeIdx*nodeEq+0)-y(cn2idx*nodeEq+0))*cond2-(y(nodeIdx*nodeEq+0)-y(cn3idx*nodeEq+0))*cond3)/el;
        else
        yprime(nodeIdx*nodeEq+0,1)=((y(cn1idx*nodeEq+0)-y(nodeIdx*nodeEq+0))*cond1-(y(nodeIdx*nodeEq+0)-y(cn2idx*nodeEq+0))*cond2-(y(nodeIdx*nodeEq+0)-y(cn3idx*nodeEq+0))*cond3)/el;    
        end
    else            %venous side
       if (cn1==-1)
        yprime(nodeIdx*nodeEq+0,1)= ((pOut-y(nodeIdx*nodeEq+0))*cond1-(y(nodeIdx*nodeEq+0)-y(cn2idx*nodeEq+0))*cond2-(y(nodeIdx*nodeEq+0)-y(cn3idx*nodeEq+0))*cond3)/el;
        else
        yprime(nodeIdx*nodeEq+0,1)= ((y(cn1idx*nodeEq+0)-y(nodeIdx*nodeEq+0))*cond1-(y(nodeIdx*nodeEq+0)-y(cn2idx*nodeEq+0))*cond2-(y(nodeIdx*nodeEq+0)-y(cn3idx*nodeEq+0))*cond3)/el;    
        end        
    end
    

    
 end
 
 %for kk=1:1:vesselsN
 %% Vessels equations    
 %end

end

end