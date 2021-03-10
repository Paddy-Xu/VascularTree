%Sniffers,buzzers,toggles and blinkers: dynamics of regulatory and
%signaling patways in the cell
%John J Tyson e t.c.
%Negative-feedback oscillator

function [t,y]=nephronTree(nodesMatrix,vesselsMatrix,nephronsMatrix)

%Vessel tree parameters
pIn=14; %artery border condition
el=3000;% Elasticity

vesselsN=7;
vesselEq=0;

nodesN=3;
nodeEq=1;

nephronsN=4;
nephronEq=6;


for j=1:1:vesselsN;
    length = vesselsMatrix(j,3);
    diameter = vesselsMatrix(j,4);
    visc=(1+(6*exp(-0.085*diameter)+2.2-2.44*exp(-0.06*(diameter^0.645)))*((diameter/(diameter-1.1))^2))*((diameter/(diameter-1.1))^2)/1000;
    resistance(j)=(128*visc*length)/(pi*(diameter^4));
end


% single nephron variable staff =========================
%this parameters we do change
T=(13.5);	%Henle's loop delay etc. [s]%
alpha=(20.0);	%Transglomerular feedback amplification [1]%
Pref=(1.3);
Ptau=(.1);
Pscale=(.035);
PinNeph=0; %Arterial blood preasure [kPa]%

%this are usually constants
Ctub=(3.0);	%Elasticity [nL/kPa]%
Ha=(0.5);	%Arterial hematochrite (?) [1]%
Pe=(1.3);	%Efferent arterial pressure [kPa]%
Pd=(0.6);	%CDistale tubulus hydrostatic pressure [kPa]%
Fhen0=(0.2);	%Henle's loop equilibrium flow [nL/s]%
Freab=(0.3);	%Proximate tubulus reabsorption [nL/s]%
Rhen=(5.3);	%Henle's loop resistance [kPa*s/nL]%

Ra0=(2.3);	%Afferent arterial equilibrium resistance [kPa*s/nL]%

Re=(1.9);	%Efferent arterial resistance [kPa*s/nL]%
omega=(20.0);	%Damped oscillator parameter [kPa*s^2]%
d=(0.04);	%! DISCREPANCY - Damped oscillator parameter [1/s]%
beta=(0.67);	%Non-variable fraction of efferent arterial [1]%
psimin=(0.20);	%Lower activation limit [1]%
psimax=(0.44);	%Upper activation limit [1]%
psieq=(0.38);	%Equilibrium activation [1]%
Cvasc=(3.0);	%Elasticity [nL/kPa]%




% Parameters used only in cubic equation, therefore units [L] and [g] are OK%
Ca=(54.0);	%Afferent plasma protein concentration [g/L]%
a=(22.0e-3);	%Protein concentration parameter [kPa*L/g]%
b=(0.39e-3);	%Protein concentration parameter [kPa*L^2/g^2]%


tspan = [0 1000];

for i=1:1:(nodesN*nodeEq)
    yzero(i)=1;
end
for i=nodesN*nodeEq+1:nephronEq:(nodesN*nodeEq+nephronsN*nephronEq)
    yzero(i)=2.07828 - (i-nodesN*nodeEq-1)/nephronEq;
    yzero(i+1)=1.1891- (i-nodesN*nodeEq-1)/nephronEq;
    yzero(i+2)= -0.186993- (i-nodesN*nodeEq-1)/nephronEq;
    yzero(i+3)= 1.44853- (i-nodesN*nodeEq-1)/nephronEq;
    yzero(i+4)= 1.25917- (i-nodesN*nodeEq-1)/nephronEq;
    yzero(i+5)= 1.0752- (i-nodesN*nodeEq-1)/nephronEq;
end


psiMatrix=preproc(nodesMatrix,vesselsMatrix,nephronsMatrix);


sol=ode23s(@diff_equ, tspan, yzero);
t = linspace(0, 1000, 15000);
y = deval(sol,t);

plot(t,y(13,:))
figure
plot(t,y(1,:),t,y(7,:))



% Pel & Pact yields [kPa]%
    function result=Pel(x)
        result=(1.6*((x)-1)+6e-3*exp(10*((x)-0.8)));
    end

    function result=Pact(x)
        result=(6.3+7.2*(x)+4.7/(1+exp(13*(0.4-(x)))));
    end

    function [x,y]=Impow13(x,y)
        impowv=atan2(y,x)/3.0;
        impowr=hypot(x,y)^(1.0/3.0);
        y=sin(impowv)*impowr;
        x=cos(impowv)*impowr;
    end

    function result=Cubic(s,p,q,r)  % This function is mostly CRAP. 
        %Works only in our case (but taken from original works and used in a lot of papers). 
        %Better rewrite it when i have time
        i=2;
        if (s~=0.0)
            p=p/s;
            q=q/s;
            r=r/s;
            aa=(3.0*q-p*p)/3.0;
            bb=(2.0*p*p*p-9*p*q+27.0*r)/27.0;
            dummy=bb*bb/4.0+aa*aa*aa/27.0;
            if (dummy>=0.0)
                A=(-bb)/2.0+sqrt(dummy);
                Aim=0.0;
                B=(-bb)/2.0-sqrt(dummy);
                Bim=0.0;
                if (A<0.0)
                    A=0.0-(0.0-A)^(1.0/3.0);
                else
                    A=A^(1.0/3.0);
                end
                if (B<0.0)
                    B=0.0-(0.0-B)^(1.0/3.0);
                else
                    B=B^(1.0/3.0);
                end
            else
                A=(-bb)/2.0;
                Aim=sqrt(0.0-dummy);
                B=(-bb)/2.0;
                Bim=0.0-sqrt(0.0-dummy);
                [A,Aim]=Impow13(A,Aim);
                [B,Bim]=Impow13(B,Bim);
            end
            %Commetned imres part since it is not used at all
            %        if (imres~=0)
            %            imres(0)=Aim+Bim;
            %            imres(1)=(-(Aim+Bim))/2.0+sqrt(3.0)*(A-B)/2.0;
            %            imres(i)=(-(Aim+Bim))/2.0-sqrt(3.0)*(A-B)/2.0;
            %        end
            res(1)=A+B-p/3;
            res(2)=(-(A+B))/2.0-sqrt(3.0)*(Aim-Bim)/2.0-p/3;
            res(3)=(-(A+B))/2.0+sqrt(3.0)*(Aim-Bim)/2.0-p/3;
            result=res(1);
        end
    end

    function yPrime = diff_equ(t,y)
        yPrime=zeros(nephronsN*nephronEq+nodesN*nodeEq,1);
        previousEq=0;
        
        for k=1:1:nephronsN
            psi(k)=psimax - (psimax - psimin)/ (1.0 + (psieq-psimin)/(psimax-psieq)*exp(alpha*(3*y(previousEq+(k-1)*nephronEq+6)/T/Fhen0 - 1.0)));
        end
        
        for k=1:1:nephronsN
            
            
            cnidx = find ( (nodesMatrix(:,9) == nephronsMatrix(k,1) | nodesMatrix(:,10) == nephronsMatrix(k,1) ), 1, 'first');
            
            if(nodesN>0)
                PinNeph= y(nephronsN*nephronEq+(cnidx-1)*nodeEq+1,1); % Afferent arteriols in the matrix are not used, it means all of them are the same! WRONG!
                cvidx = find ( (vesselsMatrix(:,1) == nephronsMatrix(k,1)), 1, 'first');
                Ra0=resistance(cvidx)*1000;
                % convert to kPa/s;
            end
            
            psiLocal=0;
            for kk=1:1:nephronsN
                psiLocal= psiLocal+psi(kk)*psiMatrix(k,kk);
                
            end
            
            
            pEq=Pel(y(previousEq+(k-1)*nephronEq+2))+psiLocal*Pact(y(previousEq+(k-1)*nephronEq+2));
            rA=Ra0*(beta+(1.0-beta)/((y(previousEq+(k-1)*nephronEq+2))^4));
            r=Re/rA;
            
            tmpA=b+r*b*Ha;
            tmpB=a+r*b*Ca*(1-Ha)+r*a*Ha;
            tmpC=y(previousEq+(k-1)*nephronEq+1)-Pe+r*a*Ca*(1.0-Ha)+r*(y(previousEq+(k-1)*nephronEq+1)-PinNeph)*Ha;
            tmpD = (y(previousEq+(k-1)*nephronEq+1)-PinNeph)*r*Ca*(1.0-Ha);
            
            ce=Cubic(tmpA, tmpB, tmpC, tmpD);
            
            pG   = b*(ce^2) + a*ce + y(previousEq+(k-1)*nephronEq+1);
            pAv  = (PinNeph - (PinNeph-pG) * beta * Ra0/rA + pG)/2;
            psFun=(1.0-Ha)*(1.0-Ca/ce)*(PinNeph-pG)/rA;
            G=(PinNeph-pG)/rA; %flow in afferent areteriole
            
            if(nodesN>0)
                yPrime(nephronsN*nephronEq+(cnidx-1)*nodeEq+1,1)=yPrime(nephronsN*nephronEq+(cnidx-1)*nodeEq+1,1)-G; %important for tree
            end
            
            yPrime(previousEq+(k-1)*nephronEq+1,1)=1.0/Ctub*( psFun - Freab -(y(previousEq+(k-1)*nephronEq+1)-Pd)/Rhen);
            yPrime(previousEq+(k-1)*nephronEq+2,1)=y(previousEq+(k-1)*nephronEq+3);
            yPrime(previousEq+(k-1)*nephronEq+3,1)=(pAv-pEq)/omega - d*y(previousEq+(k-1)*nephronEq+3);
            yPrime(previousEq+(k-1)*nephronEq+4,1)=(y(previousEq+(k-1)*nephronEq+1)-Pd)/Rhen - 3*y(previousEq+(k-1)*nephronEq+4)/T;
            yPrime(previousEq+(k-1)*nephronEq+5,1)=3*(y(previousEq+(k-1)*nephronEq+4)-y(previousEq+(k-1)*nephronEq+5))/T;
            yPrime(previousEq+(k-1)*nephronEq+6,1)=3*(y(previousEq+(k-1)*nephronEq+5)-y(previousEq+(k-1)*nephronEq+6))/T;
        end
        
        previousEq=nephronsN*nephronEq;
        for k=1:1:nodesN
            cn1idx = find (nodesMatrix(:,1) == nodesMatrix(k,3), 1, 'first');
            cn2idx = find (nodesMatrix(:,1) == nodesMatrix(k,4), 1, 'first');
            cn3idx = find (nodesMatrix(:,1) == nodesMatrix(k,5), 1, 'first');
            
            cv1idx = find (vesselsMatrix(:,1) == nodesMatrix(k,6), 1, 'first');
            cv2idx = find (vesselsMatrix(:,1) == nodesMatrix(k,7),1, 'first');
            cv3idx = find (vesselsMatrix(:,1) == nodesMatrix(k,8),1, 'first');
            
            
            if (nodesMatrix (k,3)==-1)
                cp1=pIn;
            else
                cp1=y(previousEq+(cn1idx-1)*nodeEq+1);
            end
            if (nodesMatrix (k,4)==-1)
                cp2=- y(previousEq+(k-1)*nodeEq+1);
            else
                cp2=y(previousEq+(cn2idx-1)*nodeEq+1);
            end
            if (nodesMatrix (k,5)==-1)
                cp3=- y(previousEq+(k-1)*nodeEq+1);
            else
                cp3= y(previousEq+(cn3idx-1)*nodeEq+1);
            end
            
            
            
            yPrime (previousEq+(k-1)*nodeEq+1,1) =yPrime (previousEq+(k-1)*nodeEq+1,1) +  ((cp1 - y(previousEq+(k-1)*nodeEq+1)) /resistance(cv1idx) - (y(previousEq+(k-1)*nodeEq+1) - cp2)/resistance(cv2idx) - (y(previousEq+(k-1)*nodeEq+1) - cp3)/resistance(cv3idx))/el;
            
        end
    end


end

