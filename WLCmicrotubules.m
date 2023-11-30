function [wlcseries]=WLCmicrotubules(forces,KT,A,l,steptot)

% fcand=[0,0.1,1,10,100];
fcand=forces;
DNAseriestot=[];

for ff=1:1:size(fcand,2)
%Parameters
ff
f=fcand(ff); %in pN
probmax=exp(f*l*1/KT);

% DNA streching modeling
DNAt=[];
DNAt(1,:)=[0 0 0];
inirnd=rand(1,3)*2-1;
DNAt(2,:)=inirnd./sqrt(sum(inirnd.^2));
DNAseries=DNAt;
indx=3;
for tt=1:1:10000000  % 100 steps are used which means 10 amino acid per segment
    dirtemp=(2*rand(3,1)-1);
    direction=dirtemp/sqrt(sum(dirtemp.^2));
    costheta=direction(3);
    phi2=2*(1-DNAt(indx-1,:)*direction);
    prob=exp(f*l*costheta/KT-A/2/l*phi2);
    y=rand*probmax;
    if y<prob
        DNAt=cat(1,DNAt,direction');
        DNAseries=cat(1,DNAseries,DNAseries(indx-1,:)+direction');
        indx=indx+1;
    else
        continue
    end
    
    if indx>steptot
        break
    end
end
DNAseriestot=cat(3,DNAseriestot,DNAseries);

end

wlcseries=DNAseriestot;

    