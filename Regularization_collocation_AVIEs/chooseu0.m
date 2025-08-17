function U0=chooseu0(y_delta,u_exact,k,delta,m,h,c,d,x)
Deltat=sqrt(delta+h^d);
while Deltat>1
    Deltat=0.5*Deltat;
end
if x==0
    l=fix(Deltat/h);v=(Deltat-l*h)/h;
    if v==0
        l=0;v=1;
    end
    if size(c,2)==1
        temp=y_delta(l*m+1)/(abs(k(0,0))*Deltat);
    elseif size(c,2)==2
        l1=@(s)((s-c(2))/(c(1)-c(2)));l2=@(s)((s-c(1))/(c(2)-c(1)));
        temp=(l1(v)*y_delta(l*m+1)+l2(v)*y_delta(l*m+2))/(abs(k(0,0))*Deltat);
    end
    if temp>0
        U0=sqrt(temp);
    else
        U0=0;
    end
elseif x==1
    U0=u_exact(0);
else
    %c=0.01;
    U0=u_exact(0)+rand*(sqrt(delta)+h^(d/2));
end
end