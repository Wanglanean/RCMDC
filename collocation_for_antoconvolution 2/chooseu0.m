function U0=chooseu0(u_exact,delta,h,d,x)
if x==0
    U0=u_exact(0)+rand*(sqrt(delta)+h^(d/2));
elseif x==1
    U0=u_exact(0);
end
end