function [Uh_adp,Uh_opt,y_delta,alpha,error_adp,error_opt,alpha_optmin]=AVIEmain2(delta,h,fun_idx,u_exact,k,m)
%%%%collocation points
if m==1
    c=1/2;
elseif m==2
    c=[(3-sqrt(3))/6,(3+sqrt(3))/6];
end
d=m;

%%%%mesh
N=round(1/h);h=1/N;
T=zeros(N*m,1);
for n=0:N-1
    for j=1:m
        T(m*n+j)=n*h+c(j)*h;
    end
end

%%%%compute y_delta
y=computey(u_exact,k,T);
u=u_exact(T);
delta_vec=delta*(-1+2*rand(N*m,1));
y_delta=y+delta_vec;
%%%%choose variable regularization parameters
c_min=1;c_max=2;
alpha_min=c_min*sqrt(delta);
alpha_max=c_max*sqrt(delta);
alpha=linspace(alpha_min,alpha_max,N*m);
%%% Solution with variable regularization parameter
x=0;
if (fun_idx == 1 || fun_idx == 3 || fun_idx == 4)
    if m==1
        Uh_adp = collo_solve1(y_delta, u_exact, k, N, c, d, alpha, delta, x);
    elseif m==2
        Uh_adp = collo_solve2(y_delta, u_exact, k, N, c, d, alpha, delta, x);
    end
elseif fun_idx == 2
    if m==1
        Uh_adp = collo_solve1k(y_delta, u_exact, k, N, c, d, alpha, delta, x);
    elseif m==2
        Uh_adp = collo_solve2k(y_delta, u_exact, k, N, c, d, alpha, delta, x);
    end
end

%%% Solution with fixed regularization parameter
alpha_opt=linspace(c_min*sqrt(delta),c_min*sqrt(delta),N*m);alpha_optmin=c_min*sqrt(delta);
if (fun_idx == 1 || fun_idx == 3 || fun_idx == 4)
     if m==1
          Uh_opt = collo_solve1(y_delta, u_exact, k, N, c, d, alpha_opt, delta, x);
     elseif m==2
          Uh_opt = collo_solve2(y_delta,u_exact, k, N, c, d, alpha, delta, x);
     end      
elseif fun_idx==2
     if m==1
          Uh_opt = collo_solve1k(y_delta, u_exact, k, N, c, d, alpha_opt, delta, x);
     elseif m==2
         Uh_opt = collo_solve2k(y_delta, u_exact, k, N, c, d, alpha_opt, delta, x);
     end    
end
%%%%% Compute error norms for fixed regularization
[erroropt2, erroroptinf] = normcompute(u_exact, Uh_opt, h, c, fun_idx);  
erroroptgrid = norm(u - Uh_opt, inf);
error_opt=[erroropt2,erroroptinf,erroroptgrid];
%%%%% Compute error norms for variable regularization
[erroradp2, erroradpinf] = normcompute(u_exact, Uh_adp, h, c, fun_idx);  
erroradpgrid = norm(u - Uh_adp, inf);
error_adp=[erroradp2,erroradpinf,erroradpgrid];
end