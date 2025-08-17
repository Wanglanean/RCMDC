function [Uh_adp,Uh_opt,y_delta,alpha,error_adp,error_opt,alpha_optmin]=AVIEmain_Ex4(delta,h,noise,fun_idx,m)
delta1=delta-noise;
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

%%%%%%%%%%%%
%%%piecewise constant function and kernel function
u_exact=@(t) (t>=0 & t<=0.5).*0.5+(t>0.5 & t<=0.8).*0.25+(t>0.8 & t<=1).*0.75;
k=@(t,s) 1;

%%%%compute y_delta
y=computey(u_exact,k,T);
u=u_exact(T);
u0=u_exact(0);
delta_vec=zeros(N*m,1);
y_delta=y+delta_vec;
%%%%choose variable regularization parameters
c_min=0.1;c_max=1;
idx = (T >= 0.45 & T <= 0.55) | (T >= 0.75 & T <= 0.85); % t=0.5 或 t=0.8 附近
delta_vec(idx) =(delta1/2+ delta1/2 * rand(sum(idx), 1) ).* sign(u(idx) - u0);
c_mid=(c_min+c_max)/3.5;
alpha=linspace(c_mid*sqrt(delta),c_mid*sqrt(delta),N*m)';
%%%%%%%%Idealized estimation
%alpha(idx)=delta_vec(idx)./(u(idx)-u0);
%%%%%%%%Practical estimation
x=1;
if m==1
    u_priori=collo_solve1(y_delta,u_exact,k,N,c,d,alpha,delta,x);
elseif m==2
    u_priori=collo_solve2(y_delta,u_exact,k,N,c,d,alpha,delta,x);
end
alpha(idx)=delta_vec(idx)./(u_priori(idx)-u0);
alpha_min=c_min*sqrt(delta);
alpha_max=c_max*sqrt(delta);
alpha=max(alpha,alpha_min);
alpha=min(alpha,alpha_max);
y_delta=y+alpha.*(u-u0)+noise*(-1+2*rand(N*m,1));
%%% Solution with variable regularization parameter
if m==1
    Uh_adp=collo_solve1(y_delta,u_exact,k,N,c,d,alpha,delta,x);
elseif m==2
    Uh_adp=collo_solve2(y_delta,u_exact,k,N,c,d,alpha,delta,x);
end
%%% Solution with fixed regularization parameter
nopt=20;
alpha_opt = zeros(N*m, nopt);
for i = 1:N*m
    alpha_opt(i,:) = linspace(c_min*sqrt(delta),c_max*sqrt(delta),nopt);  
end 
error_opt_2 = zeros(1, nopt); error_opt_inf = zeros(1, nopt); error_opt_grid = zeros(1, nopt); 
Uh_optvec = zeros(N*m, nopt);
for i = 1:nopt
    %%%%Solve using fixed regularization
    if m==1
        Uh_optvec(:, i) = collo_solve1(y_delta, u_exact,k, N, c, d, alpha_opt(:, i), delta, x);
    elseif m==2
        Uh_optvec(:, i) = collo_solve2(y_delta, u_exact,k,N, c, d, alpha_opt(:, i), delta, x);
    end      
    %%%%Compute error norms for fixed regularization
    [error_opt_2(i), error_opt_inf(i)] = normcompute(u_exact, Uh_optvec(:, i), h, c, fun_idx);  
    error_opt_grid(i) = norm(u - Uh_optvec(:, i), inf);
end
%%%%%%%%%%%%%%%%%%%%%
%%%%% Compute error norms for fixed regularization
[min_erroropt2, min_erroropt2_index] = min(error_opt_2);  % Find minimum L2 error and its index
error_opt=[min_erroropt2,error_opt_inf(min_erroropt2_index),error_opt_grid(min_erroropt2_index)];
alpha_optmin=alpha_opt(1,min_erroropt2_index);
Uh_opt=Uh_optvec(:, min_erroropt2_index);
%%%%% Compute error norms for variable regularization
[erroradp2, erroradpinf] = normcompute(u_exact, Uh_adp, h, c, fun_idx);  
erroradpgrid = norm(u - Uh_adp, inf);
error_adp=[erroradp2,erroradpinf,erroradpgrid];

