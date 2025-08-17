function [Uh_adp,Uh_opt,y_delta,alpha,error_adp,error_opt,alpha_optmin]=AVIEmain_EX123(delta,h,noise,fun_idx,u_exact,k,m)
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

%%%%compute y_delta
y=computey(u_exact,k,T);
u=u_exact(T);
u0=u_exact(0);
delta_vec=(delta1/2+(delta1/2)*rand(N*m,1)).*(u0-u)*0.5;
y_delta=y-delta_vec;
c_min=0.01;c_max=1;
x=1;alpha=linspace(1*sqrt(delta),1*sqrt(delta),N*m);
%%%%choose vaiable regularization parameters
%%%%%%%Practical estimation
if fun_idx==3
    alpha=linspace(0.5*sqrt(delta),0.5*sqrt(delta),N*m);
end
if (fun_idx == 1) || (fun_idx == 3)
    if m==1
        u_priori = collo_solve1(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    elseif m==2
        u_priori = collo_solve2(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    end
elseif fun_idx == 2
    if m==1
        u_priori = collo_solve1k(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    elseif m==2
        u_priori = collo_solve2k(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    end
end
alpha=-delta_vec./(u_priori-u0);
%%%%%%%Idealized estimation
%alpha=-delta_vec./(u-u0);
%y_delta=y+alpha.*(u_priori-u0)+1*noise*(-1+2*rand(N*m,1));
alpha_min=c_min*sqrt(delta);
alpha_max=c_max*sqrt(delta);
alpha=max(alpha,alpha_min);
alpha=min(alpha,alpha_max);
%%% Solution with variable regularization parameter
if (fun_idx == 1) || (fun_idx == 3)
    if m==1
        Uh_adp = collo_solve1(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    elseif m==2
        Uh_adp = collo_solve2(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    end
elseif fun_idx == 2
    if m==1
        Uh_adp = collo_solve1k(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    elseif m==2
        Uh_adp = collo_solve2k(y_delta, u_exact,k, N, c, d, alpha, delta, x);
    end
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
    if (fun_idx == 1) || (fun_idx == 3)
        if m==1
            Uh_optvec(:, i) = collo_solve1(y_delta, u_exact,k, N, c, d, alpha_opt(:, i), delta, x);
        elseif m==2
            Uh_optvec(:, i) = collo_solve2(y_delta, u_exact,k, N, c, d, alpha_opt(:, i), delta, x);
        end      
    elseif fun_idx==2
        if m==1
            Uh_optvec(:, i) = collo_solve1k(y_delta, u_exact,k, N, c, d, alpha_opt(:, i), delta, x);
        elseif m==2
            Uh_optvec(:, i) = collo_solve2k(y_delta, u_exact,k, N, c, d, alpha_opt(:, i), delta, x);
        end    
    end
    %%%%Compute error norms for fixed regularization
    [error_opt_2(i), error_opt_inf(i)] = normcompute(u_exact, Uh_optvec(:, i), h, c, fun_idx);  
    error_opt_grid(i) = norm(u - Uh_optvec(:, i), inf);
end


%%%%% Compute error norms for fixed regularization
[min_erroropt2, min_erroropt2_index] = min(error_opt_2);  % Find minimum L2 error and its index
error_opt=[min_erroropt2,error_opt_inf(min_erroropt2_index),error_opt_grid(min_erroropt2_index)];
alpha_optmin=alpha_opt(1,min_erroropt2_index);
Uh_opt=Uh_optvec(:, min_erroropt2_index);
%%%%% Compute error norms for variable regularization
[erroradp2, erroradpinf] = normcompute(u_exact, Uh_adp, h, c, fun_idx);  
erroradpgrid = norm(u - Uh_adp, inf);
error_adp=[erroradp2,erroradpinf,erroradpgrid];

