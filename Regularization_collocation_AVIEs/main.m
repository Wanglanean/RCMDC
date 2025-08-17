clc,clear;
c=1/2;
%c=[(3-sqrt(3))/6,(3+sqrt(3))/6];
m=size(c,2);
delta=[0.01,0.05,0.1,0.5,1];
h=[0.0001,0.0005,0.001,0.005,0.01];
%%%exact function
u_exact1=@(t) 2+cos(4*pi*t);
u_exact2=@(t) t.^2-2*t+2;
u_exact3=@(t) (t>=0 & t<=0.25).*(-2*t+1)+(t>0.25 & t<=0.5).*(2*t)+(t>0.5 & t<=0.75).*(-2*t+2)+(t>0.75 & t<=1).*(2*t-1);
u_exact4=@(t) (t>=0 & t<=0.5).*0.5+(t>0.5 & t<=0.8).*0.25+(t>0.8 & t<=1).*0.75;
%%%%kernel function
k1=@(t,s) 1;
k2=@(t,s) t-2*s+1;
fuc_idx1=1;fuc_idx2=2;fuc_idx3=3;fuc_idx4=4;
%%%%Examplel
Nmax=1/h(1)*m;
noise=[0.00001,0.0005,0.001,0.005,0.01];
N=[1/h(1)*m,1/h(2)*m,1/h(3)*m,1/h(4)*m,1/h(5)*m];
Uh1_adp=zeros(Nmax,5);Uh1_opt=zeros(Nmax,5);alpha1=zeros(Nmax,5);y_delta1=zeros(Nmax,5);
error1_adp=zeros(5,3);error1_opt=zeros(5,3);alpha1_optmin=zeros(5,1);
for i=1:5
    [Uh1_adp(1:N(i),i),Uh1_opt(1:N(i),i),y_delta1(1:N(i),i),alpha1(1:N(i),i),error1_adp(i,:),error1_opt(i,:),alpha1_optmin(i)]...
        =AVIEmain_EX123(delta(i),h(i),noise(i),fuc_idx1,u_exact1,k1,m);
end
%%%%Example2
Uh2_adp=zeros(Nmax,5);Uh2_opt=zeros(Nmax,5);alpha2=zeros(Nmax,5);y_delta2=zeros(Nmax,5);
error2_adp=zeros(5,3);error2_opt=zeros(5,3);alpha2_optmin=zeros(5,1);
for i=1:5
    [Uh2_adp(1:N(i),i),Uh2_opt(1:N(i),i),y_delta2(1:N(i),i),alpha2(1:N(i),i),error2_adp(i,:),error2_opt(i,:),alpha2_optmin(i)]...
        =AVIEmain_EX123(delta(i),h(i),noise(i),fuc_idx2,u_exact2,k2,m);
end
%%%%Example3
Uh3_adp=zeros(Nmax,5);Uh3_opt=zeros(Nmax,5);alpha3=zeros(Nmax,5);y_delta3=zeros(Nmax,5);
error3_adp=zeros(5,3);error3_opt=zeros(5,3);alpha3_optmin=zeros(5,1);
for i=1:5
    [Uh3_adp(1:N(i),i),Uh3_opt(1:N(i),i),y_delta3(1:N(i),i),alpha3(1:N(i),i),error3_adp(i,:),error3_opt(i,:),alpha3_optmin(i)]...
        =AVIEmain_EX123(delta(i),h(i),noise(i),fuc_idx3,u_exact3,k1,m);
end
%%%%Example4
Uh4_adp=zeros(Nmax,5);Uh4_opt=zeros(Nmax,5);alpha4=zeros(Nmax,5);y_delta4=zeros(Nmax,5);
error4_adp=zeros(5,3);error4_opt=zeros(5,3);alpha4_optmin=zeros(5,1);
for i=1:5
    [Uh4_adp(1:N(i),i),Uh4_opt(1:N(i),i),y_delta4(1:N(i),i),alpha4(1:N(i),i),error4_adp(i,:),error4_opt(i,:),alpha4_optmin(i)]...
        =AVIEmain_Ex4(delta(i),h(i),noise(i),fuc_idx4,m);
end
%%%%%%%figure
set(groot, 'DefaultAxesFontSize', 12);  
set(groot, 'DefaultTextFontSize', 12);  
set(groot, 'DefaultAxesLineWidth', 0.8);
[T1,T2,T3,T4,T5]=generateT(h,c);N=round(1./h)*m;
figure(1)
plot(T1, u_exact1(T1), "k-", T1, Uh1_adp(1:N(1),1), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)
print('Ex1m2delta1alpha', '-dpng', '-r600'); 
figure(2)
plot(T3, u_exact1(T3), "k-", T3, Uh1_adp(1:N(3),3), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)
print('Ex1m2delta2alpha', '-dpng', '-r600'); 
figure(3)
plot(T5, u_exact1(T5), "k-", T5, Uh1_adp(1:N(5),5), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)  
print('Ex1m2delta3alpha', '-dpng', '-r600'); 
figure(4)
plot(T1, u_exact1(T1), "k-", T1, Uh1_opt(1:N(1),1), "k--") 
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)
print('Ex1m2delta1alpha_fix', '-dpng', '-r600');
%%%%%%%
figure(5)
plot(T1, u_exact2(T1), "k-", T1, Uh2_adp(1:N(1),1), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)
print('Ex2m2delta1alpha', '-dpng', '-r600'); 
figure(6)
plot(T3, u_exact2(T3), "k-", T3, Uh2_adp(1:N(3),3), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)
print('Ex2m2delta2alpha', '-dpng', '-r600'); 
figure(7)
plot(T5, u_exact2(T5), "k-", T5, Uh2_adp(1:N(5),5), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12)  
print('Ex2m2delta3alpha', '-dpng', '-r600'); 
figure(8)
plot(T1, u_exact2(T1), "k-", T1, Uh2_opt(1:N(1),1), "k--") 
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12) 
print('Ex2m2delta1alpha_fix', '-dpng', '-r600');
%%%%%%%%
figure(9)
plot(T1, u_exact3(T1), "k-", T1, Uh3_adp(1:N(1),1), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest')
print('Ex3m2delta1alpha', '-dpng', '-r600'); 
figure(10)
plot(T3, u_exact3(T3), "k-", T3, Uh3_adp(1:N(3),3), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest')
print('Ex3m2delta2alpha', '-dpng', '-r600'); 
figure(11)
plot(T5, u_exact3(T5), "k-", T5, Uh3_adp(1:N(5),5), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest')  
print('Ex3m2delta3alpha', '-dpng', '-r600'); 
figure(12)
plot(T1, u_exact3(T1), "k-", T1, Uh3_opt(1:N(1),1), "k--") 
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest') 
print('Ex3m2delta1alpha_fix', '-dpng', '-r600');
%%%%%%%%
figure(13)
plot(T1, u_exact4(T1), "k-", T1, Uh4_adp(1:N(1),1), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest')
print('Ex4m2delta1alpha', '-dpng', '-r600'); 
figure(14)
plot(T3, u_exact4(T3), "k-", T3, Uh4_adp(1:N(3),3), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest')
print('Ex4m2delta2alpha', '-dpng', '-r600'); 
figure(15)
plot(T5, u_exact4(T5), "k-", T5, Uh4_adp(1:N(5),5), "k--")
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest')  
print('Ex4m2delta3alpha', '-dpng', '-r600'); 
figure(16)
plot(T1, u_exact4(T1), "k-", T1, Uh4_opt(1:N(1),1), "k--") 
legend('$$u^{\dagger}$$','$$u_h^{\alpha,\delta}$$','Interpreter', 'latex', 'FontSize', 12,'Location', 'northwest') 
print('Ex4m2delta1alpha_fix', '-dpng', '-r600');
%%%convergence rate
r1_adp = showratenew(delta, error1_adp);
r1_opt = showratenew(delta, error1_opt);
r2_adp = showratenew(delta, error2_adp);
r2_opt = showratenew(delta, error2_opt);
r3_adp = showratenew(delta, error3_adp);
r3_opt = showratenew(delta, error3_opt);
r4_adp = showratenew(delta, error4_adp);
r4_opt = showratenew(delta, error4_opt);