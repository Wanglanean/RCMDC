clc,clear
syms t s t n c c1 c2 h k l v
% k=@(t,s) 1;
% k=@(t,s) t-2*s+1;
%%%m=1
% a0=int(k(c*h,v*h),v,0,c)
% a=int(k(n*h+c*h,l*h+v*h),v,0,c)
% b=int(k(n*h+c*h,l*h+v*h),v,c,1)
% an=int(k(n*h+c*h,n*h+v*h),v,0,c)
%%%%m=2
% l1=@(s)((s-c2)/(c1-c2));l2=@(s)((s-c1)/(c2-c1));
%%%系数a
%n=0
% a111=int(k(c1*h,s*h)*l1(c1-s)*l1,s,0,c1)
% a112=int(k(c1*h,s*h)*l1(c1-s)*l2,s,0,c1)
% a121=int(k(c1*h,s*h)*l2(c1-s)*l1,s,0,c1)
% a122=int(k(c1*h,s*h)*l2(c1-s)*l2,s,0,c1)
% a211=int(k(c2*h,s*h)*l1(c2-s)*l1,s,0,c2)
% a212=int(k(c2*h,s*h)*l1(c2-s)*l2,s,0,c2)
% a221=int(k(c2*h,s*h)*l2(c2-s)*l1,s,0,c2)
% a222=int(k(c2*h,s*h)*l2(c2-s)*l2,s,0,c2)
%%n>0
% a111=int(k(n*h+c1*h,l*h+s*h)*l1(c1-s)*l1,s,0,c1)
% a112=int(k(n*h+c1*h,l*h+s*h)*l1(c1-s)*l2,s,0,c1)
% a121=int(k(n*h+c1*h,l*h+s*h)*l2(c1-s)*l1,s,0,c1)
% a122=int(k(n*h+c1*h,l*h+s*h)*l2(c1-s)*l2,s,0,c1)
% a211=int(k(n*h+c2*h,l*h+s*h)*l1(c2-s)*l1,s,0,c2)
% a212=int(k(n*h+c2*h,l*h+s*h)*l1(c2-s)*l2,s,0,c2)
% a221=int(k(n*h+c2*h,l*h+s*h)*l2(c2-s)*l1,s,0,c2)
% a222=int(k(n*h+c2*h,l*h+s*h)*l2(c2-s)*l2,s,0,c2)
%%l=n
% a111=int(k(n*h+c1*h,n*h+s*h)*l1(c1-s)*l1,s,0,c1)
% a112=int(k(n*h+c1*h,n*h+s*h)*l1(c1-s)*l2,s,0,c1)
% a121=int(k(n*h+c1*h,n*h+s*h)*l2(c1-s)*l1,s,0,c1)
% a122=int(k(n*h+c1*h,n*h+s*h)*l2(c1-s)*l2,s,0,c1)
% a211=int(k(n*h+c2*h,n*h+s*h)*l1(c2-s)*l1,s,0,c2)
% a212=int(k(n*h+c2*h,n*h+s*h)*l1(c2-s)*l2,s,0,c2)
% a221=int(k(n*h+c2*h,n*h+s*h)*l2(c2-s)*l1,s,0,c2)
% a222=int(k(n*h+c2*h,n*h+s*h)*l2(c2-s)*l2,s,0,c2)
%%%%系数b
%%n=0
% b111=int(l1(c1-s+1)*l1,s,c1,1)
% b112=int(l1(c1-s+1)*l2,s,c1,1)
% b121=int(l2(c1-s+1)*l1,s,c1,1)
% b122=int(l2(c1-s+1)*l2,s,c1,1)
% b211=int(l1(c2-s+1)*l1,s,c2,1)
% b212=int(l1(c2-s+1)*l2,s,c2,1)
% b221=int(l2(c2-s+1)*l1,s,c2,1)
% b222=int(l2(c2-s+1)*l2,s,c2,1)

% b111=int(k(n*h+c1*h,l*h+s*h)*l1(c1-s+1)*l1,s,c1,1)
% b112=int(k(n*h+c1*h,l*h+s*h)*l1(c1-s+1)*l2,s,c1,1)
% b121=int(k(n*h+c1*h,l*h+s*h)*l2(c1-s+1)*l1,s,c1,1)
% b122=int(k(n*h+c1*h,l*h+s*h)*l2(c1-s+1)*l2,s,c1,1)
% b211=int(k(n*h+c2*h,l*h+s*h)*l1(c2-s+1)*l1,s,c2,1)
% b212=int(k(n*h+c2*h,l*h+s*h)*l1(c2-s+1)*l2,s,c2,1)
% b221=int(k(n*h+c2*h,l*h+s*h)*l2(c2-s+1)*l1,s,c2,1)
% b222=int(k(n*h+c2*h,l*h+s*h)*l2(c2-s+1)*l2,s,c2,1)

% syms Uh1 Uh2 Uh
% %m=1
% u_exact1=@(t) 2+cos(4*pi*t);
% u_exact2=@(t) t.^2-2*t+2;
% u_exact3=@(t) (t>=0 & t<=0.5).*0.5+(t>0.5 & t<=0.8).*0.25+(t>0.8 & t<=1).*0.75;
% error1=int((Uh-u_exact1(n*h+s*h))^2,s,0,1)*h
% error2=int((Uh-u_exact2(n*h+s*h))^2,s,0,1)*h
% error5=int((Uh-u_exact3(n*h+s*h))^2,s,0,1)*h
% %%m=2
% error3=int((Uh1*l1(s)+Uh2*l2(s)-u_exact1(n*h+s*h))^2,s,0,1)*h
% error4=int((Uh1*l1(s)+Uh2*l2(s)-u_exact2(n*h+s*h))^2,s,0,1)*h
% error6=int((Uh1*l1(s)+Uh2*l2(s)-u_exact3(n*h+s*h))^2,s,0,1)*h
