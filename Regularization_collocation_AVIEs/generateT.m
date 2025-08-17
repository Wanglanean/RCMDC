function [T1,T2,T3,T4,T5]=generateT(h,c)
N=round(1./h);
m=size(c,2);
T1=zeros(round(N(1)),1);T2=zeros(round(N(2)),1);T3=zeros(round(N(3)),1);
T4=zeros(N(4),1);T5=zeros(N(5),1);
for n=0:N(1)-1
    for j=1:m
        T1(m*n+j)=n*h(1)+c(j)*h(1);
    end
end
for n=0:N(2)-1
    for j=1:m
        T2(m*n+j)=n*h(2)+c(j)*h(2);
    end
end
for n=0:N(3)-1
    for j=1:m
        T3(m*n+j)=n*h(3)+c(j)*h(3);
    end
end
for n=0:N(4)-1
    for j=1:m
        T4(m*n+j)=n*h(4)+c(j)*h(4);
    end
end
for n=0:N(5)-1
    for j=1:m
        T5(m*n+j)=n*h(5)+c(j)*h(5);
    end
end