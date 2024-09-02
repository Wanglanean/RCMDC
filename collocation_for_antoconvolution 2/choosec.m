function [c_min,c_max]=choosec(C9,C10,M_line,Kd,normU,m,eta,eps2,delta)
b1=(C10*eta+sqrt(C10^2*eta^2+4*eps2*C9))/(2*eps2);
b2=(sqrt(32*m^3)*M_line*Kd*normU)/sqrt(delta);
c_min=max(b1,b2)+1;
c_max=eta*c_min;
end