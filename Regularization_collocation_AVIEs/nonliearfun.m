function F=nonliearfun(u,h,alpha,U0,y1,y2,a111,a112,a121,a122,a211,a212,a221,a222)
    F(1)=y1+alpha(1)*(U0-u(1))-h*(a111*u(1)*u(1)+a112*u(1)*u(2)+a121*u(2)*u(1)+a122*u(2)*u(2));
    F(2)=y2+alpha(2)*(U0-u(2))-h*(a211*u(1)*u(1)+a212*u(1)*u(2)+a221*u(2)*u(1)+a222*u(2)*u(2));
end