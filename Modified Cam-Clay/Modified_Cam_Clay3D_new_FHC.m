%By Fang Huangcheng @PolyU
%Implicit return algorithm; Backward Euler Method
%Last update @2024/1/31
%Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
%Ref:KRISTIAN KRABBENHØFT, 2002, BASIC COMPUTATIONAL PLASTICITY
%plaspara: [kappa,niu,lenda,M]
%stress:6x1 vector,[Sx;Sy;Sz;Sxy;Syz;Szx]
%pstrain(plastic strain):6x1 vector,[PEx;PEy;PEz;PExy;PEyz;PEzx]
%dstrain(strain increment):6x1 vector,[dEx;dEy;dEz;dExy;dEyz;dEzx]
%histPara:[voidratio,pc]
function [stress,pstrain,histPara,Dct,is_plastic]=Modified_Cam_Clay3D_new(plasPara,stress,pstrain,dstrain,histPara)
kappa=plasPara(1);niu=plasPara(2);
lenda=plasPara(3);M=plasPara(4);
voidratio=histPara(1);pc=histPara(2);MM=M^2;
%===========================Elastic prediction=============================
p0=mean(stress(1:3));
dstrainv=sum(dstrain(1:3));
[dsigma,~,D_elastic]=Porous_elastic(kappa,niu,voidratio,p0,dstrain,dstrainv);
[s,p,q]=Get_stress_info(stress+dsigma);
[F,refF]=Get_yield_value(p,q,MM,pc);
%=============================Elastic output===============================
if F<refF*1e-5
    Dct=D_elastic;
    stress=stress+dsigma;
    is_plastic=false;
    histPara(1)=max(exp(dstrainv).*(1+voidratio)-1,0);
    return;
end
%======================Plastic Constitutive integral=======================
if abs(dstrainv)>1e-10
    theta=(exp(dstrainv)-1)/dstrainv*(1+voidratio)/(lenda-kappa);
else
    theta=(1+voidratio)/(lenda-kappa);
end
iternum=0;tol=1e-10;eye6=eye(6);dLamda=0;pc0=pc;
while true
    %Calculate derivative
    [dGdSig,ddGddSig]=Get_potential_derivative(s,p,pc,MM);
    dpstrain=dLamda*dGdSig;dGdSigv=sum(dGdSig(1:3));

    % Calculate residuals
    [dsigma_corrected,K,D_elastic]=Porous_elastic(kappa,niu,voidratio,p0,dstrain-dpstrain,dstrainv);
    sigma_res=dsigma_corrected-dsigma;
    pc_res=pc0*exp(-theta*dLamda*dGdSigv)-pc;  

    % Calculate tangent matrix
    tempv=pc0*exp(-theta*dLamda*dGdSigv)*theta;
    k=[eye6+dLamda*D_elastic*ddGddSig, D_elastic*dGdSig, dLamda*[K;K;K;0;0;0];
                             dGdSig.',                0,                    p;
       [2/3,2/3,2/3,0,0,0]*tempv*dLamda,  tempv*dGdSigv,      1+tempv*dLamda];
    
    % Convergence criterion
    if abs(F)<refF*tol&&norm(sigma_res)<tol*norm(dsigma); break;end

    %Update variables
    x=k\[sigma_res;-F;pc_res];
    dsigma=dsigma+x(1:6);
    dLamda=dLamda+x(7);
    pc=pc+x(8);

    %Calculate yield value
    [s,p,q]=Get_stress_info(stress+dsigma);
    [F,refF]=Get_yield_value(p,q,MM,pc);

    %count
    iternum=iternum+1;
    if iternum>10000
        error('Material nonlinear iteration cannot converge in UMat')
    end
end
%========================Plastic output====================================
is_plastic=true;
stress=stress+dsigma;
pstrain=pstrain+dpstrain;
Dct=k\[D_elastic;zeros(2,6)];%Consistent tangent matrix
Dct=Dct(1:6,1:6);
histPara(2)=pc;
histPara(1)=max(exp(dstrainv).*(1+voidratio)-1,0);
end

function [dsigma,K,D_elastic]=Porous_elastic(kappa,niu,e,p,destrain,dstrainv)
destrainv=sum(destrain(1:3));
yita=destrainv/dstrainv;
if abs(destrainv)>1e-10&&abs(dstrainv)>1e-10
    dp_negative=p.*(exp(yita/kappa*(1+e).*(1-exp(dstrainv)))-1);
    K=dp_negative./destrainv;
elseif abs(destrainv)>1e-10&&abs(dstrainv)<1e-10
    dp_negative=p.*(exp(-(1+e)/kappa*destrainv)-1);
    K=dp_negative./destrainv;
elseif abs(destrainv)<1e-10&&abs(dstrainv)<1e-10
    K=-(1+e)/kappa*p;
    dp_negative=K*destrainv;
else
    K=-exp(dstrainv).*(1+e)/kappa*p;
    dp_negative=K*destrainv;
end
G=3*(1-2*niu)/2/(1+niu)*K;
G2=2*G;H=K-G2/3;
dsigma=[H*destrainv+G2*destrain(1:3);G*destrain(4:6)];
%==========================================================================
e=exp(dstrainv).*(1+e)-1;
K=-1/kappa*(1+e)*(p+dp_negative);
G=3*(1-2*niu)/2/(1+niu)*K;
G2=2*G;H=K-G2/3;F=G2+H;
D_elastic=[F, H, H, 0, 0, 0;
           H, F, H, 0, 0, 0;
           H, H, F, 0, 0, 0;
           0, 0, 0, G, 0, 0;
           0, 0, 0, 0, G, 0;
           0, 0, 0, 0, 0, G;];
end

function [F,refF]=Get_yield_value(p,q,MM,pc)
F=q^2/MM+p*(p+pc);
refF=pc^2;
end

function [dGdSig,ddGddSig]=Get_potential_derivative(s,p,pc,MM)
dGdSig=[(2*p + pc)/3 + 3*s(1)/MM;
        (2*p + pc)/3 + 3*s(2)/MM;
        (2*p + pc)/3 + 3*s(3)/MM;
                     (6*s(4))/MM;
                     (6*s(5))/MM;
                     (6*s(6))/MM];
ddGddSig=[2/MM + 2/9, 2/9 - 1/MM, 2/9 - 1/MM,     0,     0,     0;
          2/9 - 1/MM, 2/MM + 2/9, 2/9 - 1/MM,     0,     0,     0;
          2/9 - 1/MM, 2/9 - 1/MM, 2/MM + 2/9,     0,     0,     0;
                  0,           0,           0, 6/MM,     0,     0;
                  0,           0,           0,     0, 6/MM,     0;
                  0,           0,           0,     0,     0, 6/MM];
end


function [s,sm,sbar]=Get_stress_info(s)
sm=(s(1)+s(2)+s(3))/3;
s(1)=s(1)-sm;s(2)=s(2)-sm;s(3)=s(3)-sm;
J2=0.5*(s(1).^2+s(2).^2+s(3).^2)+s(4).^2+s(5).^2+s(6).^2;
sbar=sqrt(3*J2);
end