%By Fang Huangcheng @BJTU 
%Email: valy_f@bjtu.edu.cn
%Ref:https://doi.org/10.16285/j.rsm.2010.07.041
%D_elastic:6x6 matrix
%mc_angle:2x1 vector, [friction angle; dilation angle],radian
%hardpara:nx3 vector,[Equivalent plastic strain,Cohesion,tangent];%%perfectly-plastic:[0,c,0]
%stress:6x1 vector,[Sx;Sy;Sz;Sxy;Syz;Szx]
%pstrain(plastic strain):6x1 vector,[PEx;PEy;PEz;PExy;PEyz;PEzx]
%dstrain(strain increment):6x1 vector,[dEx;dEy;dEz;dExy;dEyz;dEzx]
function [D_eq,stress,pstrain,is_plastic]=Mohr_Coulomb3D(D_elastic,mc_angle,hardpara,stress,pstrain,dstrain)
fai=mc_angle(1);pfai=mc_angle(2);theta_T=29.9/180*pi;
%--------------------------------------------------------------------------
Ks=[[3+tan(theta_T)*tan(3*theta_T),1/sqrt(3)*(tan(3*theta_T)-3*tan(theta_T))*sin(fai)]*cos(theta_T)/3;...
    [1/sqrt(3)*sin(fai)*cos(theta_T),sin(theta_T)]/cos(3*theta_T)/3];
%--------------------------------------------------------------------------
[c,H]=Hardening(hardpara,pstrain);H=max(H,0);
%--------------------------------------------------------------------------
dsigma_prediction=D_elastic*dstrain;
sigma_prediction=stress+dsigma_prediction;
[s,sigmam,sigmabar,Lode,K,dKdtheta]=Get_stress_info(sigma_prediction,fai,theta_T,Ks);
F=sigmam*sin(fai)+sigmabar*K-c*cos(fai);
refF=abs(sigmam*sin(fai))+c*cos(fai)+1e-32;
if F/refF<1e-8
    D_eq=D_elastic;
    stress=sigma_prediction;
    is_plastic=0;
    return
end
%--------------------------------------------------------------------------
NF2=F/refF;
N=min(floor(NF2/0.02)+1,5);
Aincrement2=1/N;
%--------------------------------------------------------------------------
flag=0;
for di=1:1:N
    sigma_prediction=stress+dsigma_prediction*Aincrement2;
    [s,sigmam,sigmabar,Lode,K,dKdtheta]=Get_stress_info(sigma_prediction,fai,theta_T,Ks);
    F=sigmam*sin(fai)+sigmabar*K-c*cos(fai);
    while abs(F)/refF>1e-10
        flag=flag+1;
        dF=Get_yield_derivative(fai,s,sigmabar,Lode,K,dKdtheta);
        dG=Get_potential_derivative(pfai,s,sigmabar,Lode,K,dKdtheta);
        dpstrain=F/(H+dF'*D_elastic*dG)*dG;
        dsigma_p=D_elastic*dpstrain;
        sigma_prediction=sigma_prediction-dsigma_p;
        pstrain=pstrain+dpstrain;
        
        [c,H]=Hardening(hardpara,pstrain);H=max(H,0);

        [s,sigmam,sigmabar,Lode,K,dKdtheta]=Get_stress_info(sigma_prediction,fai,theta_T,Ks);
        F=sigmam*sin(fai)+sigmabar*K-c*cos(fai);
        if flag>10000
            error('Material nonlinear iteration cannot converge in Mohr_Coulomb')
        end
    end
    stress=sigma_prediction;
end
dF=Get_yield_derivative(fai,s,sigmabar,Lode,K,dKdtheta);
dG=Get_potential_derivative(pfai,s,sigmabar,Lode,K,dKdtheta);

D_eq=D_elastic-(D_elastic*dG)*(dF'*D_elastic)/(H+dF'*D_elastic*dG);
is_plastic=1;
end

function [s,sm,sbar,Lode,K,dKdtheta]=Get_stress_info(s,fai,theta_T,Ks)
sm=(s(1)+s(2)+s(3))/3;
s(1)=s(1)-sm;s(2)=s(2)-sm;s(3)=s(3)-sm;
J2=0.5*(s(1).^2+s(2).^2+s(3).^2)+s(4).^2+s(5).^2+s(6).^2;
J3=s(1)*s(2)*s(3)+2*s(4)*s(5)*s(6)-s(1)*s(5).^2-s(2)*s(6).^2-s(3)*s(4).^2;
sbar=sqrt(J2);
if sbar~=0
    Lode=real(asin(complex(-1.5*sqrt(3)*J3/sbar^3))/3);
else
    Lode=0;sbar=1e-32;
end
if abs(Lode)>theta_T
    Ka=Ks(:,1)+Ks(:,2)*sign(Lode);
    K=Ka(1)-Ka(2)*sin(3.0*Lode);
    dKdtheta=-3.0*Ka(2)*cos(3.0*Lode);
else
    K=cos(Lode)-1.0/sqrt(3.0)*sin(fai)*sin(Lode);
    dKdtheta=-sin(Lode)-1.0/sqrt(3.0)*sin(fai)*cos(Lode);
end
end

function dF=Get_yield_derivative(fai,s,sigmabar,Lode,K,dKdtheta)
J2=sigmabar^2;
dsigmam=[1/3;1/3;1/3;0.0;0.0;0.0];
dsigmabar=[s(1)/2/sigmabar;s(2)/2/sigmabar;s(3)/2/sigmabar;s(4)/sigmabar;s(5)/sigmabar;s(6)/sigmabar;];
dJ3=[s(2)*s(3)-s(5)*s(5)+J2/3;
     s(1)*s(3)-s(6)*s(6)+J2/3;
     s(1)*s(2)-s(4)*s(4)+J2/3;
     (s(5)*s(6)-s(3)*s(4))*2;
     (s(6)*s(4)-s(1)*s(5))*2;
     (s(4)*s(5)-s(2)*s(6))*2;];
 
C1=sin(fai);
C2=(K-tan(3.0*Lode)*dKdtheta);
C3=-sqrt(3.0)/2.0/cos(3.0*Lode)/J2*dKdtheta;
dF=C1*dsigmam+C2*dsigmabar+C3*dJ3;
end

function dG=Get_potential_derivative(pfai,s,sigmabar,Lode,K,dKdtheta)
J2=sigmabar^2;
dsigmam=[1/3;1/3;1/3;0.0;0.0;0.0];
dsigmabar=[s(1)/2/sigmabar;s(2)/2/sigmabar;s(3)/2/sigmabar;s(4)/sigmabar;s(5)/sigmabar;s(6)/sigmabar;];
dJ3=[s(2)*s(3)-s(5)*s(5)+J2/3;
     s(1)*s(3)-s(6)*s(6)+J2/3;
     s(1)*s(2)-s(4)*s(4)+J2/3;
     (s(5)*s(6)-s(3)*s(4))*2;
     (s(6)*s(4)-s(1)*s(5))*2;
     (s(4)*s(5)-s(2)*s(6))*2;];
 
C1=sin(pfai);
C2=(K-tan(3.0*Lode)*dKdtheta);
C3=-sqrt(3.0)/2.0/cos(3.0*Lode)/J2*dKdtheta;
dG=C1*dsigmam+C2*dsigmabar+C3*dJ3;
end

function [sigmarlim,H]=Hardening(hardpara,pstrain)
pe=((pstrain(1)-pstrain(2)).^2+(pstrain(2)-pstrain(3)).^2+(pstrain(3)-pstrain(1)).^2)/2+(pstrain(4).^2+pstrain(5).^2+pstrain(6).^2)*0.75;
pe=sqrt(pe);
[~,loc]=min(abs(hardpara(:,1)-pe));
if hardpara(loc,1)-pe>0
    loc=loc-1;
end
H=hardpara(loc,3);sigmarlim=hardpara(loc,2)+H*(pe-hardpara(loc,1));
end