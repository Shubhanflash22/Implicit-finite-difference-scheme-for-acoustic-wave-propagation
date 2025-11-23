function numerical_anisotropy3
clear all; clc
nt=301; np=101;nphi=101;
theta= linspace(0,2*pi,nt);
p= linspace(0,pi/2,np);
phi = linspace(0, 2*pi, nphi);
%zeta= 1/sqrt(2);
zeta_coeff= [.9,.2,.5,.8,1,1.2,.9];
pi1=2;     pi2=5;

for i=1%1:5
    figure()
    zeta= zeta_coeff(i)/ sqrt(3); 
    c1=zeros(nt,np);
    c2=zeros(nt,np);
    c3=zeros(nt,np);
    c4=zeros(nt,np);
    c5=zeros(nt,np);
    %c1= Spectral_analysis_CD(2,theta,p,zeta);        plotting_polar(theta,1-c1,[pi1,pi2,1]);   plotting_surf(theta,p,c1,[pi1,pi2,1+pi2]);
    %c2= Spectral_analysis_CD(4,theta,p,zeta);        plotting_polar(theta,1-c2,[pi1,pi2,2]);   plotting_surf(theta,p,c2,[pi1,pi2,2+pi2]);
    %c3= Spectral_analysis_CD(6,theta,p,zeta);        plotting_polar(theta,1-c3,[pi1,pi2,3]);   plotting_surf(theta,p,c3,[pi1,pi2,3+pi2]);
    %c4= Spectral_analysis_CD(8,theta,p,zeta);        plotting_polar(theta,1-c4,[pi1,pi2,4]);   plotting_surf(theta,p,c4,[pi1,pi2,4+pi2]);
    %c5= Spectral_analysis_compact(6,theta,p,zeta);   plotting_polar(theta,1-c5,[pi1,pi2,5]);   plotting_surf(theta,p,c5,[pi1,pi2,5+pi2]);
end
end

function [vel]= Spectral_analysis_CD(select,theta,p,zeta)

n1=length(theta);   n2=length(p);
vel=zeros(n1,n2);
for j=1:np
    for l=1:nt
        pct= p(j)*cos(theta(l))*cos(phi(k));
        pst= p(j)*cos(theta(l))*sin(phi(k));
        psz= p(j)*sin(theta(l));
        select = pi1 + pi2;
        if (select==2)
            D2kx= -(-2+2*cos(pct)) ;
            D2kz= -(-2+2*cos(pst)) ;
            D2ky= -(-2+2*cos(psz)) ;
        elseif (select==4)
            D2kx= -(-30 + 16*2* cos(pct) - 2* cos(2*pct))/12;
            D2kz= -(-30 + 16*2* cos(pst) - 2* cos(2*pst))/12;
            D2ky= -(-30 + 16*2* cos(psz) - 2* cos(2*psz))/12;
        elseif (select==6)
            D2kx= -(-490 + 270*2*cos(pct) - 27*2* cos(2*pct) + 2*2* cos(3*pct))/180;
            D2kz= -(-490 + 270*2*cos(pst) - 27*2* cos(2*pst) + 2*2* cos(3*pst))/180;
            D2ky= -(-490 + 270*2*cos(psz) - 27*2* cos(2*psz) + 2*2* cos(3*psz))/180;
        elseif (select==8)
            D2kx= -( -14350 +8064*2*cos(pct) - 1008*2* cos(pct*2) + 128*2* cos(3*pct)  -  9*2*cos(4*pct))/5040;
            D2kz= -( -14350 +8064*2*cos(pst) - 1008*2* cos(pst*2) + 128*2* cos(3*pst)  -  9*2*cos(4*pst))/5040;
            D2ky= -( -14350 +8064*2*cos(psz) - 1008*2* cos(psz*2) + 128*2* cos(3*psz)  -  9*2*cos(4*psz))/5040;
        end
    temp=0.5*zeta*sqrt(D2kx + D2kz);
    c = 2*asin(temp)./(zeta*p(j));
    vel(:,j)=c;
    %polarplot(theta,c); hold on
end

end

function [vel]= Spectral_analysis_compact(select,theta,p,zeta)
n1=length(theta);   n2=length(p);
vel=zeros(n1,n2);
for j=1:length(p)
    pct= p(j)*cos(theta);
    pst= p(j)*sin(theta);
    if (select==6)
        a=12/11; b=3/11; c=0;
        alpha=2/11;  beta=0;
        num= a*2*(1-cos(pct)) + (b/4)*2*(1-cos(2*pct)) + (c/9)*2*(1-cos(3*pct));
        den= 1 + alpha*2*cos(pct) + beta*2*cos(2*pct);
        D2kx= num./den;
        num= a*2*(1-cos(pst)) + (b/4)*2*(1-cos(2*pst)) + (c/9)*2*(1-cos(3*pst));
        den= 1 + alpha*2*cos(pst) + beta*2*cos(2*pst);
        D2kz= num./den;
    end
    temp=zeta*sqrt(D2kx + D2kz);
    c = asin(.5*temp)*2 /(zeta*p(j) );
    vel(:,j)=c;
end
end

function plotting_polar(theta,c,fid)
subplot(fid(1),fid(2),fid(3))
for i=1:10:size(c,2)
    polarplot(theta,c(:,i));        hold on
end
%set(gca,'ThetaTickLabel',[])
%set(gca,'RTickLabel',[])
shg
end

function plotting_surf(theta,p,c,fid)
[T,R]=meshgrid(theta,p);
[X,Y]=pol2cart(T,R);
Z=c;
subplot(fid(1),fid(2),fid(3))
mesh(X,Y,Z')
xlabel('k_x \rightarrow')
ylabel('\leftarrow k_z ')
%zlabel('F(\partial_x^2)')
%set(gca,'ZTickLabel',[])
%zlim([.47 .53])
%colorbar 
%colormap(gray)
shg
end

