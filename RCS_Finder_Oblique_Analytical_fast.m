function RCS_Normalized=RCS_Finder_Oblique_Analytical_fast(OT,A,XminElements,XmaxElements,YminElements,YmaxElements,ZminElements,ZmaxElements,xmin,xmax,ymin,ymax,zmin,zmax,k0,THETA,PHI,THETA_IN,PHI_IN)

%% 3D bases
% N1=(1/ly/lz)*(yc+ly/2-y)*(zc+lz/2-z);
% N2=(1/ly/lz)*(-yc+ly/2+y)*(zc+lz/2-z);
% N3=(1/ly/lz)*(yc+ly/2-y)*(-zc+lz/2+z);
% N4=(1/ly/lz)*(-yc+ly/2+y)*(-zc+lz/2+z);
% 
% N5=(1/lx/lz)*(xc+lx/2-x)*(zc+lz/2-z);
% N6=(1/lx/lz)*(xc+lx/2-x)*(-zc+lz/2+z);
% N7=(1/lx/lz)*(-xc+lx/2+x)*(zc+lz/2-z);
% N8=(1/lx/lz)*(-xc+lx/2+x)*(-zc+lz/2+z);
% 
% N9=(1/lx/ly)*(xc+lx/2-x)*(yc+ly/2-y);
% N10=(1/lx/ly)*(-xc+lx/2+x)*(yc+ly/2-y);
% N11=(1/lx/ly)*(xc+lx/2-x)*(-yc+ly/2+y);
% N12=(1/lx/ly)*(-xc+lx/2+x)*(-yc+ly/2+y);

Lambda=2*pi/k0;
EPS0=8.8541878176e-12;
MU0=4*pi*1e-7;
ETA0=sqrt(MU0/EPS0);
C_O=2.997924580003452e+08;
omega=2*pi*C_O/Lambda;

%% z=zmax Plane , top face

NT_Zmax=0; % N theta on Z=zmax
LT_Zmax=0; % L theta on Z=zmax
NP_Zmax=0; % N phi on Z=zmax
LP_Zmax=0; % L phi on Z=zmax

for e=1:length(ZmaxElements)
    
    xmine=OT.BinBoundaries(ZmaxElements(e),1);
    ymine=OT.BinBoundaries(ZmaxElements(e),2);
    xmaxe=OT.BinBoundaries(ZmaxElements(e),4);
    ymaxe=OT.BinBoundaries(ZmaxElements(e),5);
    zmine=OT.BinBoundaries(ZmaxElements(e),3);
    zmaxe=OT.BinBoundaries(ZmaxElements(e),6);

    xc=0.5*(xmine+xmaxe);
    yc=0.5*(ymine+ymaxe);
    zc=0.5*(zmine+zmaxe);
    
    lx=OT.BinEdgeSize(ZmaxElements(e),1);
    ly=OT.BinEdgeSize(ZmaxElements(e),2);
    lz=OT.BinEdgeSize(ZmaxElements(e),3);

    % Eincident here has only z component
    Exi1=0;Exi2=0;Exi3=0;Exi4=0;
    Eyi5=0;Eyi6=0;Eyi7=0;Eyi8=0;
    
    Ezi9=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi10=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi11=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi12=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));

E1=A(OT.ConnectivityArrayEdge(ZmaxElements(e),1))-Exi1;
E2=A(OT.ConnectivityArrayEdge(ZmaxElements(e),2))-Exi2;
E3=A(OT.ConnectivityArrayEdge(ZmaxElements(e),3))-Exi3;
E4=A(OT.ConnectivityArrayEdge(ZmaxElements(e),4))-Exi4;
E5=A(OT.ConnectivityArrayEdge(ZmaxElements(e),5))-Eyi5;
E6=A(OT.ConnectivityArrayEdge(ZmaxElements(e),6))-Eyi6;
E7=A(OT.ConnectivityArrayEdge(ZmaxElements(e),7))-Eyi7;
E8=A(OT.ConnectivityArrayEdge(ZmaxElements(e),8))-Eyi8;
E9=A(OT.ConnectivityArrayEdge(ZmaxElements(e),9))-Ezi9;
E10=A(OT.ConnectivityArrayEdge(ZmaxElements(e),10))-Ezi10;
E11=A(OT.ConnectivityArrayEdge(ZmaxElements(e),11))-Ezi11;
E12=A(OT.ConnectivityArrayEdge(ZmaxElements(e),12))-Ezi12;
 
% 3D bases on Zmax plane

% N3=@(y)(1/ly/lz)*(yc+ly/2-y)*(-zc+lz/2+zmaxe);
% N4=@(y)(1/ly/lz)*(-yc+ly/2+y)*(-zc+lz/2+zmaxe);
% N6=@(x)(1/lx/lz)*(xc+lx/2-x)*(-zc+lz/2+zmaxe);
% N8=@(x)(1/lx/lz)*(-xc+lx/2+x)*(-zc+lz/2+zmaxe);
%  

% curlx=@(x) -E5*((-1/lx/lz)*(xc+lx/2-x))-E6*((1/lx/lz)*(xc+lx/2-x))-E7*((-1/lx/lz)*(-xc+lx/2+x))-E8*((1/lx/lz)*(-xc+lx/2+x))+E9*((-1/lx/ly)*(xc+lx/2-x))+E10*((-1/lx/ly)*(-xc+lx/2+x))+E11*((1/lx/ly)*(xc+lx/2-x))+E12*((1/lx/ly)*(-xc+lx/2+x));
% curly=@(y) (E1*((-1/ly/lz)*(yc+ly/2-y))+E2*((-1/ly/lz)*(-yc+ly/2+y))+E3*((1/ly/lz)*(yc+ly/2-y))+E4*((1/ly/lz)*(-yc+ly/2+y)))-(E9*((-1/lx/ly)*(yc+ly/2-y))+E10*((1/lx/ly)*(yc+ly/2-y))+E12*((1/lx/ly)*(-yc+ly/2+y))+E11*((-1/lx/ly)*(-yc+ly/2+y)));
% Jx=@(y)(1/1j/omega/MU0)*(curly(y));
% Jy=@(x)-(1/1j/omega/MU0)*(curlx(x));
% 
% Mx=@(x) (E6*N6(x)+E8*N8(x));
% My=@(y) -(E3*N3(y)+E4*N4(y));

%%%%% First Part of NT
EE=(1/1j/omega/MU0)*cos(THETA)*cos(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E3-E1)/ly/lz+(E9-E10)/lx/ly;
gg=(E4-E2)/ly/lz+(E11-E12)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1= ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));

FF2= gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));

FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2 ;

 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2;
 FF=FF1+FF2;

end
AA=EE*FF*GG;

%%%%% Second Part of NT
EE=(-1/1j/omega/MU0)*cos(THETA)*sin(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E5-E6)/lx/lz+(-E9+E11)/lx/ly;
gg=(E7-E8)/lx/lz+(-E10+E12)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1= ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));

FF2= gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));

FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2 +(ff/2)*(hh-xmine).^2;

 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2;   
end

BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (Jx(y)*cos(THETA)*cos(PHI)+Jy(x)*cos(THETA)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
NT_Zmax=NT_Zmax+TEMP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% First Part of NP
EE=(-1/1j/omega/MU0)*sin(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E3-E1)/ly/lz+(E9-E10)/lx/ly;
gg=(E4-E2)/ly/lz+(E11-E12)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1= ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2= gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));

FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2-(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2;   
end

AA=EE*FF*GG;

%%%%% Second Part of NP
EE=(-1/1j/omega/MU0)*cos(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E5-E6)/lx/lz+(-E9+E11)/lx/ly;
gg=(E7-E8)/lx/lz+(-E10+E12)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1= ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));

FF2= gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));

FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;

 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2;   
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (-Jx(y)*sin(PHI)+Jy(x)*cos(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
NP_Zmax=NP_Zmax+TEMP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% First Part of LT
EE=cos(THETA)*cos(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=E6*(-zc+lz/2+zmaxe)/lx/lz;
gg=E8*(-zc+lz/2+zmaxe)/lx/lz;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1= ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2= gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));

FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2 -(gg/2)*(kk+xmine).^2;
 FF=FF1+FF2;   
end
AA=EE*FF*GG;

%%%%% Second Part of LT
EE=cos(THETA)*sin(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=-E3*(-zc+lz/2+zmaxe)/ly/lz;
gg=-E4*(-zc+lz/2+zmaxe)/ly/lz;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1= ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));

FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2+ (ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2-(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2;   
end
BB=EE*FF*GG;

TEMP=AA+BB;
% TEMP=quad2d(@(x,y) (Mx(x)*cos(THETA)*cos(PHI)+My(y)*cos(THETA)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
LT_Zmax=LT_Zmax+TEMP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% First Part of LP
EE=-sin(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=E6*(-zc+lz/2+zmaxe)/lx/lz;
gg=E8*(-zc+lz/2+zmaxe)/lx/lz;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2 +(ff/2)*(hh-xmine).^2;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2  ;
 FF=FF1+FF2;   
end
AA=EE*FF*GG;

%%%%% Second Part of LP
EE=cos(PHI)*exp(1j*k0*zmaxe*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=-E3*(-zc+lz/2+zmaxe)/ly/lz;
gg=-E4*(-zc+lz/2+zmaxe)/ly/lz;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 + (ff/2)*(hh-ymine).^2;
 FF2=(gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2;
 FF=FF1+FF2;   
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (-Mx(x)*sin(PHI)+My(y)*cos(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
LP_Zmax=LP_Zmax+TEMP;


end


%% z=zmin Plane , bottom face

NT_Zmin=0; % N theta on Z=zmin
LT_Zmin=0; % L theta on Z=zmin
NP_Zmin=0; % N phi on Z=zmin
LP_Zmin=0; % L phi on Z=zmin

for e=1:length(ZminElements)
    
    xmine=OT.BinBoundaries(ZminElements(e),1);
    ymine=OT.BinBoundaries(ZminElements(e),2);
    xmaxe=OT.BinBoundaries(ZminElements(e),4);
    ymaxe=OT.BinBoundaries(ZminElements(e),5);
    zmine=OT.BinBoundaries(ZminElements(e),3);
    zmaxe=OT.BinBoundaries(ZminElements(e),6);

    xc=0.5*(xmine+xmaxe);
    yc=0.5*(ymine+ymaxe);
    zc=0.5*(zmine+zmaxe);
    
    lx=OT.BinEdgeSize(ZminElements(e),1);
    ly=OT.BinEdgeSize(ZminElements(e),2);
    lz=OT.BinEdgeSize(ZminElements(e),3);

    Exi1=0;Exi2=0;Exi3=0;Exi4=0;
    Eyi5=0;Eyi6=0;Eyi7=0;Eyi8=0;
    
    Ezi9=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi10=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi11=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi12=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));

E1=A(OT.ConnectivityArrayEdge(ZminElements(e),1))-Exi1;
E2=A(OT.ConnectivityArrayEdge(ZminElements(e),2))-Exi2;
E3=A(OT.ConnectivityArrayEdge(ZminElements(e),3))-Exi3;
E4=A(OT.ConnectivityArrayEdge(ZminElements(e),4))-Exi4;
E5=A(OT.ConnectivityArrayEdge(ZminElements(e),5))-Eyi5;
E6=A(OT.ConnectivityArrayEdge(ZminElements(e),6))-Eyi6;
E7=A(OT.ConnectivityArrayEdge(ZminElements(e),7))-Eyi7;
E8=A(OT.ConnectivityArrayEdge(ZminElements(e),8))-Eyi8;
E9=A(OT.ConnectivityArrayEdge(ZminElements(e),9))-Ezi9;
E10=A(OT.ConnectivityArrayEdge(ZminElements(e),10))-Ezi10;
E11=A(OT.ConnectivityArrayEdge(ZminElements(e),11))-Ezi11;
E12=A(OT.ConnectivityArrayEdge(ZminElements(e),12))-Ezi12;
 
% 3D bases on Zmin plane

% N1=@(y)(1/ly/lz)*(yc+ly/2-y)*(zc+lz/2-zmine);
% N2=@(y)(1/ly/lz)*(-yc+ly/2+y)*(zc+lz/2-zmine);
% N5=@(x)(1/lx/lz)*(xc+lx/2-x)*(zc+lz/2-zmine);
% N7=@(x)(1/lx/lz)*(-xc+lx/2+x)*(zc+lz/2-zmine);
% 
% curlx=@(x) -E5*((-1/lx/lz)*(xc+lx/2-x))-E6*((1/lx/lz)*(xc+lx/2-x))-E7*((-1/lx/lz)*(-xc+lx/2+x))-E8*((1/lx/lz)*(-xc+lx/2+x))+E9*((-1/lx/ly)*(xc+lx/2-x))+E10*((-1/lx/ly)*(-xc+lx/2+x))+E11*((1/lx/ly)*(xc+lx/2-x))+E12*((1/lx/ly)*(-xc+lx/2+x));
% curly=@(y) (E1*((-1/ly/lz)*(yc+ly/2-y))+E2*((-1/ly/lz)*(-yc+ly/2+y))+E3*((1/ly/lz)*(yc+ly/2-y))+E4*((1/ly/lz)*(-yc+ly/2+y)))-(E9*((-1/lx/ly)*(yc+ly/2-y))+E10*((1/lx/ly)*(yc+ly/2-y))+E12*((1/lx/ly)*(-yc+ly/2+y))+E11*((-1/lx/ly)*(-yc+ly/2+y)));
% Jx=@(y)-(1/1j/omega/MU0)*(curly(y));
% Jy=@(x)(1/1j/omega/MU0)*(curlx(x));
% 
% Mx=@(x) -(E5*N5(x)+E7*N7(x));
% My=@(y) (E1*N1(y)+E2*N2(y));

%%%%% First Part of NT
EE=(-1/1j/omega/MU0)*cos(THETA)*cos(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E3-E1)/ly/lz+(E9-E10)/lx/ly;
gg=(E4-E2)/ly/lz+(E11-E12)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2+(ff/2)*(hh-ymine).^2 ;

 FF2= (gg/2)*(kk+ymaxe).^2-(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2;   
end
AA=EE*FF*GG;

%%%%% Second Part of NT
EE=(1/1j/omega/MU0)*cos(THETA)*sin(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E5-E6)/lx/lz+(-E9+E11)/lx/ly;
gg=(E7-E8)/lx/lz+(-E10+E12)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2 +(ff/2)*(hh-xmine).^2;
 FF2= (gg/2)*(kk+xmaxe).^2 - (gg/2)*(kk+xmine).^2;
 FF=FF1+FF2; 
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (Jx(y)*cos(THETA)*cos(PHI)+Jy(x)*cos(THETA)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
NT_Zmin=NT_Zmin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% First Part of NP
EE=(1/1j/omega/MU0)*sin(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E3-E1)/ly/lz+(E9-E10)/lx/ly;
gg=(E4-E2)/ly/lz+(E11-E12)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1= ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2= gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2 ;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2;
 FF=FF1+FF2;   
end
AA=EE*FF*GG;

%%%%% Second Part of NP
EE=(1/1j/omega/MU0)*cos(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E5-E6)/lx/lz+(-E9+E11)/lx/ly;
gg=(E7-E8)/lx/lz+(-E10+E12)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2 +(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (-Jx(y)*sin(PHI)+Jy(x)*cos(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
NP_Zmin=NP_Zmin+TEMP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% First Part of LT
EE=cos(THETA)*cos(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=-E5*(zc+lz/2-zmine)/lx/lz;
gg=-E7*(zc+lz/2-zmine)/lx/lz;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2  ;
 FF=FF1+FF2; 
end
AA=EE*FF*GG;

%%%%% Second Part of LT
EE=cos(THETA)*sin(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=E1*(zc+lz/2-zmine)/ly/lz;
gg=E2*(zc+lz/2-zmine)/ly/lz;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2+(ff/2)*(hh-ymine).^2 ;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2;
 FF=FF1+FF2; 
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (Mx(x)*cos(THETA)*cos(PHI)+My(y)*cos(THETA)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
LT_Zmin=LT_Zmin+TEMP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% First Part of LP
EE=-sin(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=-E5*(zc+lz/2-zmine)/lx/lz;
gg=-E7*(zc+lz/2-zmine)/lx/lz;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1= ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2 +(ff/2)*(hh-xmine).^2;
 FF2= (gg/2)*(kk+xmaxe).^2 -(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end
AA=EE*FF*GG;

%%%%% Second Part of LP
EE=cos(PHI)*exp(1j*k0*zmine*cos(THETA));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=E1*(zc+lz/2-zmine)/ly/lz;
gg=E2*(zc+lz/2-zmine)/ly/lz;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2+(ff/2)*(hh-ymine).^2 ;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2;
 FF=FF1+FF2;  
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,y) (-Mx(x)*sin(PHI)+My(y)*cos(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
LP_Zmin=LP_Zmin+TEMP;


end



%% y=ymin plane front face

NT_Ymin=0; % N theta on Y=ymin
LT_Ymin=0; % L theta on Y=ymin
NP_Ymin=0; % N phi on Y=ymin
LP_Ymin=0; % L phi on Y=ymin



for e=1:length(YminElements)
    xmine=OT.BinBoundaries(YminElements(e),1);
    zmine=OT.BinBoundaries(YminElements(e),3);
    xmaxe=OT.BinBoundaries(YminElements(e),4);
    zmaxe=OT.BinBoundaries(YminElements(e),6);
    ymine=OT.BinBoundaries(YminElements(e),2);
    ymaxe=OT.BinBoundaries(YminElements(e),5);

    xc=0.5*(xmine+xmaxe);
    yc=0.5*(ymine+ymaxe);
    zc=0.5*(zmine+zmaxe);
    lx=OT.BinEdgeSize(YminElements(e),1);
    ly=OT.BinEdgeSize(YminElements(e),2);
    lz=OT.BinEdgeSize(YminElements(e),3);

    Exi1=0;Exi2=0;Exi3=0;Exi4=0;
    Eyi5=0;Eyi6=0;Eyi7=0;Eyi8=0;
    
    Ezi9=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi10=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi11=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi12=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));

E1=A(OT.ConnectivityArrayEdge(YminElements(e),1))-Exi1;
E2=A(OT.ConnectivityArrayEdge(YminElements(e),2))-Exi2;
E3=A(OT.ConnectivityArrayEdge(YminElements(e),3))-Exi3;
E4=A(OT.ConnectivityArrayEdge(YminElements(e),4))-Exi4;
E5=A(OT.ConnectivityArrayEdge(YminElements(e),5))-Eyi5;
E6=A(OT.ConnectivityArrayEdge(YminElements(e),6))-Eyi6;
E7=A(OT.ConnectivityArrayEdge(YminElements(e),7))-Eyi7;
E8=A(OT.ConnectivityArrayEdge(YminElements(e),8))-Eyi8;
E9=A(OT.ConnectivityArrayEdge(YminElements(e),9))-Ezi9;
E10=A(OT.ConnectivityArrayEdge(YminElements(e),10))-Ezi10;
E11=A(OT.ConnectivityArrayEdge(YminElements(e),11))-Ezi11;
E12=A(OT.ConnectivityArrayEdge(YminElements(e),12))-Ezi12;
 
% 3D bases on ymin plane

% N1=@(z)(1/ly/lz)*(yc+ly/2-ymine)*(zc+lz/2-z);
% N3=@(z)(1/ly/lz)*(yc+ly/2-ymine)*(-zc+lz/2+z);
% N9=@(x)(1/lx/ly)*(xc+lx/2-x)*(yc+ly/2-ymine);
% N10=@(x)(1/lx/ly)*(-xc+lx/2+x)*(yc+ly/2-ymine);
% 
% 
% curlx=@(x) E9*((-1/lx/ly)*(xc+lx/2-x))+E10*((-1/lx/ly)*(-xc+lx/2+x))+E11*((1/lx/ly)*(xc+lx/2-x))+E12*((1/lx/ly)*(-xc+lx/2+x))-E5*((-1/lx/lz)*(xc+lx/2-x))-E6*((1/lx/lz)*(xc+lx/2-x))-E7*((-1/lx/lz)*(-xc+lx/2+x))-E8*((1/lx/lz)*(-xc+lx/2+x));
% curlz=@(z) E5*((-1/lx/lz)*(zc+lz/2-z))+E6*((-1/lx/lz)*(-zc+lz/2+z))+E7*((1/lx/lz)*(zc+lz/2-z))+E8*((1/lx/lz)*(-zc+lz/2+z))-(E1*((-1/ly/lz)*(zc+lz/2-z))+E2*((1/ly/lz)*(zc+lz/2-z))+E3*((-1/ly/lz)*(-zc+lz/2+z))+E4*((1/ly/lz)*(-zc+lz/2+z)));
% Jx=@(z) (1/1j/omega/MU0)*curlz(z);
% Jz=@(x) -(1/1j/omega/MU0)*curlx(x);
% 
% Mx=@(x) (E9*N9(x)+E10*N10(x));
% Mz=@(z) -(E1*N1(z)+E3*N3(z));

%%%%%%%%%%%%%%%% First part of NT
EE=(1/1j/omega/MU0)*cos(THETA)*cos(PHI)*exp(1j*k0*ymine*sin(THETA)*sin(PHI));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2  ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end
AA=EE*FF*GG;
%%%%%%%%%%%%%%%% Second part of NT
EE=(1/1j/omega/MU0)*sin(THETA)*exp(1j*k0*ymine*sin(THETA)*sin(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=(E5-E6)/lx/lz+(-E9+E11)/lx/ly;
gg=(E7-E8)/lx/lz+(-E10+E12)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end
BB=EE*FF*GG;

TEMP=AA+BB;


% TEMP=quad2d(@(x,z) (Jx(z)*cos(THETA)*cos(PHI)-Jz(x)*sin(THETA)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymine*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
NT_Ymin=NT_Ymin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=(-1/1j/omega/MU0)*sin(PHI)*exp(1j*k0*ymine*sin(THETA)*sin(PHI));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end
TEMP=EE*FF*GG;


% TEMP=quad2d(@(x,z) (-Jx(z)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymine*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
NP_Ymin=NP_Ymin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%First part of LT

EE=cos(THETA)*cos(PHI)*exp(1j*k0*ymine*sin(THETA)*sin(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E9*(yc+ly/2-ymine)/lx/ly;
gg=E10*(yc+ly/2-ymine)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end
AA=EE*FF*GG;

%%%%%%%%%%%%%%%% Second part of LT

EE=sin(THETA)*exp(1j*k0*ymine*sin(THETA)*sin(PHI));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG=(1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=E1*(yc+ly/2-ymine)/ly/lz;
gg=E3*(yc+ly/2-ymine)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,z) (Mx(x)*cos(THETA)*cos(PHI)-Mz(z)*sin(THETA)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymine*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
LT_Ymin=LT_Ymin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=-sin(PHI)*exp(1j*k0*ymine*sin(THETA)*sin(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG= (1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E9*(yc+ly/2-ymine)/lx/ly;
gg=E10*(yc+ly/2-ymine)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2 -(gg/2)*(kk+xmine).^2;
 FF=FF1+FF2; 
end
TEMP=EE*FF*GG;

% TEMP=quad2d(@(x,z) (-Mx(x)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymine*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
LP_Ymin=LP_Ymin+TEMP;

end



%% y=ymax plane back face

NT_Ymax=0; % N theta on Y=ymax
LT_Ymax=0; % L theta on Y=ymax
NP_Ymax=0; % N phi on Y=ymax
LP_Ymax=0; % L phi on Y=ymax



for e=1:length(YmaxElements)
    
    xmine=OT.BinBoundaries(YmaxElements(e),1);
    zmine=OT.BinBoundaries(YmaxElements(e),3);
    xmaxe=OT.BinBoundaries(YmaxElements(e),4);
    zmaxe=OT.BinBoundaries(YmaxElements(e),6);
    ymine=OT.BinBoundaries(YmaxElements(e),2);
    ymaxe=OT.BinBoundaries(YmaxElements(e),5);

    xc=0.5*(xmine+xmaxe);
    yc=0.5*(ymine+ymaxe);
    zc=0.5*(zmine+zmaxe);
    
    lx=OT.BinEdgeSize(YmaxElements(e),1);
    ly=OT.BinEdgeSize(YmaxElements(e),2);
    lz=OT.BinEdgeSize(YmaxElements(e),3);

    Exi1=0;Exi2=0;Exi3=0;Exi4=0;
    Eyi5=0;Eyi6=0;Eyi7=0;Eyi8=0;
    
    Ezi9=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi10=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi11=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi12=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));

E1=A(OT.ConnectivityArrayEdge(YmaxElements(e),1))-Exi1;
E2=A(OT.ConnectivityArrayEdge(YmaxElements(e),2))-Exi2;
E3=A(OT.ConnectivityArrayEdge(YmaxElements(e),3))-Exi3;
E4=A(OT.ConnectivityArrayEdge(YmaxElements(e),4))-Exi4;
E5=A(OT.ConnectivityArrayEdge(YmaxElements(e),5))-Eyi5;
E6=A(OT.ConnectivityArrayEdge(YmaxElements(e),6))-Eyi6;
E7=A(OT.ConnectivityArrayEdge(YmaxElements(e),7))-Eyi7;
E8=A(OT.ConnectivityArrayEdge(YmaxElements(e),8))-Eyi8;
E9=A(OT.ConnectivityArrayEdge(YmaxElements(e),9))-Ezi9;
E10=A(OT.ConnectivityArrayEdge(YmaxElements(e),10))-Ezi10;
E11=A(OT.ConnectivityArrayEdge(YmaxElements(e),11))-Ezi11;
E12=A(OT.ConnectivityArrayEdge(YmaxElements(e),12))-Ezi12;
 
 
% 3D bases on ymin plane

% N2=@(z)(1/ly/lz)*(-yc+ly/2+ymaxe)*(zc+lz/2-z);
% N4=@(z)(1/ly/lz)*(-yc+ly/2+ymaxe)*(-zc+lz/2+z);
% N11=@(x)(1/lx/ly)*(xc+lx/2-x)*(-yc+ly/2+ymaxe);
% N12=@(x)(1/lx/ly)*(-xc+lx/2+x)*(-yc+ly/2+ymaxe);
% 
% curlx=@(x) E9*((-1/lx/ly)*(xc+lx/2-x))+E10*((-1/lx/ly)*(-xc+lx/2+x))+E11*((1/lx/ly)*(xc+lx/2-x))+E12*((1/lx/ly)*(-xc+lx/2+x))-E5*((-1/lx/lz)*(xc+lx/2-x))-E6*((1/lx/lz)*(xc+lx/2-x))-E7*((-1/lx/lz)*(-xc+lx/2+x))-E8*((1/lx/lz)*(-xc+lx/2+x));
% curlz=@(z) E5*((-1/lx/lz)*(zc+lz/2-z))+E6*((-1/lx/lz)*(-zc+lz/2+z))+E7*((1/lx/lz)*(zc+lz/2-z))+E8*((1/lx/lz)*(-zc+lz/2+z))-(E1*((-1/ly/lz)*(zc+lz/2-z))+E2*((1/ly/lz)*(zc+lz/2-z))+E3*((-1/ly/lz)*(-zc+lz/2+z))+E4*((1/ly/lz)*(-zc+lz/2+z)));
% Jx=@(z) -(1/1j/omega/MU0)*curlz(z);
% Jz=@(x) (1/1j/omega/MU0)*curlx(x);
% 
% Mx=@(x) -(E11*N11(x)+E12*N12(x));
% Mz=@(z) (E2*N2(z)+E4*N4(z));

%%%%%%%%%%%%%%%% First part of NT
EE=(-1/1j/omega/MU0)*cos(THETA)*cos(PHI)*exp(1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG=(1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1=-(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2  ;
 FF2=(gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end
AA=EE*FF*GG;
%%%%%%%%%%%%%%%% Second part of NT
EE=(-1/1j/omega/MU0)*sin(THETA)*exp(1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=(E5-E6)/lx/lz+(-E9+E11)/lx/ly;
gg=(E7-E8)/lx/lz+(-E10+E12)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end
BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,z) (Jx(z)*cos(THETA)*cos(PHI)-Jz(x)*sin(THETA)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymaxe*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
NT_Ymax=NT_Ymax+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=(1/1j/omega/MU0)*sin(PHI)*exp(1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG=(1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1=-(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2=(gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end
TEMP=EE*FF*GG;


% TEMP=quad2d(@(x,z) (-Jx(z)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymaxe*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
NP_Ymax=NP_Ymax+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%First part of LT

EE=-cos(THETA)*cos(PHI)*exp(1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E11*(-yc+ly/2+ymaxe)/lx/ly;
gg=E12*(-yc+ly/2+ymaxe)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2+(ff/2)*(hh-xmine).^2 ;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end
AA=EE*FF*GG;

%%%%%%%%%%%%%%%% Second part of LT

EE=-sin(THETA)*exp(1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*xmaxe)-(1/aa)*exp(aa*xmine);
else
    GG=lx;
end

ff=E2*(-yc+ly/2+ymaxe)/ly/lz;
gg=E4*(-yc+ly/2+ymaxe)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2 -(gg/2)*(kk+zmine).^2;
 FF=FF1+FF2; 
end

BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(x,z) (Mx(x)*cos(THETA)*cos(PHI)-Mz(z)*sin(THETA)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymaxe*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
LT_Ymax=LT_Ymax+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=sin(PHI)*exp(1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG= (1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E11*(-yc+ly/2+ymaxe)/lx/ly;
gg=E12*(-yc+ly/2+ymaxe)/lx/ly;
hh=xc+lx/2;
kk=-xc+lx/2;
bb=1j*k0*sin(THETA)*cos(PHI);

if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
FF1=ff*(((hh-xmaxe)/bb).*exp(bb*xmaxe)+(1/bb/bb).*exp(bb*xmaxe))-ff*(((hh-xmine)/bb).*exp(bb*xmine)+(1/bb/bb).*exp(bb*xmine));
FF2=gg*(((kk+xmaxe)/bb).*exp(bb*xmaxe)-(1/bb/bb).*exp(bb*xmaxe))-gg*(((kk+xmine)/bb).*exp(bb*xmine)-(1/bb/bb).*exp(bb*xmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-xmaxe).^2 +(ff/2)*(hh-xmine).^2;
 FF2= (gg/2)*(kk+xmaxe).^2-(gg/2)*(kk+xmine).^2 ;
 FF=FF1+FF2; 
end

TEMP=EE*FF*GG;

% TEMP=quad2d(@(x,z) (-Mx(x)*sin(PHI)).*exp(1j*k0*(x*sin(THETA)*cos(PHI)+ymaxe*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
LP_Ymax=LP_Ymax+TEMP;

end


 
%% x=xmin plane left face

NT_Xmin=0; % N theta on X=xmin
LT_Xmin=0; % L theta on X=xmin
NP_Xmin=0; % N phi on X=xmin
LP_Xmin=0; % L phi on X=xmin

for e=1:length(XminElements)
    
    ymine=OT.BinBoundaries(XminElements(e),2);
    zmine=OT.BinBoundaries(XminElements(e),3);
    ymaxe=OT.BinBoundaries(XminElements(e),5);
    zmaxe=OT.BinBoundaries(XminElements(e),6);
    xmine=OT.BinBoundaries(XminElements(e),1);
    xmaxe=OT.BinBoundaries(XminElements(e),4);

    xc=0.5*(xmine+xmaxe);
    yc=0.5*(ymine+ymaxe);
    zc=0.5*(zmine+zmaxe);
    
    lx=OT.BinEdgeSize(XminElements(e),1);
    ly=OT.BinEdgeSize(XminElements(e),2);
    lz=OT.BinEdgeSize(XminElements(e),3);

    Exi1=0;Exi2=0;Exi3=0;Exi4=0;
    Eyi5=0;Eyi6=0;Eyi7=0;Eyi8=0;
    
    Ezi9=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi10=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi11=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi12=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));

E1=A(OT.ConnectivityArrayEdge(XminElements(e),1))-Exi1;
E2=A(OT.ConnectivityArrayEdge(XminElements(e),2))-Exi2;
E3=A(OT.ConnectivityArrayEdge(XminElements(e),3))-Exi3;
E4=A(OT.ConnectivityArrayEdge(XminElements(e),4))-Exi4;
E5=A(OT.ConnectivityArrayEdge(XminElements(e),5))-Eyi5;
E6=A(OT.ConnectivityArrayEdge(XminElements(e),6))-Eyi6;
E7=A(OT.ConnectivityArrayEdge(XminElements(e),7))-Eyi7;
E8=A(OT.ConnectivityArrayEdge(XminElements(e),8))-Eyi8;
E9=A(OT.ConnectivityArrayEdge(XminElements(e),9))-Ezi9;
E10=A(OT.ConnectivityArrayEdge(XminElements(e),10))-Ezi10;
E11=A(OT.ConnectivityArrayEdge(XminElements(e),11))-Ezi11;
E12=A(OT.ConnectivityArrayEdge(XminElements(e),12))-Ezi12;
 
% 3D bases on ymin plane

% N5=@(z)(1/lx/lz)*(xc+lx/2-xmine)*(zc+lz/2-z);
% N6=@(z)(1/lx/lz)*(xc+lx/2-xmine)*(-zc+lz/2+z);
% N9=@(y)(1/lx/ly)*(xc+lx/2-xmine).*(yc+ly/2-y);
% N11=@(y)(1/lx/ly)*(xc+lx/2-xmine).*(-yc+ly/2+y);
% 
% 
% 
% curlz=@(z) E5*((-1/lx/lz)*(zc+lz/2-z))+E6*((-1/lx/lz)*(-zc+lz/2+z))+E7*((1/lx/lz)*(zc+lz/2-z))+E8*((1/lx/lz)*(-zc+lz/2+z))-(E1*((-1/ly/lz)*(zc+lz/2-z))+E2*((1/ly/lz)*(zc+lz/2-z))+E3*((-1/ly/lz)*(-zc+lz/2+z))+E4*((1/ly/lz)*(-zc+lz/2+z)));
% curly=@(y) (E1*((-1/ly/lz)*(yc+ly/2-y))+E2*((-1/ly/lz)*(-yc+ly/2+y))+E3*((1/ly/lz)*(yc+ly/2-y))+E4*((1/ly/lz)*(-yc+ly/2+y)))-(E9*((-1/lx/ly)*(yc+ly/2-y))+E10*((1/lx/ly)*(yc+ly/2-y))+E12*((1/lx/ly)*(-yc+ly/2+y))+E11*((-1/lx/ly)*(-yc+ly/2+y)));
% Jy=@(z) -(1/1j/omega/MU0)*curlz(z);
% Jz=@(y) (1/1j/omega/MU0)*(curly(y));

% My=@(y) -(E9*N9(y)+E11*N11(y));
% Mz=@(z) (E5*N5(z)+E6*N6(z));

%%%%%%% First Part of NT

EE=(-1/1j/omega/MU0)*cos(THETA)*sin(PHI)*exp(1j*k0*xmine*sin(THETA)*cos(PHI));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG=(1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2 -(gg/2)*(kk+zmine).^2;
 FF=FF1+FF2; 
end

AA=EE*FF*GG;


%%%%%%%%%%%%%%%%%%% Second part of NT

EE=(-1/1j/omega/MU0)*sin(THETA)*exp(1j*k0*xmine*sin(THETA)*cos(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=(E3-E1)/ly/lz+(E9-E10)/lx/ly;
gg=(E4-E2)/ly/lz+(E11-E12)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2-(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2; 
end

BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(y,z) (Jy(z)*cos(THETA)*sin(PHI)-Jz(y)*sin(THETA)).*exp(1j*k0*(xmine*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
NT_Xmin=NT_Xmin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=(-1/1j/omega/MU0)*cos(PHI)*exp(1j*k0*xmine*sin(THETA)*cos(PHI));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1=-(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2=(gg/2)*(kk+zmaxe).^2+(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end

TEMP=EE*FF*GG;

% TEMP=quad2d(@(y,z) (Jy(z)*cos(PHI)).*exp(1j*k0*(xmine*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
NP_Xmin=NP_Xmin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% First Part of LT
EE=-cos(THETA)*sin(PHI)*exp(1j*k0*xmine*sin(THETA)*cos(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG= (1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E9*(xc+lx/2-xmine)/lx/ly;
gg=E11*(xc+lx/2-xmine)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2; 
end

AA=EE*FF*GG;

%%%%%%%%%%%%%%%%%%%%% Second Part of LT
EE=-sin(THETA)*exp(1j*k0*xmine*sin(THETA)*cos(PHI));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=E5*(xc+lx/2-xmine)/lx/lz;
gg=E6*(xc+lx/2-xmine)/lx/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end

BB=EE*FF*GG;

TEMP=AA+BB;


% TEMP=quad2d(@(y,z) (My(y)*cos(THETA)*sin(PHI)-Mz(z)*sin(THETA)).*exp(1j*k0*(xmine*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
LT_Xmin=LT_Xmin+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=-cos(PHI)*exp(1j*k0*xmine*sin(THETA)*cos(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E9*(xc+lx/2-xmine)/lx/ly;
gg=E11*(xc+lx/2-xmine)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2; 
end

TEMP=EE*FF*GG;

% TEMP=quad2d(@(y,z) (My(y)*cos(PHI)).*exp(1j*k0*(xmine*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
LP_Xmin=LP_Xmin+TEMP;




end



%% x=xmax plane right face

NT_Xmax=0; % N theta on X=xmax
LT_Xmax=0; % L theta on X=xmax
NP_Xmax=0; % N phi on X=xmax
LP_Xmax=0; % L phi on X=xmax

for e=1:length(XmaxElements)
    
    ymine=OT.BinBoundaries(XmaxElements(e),2);
    zmine=OT.BinBoundaries(XmaxElements(e),3);
    ymaxe=OT.BinBoundaries(XmaxElements(e),5);
    zmaxe=OT.BinBoundaries(XmaxElements(e),6);
    xmine=OT.BinBoundaries(XmaxElements(e),1);
    xmaxe=OT.BinBoundaries(XmaxElements(e),4);

    xc=0.5*(xmine+xmaxe);
    yc=0.5*(ymine+ymaxe);
    zc=0.5*(zmine+zmaxe);
    
    lx=OT.BinEdgeSize(XmaxElements(e),1);
    ly=OT.BinEdgeSize(XmaxElements(e),2);
    lz=OT.BinEdgeSize(XmaxElements(e),3);

    Exi1=0;Exi2=0;Exi3=0;Exi4=0;
    Eyi5=0;Eyi6=0;Eyi7=0;Eyi8=0;
    
    Ezi9=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi10=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymine*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi11=exp(-1j*k0*(xmine*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));
    Ezi12=exp(-1j*k0*(xmaxe*sin(THETA_IN)*cos(PHI_IN)+ymaxe*sin(THETA_IN)*sin(PHI_IN)+zc*cos(THETA_IN)));

E1=A(OT.ConnectivityArrayEdge(XmaxElements(e),1))-Exi1;
E2=A(OT.ConnectivityArrayEdge(XmaxElements(e),2))-Exi2;
E3=A(OT.ConnectivityArrayEdge(XmaxElements(e),3))-Exi3;
E4=A(OT.ConnectivityArrayEdge(XmaxElements(e),4))-Exi4;
E5=A(OT.ConnectivityArrayEdge(XmaxElements(e),5))-Eyi5;
E6=A(OT.ConnectivityArrayEdge(XmaxElements(e),6))-Eyi6;
E7=A(OT.ConnectivityArrayEdge(XmaxElements(e),7))-Eyi7;
E8=A(OT.ConnectivityArrayEdge(XmaxElements(e),8))-Eyi8;
E9=A(OT.ConnectivityArrayEdge(XmaxElements(e),9))-Ezi9;
E10=A(OT.ConnectivityArrayEdge(XmaxElements(e),10))-Ezi10;
E11=A(OT.ConnectivityArrayEdge(XmaxElements(e),11))-Ezi11;
E12=A(OT.ConnectivityArrayEdge(XmaxElements(e),12))-Ezi12;
 
% 3D bases on ymin plane

% N7=@(z)(1/lx/lz)*(-xc+lx/2+xmaxe)*(zc+lz/2-z);
% N8=@(z)(1/lx/lz)*(-xc+lx/2+xmaxe)*(-zc+lz/2+z);
% N10=@(y)(1/lx/ly)*(-xc+lx/2+xmaxe).*(yc+ly/2-y);
% N12=@(y)(1/lx/ly)*(-xc+lx/2+xmaxe).*(-yc+ly/2+y);
% 
% 
% 
% curlz=@(z) E5*((-1/lx/lz)*(zc+lz/2-z))+E6*((-1/lx/lz)*(-zc+lz/2+z))+E7*((1/lx/lz)*(zc+lz/2-z))+E8*((1/lx/lz)*(-zc+lz/2+z))-(E1*((-1/ly/lz)*(zc+lz/2-z))+E2*((1/ly/lz)*(zc+lz/2-z))+E3*((-1/ly/lz)*(-zc+lz/2+z))+E4*((1/ly/lz)*(-zc+lz/2+z)));
% curly=@(y) (E1*((-1/ly/lz)*(yc+ly/2-y))+E2*((-1/ly/lz)*(-yc+ly/2+y))+E3*((1/ly/lz)*(yc+ly/2-y))+E4*((1/ly/lz)*(-yc+ly/2+y)))-(E9*((-1/lx/ly)*(yc+ly/2-y))+E10*((1/lx/ly)*(yc+ly/2-y))+E12*((1/lx/ly)*(-yc+ly/2+y))+E11*((-1/lx/ly)*(-yc+ly/2+y)));
% Jy=@(z) (1/1j/omega/MU0)*curlz(z);
% Jz=@(y) -(1/1j/omega/MU0)*(curly(y));
% 
% My=@(y) (E10*N10(y)+E12*N12(y));
% Mz=@(z) -(E7*N7(z)+E8*N8(z));

%%%%%%%%%%%%%%%%%%%%%%%%%%% First Part of NT
EE=(1/1j/omega/MU0)*cos(THETA)*sin(PHI)*exp(1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG=(1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end

AA=EE*FF*GG;

%%%%%%%%%%%%%%%%%%% Second part of NT

EE=(1/1j/omega/MU0)*sin(THETA)*exp(1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=(E3-E1)/ly/lz+(E9-E10)/lx/ly;
gg=(E4-E2)/ly/lz+(E11-E12)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2; 
end

BB=EE*FF*GG;

TEMP=AA+BB;

% TEMP=quad2d(@(y,z) (Jy(z)*cos(THETA)*sin(PHI)-Jz(y)*sin(THETA)).*exp(1j*k0*(xmaxe*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
NT_Xmax=NT_Xmax+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=(1/1j/omega/MU0)*cos(PHI)*exp(1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG= (1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=(E7-E5)/lx/lz+(E1-E2)/ly/lz;
gg=(E8-E6)/lx/lz+(E3-E4)/ly/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end

TEMP=EE*FF*GG;

% TEMP=quad2d(@(y,z) (Jy(z)*cos(PHI)).*exp(1j*k0*(xmaxe*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
NP_Xmax=NP_Xmax+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% First Part of LT
EE=cos(THETA)*sin(PHI)*exp(1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E10*(-xc+lx/2+xmaxe)/lx/ly;
gg=E12*(-xc+lx/2+xmaxe)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2; 
end

AA=EE*FF*GG;

%%%%%%%%%%%%%%%%%%%%% Second Part of LT
EE=sin(THETA)*exp(1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
GG=(1/aa)*exp(aa*ymaxe)-(1/aa)*exp(aa*ymine);
else
    GG=ly;
end

ff=E7*(-xc+lx/2+xmaxe)/lx/lz;
gg=E8*(-xc+lx/2+xmaxe)/lx/lz;
hh=zc+lz/2;
kk=-zc+lz/2;
bb=1j*k0*cos(THETA);

if cosd(THETA*180/pi)~=0
FF1= ff*(((hh-zmaxe)/bb).*exp(bb*zmaxe)+(1/bb/bb).*exp(bb*zmaxe))-ff*(((hh-zmine)/bb).*exp(bb*zmine)+(1/bb/bb).*exp(bb*zmine));
FF2=gg*(((kk+zmaxe)/bb).*exp(bb*zmaxe)-(1/bb/bb).*exp(bb*zmaxe))-gg*(((kk+zmine)/bb).*exp(bb*zmine)-(1/bb/bb).*exp(bb*zmine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-zmaxe).^2+(ff/2)*(hh-zmine).^2 ;
 FF2= (gg/2)*(kk+zmaxe).^2-(gg/2)*(kk+zmine).^2 ;
 FF=FF1+FF2; 
end

BB=EE*FF*GG;

TEMP=AA+BB;


% TEMP=quad2d(@(y,z) (My(y)*cos(THETA)*sin(PHI)-Mz(z)*sin(THETA)).*exp(1j*k0*(xmaxe*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
LT_Xmax=LT_Xmax+TEMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE=cos(PHI)*exp(1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
GG=(1/aa)*exp(aa*zmaxe)-(1/aa)*exp(aa*zmine);
else
    GG=lz;
end

ff=E10*(-xc+lx/2+xmaxe)/lx/ly;
gg=E12*(-xc+lx/2+xmaxe)/lx/ly;
hh=yc+ly/2;
kk=-yc+ly/2;
bb=1j*k0*sin(THETA)*sin(PHI);

if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
FF1=ff*(((hh-ymaxe)/bb).*exp(bb*ymaxe)+(1/bb/bb).*exp(bb*ymaxe))-ff*(((hh-ymine)/bb).*exp(bb*ymine)+(1/bb/bb).*exp(bb*ymine));
FF2=gg*(((kk+ymaxe)/bb).*exp(bb*ymaxe)-(1/bb/bb).*exp(bb*ymaxe))-gg*(((kk+ymine)/bb).*exp(bb*ymine)-(1/bb/bb).*exp(bb*ymine));
FF=FF1+FF2;
else
 FF1= -(ff/2)*(hh-ymaxe).^2 +(ff/2)*(hh-ymine).^2;
 FF2= (gg/2)*(kk+ymaxe).^2 -(gg/2)*(kk+ymine).^2 ;
 FF=FF1+FF2; 
end

TEMP=EE*FF*GG;

% TEMP=quad2d(@(y,z) (My(y)*cos(PHI)).*exp(1j*k0*(xmaxe*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
LP_Xmax=LP_Xmax+TEMP;

end

NT=NT_Zmin+NT_Zmax+NT_Xmin+NT_Xmax+NT_Ymin+NT_Ymax;
NP=NP_Zmin+NP_Zmax+NP_Xmin+NP_Xmax+NP_Ymin+NP_Ymax;
LT=LT_Zmin+LT_Zmax+LT_Xmin+LT_Xmax+LT_Ymin+LT_Ymax;
LP=LP_Zmin+LP_Zmax+LP_Xmin+LP_Xmax+LP_Ymin+LP_Ymax;

Er=0;
ET=-(LP+ETA0*NT);
EP=(LT-ETA0*NP);

ES=[Er ET EP]; % terme green function dar in rabete zarb nashode , meghdar vaghei ES bayad jk*exp(-jkr)/4pi*r dar in vector zarb shavad

RCS_Normalized=(k0^2/4/pi)*norm(ES).^2;


% RCS_Normalized=RCS/radius/radius; % the rcs should be divided by the area of face of cube


 
end
