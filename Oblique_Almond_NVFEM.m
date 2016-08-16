function [RCS_Normalized2,ITER,SOLVETIME]=Oblique_Almond_NVFEM(THETA,PHI,OT,xmin,xmax,ymin,ymax,zmin,zmax,MASS_final,k0)

Nn=length(MASS_final);

%% ========================== impose 1st order ABC (Absorbing Boundary Condition) ========================

ZminElements=find(OT.BinBoundaries(:,3)<=zmin+.000001 & OT.BinBoundaries(:,3)>=zmin-.000001  );
ZmaxElements=find(OT.BinBoundaries(:,6)<=zmax+.000001 & OT.BinBoundaries(:,6)>=zmax-.000001 );
YminElements=find(OT.BinBoundaries(:,2)<=ymin+.000001 & OT.BinBoundaries(:,2)>=ymin-.000001 );
YmaxElements=find(OT.BinBoundaries(:,5)<=ymax+.000001 & OT.BinBoundaries(:,5)>=ymax-.000001 );
XminElements=find(OT.BinBoundaries(:,1)<=xmin+.000001 & OT.BinBoundaries(:,1)>=xmin-.000001 );
XmaxElements=find(OT.BinBoundaries(:,4)<=xmax+.000001 & OT.BinBoundaries(:,4)>=xmax-.000001  );


ABC_M_loc=(1/6)*[2 1 0 0 ; 1 2 0 0 ; 0 0 2 1 ; 0 0 1 2]; 

Coeff_Zmin = 1j*k0*(OT.BinEdgeSize(ZminElements,1).*OT.BinEdgeSize(ZminElements,2)).';
ABC_M_loc_Zmin = repmat(ABC_M_loc(:),1,length(ZminElements)).*repmat(Coeff_Zmin,16,1);
row_Zmin = OT.ConnectivityArrayEdge(ZminElements,[1 2 5 7 1 2 5 7 1 2 5 7 1 2 5 7]).';
col_Zmin = OT.ConnectivityArrayEdge(ZminElements,[1 1 1 1 2 2 2 2 5 5 5 5 7 7 7 7]).';

Coeff_Zmax = 1j*k0*(OT.BinEdgeSize(ZmaxElements,1).*OT.BinEdgeSize(ZmaxElements,2)).';
ABC_M_loc_Zmax = repmat(ABC_M_loc(:),1,length(ZmaxElements)).*repmat(Coeff_Zmax,16,1);
row_Zmax = OT.ConnectivityArrayEdge(ZmaxElements,[3 4 6 8 3 4 6 8 3 4 6 8 3 4 6 8]).';
col_Zmax = OT.ConnectivityArrayEdge(ZmaxElements,[3 3 3 3 4 4 4 4 6 6 6 6 8 8 8 8]).';

Coeff_Ymin = 1j*k0*(OT.BinEdgeSize(YminElements,1).*OT.BinEdgeSize(YminElements,3)).';
ABC_M_loc_Ymin = repmat(ABC_M_loc(:),1,length(YminElements)).*repmat(Coeff_Ymin,16,1);
row_Ymin = OT.ConnectivityArrayEdge(YminElements,[1 3 9 10 1 3 9 10 1 3 9 10 1 3 9 10]).';
col_Ymin = OT.ConnectivityArrayEdge(YminElements,[1 1 1 1 3 3 3 3 9 9 9 9 10 10 10 10]).';

Coeff_Ymax = 1j*k0*(OT.BinEdgeSize(YmaxElements,1).*OT.BinEdgeSize(YmaxElements,3)).';
ABC_M_loc_Ymax = repmat(ABC_M_loc(:),1,length(YmaxElements)).*repmat(Coeff_Ymax,16,1);
row_Ymax = OT.ConnectivityArrayEdge(YmaxElements,[2 4 11 12 2 4 11 12 2 4 11 12 2 4 11 12]).';
col_Ymax = OT.ConnectivityArrayEdge(YmaxElements,[2 2 2 2 4 4 4 4 11 11 11 11 12 12 12 12]).';


Coeff_Xmin = 1j*k0*(OT.BinEdgeSize(XminElements,2).*OT.BinEdgeSize(XminElements,3)).';
ABC_M_loc_Xmin = repmat(ABC_M_loc(:),1,length(XminElements)).*repmat(Coeff_Xmin,16,1);
row_Xmin = OT.ConnectivityArrayEdge(XminElements,[5 6 9 11 5 6 9 11 5 6 9 11 5 6 9 11]).';
col_Xmin = OT.ConnectivityArrayEdge(XminElements,[5 5 5 5 6 6 6 6 9 9 9 9 11 11 11 11]).';


Coeff_Xmax = 1j*k0*(OT.BinEdgeSize(XmaxElements,2).*OT.BinEdgeSize(XmaxElements,3)).';
ABC_M_loc_Xmax = repmat(ABC_M_loc(:),1,length(XmaxElements)).*repmat(Coeff_Xmax,16,1);
row_Xmax = OT.ConnectivityArrayEdge(XmaxElements,[7 8 10 12 7 8 10 12 7 8 10 12 7 8 10 12]).';
col_Xmax = OT.ConnectivityArrayEdge(XmaxElements,[7 7 7 7 8 8 8 8 10 10 10 10 12 12 12 12]).';

row = [row_Xmin row_Xmax row_Ymin row_Ymax row_Zmin row_Zmax];
col = [col_Xmin col_Xmax col_Ymin col_Ymax col_Zmin col_Zmax];
ABC_M_loc = [ABC_M_loc_Xmin ABC_M_loc_Xmax ABC_M_loc_Ymin ABC_M_loc_Ymax ABC_M_loc_Zmin ABC_M_loc_Zmax];

ABC_M = sparse(row,col,ABC_M_loc);

if size(ABC_M,1)~=size(MASS_final,1)
    ABC_M(size(MASS_final,1),size(MASS_final,1))=0;
end

MASS_final=MASS_final+ABC_M;
ABC_M=[];




%% impose RHS for ABC surface

% THETA=pi/2;
% incident field polarization is in az direction and k is in x-y plane

 b=zeros(Nn,1); % RHS vector

% %  zmin plane
 for e=1:length(ZminElements)
    
    lz=OT.BinEdgeSize(ZminElements(e),3);
    lx=OT.BinEdgeSize(ZminElements(e),1);
    ly=OT.BinEdgeSize(ZminElements(e),2);

    ymaxe=OT.BinBoundaries(ZminElements(e),5);
    ymine=OT.BinBoundaries(ZminElements(e),2);
    xmaxe=OT.BinBoundaries(ZminElements(e),4);
    xmine=OT.BinBoundaries(ZminElements(e),1);
    zmaxe=OT.BinBoundaries(ZminElements(e),6);
    zmine=OT.BinBoundaries(ZminElements(e),3);
    
    zc=0.5*(zmaxe+zmine);
    xc=0.5*(xmaxe+xmine);
    yc=0.5*(ymaxe+ymine);

%     temp1=-quad2d(@(x,y) (1/ly/lz).*(yc+ly/2-y).*(zc+lz/2-zmine).*(sin(THETA)*cos(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
%     temp2=-quad2d(@(x,y) (1/ly/lz).*(-yc+ly/2+y).*(zc+lz/2-zmine).*(sin(THETA)*cos(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
%     temp5=-quad2d(@(x,y) (1/lx/lz).*(xc+lx/2-x).*(zc+lz/2-zmine).*(sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
%     temp7=-quad2d(@(x,y) (1/lx/lz).*(-xc+lx/2+x).*(zc+lz/2-zmine).*(sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmine*cos(THETA))),xmine,xmaxe,ymine,ymaxe);

AA=(-1/ly/lz)*(zc+lz/2-zmine)*(sin(THETA)*cos(PHI))*exp(-1j*k0*zmine*cos(THETA));

aa=-1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
CC_temp=@(x) (1/aa).*exp(aa*x);
CC=CC_temp(xmaxe)-CC_temp(xmine);
else
    CC=lx;
end

bb=yc+ly/2; cc=-1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
BB_temp=@(y) ((bb-y)./cc).*exp(cc*y)+(1/cc/cc).*exp(cc*y);
BB=BB_temp(ymaxe)-BB_temp(ymine);
else
  BB_temp=@(y) (-1/2)*(bb-y).^2;
  BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp1=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb=-yc+ly/2;
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
    BB_temp=@(y) ((bb+y)./cc).*exp(cc*y)-(1/cc/cc).*exp(cc*y);
    BB=BB_temp(ymaxe)-BB_temp(ymine);
else
    BB_temp=@(y) (1/2)*(bb+y).^2;
    BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp2=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA=(-1/lx/lz)*(zc+lz/2-zmine)*(sin(THETA)*sin(PHI))*exp(-1j*k0*zmine*cos(THETA));

aa=-1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
CC_temp=@(y) (1/aa).*exp(aa*y);
CC=CC_temp(ymaxe)-CC_temp(ymine);
else
    CC=ly;
end

bb=xc+lx/2; cc=-1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb-x)./cc).*exp(cc*x)+(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (-1/2)*(bb-x).^2; 
BB=BB_temp(xmaxe)-BB_temp(xmine);
end   
temp5=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bb=-xc+lx/2;
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb+x)./cc).*exp(cc*x)-(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (1/2)*(bb+x).^2;
BB=BB_temp(xmaxe)-BB_temp(xmine);
end

temp7=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    b_loc_zmin=1j*k0*[temp1; temp2;temp5;temp7 ];        
    b(OT.ConnectivityArrayEdge(ZminElements(e),[1 2 5 7]),1)=b(OT.ConnectivityArrayEdge(ZminElements(e),[1 2 5 7 ]),1)+b_loc_zmin(1:4);
 end

 % %  zmax plane
 for e=1:length(ZmaxElements)
    
    lz=OT.BinEdgeSize(ZmaxElements(e),3);
    lx=OT.BinEdgeSize(ZmaxElements(e),1);
    ly=OT.BinEdgeSize(ZmaxElements(e),2);

    ymaxe=OT.BinBoundaries(ZmaxElements(e),5);
    ymine=OT.BinBoundaries(ZmaxElements(e),2);
    xmaxe=OT.BinBoundaries(ZmaxElements(e),4);
    xmine=OT.BinBoundaries(ZmaxElements(e),1);
    zmaxe=OT.BinBoundaries(ZmaxElements(e),6);
    zmine=OT.BinBoundaries(ZmaxElements(e),3);
    
    zc=0.5*(zmaxe+zmine);
    xc=0.5*(xmaxe+xmine);
    yc=0.5*(ymaxe+ymine);

%     temp3=quad2d(@(x,y) (1/ly/lz).*(yc+ly/2-y).*(-zc+lz/2+zmaxe).*(sin(THETA)*cos(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
%     temp4=quad2d(@(x,y) (1/ly/lz).*(-yc+ly/2+y).*(-zc+lz/2+zmaxe).*(sin(THETA)*cos(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
%     temp6=quad2d(@(x,y) (1/lx/lz).*(xc+lx/2-x).*(-zc+lz/2+zmaxe).*(sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);
%     temp8=quad2d(@(x,y) (1/lx/lz).*(-xc+lx/2+x).*(-zc+lz/2+zmaxe).*(sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+zmaxe*cos(THETA))),xmine,xmaxe,ymine,ymaxe);

AA=(1/ly/lz)*(-zc+lz/2+zmaxe)*(sin(THETA)*cos(PHI))*exp(-1j*k0*zmaxe*cos(THETA));

aa=-1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
CC_temp=@(x) (1/aa).*exp(aa*x);
CC=CC_temp(xmaxe)-CC_temp(xmine);
else
    CC=lx;
end

bb=yc+ly/2; cc=-1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
BB_temp=@(y) ((bb-y)./cc).*exp(cc*y)+(1/cc/cc).*exp(cc*y);
BB=BB_temp(ymaxe)-BB_temp(ymine);
else
  BB_temp=@(y) (-1/2)*(bb-y).^2;
  BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp3=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bb=-yc+ly/2;
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
    BB_temp=@(y) ((bb+y)./cc).*exp(cc*y)-(1/cc/cc).*exp(cc*y);
    BB=BB_temp(ymaxe)-BB_temp(ymine);
else
    BB_temp=@(y) (1/2)*(bb+y).^2;
    BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp4=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA=(1/lx/lz)*(-zc+lz/2+zmaxe)*(sin(THETA)*sin(PHI))*exp(-1j*k0*zmaxe*cos(THETA));

aa=-1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
CC_temp=@(y) (1/aa).*exp(aa*y);
CC=CC_temp(ymaxe)-CC_temp(ymine);
else
    CC=ly;
end

bb=xc+lx/2; cc=-1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb-x)./cc).*exp(cc*x)+(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (-1/2)*(bb-x).^2; 
BB=BB_temp(xmaxe)-BB_temp(xmine);
end   
temp6=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bb=-xc+lx/2;
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb+x)./cc).*exp(cc*x)-(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (1/2)*(bb+x).^2;
BB=BB_temp(xmaxe)-BB_temp(xmine);
end

temp8=AA*BB*CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b_loc_zmax=1j*k0*[temp3;temp4;temp6; temp8 ];        
    b(OT.ConnectivityArrayEdge(ZmaxElements(e),[3 4 6 8]),1)=b(OT.ConnectivityArrayEdge(ZmaxElements(e),[3 4 6 8 ]),1)+b_loc_zmax(1:4);
 end
 

 
 
% % front face

for e=1:length(YminElements)
    
    ly=OT.BinEdgeSize(YminElements(e),2);
    lz=OT.BinEdgeSize(YminElements(e),3);
    lx=OT.BinEdgeSize(YminElements(e),1);

    zmaxe=OT.BinBoundaries(YminElements(e),6);
    zmine=OT.BinBoundaries(YminElements(e),3);
    xmaxe=OT.BinBoundaries(YminElements(e),4);
    xmine=OT.BinBoundaries(YminElements(e),1);
    ymaxe=OT.BinBoundaries(YminElements(e),5);
    ymine=OT.BinBoundaries(YminElements(e),2);
    
    zc=0.5*(zmaxe+zmine);
    yc=0.5*(ymaxe+ymine);
    xc=0.5*(xmaxe+xmine);

% temp9=quad2d(@(x,z) (1/lx/ly).*(xc+lx/2-x).*(yc+ly/2-ymine).*(1+sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+ymine*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
% temp10=quad2d(@(x,z) (1/lx/ly).*(-xc+lx/2+x).*(yc+ly/2-ymine).*(1+sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+ymine*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);

AA=(1/lx/ly)*(yc+ly/2-ymine)*(1+sin(THETA)*sin(PHI))*exp(-1j*k0*ymine*sin(THETA)*sin(PHI));

aa=-1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
    CC_temp=@(z) (1/aa)*exp(aa*z);
    CC=CC_temp(zmaxe)-CC_temp(zmine);
else
    CC=lz;
end

bb=xc+lx/2; cc=-1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb-x)./cc).*exp(cc*x)+(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (-1/2)*(bb-x).^2; 
BB=BB_temp(xmaxe)-BB_temp(xmine);
end

temp9=AA*BB*CC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb=-xc+lx/2;
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb+x)./cc).*exp(cc*x)-(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (1/2)*(bb+x).^2; 
BB=BB_temp(xmaxe)-BB_temp(xmine);
end

temp10=AA*BB*CC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_loc_ymin=1j*k0*[temp9 ;temp10];        
b(OT.ConnectivityArrayEdge(YminElements(e),[ 9 10 ]),1)=b(OT.ConnectivityArrayEdge(YminElements(e),[ 9 10 ]),1)+b_loc_ymin(1:2);
end


% % % Back face - ymax


for e=1:length(YmaxElements)
    
    ly=OT.BinEdgeSize(YmaxElements(e),2);
    lz=OT.BinEdgeSize(YmaxElements(e),3);
    lx=OT.BinEdgeSize(YmaxElements(e),1);

    zmaxe=OT.BinBoundaries(YmaxElements(e),6);
    zmine=OT.BinBoundaries(YmaxElements(e),3);
    xmaxe=OT.BinBoundaries(YmaxElements(e),4);
    xmine=OT.BinBoundaries(YmaxElements(e),1);
    ymaxe=OT.BinBoundaries(YmaxElements(e),5);
    ymine=OT.BinBoundaries(YmaxElements(e),2);
    
    zc=0.5*(zmaxe+zmine);
    yc=0.5*(ymaxe+ymine);
    xc=0.5*(xmaxe+xmine);

% temp11=quad2d(@(x,z) (1/lx/ly).*(xc+lx/2-x).*(ymaxe-yc+ly/2).*(1-sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+ymaxe*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
% temp12=quad2d(@(x,z) (1/lx/ly).*(-xc+lx/2+x).*(ymaxe-yc+ly/2).*(1-sin(THETA)*sin(PHI)).*exp(-1j*k0*(x*sin(THETA)*cos(PHI)+ymaxe*sin(THETA)*sin(PHI)+z*cos(THETA))),xmine,xmaxe,zmine,zmaxe);
    
AA=(1/lx/ly)*(-yc+ly/2+ymaxe)*(1-sin(THETA)*sin(PHI))*exp(-1j*k0*ymaxe*sin(THETA)*sin(PHI));

aa=-1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
    CC_temp=@(z) (1/aa)*exp(aa*z);
    CC=CC_temp(zmaxe)-CC_temp(zmine);
else
    CC=lz;
end

bb=xc+lx/2; cc=-1j*k0*sin(THETA)*cos(PHI);
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb-x)./cc).*exp(cc*x)+(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (-1/2)*(bb-x).^2; 
BB=BB_temp(xmaxe)-BB_temp(xmine);
end

temp11=AA*BB*CC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb=-xc+lx/2;
if sind(THETA*180/pi)*cosd(PHI*180/pi)~=0
BB_temp=@(x) ((bb+x)./cc).*exp(cc*x)-(1/cc/cc).*exp(cc*x);
BB=BB_temp(xmaxe)-BB_temp(xmine);
else
BB_temp=@(x) (1/2)*(bb+x).^2; 
BB=BB_temp(xmaxe)-BB_temp(xmine);
end

temp12=AA*BB*CC;


    b_loc_ymax=1j*k0*[ temp11 ; temp12 ];        
    b(OT.ConnectivityArrayEdge(YmaxElements(e),[ 11 12]),1)=b(OT.ConnectivityArrayEdge(YmaxElements(e),[ 11 12 ]),1)+b_loc_ymax(1:2);
end


% left face - xmin

for e=1:length(XminElements)
    
    lz=OT.BinEdgeSize(XminElements(e),3);
    lx=OT.BinEdgeSize(XminElements(e),1);
    ly=OT.BinEdgeSize(XminElements(e),2);

    zmaxe=OT.BinBoundaries(XminElements(e),6);
    zmine=OT.BinBoundaries(XminElements(e),3);
    ymaxe=OT.BinBoundaries(XminElements(e),5);
    ymine=OT.BinBoundaries(XminElements(e),2);
    xmaxe=OT.BinBoundaries(XminElements(e),4);
    xmine=OT.BinBoundaries(XminElements(e),1);
    
    zc=0.5*(zmaxe+zmine);
    xc=0.5*(xmaxe+xmine);
    yc=0.5*(ymaxe+ymine);


%         temp9=quad2d(@(y,z) (1/ly/lx).*(yc+ly/2-y).*(xc+lx/2-xmine).*(1+sin(THETA)*cos(PHI)).*exp(-1j*k0*(xmine*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
%         temp11=quad2d(@(y,z) (1/ly/lx).*(-yc+ly/2+y).*(xc+lx/2-xmine).*(1+sin(THETA)*cos(PHI)).*exp(-1j*k0*(xmine*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);

AA=(1/lx/ly)*(xc+lx/2-xmine)*(1+sin(THETA)*cos(PHI))*exp(-1j*k0*xmine*sin(THETA)*cos(PHI));

aa=-1j*k0*cos(THETA);
if cosd(THETA*180/pi)~=0
    CC_temp=@(z) (1/aa)*exp(aa*z);
    CC=CC_temp(zmaxe)-CC_temp(zmine);
else
    CC=lz;
end

bb=yc+ly/2; cc=-1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
BB_temp=@(y) ((bb-y)./cc).*exp(cc*y)+(1/cc/cc).*exp(cc*y);
BB=BB_temp(ymaxe)-BB_temp(ymine);
else
BB_temp=@(y) (-1/2)*(bb-y).^2; 
BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp9=AA*BB*CC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb=-yc+ly/2;
if  sind(THETA*180/pi)*sind(PHI*180/pi)~=0
BB_temp=@(y) ((bb+y)./cc).*exp(cc*y)-(1/cc/cc).*exp(cc*y);
BB=BB_temp(ymaxe)-BB_temp(ymine);
else
BB_temp=@(y) (1/2)*(bb+y).^2; 
BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp11=AA*BB*CC;


    b_loc_xmin=1j*k0*[ temp9 ;temp11]; 

    b(OT.ConnectivityArrayEdge(XminElements(e),[ 9 11]),1)=b(OT.ConnectivityArrayEdge(XminElements(e),[ 9 11]),1)+b_loc_xmin(1:2);
end
% 


%Right face - xmax

for e=1:length(XmaxElements)
    
    lz=OT.BinEdgeSize(XmaxElements(e),3);
    lx=OT.BinEdgeSize(XmaxElements(e),1);
    ly=OT.BinEdgeSize(XmaxElements(e),2);

    zmaxe=OT.BinBoundaries(XmaxElements(e),6);
    zmine=OT.BinBoundaries(XmaxElements(e),3);
    ymaxe=OT.BinBoundaries(XmaxElements(e),5);
    ymine=OT.BinBoundaries(XmaxElements(e),2);
    xmaxe=OT.BinBoundaries(XmaxElements(e),4);
    xmine=OT.BinBoundaries(XmaxElements(e),1);
    
    zc=0.5*(zmaxe+zmine);
    xc=0.5*(xmaxe+xmine);
    yc=0.5*(ymaxe+ymine);


%         temp10=quad2d(@(y,z) (1/ly/lx).*(yc+ly/2-y).*(-xc+lx/2+xmaxe).*(1-sin(THETA)*cos(PHI)).*exp(-1j*k0*(xmaxe*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);
%         temp12=quad2d(@(y,z) (1/ly/lx).*(-yc+ly/2+y).*(-xc+lx/2+xmaxe).*(1-sin(THETA)*cos(PHI)).*exp(-1j*k0*(xmaxe*sin(THETA)*cos(PHI)+y*sin(THETA)*sin(PHI)+z*cos(THETA))),ymine,ymaxe,zmine,zmaxe);

AA=(1/lx/ly)*(-xc+lx/2+xmaxe)*(1-sin(THETA)*cos(PHI))*exp(-1j*k0*xmaxe*sin(THETA)*cos(PHI));

aa=-1j*k0*cos(THETA);
if  cos(THETA*180/pi)~=0
    CC_temp=@(z) (1/aa)*exp(aa*z);
    CC=CC_temp(zmaxe)-CC_temp(zmine);
else
    CC=lz;
end

bb=yc+ly/2; cc=-1j*k0*sin(THETA)*sin(PHI);
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
BB_temp=@(y) ((bb-y)./cc).*exp(cc*y)+(1/cc/cc).*exp(cc*y);
BB=BB_temp(ymaxe)-BB_temp(ymine);
else
BB_temp=@(y) (-1/2)*(bb-y).^2; 
BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp10=AA*BB*CC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bb=-yc+ly/2;
if sind(THETA*180/pi)*sind(PHI*180/pi)~=0
BB_temp=@(y) ((bb+y)./cc).*exp(cc*y)-(1/cc/cc).*exp(cc*y);
BB=BB_temp(ymaxe)-BB_temp(ymine);
else
BB_temp=@(y) (1/2)*(bb+y).^2; 
BB=BB_temp(ymaxe)-BB_temp(ymine);
end

temp12=AA*BB*CC;


    b_loc_xmax=1j*k0*[ temp10 ;temp12]; 

    b(OT.ConnectivityArrayEdge(XmaxElements(e),[ 10 12]),1)=b(OT.ConnectivityArrayEdge(XmaxElements(e),[ 10 12]),1)+b_loc_xmax(1:2);
end
% 




%% ================================ impose linear constraints on Hanging nodes ================================%%

if ~isempty(OT.HangingEdge)

tic

P=speye([Nn Nn]);
C=sparse(Nn,1);


SIZEP=size(P,1);

P(OT.HangingEdgeFace(:,1),OT.HangingEdgeFace(:,1))=0;


PP=sparse(OT.HangingEdgeFace(:,1),OT.HangingEdgeFace(:,2) ,0.5);
PP(size(P,1),size(P,2))=0;
P=P+PP;

PP=sparse(OT.HangingEdgeFace(:,1),OT.HangingEdgeFace(:,3) ,0.5);
PP(size(P,1),size(P,2))=0;
P=P+PP;

 b=P'*(b-MASS_final*C)+C;  % agar exc gheire sefr bood bayad
    MASS_final=P'*MASS_final*P+(speye(Nn,Nn)-P);
    
 P=speye([Nn Nn]);
C=sparse(Nn,1);


SIZEP=size(P,1);   

P(OT.HangingEdge(1,:),OT.HangingEdge(1,:))=0;


PP=sparse(OT.HangingEdge(1,:),OT.HangingEdge(2,:) ,1);
PP(size(P,1),size(P,2))=0;
P=P+PP;


   
    b=P'*(b-MASS_final*C)+C;  % agar exc gheire sefr bood bayad
    MASS_final=P'*MASS_final*P+(speye(Nn,Nn)-P);

    
toc
end

%% ================================ impose Dirichlet boundary conditions on PEC elements ================================%%

PEC_elements=find(OT.alpha==0);
PEC_Edge=unique(OT.ConnectivityArrayEdge(PEC_elements,:));


P=speye(Nn,Nn);

P(PEC_Edge,PEC_Edge)=0;

C=sparse(Nn,1);

 b=P'*(b-MASS_final*C)+C;  
 MASS_final=P'*MASS_final*P+(speye(Nn,Nn)-P);
 
%% ================================ Solving the system of equations ================================%%


tempEdge=find(sum(abs(MASS_final)) == 0);
% 
C=sparse(Nn,1);
P=speye(Nn,Nn);
P(tempEdge,tempEdge)=0;
b=P'*(b-MASS_final*C)+C;
MASS_final=P'*MASS_final*P+(speye(Nn,Nn)-P);
 
 
if Nn-length(tempEdge)<60000
    tic
A=MASS_final\b;
 ITER=1;

SOLVETIME=toc
else
    tic
 [A,~,~,ITER,~]=MYpcgComplex(MASS_final,b,1e-6,100000,triu(MASS_final,0),tril(MASS_final,0)); %% Preconditioned Conjugate Gradient
SOLVETIME=toc
end

DOF=Nn-length(tempEdge);



%% ================================ monostatic RCS Computation ================================%%

THETA_IN=THETA;
PHI_IN=PHI;

THETA=pi-THETA;
PHI=-pi+PHI;




RCS_Normalized2=RCS_Finder_Oblique_Analytical_fast(OT,A,XminElements,XmaxElements,YminElements,YmaxElements,ZminElements,ZmaxElements,xmin,xmax,ymin,ymax,zmin,zmax,k0,THETA,PHI,THETA_IN,PHI_IN);
str = sprintf('RCS_Normalized is  %d  by %d DOFs',  RCS_Normalized2,DOF);
disp(str)




%% plot fields
% z0=0.001;
%  plot_field(xmin,ymin,xmax,ymax,z0,A,OT) %plot in x-y plane
%  
%  y0=ymax-0.001;
% plot_fieldxz(xmin,zmin,xmax,zmax,y0,A,OT) %plot in x-z plane
% % % 
%  x0=xmax-0.001;
% plot_fieldyz(ymin,zmin,ymax,zmax,x0,A,OT)




end