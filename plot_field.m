
function plot_field(xmin,ymin,xmax,ymax,z0,A,OT)



x=xmin:(xmax-xmin)/50/pi:xmax;
y=ymin:(xmax-xmin)/50/pi:ymax;

[x,y]=ndgrid(x,y);

ESample_x=zeros(size(x));
ESample_y=zeros(size(x));
ESample_z=zeros(size(x));



for iX = 1:size(x, 1);
    for iY = 1:size(y, 2);
x0=x(iX,iY);
y0=y(iX,iY);

        Sample_point=[x0 y0 z0];

temp=find(OT.BinBoundaries(:,1)<=x0 & OT.BinBoundaries(:,2)<=y0 & OT.BinBoundaries(:,3)<=z0 & OT.BinBoundaries(:,4)>=x0 & OT.BinBoundaries(:,5)>=y0 & OT.BinBoundaries(:,6)>=z0 );

E1=A(OT.ConnectivityArrayEdge(temp,1));
E2=A(OT.ConnectivityArrayEdge(temp,2));
E3=A(OT.ConnectivityArrayEdge(temp,3));
E4=A(OT.ConnectivityArrayEdge(temp,4));
E5=A(OT.ConnectivityArrayEdge(temp,5));
E6=A(OT.ConnectivityArrayEdge(temp,6));
E7=A(OT.ConnectivityArrayEdge(temp,7));
E8=A(OT.ConnectivityArrayEdge(temp,8));
E9=A(OT.ConnectivityArrayEdge(temp,9));
E10=A(OT.ConnectivityArrayEdge(temp,10));
E11=A(OT.ConnectivityArrayEdge(temp,11));
E12=A(OT.ConnectivityArrayEdge(temp,12));

% Center point of sample element
xc=0.5*(OT.BinBoundaries(temp,1)+OT.BinBoundaries(temp,4)) ;
yc=0.5*(OT.BinBoundaries(temp,2)+OT.BinBoundaries(temp,5)) ;
zc=0.5*(OT.BinBoundaries(temp,3)+OT.BinBoundaries(temp,6)) ;

% Edge size of sample element
lx=OT.BinEdgeSize(temp,1);
ly=OT.BinEdgeSize(temp,2);
lz=OT.BinEdgeSize(temp,3);

% 3D bases
N1=(1/ly/lz)*(yc+ly/2-y0)*(zc+lz/2-z0);
N2=(1/ly/lz)*(-yc+ly/2+y0)*(zc+lz/2-z0);
N3=(1/ly/lz)*(yc+ly/2-y0)*(-zc+lz/2+z0);
N4=(1/ly/lz)*(-yc+ly/2+y0)*(-zc+lz/2+z0);

N5=(1/lx/lz)*(xc+lx/2-x0)*(zc+lz/2-z0);
N6=(1/lx/lz)*(xc+lx/2-x0)*(-zc+lz/2+z0);
N7=(1/lx/lz)*(-xc+lx/2+x0)*(zc+lz/2-z0);
N8=(1/lx/lz)*(-xc+lx/2+x0)*(-zc+lz/2+z0);

N9=(1/lx/ly)*(xc+lx/2-x0)*(yc+ly/2-y0);
N10=(1/lx/ly)*(-xc+lx/2+x0)*(yc+ly/2-y0);
N11=(1/lx/ly)*(xc+lx/2-x0)*(-yc+ly/2+y0);
N12=(1/lx/ly)*(-xc+lx/2+x0)*(-yc+ly/2+y0);

ESample_x(iX,iY)=E1*N1+E2*N2+E3*N3+E4*N4;
ESample_y(iX,iY)=E5*N5+E6*N6+E7*N7+E8*N8;
ESample_z(iX,iY)=E9*N9+E10*N10+E11*N11+E12*N12;
    end
end
    set(figure, 'color', 'white');
imagesc(x(:,1), y(1,:), (abs(ESample_x))'); 
title('|E_x| x-y plane (V/m)')

set(figure, 'color', 'white');
imagesc(x(:,1), y(1,:), (abs(ESample_y))');
title('|E_y| x-y plane (V/m)')

set(figure, 'color', 'white');
imagesc(x(:,1), y(1,:), (abs(ESample_z))');
title('|E_z| x-y plane (V/m)')
end
