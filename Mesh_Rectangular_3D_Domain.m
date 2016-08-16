function [no2co,el2no,el2ed,ed2no,el2fc,fc2no,BinBoundaries] = Mesh_Rectangular_3D_Domain(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,Nx,Ny,Nz)

%%%%%%% inputs 
% Xmin,Xmax,Ymin,Ymax,Zmin,Zmax are the boundaries of the domain
% Nx,Ny,Nz are the number of elements in x,y and z direction.

%%%%%%% Outputs:
% no2co: node to cooridinate is 3Xn matrice, and it shows the location of the nodes
% el2no: element to node is connectivity array for nodes.
% el2ed: element to edge is connectivity array for edges.
% el2fc: element to face is connectivity array for faces.
% fc2no: face to node
% BinBoundaries shows the boundaries of each element . for example mth row
% corresponds to [xmin ymin zmin xmax ymax zmax] of element number m. 





hx=(Xmax-Xmin)/Nx; %Edge size in x-direction
hy=(Ymax-Ymin)/Ny;
hz=(Zmax-Zmin)/Nz;



%% ConnectivityArray
[x,y,z]=meshgrid( Xmin:hx:Xmax , Ymin:hy:Ymax , Zmin:hz:Zmax );
nod_num = reshape(1:size(x,1)*size(x,2)*size(x,3),size(x,1),size(x,2),size(x,3));
Nn1 = nod_num(1:end-1,1:end-1,1:end-1);
Nn2 = nod_num(1:end-1,2:end  ,1:end-1);
Nn3 = nod_num(2:end  ,1:end-1,1:end-1);
Nn4 = nod_num(2:end  ,2:end  ,1:end-1);
Nn5 = nod_num(1:end-1,1:end-1,2:end  );
Nn6 = nod_num(1:end-1,2:end  ,2:end  );
Nn7 = nod_num(2:end  ,1:end-1,2:end  );
Nn8 = nod_num(2:end  ,2:end  ,2:end  );
el2no_hex = [Nn1(:) Nn2(:) Nn4(:) Nn3(:) Nn5(:) Nn6(:) Nn8(:) Nn7(:)].';

temp = [nod_num(:) x(:) y(:) z(:)];
no2co = sortrows(temp,1).';
no2co = no2co(2:4,:);
el2no=el2no_hex;

%% BinBoundaries
xxmin=no2co(1,el2no(1,:))';
xxmax=no2co(1,el2no(2,:))';
yymin=no2co(2,el2no(1,:))';
yymax=no2co(2,el2no(3,:))';
zzmin=no2co(3,el2no(1,:))';
zzmax=no2co(3,el2no(5,:))';

BinBoundaries=[xxmin yymin zzmin xxmax yymax zzmax];

%% plot the mesh

pos=BinBoundaries;

for e=1:size(pos,1)
    line([pos(e,1) pos(e,4)],[pos(e,2) pos(e,2)],[pos(e,3) pos(e,3)])
    line([pos(e,4) pos(e,4)],[pos(e,2) pos(e,5)],[pos(e,3) pos(e,3)])
    line([pos(e,4) pos(e,1)],[pos(e,5) pos(e,5)],[pos(e,3) pos(e,3)])
    line([pos(e,1) pos(e,1)],[pos(e,5) pos(e,2)],[pos(e,3) pos(e,3)])

    line([pos(e,1) pos(e,4)],[pos(e,2) pos(e,2)],[pos(e,6) pos(e,6)])
    line([pos(e,4) pos(e,4)],[pos(e,2) pos(e,5)],[pos(e,6) pos(e,6)])
    line([pos(e,4) pos(e,1)],[pos(e,5) pos(e,5)],[pos(e,6) pos(e,6)])
    line([pos(e,1) pos(e,1)],[pos(e,5) pos(e,2)],[pos(e,6) pos(e,6)])

    line([pos(e,1) pos(e,1)],[pos(e,2) pos(e,2)],[pos(e,3) pos(e,6)])
    line([pos(e,4) pos(e,4)],[pos(e,2) pos(e,2)],[pos(e,3) pos(e,6)])
    line([pos(e,1) pos(e,1)],[pos(e,5) pos(e,5)],[pos(e,3) pos(e,6)])
    line([pos(e,4) pos(e,4)],[pos(e,5) pos(e,5)],[pos(e,3) pos(e,6)])

        text(0.5*(pos(e,1)+pos(e,4)),0.5*(pos(e,2)+pos(e,5)),0.5*(pos(e,3)+pos(e,6)),['',num2str(e)],'FontSize',6)
end

%% BinEdge

BinEdgeSize=zeros(size(BinBoundaries,1),3);
BinEdgeSize(:,1)=hx; BinEdgeSize(:,2)=hy; BinEdgeSize(:,3)=hz; 
 
%% ConnectivtyArray Edge


n1 = el2no([1 4 5 8 1 5 2 6 1 2 4 3],:);
n2 = el2no([2 3 6 7 4 8 3 7 5 6 8 7 ],:);
[ed2no,~,el2ed]=unique([n1(:) n2(:)],'rows');
el2ed = reshape(el2ed , 12 , size(el2no,2) );

%%  ConnectivityArrayFace

            f1 = el2no([ 5 1 1 2 4  1  ],:);
            f2 = el2no([ 6 2 4 3 3 2  ],:);
            f3 = el2no([ 8 4 5 6 8 5  ],:);
            f4 = el2no([ 7 3 8 7 7 6 ],:);
            
temp=[f1(:) f2(:) f3(:) f4(:)];
tempFC=sort(temp.').';
[fc2no,b,el2fc]=unique(tempFC,'rows');
fc2no=temp(b,:).';
Nf = size(fc2no,2);
el2fc = reshape(el2fc , 6, size(el2no,2) ); % faces constructing each element


return;