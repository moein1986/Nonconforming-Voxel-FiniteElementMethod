
%% Non-conforming Voxel FEM is implemented and it is applied to NASA Almond test-case to compute its Radar Cross Section  

%% ================================ Initializing ================================%%

clear all
close all
clc

Freq=1.19e9; % Frequency
C_0=2.997924580003452e+08; % speed of light
Lambda=C_0/Freq; % wavelength
k0=2*pi/Lambda; %wave number


%% ================================ plot NASA Almond ================================%%

% x<0 Part


d=9.936; % 9.936 inch  almond length inch

a=.193333;
b=.064444;

t=-.41667:.41667/500:0;
psi=-pi:2*pi/500:pi;

x=zeros(length(t),length(psi));
y=zeros(length(t),length(psi));
z=zeros(length(t),length(psi));

for i=1:length(t)
    
    for j=1:length(psi)
        
        x(i,j)=d*t(i);
        y(i,j)=a*d*(sqrt(1-(t(i)/0.41667)^2))*cos(psi(j));
        z(i,j)=b*d*(sqrt(1-(t(i)/0.41667)^2))*sin(psi(j));
        
    end
end

% inch to meter
x=x*(2.54/100);
y=y*(2.54/100);
z=z*(2.54/100);

C=.5*ones(length(z));
mesh(x,y,z,C)
xlabel X; ylabel Y; zlabel Z;
hold on

xmin_Almond=min(x(:));
%% x>0 Part


a=4.83345;
b=1.61115;

t=0:.58333/500:.58333;
psi=-pi:2*pi/500:pi;

x=zeros(length(t),length(psi));
y=zeros(length(t),length(psi));
z=zeros(length(t),length(psi));

for i=1:length(t)
    
    for j=1:length(psi)
        
        x(i,j)=d*t(i);
        y(i,j)=a*d*(sqrt(1-(t(i)/2.083350)^2)-.96)*cos(psi(j));
        z(i,j)=b*d*(sqrt(1-(t(i)/2.083350)^2)-.96)*sin(psi(j));
        
    end
end

C=.5*ones(length(z));

% inch to meter
x=x*(2.54/100);
y=y*(2.54/100);
z=z*(2.54/100);

mesh(x,y,z,C)
xlabel X; ylabel Y; zlabel Z;
hold on
axis equal


xmax_Almond=max(x(:));
ymax_Almond=max(y(:));
zmax_Almond=max(z(:));

%% ================================ 3D Mesh Generation ================================%%


Density=16; % Mesh Density (Grid Size)
D_ABC=4; % distance from scaterer to ABC (Absorbing Boundary Condition) surface is Lambda/D_ABC


hx=Lambda/Density; hy=Lambda/Density; hz=Lambda/Density;


Ny=ceil((ymax_Almond+Lambda/D_ABC)/hy);
ymax=Ny*hy; ymin=-ymax;

Nz=ceil((zmax_Almond+Lambda/D_ABC)/hz);
zmax=Nz*hz; zmin=-zmax;


Nxmax=ceil((xmax_Almond+Lambda/3)/hx);
xmax=Nxmax*hx; 

Nxmin=ceil((xmin_Almond-Lambda/4)/hx);
xmin=Nxmin*hx;


pts=[xmin ymin zmin  ; xmin ymax zmin  ; xmax ymin zmin  ; xmin ymin zmax  ];


d_almond=9.936*(2.54/100); % 9.936 inch almond length

tic

OT =OcTreeMeshArbitaryGrid(pts,d_almond,hx,hy,hz,'maxDepth',2);

toc






LeafElements=find(OT.BinIsLeaf==1);
BoundaryLeaf=LeafElements(find(OT.alpha(LeafElements)==0)); % find non_uniform elements with real indexing (before elimination)


%% Extract Leaf Elements
FinalElements=find(OT.BinIsLeaf==0);
OT.ConnectivityArray(FinalElements,:)=[];
OT.ConnectivityArrayEdge(FinalElements,:)=[];

OT.BinEdgeSize(FinalElements,:)=[];
OT.alpha(FinalElements,:)=[];
OT.BinBoundaries(FinalElements,:)=[];

OT.HangingNodes(:, find(sum(OT.HangingNodes) == 0)) = [];
OT.HangingNodesFace(:, find(sum(OT.HangingNodesFace) == 0)) = [];
OT.HangingEdge(:,find(sum(OT.HangingEdge) == 0)) = [];
OT.HangingEdgeFace(find(sum(OT.HangingEdgeFace,2) == 0),:) = [];


%% ================================ Visualize 3D mesh ================================%%

figure
pos=find(OT.alpha==0); % plot the PEC object
for i=1:length(pos)
    
    voxel([OT.BinBoundaries(pos(i),1) OT.BinBoundaries(pos(i),2) OT.BinBoundaries(pos(i),3)],[OT.BinEdgeSize(pos(i),1) OT.BinEdgeSize(pos(i),2) OT.BinEdgeSize(pos(i),3)],'r',.2)
    
    
end
view (35,45)
axis equal
% 
% figure
% pos=find(OT.alpha==1); %plot the non-PEC parts
% for i=1:length(pos)
%     
%     voxel([OT.BinBoundaries(pos(i),1) OT.BinBoundaries(pos(i),2) OT.BinBoundaries(pos(i),3)],[OT.BinEdgeSize(pos(i),1) OT.BinEdgeSize(pos(i),2) OT.BinEdgeSize(pos(i),3)],'g',.01)
%     
%     
% end
%        

Nn=max(OT.ConnectivityArrayEdge(:));  %Number of Edges


%% assign material to non_uniform elements

tic
temp=find(OT.alpha==0 );

for i=1:length(temp)
     OT.alpha=MaterialFinder_CenterPoint(OT.BinBoundaries(temp(i),:),temp(i),OT.alpha,1,2.25,d_almond);
end
toc

OT.alpha(find(OT.alpha==2.25))=0; % Assign zero to PEC materials



%% ================================ Assembly (Vectorized implementation) ================================%%
tic
Nn=(max(OT.ConnectivityArrayEdge(:)));

K1=[2 -2 1 -1 ; -2 2 -1 1 ; 1 -1 2 -2 ; -1 1 -2 2];
K2=[2 1 -2 -1 ; 1 2 -1 -2 ; -2 -1 2 1 ; -1 -2 1 2];
K3=[2 1 -2 -1 ; -2 -1 2 1; 1 2 -1 -2; -1 -2 1 2 ];
K4=[4 2 2 1 ; 2 4 1 2 ; 2 1 4 2 ; 1 2 2 4];

E_loc=(1/6)*[K1+K2 -K3 -K3' ; -K3' K1+K2 -K3 ; -K3 -K3' K1+K2];
F_loc=(1/36)*[K4 zeros(4) zeros(4) ; zeros(4) K4 zeros(4) ;  zeros(4) zeros(4) K4  ];


E_final=sparse(Nn,Nn);
F_final=sparse(Nn,Nn);


OT.ConnectivityArrayEdge=OT.ConnectivityArrayEdge';

E_loc=repmat(E_loc(:),1,size(OT.ConnectivityArrayEdge,2));
F_loc=repmat(F_loc(:),1,size(OT.ConnectivityArrayEdge,2));

E_loc=E_loc.*repmat((OT.BinEdgeSize(:,1))',144,1);

F_loc=F_loc.*repmat((OT.BinEdgeSize(:,1).^3)',144,1);
F_loc=F_loc.*repmat((OT.alpha)',144,1);


row = OT.ConnectivityArrayEdge([1:12 1:12 1:12 1:12 1:12 1:12 1:12 1:12 1:12 1:12 1:12 1:12],:);
col = OT.ConnectivityArrayEdge([ones(1,12) 2*ones(1,12) 3*ones(1,12) 4*ones(1,12) 5*ones(1,12) 6*ones(1,12) 7*ones(1,12) 8*ones(1,12) 9*ones(1,12) 10*ones(1,12) 11*ones(1,12) 12*ones(1,12)],:);

E_final=sparse(row,col,E_loc);
E_loc=[];

F_final=sparse(row,col,F_loc);

% deallocate space
F_loc=[];
row=[];
col=[];

MASS_final=E_final-k0^2*F_final;

OT.ConnectivityArrayEdge=OT.ConnectivityArrayEdge';
toc

%deallocate space
E_final=[];
F_final=[];

%% ================================ monostatic RCS analysis ================================%%

THETA=pi/2; % elevation angle is zero
i=1;

for PHI=0:pi/72:pi


[RCS1,ITER,SOLVETIME]=Oblique_Almond_NVFEM(THETA,PHI,OT,xmin,xmax,ymin,ymax,zmin,zmax,MASS_final,k0);

    DATA(i,:)=[ RCS1 ITER SOLVETIME ];

    i=i+1
    
end

DATA=[(180:-2.5:0)' DATA]; % computed RCS Data 
plot(DATA(:,1),10*log10(DATA(:,2)))






