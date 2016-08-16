function [no2co,el2no,el2ed,ed2no,el2fc,fc2no,BinNeighbor] = Mesh_Rectang_3D_Domain(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,hx,hy,hz)



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

pos=BinBoundaries;
% pos=OT.BinBoundaries;

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

% Nnod = size(no2co,2);
% el2no = [el2no_hex([2 6 7 8],:) el2no_hex([2 5 6 7],:) el2no_hex([2 4 7 8],:) el2no_hex([1 2 5 7],:) el2no_hex([2 3 4 7],:)  el2no_hex([1 2 3 7],:)];
% el2no = sort(el2no);
% tetramesh(el2no',no2co');xlabel('X');ylabel('Y');zlabel('Z')
% Nel = size(el2no,2);
n1 = el2no([1 4 5 8 1 5 2 6 1 2 4 3],:);
n2 = el2no([2 3 6 7 4 8 3 7 5 6 8 7 ],:);
[ed2no,~,el2ed]=unique([n1(:) n2(:)],'rows');
el2ed = reshape(el2ed , 12 , size(el2no,2) );
% Ned = size(ed2no,1);

%%  ConnectivityArrayFace

% f1 = el2no([ 5 1 1 2 3  1  ],:);
% f2 = el2no([ 6 2 3 4 4 2  ],:);
% f3 = el2no([ 7 3 5 6 7 5  ],:);
% f4 = el2no([ 8 4 7 8 8 6 ],:);

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

%% BinNeighbors

Ne=size(el2fc,2);
BinNeighbor=zeros(Ne,18);


temp=1:Ne;
temp=repmat(temp',1,6);
A=sparse(temp,el2fc',1);

[el,fc]=find(A);
m=[fc , el];
m=sortrows(m,1);

for ii=1:size(m,1)-1
    
    if m(ii,1)==m(ii+1,1) %% two elements are sharing  face m(ii,1)
        
        temp1=find(el2fc(:,m(ii,2))==m(ii,1));
        temp2=find(el2fc(:,m(ii+1,2))==m(ii,1));
        BinNeighbor(m(ii,2),temp1)=m(ii+1,2);
        BinNeighbor(m(ii+1,2),temp2)=m(ii,2);

    end
end
    


% % Update Face Neighbors
% for ii=1:Nf
%     [aa,bb]=find(el2fc==ii);
%     if length(aa)>1
%         BinNeighbor(aa(1),bb(1))=bb(2);
%         BinNeighbor(aa(2),bb(2))=bb(1);
%     end
% end
% 
% BinNeighbor=BinNeighbor';

%% Update Edge Neighbors

for ee=1:Ne
  
    %direction 7
  if BinNeighbor(ee,2) && BinNeighbor(ee,6)
    BinNeighbor(ee,7)=BinNeighbor(BinNeighbor(ee,6),2) ;
  end
  
  %direction 8
  if BinNeighbor(ee,2) && BinNeighbor(ee,4)
    BinNeighbor(ee,8)=BinNeighbor(BinNeighbor(ee,4),2) ;
  end
  
  %direction 9
  if BinNeighbor(ee,2) && BinNeighbor(ee,5)
    BinNeighbor(ee,9)=BinNeighbor(BinNeighbor(ee,5),2) ;
  end
  
  %direction 10
  if BinNeighbor(ee,2) && BinNeighbor(ee,3)
    BinNeighbor(ee,10)=BinNeighbor(BinNeighbor(ee,3),2) ;
  end
  
  %direction 11
  if BinNeighbor(ee,3) && BinNeighbor(ee,6)
    BinNeighbor(ee,11)=BinNeighbor(BinNeighbor(ee,6),3) ;
  end
  
  %direction 12
  if BinNeighbor(ee,4) && BinNeighbor(ee,6)
    BinNeighbor(ee,12)=BinNeighbor(BinNeighbor(ee,6),4) ;
  end
  
  %direction 13
  if BinNeighbor(ee,4) && BinNeighbor(ee,5)
    BinNeighbor(ee,13)=BinNeighbor(BinNeighbor(ee,5),4) ;
  end
  
  %direction 14
  if BinNeighbor(ee,3) && BinNeighbor(ee,5)
    BinNeighbor(ee,14)=BinNeighbor(BinNeighbor(ee,5),3) ;
  end
  
  %direction 15
  if BinNeighbor(ee,1) && BinNeighbor(ee,6)
    BinNeighbor(ee,15)=BinNeighbor(BinNeighbor(ee,6),1) ;
  end
  
  %direction 16
  if BinNeighbor(ee,1) && BinNeighbor(ee,4)
    BinNeighbor(ee,16)=BinNeighbor(BinNeighbor(ee,4),1) ;
  end
  
  %direction 17
  if BinNeighbor(ee,1) && BinNeighbor(ee,5)
    BinNeighbor(ee,17)=BinNeighbor(BinNeighbor(ee,5),1) ;
  end
  
  %direction 18
  if BinNeighbor(ee,1) && BinNeighbor(ee,3)
    BinNeighbor(ee,18)=BinNeighbor(BinNeighbor(ee,3),1) ;
  end
  
  
  
end

    
    
    

return;