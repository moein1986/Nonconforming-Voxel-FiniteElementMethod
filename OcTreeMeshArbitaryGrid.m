classdef OcTreeMeshArbitaryGrid< handle
    
    %   Created by Moein Nazari
    %   1.0     - 2013-06 Initial Release
    
    
    
    
    properties
        Points;
        
        
        PointBins;
        BinCount;
        PointCount=0;
        EdgeCount=0;
        
        BinBoundaries=zeros(2.5*10^6,6);
        BinDepths;
        BinEdgeSize=zeros(2.5*10^6,3);
        BinIsLeaf=ones(2.5*10^6,1);
        BinParents = zeros(2.5*10^6,1);
        HangingNodes=zeros(5,3.8*10^6); % [  master points   HN number  ]
        HangingNodesFace=zeros(5,3.8*10^6);
        HangingEdge=zeros(2,3.8*10^6);
        HangingEdgeFace=zeros(3.8*10^6,3); %[HN Edge , 2 masterEdge]
        
        BinChildren=zeros(2.5*10^6,8);
        ConnectivityArray=zeros(2.5*10^6,8);
        ConnectivityArrayEdge=zeros(2.5*10^6,12);
        
        alpha=zeros(2.5*10^6,1);
        BinNeighbor=zeros(2.5*10^6,18); % Face neighbors [up Down back front left right]
        Properties;
    end
    
    methods
        
        function this = OcTreeMeshArbitaryGrid(pts,radius,hx,hy,hz,varargin)
            
            
            
            
            % Initialise a single bin surrounding all given points
            numPts = size(pts,1);
            this.BinBoundaries(1,:) = [min(pts,[],1) max(pts,[],1)]; %[xmin ymin zmin xmax ymax zmax]
            this.Points = pts;
            
            
            this.PointBins = ones(numPts,1);
            this.BinDepths = 0;
            this.BinParents(1) = 0;
            this.BinCount = 0;
            this.ConnectivityArray(1,:)=[1 2 3 4 5 6 7 8];
            this.ConnectivityArrayEdge(1,:)=[1 2 3 4 5 6 7 8 9 10 11 12];
            
            % Allow custom setting of Properties
            IP = inputParser;
            %             IP.addParamValue('binCapacity',ceil(numPts)/20);
            IP.addParamValue('maxDepth',inf);
            
            %             IP.addParamValue('maxSize',inf);
            %             IP.addParamValue('minSize',1000 * eps);
            %             IP.addParamValue('style','equal');
            
            
            
            IP.parse(varargin{:});
            this.Properties = IP.Results;
            
            
            
            
            % Return on empty or trivial bins
            if numPts<2, return; end
            
            % Start dividing!
            
            this.preallocateSpace;
            this.divideGrid(hx,hy,hz);
            SB= find(this.BinIsLeaf(1:this.BinCount));
            for i=1:length(SB)
                UniformFlag(this,SB(i),1,2.25,radius)
            end
            this.divide(SB,radius);
            
            %%
            
            
            this.deallocateSpace;
            %             this.BinChildren = arrayfun(@(i)find(this.BinParents==i),1:this.BinCount,'Un',0)';
            %             this.binIsLeaf = cellfun(@isempty, this.BinChildren);
            %             new=[2.5 2.5 0 ; 5 2.5 0 ;  2.5 5 0 ; 5 5 0 ;2.5 2.5 2.5 ; 5 2.5 2.5 ;  2.5 5 2.5 ; 5 5 2.5 ];
            
        end
        
        % MATLAB performs better if arrays that grow are initialised,
        % rather than grown during a loop. These two functions do just that
        % before and after the identification of new beens.
        
        function preallocateSpace(this)
            %             numPts = size(this.Points,1);
            %             numBins = numPts;
            %             if isfinite(this.Properties.binCapacity)
            %                 numBins = ceil(200*numPts/this.Properties.binCapacity
            numBins = 2.5*10e6;
            
            %             end
            this.BinDepths(numBins) = 0;
            this.BinParents(numBins) = 0;
            %             this.BinBoundaries(numBins,1) = 0;
        end
        function deallocateSpace(this)
            this.BinDepths(this.BinCount+1:end) = [];
            this.BinParents(this.BinCount+1:end) = [];
            this.BinChildren(this.BinCount+1:end,:) = [];
            this.BinBoundaries(this.BinCount+1:end,:) = [];
            this.BinNeighbor(this.BinCount+1:end,:) = [];
            this.BinIsLeaf(this.BinCount+1:end) = [];
            this.BinEdgeSize(this.BinCount+1:end,:) = [];
            this.ConnectivityArray(this.BinCount+1:end,:)=[];
            this.ConnectivityArrayEdge(this.BinCount+1:end,:)=[];
            
            this.HangingNodes(:,this.PointCount+1:end)=[];
            this.HangingNodesFace(:,this.PointCount+1:end)=[];
            
            %             this.HangingNodes(:, find(sum(this.HangingNodes) == 0)) = [];
            %             this.HangingNodesFace(:, find(sum(this.HangingNodesFace) == 0)) = [];
            
            this.HangingEdge(:,this.EdgeCount+1:end)=[];
            %             this.HangingEdge(:,find(sum(this.HangingEdge) == 0)) = [];
            
            this.HangingEdgeFace(this.EdgeCount+1:end,:)=[];
            %             this.HangingEdgeFace(find(sum(this.HangingEdgeFace,2) == 0),:) = [];
            
            
            this.alpha(this.BinCount+1:end)=[];
            
        end
        
        function divideGrid(this,hx,hy,hz)
            
            Xmin=this.Points(1,1);  Ymin=this.Points(1,2);  Zmin=this.Points(3,3);
            Xmax=this.Points(3,1);  Ymax=this.Points(2,2);  Zmax=this.Points(4,3);
            
            
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
            
            BinBoundaries1=[xxmin yymin zzmin xxmax yymax zzmax];
            
            % pos=BinBoundaries;
            % % pos=OT.BinBoundaries;
            %
            % for e=1:size(pos,1)
            %     line([pos(e,1) pos(e,4)],[pos(e,2) pos(e,2)],[pos(e,3) pos(e,3)])
            %     line([pos(e,4) pos(e,4)],[pos(e,2) pos(e,5)],[pos(e,3) pos(e,3)])
            %     line([pos(e,4) pos(e,1)],[pos(e,5) pos(e,5)],[pos(e,3) pos(e,3)])
            %     line([pos(e,1) pos(e,1)],[pos(e,5) pos(e,2)],[pos(e,3) pos(e,3)])
            %
            %     line([pos(e,1) pos(e,4)],[pos(e,2) pos(e,2)],[pos(e,6) pos(e,6)])
            %     line([pos(e,4) pos(e,4)],[pos(e,2) pos(e,5)],[pos(e,6) pos(e,6)])
            %     line([pos(e,4) pos(e,1)],[pos(e,5) pos(e,5)],[pos(e,6) pos(e,6)])
            %     line([pos(e,1) pos(e,1)],[pos(e,5) pos(e,2)],[pos(e,6) pos(e,6)])
            %
            %     line([pos(e,1) pos(e,1)],[pos(e,2) pos(e,2)],[pos(e,3) pos(e,6)])
            %     line([pos(e,4) pos(e,4)],[pos(e,2) pos(e,2)],[pos(e,3) pos(e,6)])
            %     line([pos(e,1) pos(e,1)],[pos(e,5) pos(e,5)],[pos(e,3) pos(e,6)])
            %     line([pos(e,4) pos(e,4)],[pos(e,5) pos(e,5)],[pos(e,3) pos(e,6)])
            %
            % %         text(0.5*(pos(e,1)+pos(e,4)),0.5*(pos(e,2)+pos(e,5)),0.5*(pos(e,3)+pos(e,6)),['',num2str(e)],'FontSize',6)
            % end
            
            %% BinEdge
            
            BinEdgeSize1=zeros(size(BinBoundaries1,1),3);
            BinEdgeSize1(:,1)=hx; BinEdgeSize1(:,2)=hy; BinEdgeSize1(:,3)=hz;
            
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
            
            %% BinNeighbors
            
            Ne=size(el2fc,2);
            BinNeighbor1=zeros(Ne,18);
            
            
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
                    BinNeighbor1(m(ii,2),temp1)=m(ii+1,2);
                    BinNeighbor1(m(ii+1,2),temp2)=m(ii,2);
                    
                end
            end
            
            
            
            % % Update Face Neighbors
            % for ii=1:Nf
            %     [aa,bb]=find(el2fc==ii);
            %     if length(aa)>1
            %         BinNeighbor1(aa(1),bb(1))=bb(2);
            %         BinNeighbor1(aa(2),bb(2))=bb(1);
            %     end
            % end
            %
            % BinNeighbor1=BinNeighbor1';
            
            %% Update Edge Neighbors
            
            for ee=1:Ne
                
                %direction 7
                if BinNeighbor1(ee,2) && BinNeighbor1(ee,6)
                    BinNeighbor1(ee,7)=BinNeighbor1(BinNeighbor1(ee,6),2) ;
                end
                
                %direction 8
                if BinNeighbor1(ee,2) && BinNeighbor1(ee,4)
                    BinNeighbor1(ee,8)=BinNeighbor1(BinNeighbor1(ee,4),2) ;
                end
                
                %direction 9
                if BinNeighbor1(ee,2) && BinNeighbor1(ee,5)
                    BinNeighbor1(ee,9)=BinNeighbor1(BinNeighbor1(ee,5),2) ;
                end
                
                %direction 10
                if BinNeighbor1(ee,2) && BinNeighbor1(ee,3)
                    BinNeighbor1(ee,10)=BinNeighbor1(BinNeighbor1(ee,3),2) ;
                end
                
                %direction 11
                if BinNeighbor1(ee,3) && BinNeighbor1(ee,6)
                    BinNeighbor1(ee,11)=BinNeighbor1(BinNeighbor1(ee,6),3) ;
                end
                
                %direction 12
                if BinNeighbor1(ee,4) && BinNeighbor1(ee,6)
                    BinNeighbor1(ee,12)=BinNeighbor1(BinNeighbor1(ee,6),4) ;
                end
                
                %direction 13
                if BinNeighbor1(ee,4) && BinNeighbor1(ee,5)
                    BinNeighbor1(ee,13)=BinNeighbor1(BinNeighbor1(ee,5),4) ;
                end
                
                %direction 14
                if BinNeighbor1(ee,3) && BinNeighbor1(ee,5)
                    BinNeighbor1(ee,14)=BinNeighbor1(BinNeighbor1(ee,5),3) ;
                end
                
                %direction 15
                if BinNeighbor1(ee,1) && BinNeighbor1(ee,6)
                    BinNeighbor1(ee,15)=BinNeighbor1(BinNeighbor1(ee,6),1) ;
                end
                
                %direction 16
                if BinNeighbor1(ee,1) && BinNeighbor1(ee,4)
                    BinNeighbor1(ee,16)=BinNeighbor1(BinNeighbor1(ee,4),1) ;
                end
                
                %direction 17
                if BinNeighbor1(ee,1) && BinNeighbor1(ee,5)
                    BinNeighbor1(ee,17)=BinNeighbor1(BinNeighbor1(ee,5),1) ;
                end
                
                %direction 18
                if BinNeighbor1(ee,1) && BinNeighbor1(ee,3)
                    BinNeighbor1(ee,18)=BinNeighbor1(BinNeighbor1(ee,3),1) ;
                end
                
            end  
                
            
            
            %%Inserting the Grid Data to Octree Data sturucture
            this.BinParents(1:Ne,1)=0;
            this.BinCount=Ne; % we asumme that element 1 is the parent of all grid elements
            this.EdgeCount=max(el2ed(:));
            this.PointCount=max(el2no(:));
            this.BinBoundaries(1:Ne,:)=BinBoundaries1;
            this.ConnectivityArrayEdge(1:Ne,:)=el2ed';
            this.ConnectivityArray(1:Ne,:)=el2no';
            this.BinNeighbor(1:Ne,:)=BinNeighbor1;
            this.BinEdgeSize(1:Ne,:)=BinEdgeSize1;
            
            
            
            
        end
        
        

    function divide(this, startingBins,radius)
    
    for i = 1:length(startingBins)
        binNo = startingBins(i);
        
        %                 if this.BinDepths(binNo)==this.Properties.maxDepth-1  % modify materials for finest leaf elements
        UniformFlag(this,binNo,1,2.25,radius);
        %                 end
        % checkMaterial
        %                 BinaryMask(startingBins(i))=1;
        this.BinEdgeSize(binNo,1) = this.BinBoundaries(binNo,4)-this.BinBoundaries(binNo,1);
        this.BinEdgeSize(binNo,2) = this.BinBoundaries(binNo,5)-this.BinBoundaries(binNo,2);
        this.BinEdgeSize(binNo,3) = this.BinBoundaries(binNo,6)-this.BinBoundaries(binNo,3);
        
        
        % Prevent dividing beyond the maximum depth
        %                 if this.BinParents(binNo)+1 >= this.Properties.maxDepth
        if this.BinDepths(binNo) >= this.Properties.maxDepth
            % leaf elements
            
            continue;
        end
        
        if ~this.alpha(binNo)
            if this.BinIsLeaf(binNo)
                this.divideBin(binNo,radius);
            end
            this.divide(this.BinChildren(binNo,:),radius);
            continue;
        end
        
    end
    end
    
    
    
    
    
    function divideBin(this,binNo,radius)
    % Gather the new points (a bit more efficient to copy once)
    this.BinCount = this.BinCount + 8;
    
    
    this.BinIsLeaf(binNo)=0;
    
    % Get the old corner points and the new division point
    oldMin = this.BinBoundaries(binNo,1:3);
    oldMax = this.BinBoundaries(binNo,4:6);
    
    newDiv = mean([oldMin; oldMax], 1);
    %
    
    % Build the new boundaries of our 8 subdivisions
    minMidMax = [oldMin newDiv oldMax];
    newBounds = minMidMax([...
        1 2 3 4 5 6;
        1 2 6 4 5 9;
        1 5 3 4 8 6;
        1 5 6 4 8 9;
        4 2 3 7 5 6;
        4 2 6 7 5 9;
        4 5 3 7 8 6;
        4 5 6 7 8 9]);
    
    
    newBinInds = this.BinCount-7:this.BinCount;
    this.BinEdgeSize(newBinInds,1)=0.5*this.BinEdgeSize(binNo,1);
    this.BinEdgeSize(newBinInds,2)=0.5*this.BinEdgeSize(binNo,2);
    this.BinEdgeSize(newBinInds,3)=0.5*this.BinEdgeSize(binNo,3);
    
    this.BinChildren(binNo,:)=newBinInds;
    this.BinBoundaries(newBinInds,:) = newBounds;
    this.BinDepths(newBinInds) = this.BinDepths(binNo)+1;
    this.BinParents(newBinInds) = binNo;
    
    UniformFlag(this,newBinInds(1),1,2.25,radius);
    UniformFlag(this,newBinInds(2),1,2.25,radius);
    UniformFlag(this,newBinInds(3),1,2.25,radius);
    UniformFlag(this,newBinInds(4),1,2.25,radius);            %%
    UniformFlag(this,newBinInds(5),1,2.25,radius);
    UniformFlag(this,newBinInds(6),1,2.25,radius);
    UniformFlag(this,newBinInds(7),1,2.25,radius);
    UniformFlag(this,newBinInds(8),1,2.25,radius);
    
    %update Face and Edge Neighbor Matrix
    update_neighbor(this,binNo,radius);
    
    %% plot
    %                         pts=[0 0 0 ; 0 1 0 ; 1 0 0 ; 0 0 1 ];
    %             hold on
    %                         boxH = this.plot;
    %                    cols = lines(this.BinCount);
    %                    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
    %                    for i = 1:this.BinCount
    %                         Min = this.BinBoundaries(i,1:3);
    %                         Max = this.BinBoundaries(i,4:6);
    %                        A=mean([Min; Max], 1);
    %                        text(A(1),A(2),A(3),['',num2str(i)],'FontSize',6)
    %
    %                        set(boxH(i),'Color',cols(i,:),'LineWidth', 1+this.BinDepths(i))
    %                        set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
    %
    %                        set(boxH(i),'Color',cols(i,:))
    %
    %                        doplot3(pts(this.PointBins==i,:),'.','Color',cols(i,:))
    %                    end
    %                   axis image, view(3)
    %             %
    % %
    
    %% Update Connectivity array
    % add CenterPoint
    this.PointCount=this.PointCount+1;
    this.ConnectivityArray(newBinInds(1),7)=this.PointCount;
    this.ConnectivityArray(newBinInds(2),3)=this.PointCount;
    this.ConnectivityArray(newBinInds(3),6)=this.PointCount;
    this.ConnectivityArray(newBinInds(4),2)=this.PointCount;
    this.ConnectivityArray(newBinInds(5),8)=this.PointCount;
    this.ConnectivityArray(newBinInds(6),4)=this.PointCount;
    this.ConnectivityArray(newBinInds(7),5)=this.PointCount;
    this.ConnectivityArray(newBinInds(8),1)=this.PointCount;
    
    % update corner nodes
    this.ConnectivityArray(newBinInds(1),1)=this.ConnectivityArray(this.BinParents(newBinInds(1)),1);
    this.ConnectivityArray(newBinInds(2),5)=this.ConnectivityArray(this.BinParents(newBinInds(2)),5);
    this.ConnectivityArray(newBinInds(3),4)=this.ConnectivityArray(this.BinParents(newBinInds(3)),4);
    this.ConnectivityArray(newBinInds(4),8)=this.ConnectivityArray(this.BinParents(newBinInds(4)),8);
    this.ConnectivityArray(newBinInds(5),2)=this.ConnectivityArray(this.BinParents(newBinInds(5)),2);
    this.ConnectivityArray(newBinInds(6),6)=this.ConnectivityArray(this.BinParents(newBinInds(6)),6);
    this.ConnectivityArray(newBinInds(7),3)=this.ConnectivityArray(this.BinParents(newBinInds(7)),3);
    this.ConnectivityArray(newBinInds(8),7)=this.ConnectivityArray(this.BinParents(newBinInds(8)),7);
    
    
    % update ConnectivityArray for edges inside the divided element
    
    this.ConnectivityArrayEdge(newBinInds(1),12)=this.EdgeCount+1;
    this.ConnectivityArrayEdge(newBinInds(3),10)=this.EdgeCount+1;
    this.ConnectivityArrayEdge(newBinInds(5),11)=this.EdgeCount+1;
    this.ConnectivityArrayEdge(newBinInds(7),9)=this.EdgeCount+1;
    
    
    this.ConnectivityArrayEdge(newBinInds(2),12)=this.EdgeCount+2;
    this.ConnectivityArrayEdge(newBinInds(4),10)=this.EdgeCount+2;
    this.ConnectivityArrayEdge(newBinInds(6),11)=this.EdgeCount+2;
    this.ConnectivityArrayEdge(newBinInds(8),9)=this.EdgeCount+2;
    
    this.ConnectivityArrayEdge(newBinInds(1),8)=this.EdgeCount+3;
    this.ConnectivityArrayEdge(newBinInds(2),7)=this.EdgeCount+3;
    this.ConnectivityArrayEdge(newBinInds(5),6)=this.EdgeCount+3;
    this.ConnectivityArrayEdge(newBinInds(6),5)=this.EdgeCount+3;
    
    this.ConnectivityArrayEdge(newBinInds(3),8)=this.EdgeCount+4;
    this.ConnectivityArrayEdge(newBinInds(4),7)=this.EdgeCount+4;
    this.ConnectivityArrayEdge(newBinInds(7),6)=this.EdgeCount+4;
    this.ConnectivityArrayEdge(newBinInds(8),5)=this.EdgeCount+4;
    
    this.ConnectivityArrayEdge(newBinInds(1),4)=this.EdgeCount+5;
    this.ConnectivityArrayEdge(newBinInds(2),2)=this.EdgeCount+5;
    this.ConnectivityArrayEdge(newBinInds(3),3)=this.EdgeCount+5;
    this.ConnectivityArrayEdge(newBinInds(4),1)=this.EdgeCount+5;
    
    this.ConnectivityArrayEdge(newBinInds(5),4)=this.EdgeCount+6;
    this.ConnectivityArrayEdge(newBinInds(6),2)=this.EdgeCount+6;
    this.ConnectivityArrayEdge(newBinInds(7),3)=this.EdgeCount+6;
    this.ConnectivityArrayEdge(newBinInds(8),1)=this.EdgeCount+6;
    
    this.EdgeCount=this.EdgeCount+6;
    
    
    
    % update center nodes
    
    updateConArrayCenterNodes(this,newBinInds,binNo);
    
    %             arrayfun(@(direction) updateConArrayCenterNodes(this,newBinInds,binNo,direction) , 1:18);
    
    
    
    end
    
    
    
    function updateConArrayCenterNodes(this,newBinInds,binNo)
    
    for direction=1:18
        switch direction
            case 1 % up
                if this.BinNeighbor(binNo,direction)
                    if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.PointCount=this.PointCount+1;
                        this.ConnectivityArray(newBinInds(2),7)=this.PointCount;
                        this.ConnectivityArray(newBinInds(4),6)=this.PointCount;
                        this.ConnectivityArray(newBinInds(6),8)=this.PointCount;
                        this.ConnectivityArray(newBinInds(8),5)=this.PointCount;
                        
                        this.ConnectivityArrayEdge(newBinInds(2),8)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(4),3)=this.EdgeCount+2;
                        this.ConnectivityArrayEdge(newBinInds(6),4)=this.EdgeCount+3;
                        this.ConnectivityArrayEdge(newBinInds(8),6)=this.EdgeCount+4;
                        
                        this.ConnectivityArrayEdge(newBinInds(6),6)=this.ConnectivityArrayEdge(newBinInds(2),8);
                        this.ConnectivityArrayEdge(newBinInds(2),4)=this.ConnectivityArrayEdge(newBinInds(4),3);
                        this.ConnectivityArrayEdge(newBinInds(8),3)=this.ConnectivityArrayEdge(newBinInds(6),4);
                        this.ConnectivityArrayEdge(newBinInds(4),8)=this.ConnectivityArrayEdge(newBinInds(8),6);
                        this.EdgeCount=this.EdgeCount+4;
                        
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,5) ; this.ConnectivityArray(binNo,6); this.ConnectivityArray(binNo,7); this.ConnectivityArray(binNo,8)  ];
                        %                                 this.HangingNodesFace(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,5) ; this.ConnectivityArray(binNo,6); this.ConnectivityArray(binNo,7); this.ConnectivityArray(binNo,8)  ];
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(6),6),:)=[this.ConnectivityArrayEdge(newBinInds(6),6) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),7) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),5) ];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(2),4),:)=[this.ConnectivityArrayEdge(newBinInds(2),4) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),1) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),2)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(8),3),:)=[this.ConnectivityArrayEdge(newBinInds(8),3)  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),1) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),2)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(4),8),:)=[this.ConnectivityArrayEdge(newBinInds(4),8) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),7) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),5)];
                        
                    else
                        this.ConnectivityArray(newBinInds(2),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),3);
                        this.ConnectivityArray(newBinInds(4),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),3);
                        this.ConnectivityArray(newBinInds(6),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),3);
                        this.ConnectivityArray(newBinInds(8),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),3);
                        %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),3))=0;
                        %                                 this.HangingNodesFace(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),3))=0;
                        
                        this.ConnectivityArrayEdge(newBinInds(2),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),7);
                        this.ConnectivityArrayEdge(newBinInds(4),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),1);
                        this.ConnectivityArrayEdge(newBinInds(6),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),2);
                        this.ConnectivityArrayEdge(newBinInds(8),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),5);
                        
                        this.ConnectivityArrayEdge(newBinInds(6),6)=this.ConnectivityArrayEdge(newBinInds(2),8);
                        this.ConnectivityArrayEdge(newBinInds(2),4)=this.ConnectivityArrayEdge(newBinInds(4),3);
                        this.ConnectivityArrayEdge(newBinInds(8),3)=this.ConnectivityArrayEdge(newBinInds(6),4);
                        this.ConnectivityArrayEdge(newBinInds(4),8)=this.ConnectivityArrayEdge(newBinInds(8),6);
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),7),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),1),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),2),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),5),:)=0;
                        
                    end
                else
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(2),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),5)=this.PointCount;
                    
                    this.ConnectivityArrayEdge(newBinInds(2),8)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(4),3)=this.EdgeCount+2;
                    this.ConnectivityArrayEdge(newBinInds(6),4)=this.EdgeCount+3;
                    this.ConnectivityArrayEdge(newBinInds(8),6)=this.EdgeCount+4;
                    
                    this.ConnectivityArrayEdge(newBinInds(6),6)=this.ConnectivityArrayEdge(newBinInds(2),8);
                    this.ConnectivityArrayEdge(newBinInds(2),4)=this.ConnectivityArrayEdge(newBinInds(4),3);
                    this.ConnectivityArrayEdge(newBinInds(8),3)=this.ConnectivityArrayEdge(newBinInds(6),4);
                    this.ConnectivityArrayEdge(newBinInds(4),8)=this.ConnectivityArrayEdge(newBinInds(8),6);
                    
                    this.EdgeCount=this.EdgeCount+4;
                    
                end
                
            case 2 % down
                if this.BinNeighbor(binNo,direction)
                    if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        
                        this.PointCount=this.PointCount+1;
                        this.ConnectivityArray(newBinInds(1),3)=this.PointCount;
                        this.ConnectivityArray(newBinInds(3),2)=this.PointCount;
                        this.ConnectivityArray(newBinInds(5),4)=this.PointCount;
                        this.ConnectivityArray(newBinInds(7),1)=this.PointCount;
                        %                                 this.HangingNodesFace(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1) ; this.ConnectivityArray(binNo,2) ;this.ConnectivityArray(binNo,3) ;this.ConnectivityArray(binNo,4)  ];
                        
                        this.ConnectivityArrayEdge(newBinInds(1),7)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(3),1)=this.EdgeCount+2;
                        this.ConnectivityArrayEdge(newBinInds(5),2)=this.EdgeCount+3;
                        this.ConnectivityArrayEdge(newBinInds(7),5)=this.EdgeCount+4;
                        
                        this.ConnectivityArrayEdge(newBinInds(5),5)=this.ConnectivityArrayEdge(newBinInds(1),7);
                        this.ConnectivityArrayEdge(newBinInds(1),2)=this.ConnectivityArrayEdge(newBinInds(3),1);
                        this.ConnectivityArrayEdge(newBinInds(7),1)=this.ConnectivityArrayEdge(newBinInds(5),2);
                        this.ConnectivityArrayEdge(newBinInds(3),7)=this.ConnectivityArrayEdge(newBinInds(7),5);
                        this.EdgeCount=this.EdgeCount+4;
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(5),5),:)=[this.ConnectivityArrayEdge(newBinInds(5),5) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),6) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),8) ];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(1),2),:)=[this.ConnectivityArrayEdge(newBinInds(1),2) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),3) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),4)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(7),1),:)=[this.ConnectivityArrayEdge(newBinInds(7),1)  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),3) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),4)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(3),7),:)=[this.ConnectivityArrayEdge(newBinInds(3),7) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),6) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),8)];
                        
                        
                        
                    else
                        this.ConnectivityArray(newBinInds(1),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),8);
                        this.ConnectivityArray(newBinInds(3),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),8);
                        this.ConnectivityArray(newBinInds(5),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),8);
                        this.ConnectivityArray(newBinInds(7),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),8);
                        %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),8))=0;
                        %                                 this.HangingNodesFace(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),8))=0;
                        
                        this.ConnectivityArrayEdge(newBinInds(1),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),8);
                        this.ConnectivityArrayEdge(newBinInds(3),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),3);
                        this.ConnectivityArrayEdge(newBinInds(5),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),4);
                        this.ConnectivityArrayEdge(newBinInds(7),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),6);
                        
                        this.ConnectivityArrayEdge(newBinInds(5),5)=this.ConnectivityArrayEdge(newBinInds(1),7);
                        this.ConnectivityArrayEdge(newBinInds(1),2)=this.ConnectivityArrayEdge(newBinInds(3),1);
                        this.ConnectivityArrayEdge(newBinInds(7),1)=this.ConnectivityArrayEdge(newBinInds(5),2);
                        this.ConnectivityArrayEdge(newBinInds(3),7)=this.ConnectivityArrayEdge(newBinInds(7),5);
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),8),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),3),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),4),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),6),:)=0;
                        
                    end
                else
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(3),2)=this.PointCount;
                    this.ConnectivityArray(newBinInds(5),4)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),1)=this.PointCount;
                    
                    this.ConnectivityArrayEdge(newBinInds(1),7)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(3),1)=this.EdgeCount+2;
                    this.ConnectivityArrayEdge(newBinInds(5),2)=this.EdgeCount+3;
                    this.ConnectivityArrayEdge(newBinInds(7),5)=this.EdgeCount+4;
                    
                    this.ConnectivityArrayEdge(newBinInds(5),5)=this.ConnectivityArrayEdge(newBinInds(1),7);
                    this.ConnectivityArrayEdge(newBinInds(1),2)=this.ConnectivityArrayEdge(newBinInds(3),1);
                    this.ConnectivityArrayEdge(newBinInds(7),1)=this.ConnectivityArrayEdge(newBinInds(5),2);
                    this.ConnectivityArrayEdge(newBinInds(3),7)=this.ConnectivityArrayEdge(newBinInds(7),5);
                    this.EdgeCount=this.EdgeCount+4;
                end
                
                
            case 3 % Back
                if this.BinNeighbor(binNo,direction)
                    if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.PointCount=this.PointCount+1;
                        this.ConnectivityArray(newBinInds(1),8)=this.PointCount;
                        this.ConnectivityArray(newBinInds(2),4)=this.PointCount;
                        this.ConnectivityArray(newBinInds(3),5)=this.PointCount;
                        this.ConnectivityArray(newBinInds(4),1)=this.PointCount;
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1) ; this.ConnectivityArray(binNo,5); this.ConnectivityArray(binNo,4); this.ConnectivityArray(binNo,8)  ];
                        %                                 this.HangingNodesFace(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1) ; this.ConnectivityArray(binNo,5); this.ConnectivityArray(binNo,4); this.ConnectivityArray(binNo,8)  ];
                        
                        this.ConnectivityArrayEdge(newBinInds(1),6)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(2),11)=this.EdgeCount+2;
                        this.ConnectivityArrayEdge(newBinInds(3),9)=this.EdgeCount+3;
                        this.ConnectivityArrayEdge(newBinInds(4),5)=this.EdgeCount+4;
                        
                        this.ConnectivityArrayEdge(newBinInds(2),5)=this.ConnectivityArrayEdge(newBinInds(1),6);
                        this.ConnectivityArrayEdge(newBinInds(4),9)=this.ConnectivityArrayEdge(newBinInds(2),11);
                        this.ConnectivityArrayEdge(newBinInds(1),11)=this.ConnectivityArrayEdge(newBinInds(3),9);
                        this.ConnectivityArrayEdge(newBinInds(3),6)=this.ConnectivityArrayEdge(newBinInds(4),5);
                        this.EdgeCount=this.EdgeCount+4;
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(2),5),:)=[this.ConnectivityArrayEdge(newBinInds(2),5) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),7) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),8) ];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(4),9),:)=[this.ConnectivityArrayEdge(newBinInds(4),9) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),10) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),12)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(1),11),:)=[this.ConnectivityArrayEdge(newBinInds(1),11)  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),10) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),12)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(3),6),:)=[this.ConnectivityArrayEdge(newBinInds(3),6) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),7) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),8)];
                        
                    else
                        this.ConnectivityArray(newBinInds(1),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),3);
                        this.ConnectivityArray(newBinInds(2),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),3);
                        this.ConnectivityArray(newBinInds(3),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),3);
                        this.ConnectivityArray(newBinInds(4),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),3);
                        %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),3))=0;
                        %                                 this.HangingNodesFace(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),3))=0;
                        
                        this.ConnectivityArrayEdge(newBinInds(1),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),8);
                        this.ConnectivityArrayEdge(newBinInds(2),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),12);
                        this.ConnectivityArrayEdge(newBinInds(3),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),10);
                        this.ConnectivityArrayEdge(newBinInds(4),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),7);
                        
                        this.ConnectivityArrayEdge(newBinInds(2),5)=this.ConnectivityArrayEdge(newBinInds(1),6);
                        this.ConnectivityArrayEdge(newBinInds(4),9)=this.ConnectivityArrayEdge(newBinInds(2),11);
                        this.ConnectivityArrayEdge(newBinInds(1),11)=this.ConnectivityArrayEdge(newBinInds(3),9);
                        this.ConnectivityArrayEdge(newBinInds(3),6)=this.ConnectivityArrayEdge(newBinInds(4),5);
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),8),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),12),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),10),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),7),:)=0;
                        
                        
                    end
                else
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(2),4)=this.PointCount;
                    this.ConnectivityArray(newBinInds(3),5)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),1)=this.PointCount;
                    
                    this.ConnectivityArrayEdge(newBinInds(1),6)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(2),11)=this.EdgeCount+2;
                    this.ConnectivityArrayEdge(newBinInds(3),9)=this.EdgeCount+3;
                    this.ConnectivityArrayEdge(newBinInds(4),5)=this.EdgeCount+4;
                    
                    this.ConnectivityArrayEdge(newBinInds(2),5)=this.ConnectivityArrayEdge(newBinInds(1),6);
                    this.ConnectivityArrayEdge(newBinInds(4),9)=this.ConnectivityArrayEdge(newBinInds(2),11);
                    this.ConnectivityArrayEdge(newBinInds(1),11)=this.ConnectivityArrayEdge(newBinInds(3),9);
                    this.ConnectivityArrayEdge(newBinInds(3),6)=this.ConnectivityArrayEdge(newBinInds(4),5);
                    this.EdgeCount=this.EdgeCount+4;
                end
                
            case 4 % Front
                if this.BinNeighbor(binNo,direction)
                    if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.PointCount=this.PointCount+1;
                        this.ConnectivityArray(newBinInds(5),7)=this.PointCount;
                        this.ConnectivityArray(newBinInds(6),3)=this.PointCount;
                        this.ConnectivityArray(newBinInds(7),6)=this.PointCount;
                        this.ConnectivityArray(newBinInds(8),2)=this.PointCount;
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,2) ; this.ConnectivityArray(binNo,6); this.ConnectivityArray(binNo,3); this.ConnectivityArray(binNo,7)  ];
                        %                                 this.HangingNodesFace(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,2) ; this.ConnectivityArray(binNo,6); this.ConnectivityArray(binNo,3); this.ConnectivityArray(binNo,7)  ];
                        
                        this.ConnectivityArrayEdge(newBinInds(5),12)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(6),7)=this.EdgeCount+2;
                        this.ConnectivityArrayEdge(newBinInds(7),8)=this.EdgeCount+3;
                        this.ConnectivityArrayEdge(newBinInds(8),10)=this.EdgeCount+4;
                        
                        this.ConnectivityArrayEdge(newBinInds(7),10)=this.ConnectivityArrayEdge(newBinInds(5),12);
                        this.ConnectivityArrayEdge(newBinInds(5),8)=this.ConnectivityArrayEdge(newBinInds(6),7);
                        this.ConnectivityArrayEdge(newBinInds(8),7)=this.ConnectivityArrayEdge(newBinInds(7),8);
                        this.ConnectivityArrayEdge(newBinInds(6),12)=this.ConnectivityArrayEdge(newBinInds(8),10);
                        this.EdgeCount=this.EdgeCount+4;
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(7),10),:)=[this.ConnectivityArrayEdge(newBinInds(7),10) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),9) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),11) ];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(5),8),:)=[this.ConnectivityArrayEdge(newBinInds(5),8) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),5) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),6)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(8),7),:)=[this.ConnectivityArrayEdge(newBinInds(8),7)  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),5) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),6)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(6),12),:)=[this.ConnectivityArrayEdge(newBinInds(6),12) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),9) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),11)];
                        
                        
                    else
                        this.ConnectivityArray(newBinInds(5),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),1);
                        this.ConnectivityArray(newBinInds(6),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),1);
                        this.ConnectivityArray(newBinInds(7),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),1);
                        this.ConnectivityArray(newBinInds(8),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),1);
                        %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),1))=0;
                        %                                 this.HangingNodesFace(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),1))=0;
                        
                        this.ConnectivityArrayEdge(newBinInds(5),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),11);
                        this.ConnectivityArrayEdge(newBinInds(6),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),5);
                        this.ConnectivityArrayEdge(newBinInds(7),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),6);
                        this.ConnectivityArrayEdge(newBinInds(8),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),9);
                        
                        this.ConnectivityArrayEdge(newBinInds(7),10)=this.ConnectivityArrayEdge(newBinInds(5),12);
                        this.ConnectivityArrayEdge(newBinInds(5),8)=this.ConnectivityArrayEdge(newBinInds(6),7);
                        this.ConnectivityArrayEdge(newBinInds(8),7)=this.ConnectivityArrayEdge(newBinInds(7),8);
                        this.ConnectivityArrayEdge(newBinInds(6),12)=this.ConnectivityArrayEdge(newBinInds(8),10);
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),11),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),5),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),6),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),9),:)=0;
                        
                    end
                else
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(5),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),2)=this.PointCount;
                    
                    this.ConnectivityArrayEdge(newBinInds(5),12)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(6),7)=this.EdgeCount+2;
                    this.ConnectivityArrayEdge(newBinInds(7),8)=this.EdgeCount+3;
                    this.ConnectivityArrayEdge(newBinInds(8),10)=this.EdgeCount+4;
                    
                    this.ConnectivityArrayEdge(newBinInds(7),10)=this.ConnectivityArrayEdge(newBinInds(5),12);
                    this.ConnectivityArrayEdge(newBinInds(5),8)=this.ConnectivityArrayEdge(newBinInds(6),7);
                    this.ConnectivityArrayEdge(newBinInds(8),7)=this.ConnectivityArrayEdge(newBinInds(7),8);
                    this.ConnectivityArrayEdge(newBinInds(6),12)=this.ConnectivityArrayEdge(newBinInds(8),10);
                    
                    this.EdgeCount=this.EdgeCount+4;
                end
                
            case 5 %Left
                if this.BinNeighbor(binNo,direction)
                    if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.PointCount=this.PointCount+1;
                        this.ConnectivityArray(newBinInds(3),7)=this.PointCount;
                        this.ConnectivityArray(newBinInds(4),3)=this.PointCount;
                        this.ConnectivityArray(newBinInds(7),8)=this.PointCount;
                        this.ConnectivityArray(newBinInds(8),4)=this.PointCount;
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,4); this.ConnectivityArray(binNo,7); this.ConnectivityArray(binNo,8)  ];
                        %                                 this.HangingNodesFace(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,4); this.ConnectivityArray(binNo,7); this.ConnectivityArray(binNo,8)  ];
                        
                        this.ConnectivityArrayEdge(newBinInds(3),12)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(4),2)=this.EdgeCount+2;
                        this.ConnectivityArrayEdge(newBinInds(7),4)=this.EdgeCount+3;
                        this.ConnectivityArrayEdge(newBinInds(8),11)=this.EdgeCount+4;
                        
                        this.ConnectivityArrayEdge(newBinInds(7),11)=this.ConnectivityArrayEdge(newBinInds(3),12);
                        this.ConnectivityArrayEdge(newBinInds(3),4)=this.ConnectivityArrayEdge(newBinInds(4),2);
                        this.ConnectivityArrayEdge(newBinInds(8),2)=this.ConnectivityArrayEdge(newBinInds(7),4);
                        this.ConnectivityArrayEdge(newBinInds(4),12)=this.ConnectivityArrayEdge(newBinInds(8),11);
                        this.EdgeCount=this.EdgeCount+4;
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(7),11),:)=[this.ConnectivityArrayEdge(newBinInds(7),11) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),9) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),10) ];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(3),4),:)=[this.ConnectivityArrayEdge(newBinInds(3),4) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),1) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),3)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(8),2),:)=[this.ConnectivityArrayEdge(newBinInds(8),2)  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),1) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),3)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(4),12),:)=[this.ConnectivityArrayEdge(newBinInds(4),12) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),9) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),10)];
                        
                        
                    else
                        this.ConnectivityArray(newBinInds(3),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),6);
                        this.ConnectivityArray(newBinInds(4),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),6);
                        this.ConnectivityArray(newBinInds(7),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),6);
                        this.ConnectivityArray(newBinInds(8),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),6);
                        %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),6))=0;
                        %                                 this.HangingNodesFace(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),6))=0;
                        
                        this.ConnectivityArrayEdge(newBinInds(3),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),10);
                        this.ConnectivityArrayEdge(newBinInds(4),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),1);
                        this.ConnectivityArrayEdge(newBinInds(7),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),3);
                        this.ConnectivityArrayEdge(newBinInds(8),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),9);
                        
                        this.ConnectivityArrayEdge(newBinInds(7),11)=this.ConnectivityArrayEdge(newBinInds(3),12);
                        this.ConnectivityArrayEdge(newBinInds(3),4)=this.ConnectivityArrayEdge(newBinInds(4),2);
                        this.ConnectivityArrayEdge(newBinInds(8),2)=this.ConnectivityArrayEdge(newBinInds(7),4);
                        this.ConnectivityArrayEdge(newBinInds(4),12)=this.ConnectivityArrayEdge(newBinInds(8),11);
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),10),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),1),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),3),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),9),:)=0;
                        
                    end
                else
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(3),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),4)=this.PointCount;
                    
                    this.ConnectivityArrayEdge(newBinInds(3),12)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(4),2)=this.EdgeCount+2;
                    this.ConnectivityArrayEdge(newBinInds(7),4)=this.EdgeCount+3;
                    this.ConnectivityArrayEdge(newBinInds(8),11)=this.EdgeCount+4;
                    
                    this.ConnectivityArrayEdge(newBinInds(7),11)=this.ConnectivityArrayEdge(newBinInds(3),12);
                    this.ConnectivityArrayEdge(newBinInds(3),4)=this.ConnectivityArrayEdge(newBinInds(4),2);
                    this.ConnectivityArrayEdge(newBinInds(8),2)=this.ConnectivityArrayEdge(newBinInds(7),4);
                    this.ConnectivityArrayEdge(newBinInds(4),12)=this.ConnectivityArrayEdge(newBinInds(8),11);
                    this.EdgeCount=this.EdgeCount+4;
                    
                end
                
            case 6 % Right
                if this.BinNeighbor(binNo,direction)
                    if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.PointCount=this.PointCount+1;
                        this.ConnectivityArray(newBinInds(1),6)=this.PointCount;
                        this.ConnectivityArray(newBinInds(2),2)=this.PointCount;
                        this.ConnectivityArray(newBinInds(5),5)=this.PointCount;
                        this.ConnectivityArray(newBinInds(6),1)=this.PointCount;
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,2) ; this.ConnectivityArray(binNo,6); this.ConnectivityArray(binNo,3); this.ConnectivityArray(binNo,7)  ];
                        %                                 this.HangingNodesFace(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1) ; this.ConnectivityArray(binNo,2); this.ConnectivityArray(binNo,5); this.ConnectivityArray(binNo,6)  ];
                        
                        this.ConnectivityArrayEdge(newBinInds(1),10)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(2),1)=this.EdgeCount+2;
                        this.ConnectivityArrayEdge(newBinInds(5),3)=this.EdgeCount+3;
                        this.ConnectivityArrayEdge(newBinInds(6),9)=this.EdgeCount+4;
                        
                        this.ConnectivityArrayEdge(newBinInds(5),9)=this.ConnectivityArrayEdge(newBinInds(1),10);
                        this.ConnectivityArrayEdge(newBinInds(1),3)=this.ConnectivityArrayEdge(newBinInds(2),1);
                        this.ConnectivityArrayEdge(newBinInds(6),1)=this.ConnectivityArrayEdge(newBinInds(5),3);
                        this.ConnectivityArrayEdge(newBinInds(2),10)=this.ConnectivityArrayEdge(newBinInds(6),9);
                        this.EdgeCount=this.EdgeCount+4;
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(5),9),:)=[this.ConnectivityArrayEdge(newBinInds(5),9) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),11) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),12) ];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(1),3),:)=[this.ConnectivityArrayEdge(newBinInds(1),3) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),2) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),4)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(6),1),:)=[this.ConnectivityArrayEdge(newBinInds(6),1)  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),2) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),4)];
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(newBinInds(2),10),:)=[this.ConnectivityArrayEdge(newBinInds(2),10) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),11) this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),12)];
                        
                    else
                        this.ConnectivityArray(newBinInds(1),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),3);
                        this.ConnectivityArray(newBinInds(2),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),3);
                        this.ConnectivityArray(newBinInds(5),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),3);
                        this.ConnectivityArray(newBinInds(6),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),3);
                        %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),3))=0;
                        %                                 this.HangingNodesFace(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),3))=0;
                        
                        this.ConnectivityArrayEdge(newBinInds(1),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),12);
                        this.ConnectivityArrayEdge(newBinInds(2),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),2);
                        this.ConnectivityArrayEdge(newBinInds(5),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),4);
                        this.ConnectivityArrayEdge(newBinInds(6),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),11);
                        
                        this.ConnectivityArrayEdge(newBinInds(5),9)=this.ConnectivityArrayEdge(newBinInds(1),10);
                        this.ConnectivityArrayEdge(newBinInds(1),3)=this.ConnectivityArrayEdge(newBinInds(2),1);
                        this.ConnectivityArrayEdge(newBinInds(6),1)=this.ConnectivityArrayEdge(newBinInds(5),3);
                        this.ConnectivityArrayEdge(newBinInds(2),10)=this.ConnectivityArrayEdge(newBinInds(6),9);
                        
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),12),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),2),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),4),:)=0;
                        this.HangingEdgeFace(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),11),:)=0;
                        
                    end
                else
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(2),2)=this.PointCount;
                    this.ConnectivityArray(newBinInds(5),5)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),1)=this.PointCount;
                    
                    this.ConnectivityArrayEdge(newBinInds(1),10)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(2),1)=this.EdgeCount+2;
                    this.ConnectivityArrayEdge(newBinInds(5),3)=this.EdgeCount+3;
                    this.ConnectivityArrayEdge(newBinInds(6),9)=this.EdgeCount+4;
                    
                    this.ConnectivityArrayEdge(newBinInds(5),9)=this.ConnectivityArrayEdge(newBinInds(1),10);
                    this.ConnectivityArrayEdge(newBinInds(1),3)=this.ConnectivityArrayEdge(newBinInds(2),1);
                    this.ConnectivityArrayEdge(newBinInds(6),1)=this.ConnectivityArrayEdge(newBinInds(5),3);
                    this.ConnectivityArrayEdge(newBinInds(2),10)=this.ConnectivityArrayEdge(newBinInds(6),9);
                    this.EdgeCount=this.EdgeCount+4;
                end
                
            case 7% edge 1-2
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(1),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),6),5);
                            this.ConnectivityArray(newBinInds(5),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),6),5);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),6),5))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),2),3);
                            this.ConnectivityArrayEdge(newBinInds(5),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),6),3);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=0;
                            
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(1),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),3);
                            this.ConnectivityArray(newBinInds(5),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),3);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),3))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),3),2);
                            this.ConnectivityArrayEdge(newBinInds(5),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),7),2);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=0;
                            
                            continue;
                        end
                    end
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),2)=this.PointCount;
                    this.ConnectivityArray(newBinInds(5),1)=this.PointCount;
                    
                    if ~this.BinNeighbor(binNo,2) && ~this.BinNeighbor(binNo,6)
                        this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(binNo,1);
                        this.ConnectivityArrayEdge(newBinInds(5),1)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        this.ConnectivityArrayEdge(newBinInds(1),1)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(5),1)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=[this.ConnectivityArrayEdge(newBinInds(1),1) this.ConnectivityArrayEdge(binNo,1)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=[this.ConnectivityArrayEdge(newBinInds(5),1) this.ConnectivityArrayEdge(binNo,1) ];
                        
                        
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1); this.ConnectivityArray(binNo,2); 0 ; 0 ];
                        continue;
                    end
                else
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(1),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),7);
                        this.ConnectivityArray(newBinInds(5),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),7);
                        
                        this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),4);
                        this.ConnectivityArrayEdge(newBinInds(5),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),4);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            %                                 this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),4),7))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=0;
                            
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(1),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),6),5);
                            this.ConnectivityArray(newBinInds(5),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),6),5);
                            
                            this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),2),3);
                            this.ConnectivityArrayEdge(newBinInds(5),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),6),3);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                     this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),6),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(1),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),3);
                            this.ConnectivityArray(newBinInds(5),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),3);
                            this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),3),2);
                            this.ConnectivityArrayEdge(newBinInds(5),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),7),2);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                     this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),3))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=0;
                                
                                
                            end
                            continue;
                        end
                    end
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),2)=this.PointCount;
                    this.ConnectivityArray(newBinInds(5),1)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1); this.ConnectivityArray(binNo,2); 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(1),1)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(5),1)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),1))=[this.ConnectivityArrayEdge(newBinInds(1),1) this.ConnectivityArrayEdge(binNo,1)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),1))=[this.ConnectivityArrayEdge(newBinInds(5),1) this.ConnectivityArrayEdge(binNo,1) ];
                    
                    
                    continue;
                    
                    
                end
                
                
                
            case 8 % edge 2-3
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(5),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),6);
                            this.ConnectivityArray(newBinInds(7),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),6);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),2))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(5),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),6),8);
                            this.ConnectivityArrayEdge(newBinInds(7),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),8),8);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=0;
                            
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(5),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),4);
                            this.ConnectivityArray(newBinInds(7),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),4);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),2))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(5),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),1),5);
                            this.ConnectivityArrayEdge(newBinInds(7),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),3),5);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=0;
                            
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(5),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),2)=this.PointCount;
                    if ~this.BinNeighbor(binNo,2) && ~this.BinNeighbor(binNo,4)
                        this.ConnectivityArrayEdge(newBinInds(5),7)=this.ConnectivityArrayEdge(binNo,7);
                        this.ConnectivityArrayEdge(newBinInds(7),7)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        this.ConnectivityArrayEdge(newBinInds(5),7)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(7),7)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=[this.ConnectivityArrayEdge(newBinInds(5),7) this.ConnectivityArrayEdge(binNo,7)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=[this.ConnectivityArrayEdge(newBinInds(7),7) this.ConnectivityArrayEdge(binNo,7) ];
                        
                        
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,2); this.ConnectivityArray(binNo,3) ; 0 ; 0 ];
                        continue;
                    end
                else
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(5),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),8);
                        this.ConnectivityArray(newBinInds(7),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),8);
                        
                        this.ConnectivityArrayEdge(newBinInds(5),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),6);
                        this.ConnectivityArrayEdge(newBinInds(7),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),6);
                        
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            %                                 this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),2))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=0;
                            
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(5),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),6);
                            this.ConnectivityArray(newBinInds(7),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),6);
                            
                            this.ConnectivityArrayEdge(newBinInds(5),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),6),8);
                            this.ConnectivityArrayEdge(newBinInds(7),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),8),8);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,4)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),2))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(5),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),4);
                            this.ConnectivityArray(newBinInds(7),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),4);
                            
                            this.ConnectivityArrayEdge(newBinInds(5),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),1),5);
                            this.ConnectivityArrayEdge(newBinInds(7),7)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),3),5);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),2))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(5),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),2)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,2); this.ConnectivityArray(binNo,3) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(5),7)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(7),7)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),7))=[this.ConnectivityArrayEdge(newBinInds(5),7) this.ConnectivityArrayEdge(binNo,7)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),7))=[this.ConnectivityArrayEdge(newBinInds(7),7) this.ConnectivityArrayEdge(binNo,7) ];
                    
                    
                    continue;
                end
                
                
                
                
                
                
            case 9 % edge 3-4
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(3),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),8);
                            this.ConnectivityArray(newBinInds(7),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),8);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),4))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),4),4);
                            this.ConnectivityArrayEdge(newBinInds(7),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),8),4);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(3),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),2);
                            this.ConnectivityArray(newBinInds(7),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),2);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),4))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),1),1);
                            this.ConnectivityArrayEdge(newBinInds(7),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),5),1);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=0;
                            
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(3),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),4)=this.PointCount;
                    if ~this.BinNeighbor(binNo,2) && ~this.BinNeighbor(binNo,5)
                        this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(binNo,2);
                        this.ConnectivityArrayEdge(newBinInds(7),2)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,4) ; 0 ; 0 ];
                        
                        this.ConnectivityArrayEdge(newBinInds(3),2)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(7),2)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=[this.ConnectivityArrayEdge(newBinInds(3),2) this.ConnectivityArrayEdge(binNo,2)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=[this.ConnectivityArrayEdge(newBinInds(7),2) this.ConnectivityArrayEdge(binNo,2) ];
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(3),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),6);
                        this.ConnectivityArray(newBinInds(7),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),6);
                        
                        this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),3);
                        this.ConnectivityArrayEdge(newBinInds(7),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),3);
                        
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),4))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(3),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),8);
                            this.ConnectivityArray(newBinInds(7),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),8),8);
                            
                            this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),4),4);
                            this.ConnectivityArrayEdge(newBinInds(7),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),8),4);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,5)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),4))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(3),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),2);
                            this.ConnectivityArray(newBinInds(7),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),2);
                            
                            this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),1),1);
                            this.ConnectivityArrayEdge(newBinInds(7),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),5),1);
                            
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(7),4))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(3),3)=this.PointCount;
                    this.ConnectivityArray(newBinInds(7),4)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,4) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(3),2)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(7),2)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),2))=[this.ConnectivityArrayEdge(newBinInds(3),2) this.ConnectivityArrayEdge(binNo,2)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),2))=[this.ConnectivityArrayEdge(newBinInds(7),2) this.ConnectivityArrayEdge(binNo,2) ];
                    
                    continue;
                end
                
                
                
            case 10 % edge 4-1
                
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(1),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),2),8);
                            this.ConnectivityArray(newBinInds(3),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),2),8);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(3),1))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),2),6);
                            this.ConnectivityArrayEdge(newBinInds(3),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),4),6);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=0;
                            
                            
                            continue;
                        end
                    end
                    
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(1),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),3);
                            this.ConnectivityArray(newBinInds(3),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),3);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(3),1))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),5),7);
                            this.ConnectivityArrayEdge(newBinInds(3),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),7),7);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=0;
                            
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),4)=this.PointCount;
                    this.ConnectivityArray(newBinInds(3),1)=this.PointCount;
                    if ~this.BinNeighbor(binNo,2) && ~this.BinNeighbor(binNo,3)
                        this.ConnectivityArrayEdge(newBinInds(1),5)=this.ConnectivityArrayEdge(binNo,5);
                        this.ConnectivityArrayEdge(newBinInds(3),5)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,4) ; this.ConnectivityArray(binNo,1) ; 0 ; 0 ];
                        this.ConnectivityArrayEdge(newBinInds(1),5)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(3),5)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=[this.ConnectivityArrayEdge(newBinInds(1),5) this.ConnectivityArrayEdge(binNo,5)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=[this.ConnectivityArrayEdge(newBinInds(3),5) this.ConnectivityArrayEdge(binNo,5) ];
                        
                        continue;
                    end
                else
                    
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(1),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),7);
                        this.ConnectivityArray(newBinInds(3),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),6),7);
                        
                        this.ConnectivityArrayEdge(newBinInds(1),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),8);
                        this.ConnectivityArrayEdge(newBinInds(3),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),8);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(3),1))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,2)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,2))
                            this.ConnectivityArray(newBinInds(1),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),2),8);
                            this.ConnectivityArray(newBinInds(3),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,2),2),8);
                            
                            this.ConnectivityArrayEdge(newBinInds(1),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),2),6);
                            this.ConnectivityArrayEdge(newBinInds(3),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,2),4),6);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,3)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(3),1))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(1),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),3);
                            this.ConnectivityArray(newBinInds(3),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),3);
                            
                            this.ConnectivityArrayEdge(newBinInds(1),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),5),7);
                            this.ConnectivityArrayEdge(newBinInds(3),5)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),7),7);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,2)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(3),1))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),4)=this.PointCount;
                    this.ConnectivityArray(newBinInds(3),1)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,4) ; this.ConnectivityArray(binNo,1) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(1),5)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(3),5)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),5))=[this.ConnectivityArrayEdge(newBinInds(1),5) this.ConnectivityArrayEdge(binNo,5)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),5))=[this.ConnectivityArrayEdge(newBinInds(3),5) this.ConnectivityArrayEdge(binNo,5) ];
                    continue;
                end
                
                
                
                
                
                
            case 11 % edge 1-5
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(1),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),6);
                            this.ConnectivityArray(newBinInds(2),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),6);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(2),1))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),5),10);
                            this.ConnectivityArrayEdge(newBinInds(2),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),6),10);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(1),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),8);
                            this.ConnectivityArray(newBinInds(2),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),8);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(2),1))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),3),11);
                            this.ConnectivityArrayEdge(newBinInds(2),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),4),11);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=0;
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),5)=this.PointCount;
                    this.ConnectivityArray(newBinInds(2),1)=this.PointCount;
                    if ~this.BinNeighbor(binNo,3) && ~this.BinNeighbor(binNo,6)
                        this.ConnectivityArrayEdge(newBinInds(1),9)=this.ConnectivityArrayEdge(binNo,9);
                        this.ConnectivityArrayEdge(newBinInds(2),9)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,5) ; this.ConnectivityArray(binNo,1) ; 0 ; 0 ];
                        
                        this.ConnectivityArrayEdge(newBinInds(1),9)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(2),9)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=[this.ConnectivityArrayEdge(newBinInds(1),9) this.ConnectivityArrayEdge(binNo,9)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=[this.ConnectivityArrayEdge(newBinInds(2),9) this.ConnectivityArrayEdge(binNo,9) ];
                        
                        continue;
                    end
                else% this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(1),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),8),3);
                        this.ConnectivityArray(newBinInds(2),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),8),3);
                        
                        this.ConnectivityArrayEdge(newBinInds(1),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),12);
                        this.ConnectivityArrayEdge(newBinInds(2),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),8),12);
                        
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(2),1))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(1),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),6);
                            this.ConnectivityArray(newBinInds(2),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),5),6);
                            
                            this.ConnectivityArrayEdge(newBinInds(1),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),5),10);
                            this.ConnectivityArrayEdge(newBinInds(2),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),6),10);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(2),1))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(1),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),8);
                            this.ConnectivityArray(newBinInds(2),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),3),8);
                            
                            this.ConnectivityArrayEdge(newBinInds(1),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),3),11);
                            this.ConnectivityArrayEdge(newBinInds(2),9)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),4),11);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,3)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(2),1))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=0;
                                
                            end
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(1),5)=this.PointCount;
                    this.ConnectivityArray(newBinInds(2),1)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,5) ; this.ConnectivityArray(binNo,1) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(1),9)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(2),9)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(1),9))=[this.ConnectivityArrayEdge(newBinInds(1),9) this.ConnectivityArrayEdge(binNo,9)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),9))=[this.ConnectivityArrayEdge(newBinInds(2),9) this.ConnectivityArrayEdge(binNo,9) ];
                    
                    continue;
                end
                
                
                
                
                
            case 12 % edge 2-6
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(5),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),5);
                            this.ConnectivityArray(newBinInds(6),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),5);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),2))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(5),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),1),9);
                            this.ConnectivityArrayEdge(newBinInds(6),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),2),9);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(5),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),7),7);
                            this.ConnectivityArray(newBinInds(6),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),7),7);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),2))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(5),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),7),12);
                            this.ConnectivityArrayEdge(newBinInds(6),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),8),12);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=0;
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(5),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),2)=this.PointCount;
                    if ~this.BinNeighbor(binNo,4) && ~this.BinNeighbor(binNo,6)
                        this.ConnectivityArrayEdge(newBinInds(5),10)=this.ConnectivityArrayEdge(binNo,10);
                        this.ConnectivityArrayEdge(newBinInds(6),10)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        this.ConnectivityArrayEdge(newBinInds(5),10)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(6),10)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=[this.ConnectivityArrayEdge(newBinInds(5),10) this.ConnectivityArrayEdge(binNo,10)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=[this.ConnectivityArrayEdge(newBinInds(6),10) this.ConnectivityArrayEdge(binNo,10) ];
                        
                        
                        %                                 this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,2) ; this.ConnectivityArray(binNo,6) ; 0 ; 0 ];
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(5),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),8);
                        this.ConnectivityArray(newBinInds(6),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),8);
                        
                        this.ConnectivityArrayEdge(newBinInds(5),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),11);
                        this.ConnectivityArrayEdge(newBinInds(6),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),11);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),2))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(5),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),5);
                            this.ConnectivityArray(newBinInds(6),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),1),5);
                            
                            this.ConnectivityArrayEdge(newBinInds(5),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),1),9);
                            this.ConnectivityArrayEdge(newBinInds(6),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),2),9);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),2))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=0;
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(5),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),7),7);
                            this.ConnectivityArray(newBinInds(6),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),7),7);
                            
                            this.ConnectivityArrayEdge(newBinInds(5),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),7),12);
                            this.ConnectivityArrayEdge(newBinInds(6),10)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),8),12);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,4)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),2))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=0;
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(5),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),2)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,2) ; this.ConnectivityArray(binNo,6) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(5),10)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(6),10)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(5),10))=[this.ConnectivityArrayEdge(newBinInds(5),10) this.ConnectivityArrayEdge(binNo,10)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),10))=[this.ConnectivityArrayEdge(newBinInds(6),10) this.ConnectivityArrayEdge(binNo,10) ];
                    continue;
                end
                
                
                
                
                
                
            case 13 % edge 3-7
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(7),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),4),4);
                            this.ConnectivityArray(newBinInds(8),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),4),4);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),3))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(7),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),3),11);
                            this.ConnectivityArrayEdge(newBinInds(8),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),4),11);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(7),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),5),6);
                            this.ConnectivityArray(newBinInds(8),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),5),6);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),3))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(7),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),5),10);
                            this.ConnectivityArrayEdge(newBinInds(8),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),6),10);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=0;
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(7),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),3)=this.PointCount;
                    if ~this.BinNeighbor(binNo,4) && ~this.BinNeighbor(binNo,5)
                        this.ConnectivityArrayEdge(newBinInds(7),12)=this.ConnectivityArrayEdge(binNo,12);
                        this.ConnectivityArrayEdge(newBinInds(8),12)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        this.ConnectivityArrayEdge(newBinInds(7),12)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(8),12)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=[this.ConnectivityArrayEdge(newBinInds(7),12) this.ConnectivityArrayEdge(binNo,12)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=[this.ConnectivityArrayEdge(newBinInds(8),12) this.ConnectivityArrayEdge(binNo,12) ];
                        
                        %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,7) ; 0 ; 0 ];
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(7),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),5);
                        this.ConnectivityArray(newBinInds(8),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),5);
                        
                        this.ConnectivityArrayEdge(newBinInds(7),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),9);
                        this.ConnectivityArrayEdge(newBinInds(8),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),9);
                        
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,4)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),3))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(7),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),4),4);
                            this.ConnectivityArray(newBinInds(8),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),4),4);
                            
                            this.ConnectivityArrayEdge(newBinInds(7),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),3),11);
                            this.ConnectivityArrayEdge(newBinInds(8),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),4),11);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,5)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),3))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=0;
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(7),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),5),6);
                            this.ConnectivityArray(newBinInds(8),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),5),6);
                            
                            this.ConnectivityArrayEdge(newBinInds(7),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),5),10);
                            this.ConnectivityArrayEdge(newBinInds(8),12)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),6),10);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,4)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),3))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=0;
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(7),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),3)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,7) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(7),12)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(8),12)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(7),12))=[this.ConnectivityArrayEdge(newBinInds(7),12) this.ConnectivityArrayEdge(binNo,12)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),12))=[this.ConnectivityArrayEdge(newBinInds(8),12) this.ConnectivityArrayEdge(binNo,12) ];
                    
                    continue;
                end
                
            case 14
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(3),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),8),3);
                            this.ConnectivityArray(newBinInds(4),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),8),3);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),4))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(3),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),7),12);
                            this.ConnectivityArrayEdge(newBinInds(4),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),8),12);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(3),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),5);
                            this.ConnectivityArray(newBinInds(4),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),5);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),4))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(3),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),1),9);
                            this.ConnectivityArrayEdge(newBinInds(4),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),2),9);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=0;
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(3),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),4)=this.PointCount;
                    if ~this.BinNeighbor(binNo,3) && ~this.BinNeighbor(binNo,5)
                        this.ConnectivityArrayEdge(newBinInds(3),11)=this.ConnectivityArrayEdge(binNo,11);
                        this.ConnectivityArrayEdge(newBinInds(4),11)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        this.ConnectivityArrayEdge(newBinInds(3),11)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(4),11)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=[this.ConnectivityArrayEdge(newBinInds(3),11) this.ConnectivityArrayEdge(binNo,11)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=[this.ConnectivityArrayEdge(newBinInds(4),11) this.ConnectivityArrayEdge(binNo,11) ];
                        
                        %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,4) ; this.ConnectivityArray(binNo,8) ; 0 ; 0 ];
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(3),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),5),6);
                        this.ConnectivityArray(newBinInds(4),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),5),6);
                        
                        this.ConnectivityArrayEdge(newBinInds(3),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),10);
                        this.ConnectivityArrayEdge(newBinInds(4),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),6),10);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,3)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),4))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(3),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),8),3);
                            this.ConnectivityArray(newBinInds(4),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),8),3);
                            
                            this.ConnectivityArrayEdge(newBinInds(3),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),7),12);
                            this.ConnectivityArrayEdge(newBinInds(4),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),8),12);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,5)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),4))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(3),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),5);
                            this.ConnectivityArray(newBinInds(4),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),1),5);
                            
                            this.ConnectivityArrayEdge(newBinInds(3),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),1),9);
                            this.ConnectivityArrayEdge(newBinInds(4),11)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),2),9);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,3)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),4))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(3),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),4)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,4) ; this.ConnectivityArray(binNo,8) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(3),11)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(4),11)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(3),11))=[this.ConnectivityArrayEdge(newBinInds(3),11) this.ConnectivityArrayEdge(binNo,11)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),11))=[this.ConnectivityArrayEdge(newBinInds(4),11) this.ConnectivityArrayEdge(binNo,11) ];
                    continue;
                end
                
                
                
                
                
                
                
                
            case 15
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(2),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),2);
                            this.ConnectivityArray(newBinInds(6),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),2);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),5))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(2),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),1),1);
                            this.ConnectivityArrayEdge(newBinInds(6),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),5),1);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=0;
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(2),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),4),7);
                            this.ConnectivityArray(newBinInds(6),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),4),7);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),5))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(2),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),4),4);
                            this.ConnectivityArrayEdge(newBinInds(6),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),8),4);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=0;
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(2),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),5)=this.PointCount;
                    if ~this.BinNeighbor(binNo,1) && ~this.BinNeighbor(binNo,6)
                        this.ConnectivityArrayEdge(newBinInds(2),3)=this.ConnectivityArrayEdge(binNo,3);
                        this.ConnectivityArrayEdge(newBinInds(6),3)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,5) ; this.ConnectivityArray(binNo,6) ; 0 ; 0 ];
                        
                        this.ConnectivityArrayEdge(newBinInds(2),3)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(6),3)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=[this.ConnectivityArrayEdge(newBinInds(2),3) this.ConnectivityArrayEdge(binNo,3)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=[this.ConnectivityArrayEdge(newBinInds(6),3) this.ConnectivityArrayEdge(binNo,3) ];
                        
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(2),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),3);
                        this.ConnectivityArray(newBinInds(6),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),3);
                        
                        this.ConnectivityArrayEdge(newBinInds(2),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),2);
                        this.ConnectivityArrayEdge(newBinInds(6),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),2);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),5))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(2),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),2);
                            this.ConnectivityArray(newBinInds(6),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),2);
                            
                            this.ConnectivityArrayEdge(newBinInds(2),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),1),1);
                            this.ConnectivityArrayEdge(newBinInds(6),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),5),1);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,6)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,6)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,6))
                            this.ConnectivityArray(newBinInds(2),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),4),7);
                            this.ConnectivityArray(newBinInds(6),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,6),4),7);
                            
                            this.ConnectivityArrayEdge(newBinInds(2),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),4),4);
                            this.ConnectivityArrayEdge(newBinInds(6),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,6),8),4);
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(6),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(2),6)=this.PointCount;
                    this.ConnectivityArray(newBinInds(6),5)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,5) ; this.ConnectivityArray(binNo,6) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(2),3)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(6),3)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),3))=[this.ConnectivityArrayEdge(newBinInds(2),3) this.ConnectivityArrayEdge(binNo,3)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),3))=[this.ConnectivityArrayEdge(newBinInds(6),3) this.ConnectivityArrayEdge(binNo,3) ];
                    
                    continue;
                end
                
                
                
                
                
            case 16
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(6),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),5),3);
                            this.ConnectivityArray(newBinInds(8),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),5),3);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),6))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(6),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),5),7);
                            this.ConnectivityArrayEdge(newBinInds(8),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),7),7);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(6),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),2),8);
                            this.ConnectivityArray(newBinInds(8),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),2),8);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),6))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(6),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),2),6);
                            this.ConnectivityArrayEdge(newBinInds(8),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),4),6);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=0;
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(6),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),6)=this.PointCount;
                    if ~this.BinNeighbor(binNo,1) && ~this.BinNeighbor(binNo,4)
                        this.ConnectivityArrayEdge(newBinInds(6),8)=this.ConnectivityArrayEdge(binNo,8);
                        this.ConnectivityArrayEdge(newBinInds(8),8)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,6) ; this.ConnectivityArray(binNo,7) ; 0 ; 0 ];
                        
                        this.ConnectivityArrayEdge(newBinInds(6),8)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(8),8)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=[this.ConnectivityArrayEdge(newBinInds(6),8) this.ConnectivityArrayEdge(binNo,8)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=[this.ConnectivityArrayEdge(newBinInds(8),8) this.ConnectivityArrayEdge(binNo,8) ];
                        
                        continue;
                    end
                else %this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(6),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),4);
                        this.ConnectivityArray(newBinInds(8),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),4);
                        
                        this.ConnectivityArrayEdge(newBinInds(6),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),5);
                        this.ConnectivityArrayEdge(newBinInds(8),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),5);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),6))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=0;
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(6),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),5),3);
                            this.ConnectivityArray(newBinInds(8),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),5),3);
                            
                            this.ConnectivityArrayEdge(newBinInds(6),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),5),7);
                            this.ConnectivityArrayEdge(newBinInds(8),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),7),7);
                            
                            
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,4)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),6))=0;
                                
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,4)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,4))
                            this.ConnectivityArray(newBinInds(6),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),2),8);
                            this.ConnectivityArray(newBinInds(8),6)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,4),2),8);
                            
                            this.ConnectivityArrayEdge(newBinInds(6),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),2),6);
                            this.ConnectivityArrayEdge(newBinInds(8),8)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,4),4),6);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),6))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(6),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),6)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,6) ; this.ConnectivityArray(binNo,7) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(6),8)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(8),8)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(6),8))=[this.ConnectivityArrayEdge(newBinInds(6),8) this.ConnectivityArrayEdge(binNo,8)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),8))=[this.ConnectivityArrayEdge(newBinInds(8),8) this.ConnectivityArrayEdge(binNo,8) ];
                    
                    continue;
                end
                
                
                
                
            case 17
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(4),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),3),3);
                            this.ConnectivityArray(newBinInds(8),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),3),3);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),8))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),3),2);
                            this.ConnectivityArrayEdge(newBinInds(8),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),7),2);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(4),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),2),6);
                            this.ConnectivityArray(newBinInds(8),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),2),6);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),8))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),2),3);
                            this.ConnectivityArrayEdge(newBinInds(8),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),6),3);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=0;
                            
                            continue;
                        end
                    end
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(4),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),8)=this.PointCount;
                    if ~this.BinNeighbor(binNo,1) && ~this.BinNeighbor(binNo,5)
                        this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(binNo,4);
                        this.ConnectivityArrayEdge(newBinInds(8),4)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,7) ; this.ConnectivityArray(binNo,8) ; 0 ; 0 ];
                        
                        this.ConnectivityArrayEdge(newBinInds(4),4)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(8),4)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=[this.ConnectivityArrayEdge(newBinInds(4),4) this.ConnectivityArrayEdge(binNo,4)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=[this.ConnectivityArrayEdge(newBinInds(8),4) this.ConnectivityArrayEdge(binNo,4) ];
                        
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(4),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),2);
                        this.ConnectivityArray(newBinInds(8),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),2);
                        
                        this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),1);
                        this.ConnectivityArrayEdge(newBinInds(8),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),1);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),8))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=0;
                            
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(4),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),3),3);
                            this.ConnectivityArray(newBinInds(8),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),3),3);
                            
                            this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),3),2);
                            this.ConnectivityArrayEdge(newBinInds(8),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),7),2);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,5)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),8))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,5)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,5))
                            this.ConnectivityArray(newBinInds(4),7)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),2),6);
                            this.ConnectivityArray(newBinInds(8),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,5),2),6);
                            
                            this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),2),3);
                            this.ConnectivityArrayEdge(newBinInds(8),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,5),6),3);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(8),8))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=0;
                                
                            end
                            continue;
                        end
                    end
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(4),7)=this.PointCount;
                    this.ConnectivityArray(newBinInds(8),8)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,7) ; this.ConnectivityArray(binNo,8) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(4),4)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(8),4)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),4))=[this.ConnectivityArrayEdge(newBinInds(4),4) this.ConnectivityArrayEdge(binNo,4)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(8),4))=[this.ConnectivityArrayEdge(newBinInds(8),4) this.ConnectivityArrayEdge(binNo,4) ];
                    
                    continue;
                end
                
                
                
                
            case 18
                
                
                if this.BinNeighbor(binNo,direction)==-1 || ~this.BinNeighbor(binNo,direction)
                    
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(2),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),4);
                            this.ConnectivityArray(newBinInds(4),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),4);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),5))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(2),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),1),5);
                            this.ConnectivityArrayEdge(newBinInds(4),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),3),5);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=0;
                            
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(2),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),6),7);
                            this.ConnectivityArray(newBinInds(4),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),6),7);
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),5))=0;
                            
                            this.ConnectivityArrayEdge(newBinInds(2),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),6),8);
                            this.ConnectivityArrayEdge(newBinInds(4),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),8),8);
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=0;
                            
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(2),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),5)=this.PointCount;
                    if ~this.BinNeighbor(binNo,1) && ~this.BinNeighbor(binNo,3)
                        this.ConnectivityArrayEdge(newBinInds(2),6)=this.ConnectivityArrayEdge(binNo,6);
                        this.ConnectivityArrayEdge(newBinInds(4),6)=this.EdgeCount+1;
                        this.EdgeCount=this.EdgeCount+1;
                        continue;
                    else
                        %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,8) ; this.ConnectivityArray(binNo,5) ; 0 ; 0 ];
                        
                        this.ConnectivityArrayEdge(newBinInds(2),6)=this.EdgeCount+1;
                        this.ConnectivityArrayEdge(newBinInds(4),6)=this.EdgeCount+2;
                        this.EdgeCount=this.EdgeCount+2;
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=[this.ConnectivityArrayEdge(newBinInds(2),6) this.ConnectivityArrayEdge(binNo,6)   ];
                        this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=[this.ConnectivityArrayEdge(newBinInds(4),6) this.ConnectivityArrayEdge(binNo,6) ];
                        
                        continue;
                    end
                else % this.BinNeighbor(binNo,direction)
                    if ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                        this.ConnectivityArray(newBinInds(2),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),5),3);
                        this.ConnectivityArray(newBinInds(4),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),5),3);
                        
                        this.ConnectivityArrayEdge(newBinInds(2),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),5),7);
                        this.ConnectivityArrayEdge(newBinInds(4),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),7),7);
                        
                        if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            %                                     this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),5))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=0;
                            this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=0;
                            
                        end
                        continue;
                    end %this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    if this.BinNeighbor(binNo,1)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,1))
                            this.ConnectivityArray(newBinInds(2),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),4);
                            this.ConnectivityArray(newBinInds(4),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,1),1),4);
                            
                            this.ConnectivityArrayEdge(newBinInds(2),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),1),5);
                            this.ConnectivityArrayEdge(newBinInds(4),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,1),3),5);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,3)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=0;
                                
                            end
                            continue;
                        end
                    end
                    if this.BinNeighbor(binNo,3)
                        if ~this.BinIsLeaf(this.BinNeighbor(binNo,3))
                            this.ConnectivityArray(newBinInds(2),8)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),6),7);
                            this.ConnectivityArray(newBinInds(4),5)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,3),6),7);
                            
                            this.ConnectivityArrayEdge(newBinInds(2),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),6),8);
                            this.ConnectivityArrayEdge(newBinInds(4),6)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,3),8),8);
                            
                            if  ~ this.BinIsLeaf(this.BinNeighbor(binNo,1)) && ~ this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                %                                         this.HangingNodes(:,this.ConnectivityArray(newBinInds(4),5))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=0;
                                this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=0;
                                
                            end
                            continue;
                        end
                    end
                    
                    this.PointCount=this.PointCount+1;
                    this.ConnectivityArray(newBinInds(2),8)=this.PointCount;
                    this.ConnectivityArray(newBinInds(4),5)=this.PointCount;
                    %                             this.HangingNodes(:,this.PointCount)=[this.PointCount ; this.ConnectivityArray(binNo,8) ; this.ConnectivityArray(binNo,5) ; 0 ; 0 ];
                    
                    this.ConnectivityArrayEdge(newBinInds(2),6)=this.EdgeCount+1;
                    this.ConnectivityArrayEdge(newBinInds(4),6)=this.EdgeCount+2;
                    this.EdgeCount=this.EdgeCount+2;
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(2),6))=[this.ConnectivityArrayEdge(newBinInds(2),6) this.ConnectivityArrayEdge(binNo,6)   ];
                    this.HangingEdge(:,this.ConnectivityArrayEdge(newBinInds(4),6))=[this.ConnectivityArrayEdge(newBinInds(4),6) this.ConnectivityArrayEdge(binNo,6) ];
                    
                    continue;
                end
                
                
                
        end
        
    end
    end
    function update_neighbor(this,binNo,radius)
    
    
    
    for direction=1:18 % [up down back front left right] , [edge directions]
        if   (this.BinNeighbor(binNo,direction))
            Level_diff=this.BinDepths(binNo)-this.BinDepths(this.BinNeighbor(binNo,direction));
            
            
            if Level_diff<1
                
                if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                    updateNeighborMatrixMinus(this,binNo,direction);
                else
                    updateNeighborMatrixMinusUpdate(this,binNo,direction);
                end                        %             this.BinNeighbor(this.BinChildren(binNo,:),direction)=Level1_neighbor(:,direction);
                
                
            else  % level diff >=2
                
                this.divideBin(this.BinNeighbor(binNo,direction),radius)
                
                updateNeighborMatrixMinus(this,binNo,direction);
                
                
                
                
                
            end
            
            
            
        else
            updateNeighborMatrixRegular(this,binNo,direction);
            
        end
        
    end
    end
    
    function updateNeighborMatrixMinusUpdate(this,binNo,direction)
    if direction==1
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,2) this.BinChildren(this.BinNeighbor(binNo,1),1) this.BinChildren(binNo,4) this.BinChildren(this.BinNeighbor(binNo,1),3) this.BinChildren(binNo,6) this.BinChildren(this.BinNeighbor(binNo,1),5) this.BinChildren(binNo,8) this.BinChildren(this.BinNeighbor(binNo,1),7)  ]';
        
        
        this.BinNeighbor(this.BinChildren(binNo,2),16)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        this.BinNeighbor(this.BinChildren(binNo,2),17)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        this.BinNeighbor(this.BinChildren(binNo,4),15)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        this.BinNeighbor(this.BinChildren(binNo,4),16)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        this.BinNeighbor(this.BinChildren(binNo,6),17)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        this.BinNeighbor(this.BinChildren(binNo,6),18)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        this.BinNeighbor(this.BinChildren(binNo,8),15)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        this.BinNeighbor(this.BinChildren(binNo,8),18)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),8)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),9)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),7)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),8)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),9)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),10)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),7)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),10)=this.BinChildren(binNo,4);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),2)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),2)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),2)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),2)=this.BinChildren(binNo,8);
        
        
    elseif direction==2
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(this.BinNeighbor(binNo,2),2) this.BinChildren(binNo,1) this.BinChildren(this.BinNeighbor(binNo,2),4) this.BinChildren(binNo,3) this.BinChildren(this.BinNeighbor(binNo,2),6) this.BinChildren(binNo,5) this.BinChildren(this.BinNeighbor(binNo,2),8) this.BinChildren(binNo,7) ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),8)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        this.BinNeighbor(this.BinChildren(binNo,1),9)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        this.BinNeighbor(this.BinChildren(binNo,3),7)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        this.BinNeighbor(this.BinChildren(binNo,3),8)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        this.BinNeighbor(this.BinChildren(binNo,5),9)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        this.BinNeighbor(this.BinChildren(binNo,5),10)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        this.BinNeighbor(this.BinChildren(binNo,7),7)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        this.BinNeighbor(this.BinChildren(binNo,7),10)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),16)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),17)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),15)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),16)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),17)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),18)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),15)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),18)=this.BinChildren(binNo,3);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),1)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),1)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),1)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),1)=this.BinChildren(binNo,7);
        
        
    elseif direction==3
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(this.BinNeighbor(binNo,3),5) this.BinChildren(this.BinNeighbor(binNo,3),6) this.BinChildren(this.BinNeighbor(binNo,3),7) this.BinChildren(this.BinNeighbor(binNo,3),8) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4) ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),14)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        this.BinNeighbor(this.BinChildren(binNo,1),18)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        this.BinNeighbor(this.BinChildren(binNo,2),10)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        this.BinNeighbor(this.BinChildren(binNo,2),14)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        this.BinNeighbor(this.BinChildren(binNo,3),11)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        this.BinNeighbor(this.BinChildren(binNo,3),18)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        this.BinNeighbor(this.BinChildren(binNo,4),10)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        this.BinNeighbor(this.BinChildren(binNo,4),11)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),13)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),16)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),8)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),13)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),12)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),16)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),8)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),12)=this.BinChildren(binNo,2);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),4)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),4)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),4)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),4)=this.BinChildren(binNo,4);
        
        
        
    elseif direction==4
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinChildren(this.BinNeighbor(binNo,4),1) this.BinChildren(this.BinNeighbor(binNo,4),2) this.BinChildren(this.BinNeighbor(binNo,4),3) this.BinChildren(this.BinNeighbor(binNo,4),4)]'  ;
        
        this.BinNeighbor(this.BinChildren(binNo,5),13)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        this.BinNeighbor(this.BinChildren(binNo,5),16)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        this.BinNeighbor(this.BinChildren(binNo,6),8)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        this.BinNeighbor(this.BinChildren(binNo,6),13)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        this.BinNeighbor(this.BinChildren(binNo,7),12)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        this.BinNeighbor(this.BinChildren(binNo,7),16)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        this.BinNeighbor(this.BinChildren(binNo,8),8)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        this.BinNeighbor(this.BinChildren(binNo,8),12)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),14)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),18)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),10)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),14)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),11)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),18)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),10)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),11)=this.BinChildren(binNo,6);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),3)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),3)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),3)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),3)=this.BinChildren(binNo,8);
        
        
    elseif direction==5
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinChildren(this.BinNeighbor(binNo,5),1) this.BinChildren(this.BinNeighbor(binNo,5),2) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinChildren(this.BinNeighbor(binNo,5),5) this.BinChildren(this.BinNeighbor(binNo,5),6)  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,3),13)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        this.BinNeighbor(this.BinChildren(binNo,3),17)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        this.BinNeighbor(this.BinChildren(binNo,4),9)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        this.BinNeighbor(this.BinChildren(binNo,4),13)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        this.BinNeighbor(this.BinChildren(binNo,7),14)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        this.BinNeighbor(this.BinChildren(binNo,7),17)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        this.BinNeighbor(this.BinChildren(binNo,8),9)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        this.BinNeighbor(this.BinChildren(binNo,8),14)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),12)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),15)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),7)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),12)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),11)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),15)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),7)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),11)=this.BinChildren(binNo,4);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),6)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),6)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),5),6)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),6),6)=this.BinChildren(binNo,8);
        
        
    elseif direction==6
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(this.BinNeighbor(binNo,6),3) this.BinChildren(this.BinNeighbor(binNo,6),4) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(this.BinNeighbor(binNo,6),7) this.BinChildren(this.BinNeighbor(binNo,6),8) this.BinChildren(binNo,5) this.BinChildren(binNo,6)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),12)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        this.BinNeighbor(this.BinChildren(binNo,1),15)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        this.BinNeighbor(this.BinChildren(binNo,2),7)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        this.BinNeighbor(this.BinChildren(binNo,2),12)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        this.BinNeighbor(this.BinChildren(binNo,5),11)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        this.BinNeighbor(this.BinChildren(binNo,5),15)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        this.BinNeighbor(this.BinChildren(binNo,6),7)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        this.BinNeighbor(this.BinChildren(binNo,6),11)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),13)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),17)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),9)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),13)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),14)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),17)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),9)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),14)=this.BinChildren(binNo,2);
        
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),5)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),5)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),7),5)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),8),5)=this.BinChildren(binNo,6);
        
        
    elseif direction==7
        %                 Output=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(this.BinNeighbor(binNo,7),4);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(this.BinNeighbor(binNo,7),8);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,5);
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,1),direction),17)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,5),direction),17)=this.BinChildren(binNo,5);
        
        
    elseif direction==8
        %                 Output=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(this.BinNeighbor(binNo,8),2);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(this.BinNeighbor(binNo,8),4);
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,5),direction),18)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,7),direction),18)=this.BinChildren(binNo,7);
        
        
    elseif direction==9
        %                 Output=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(this.BinNeighbor(binNo,9),2);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(this.BinNeighbor(binNo,9),6);
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,3),direction),15)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,7),direction),15)=this.BinChildren(binNo,7);
        
    elseif direction==10
        %                 Output=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(this.BinNeighbor(binNo,10),6);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(this.BinNeighbor(binNo,10),8);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,3) ;
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,1),direction),16)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,3),direction),16)=this.BinChildren(binNo,3);
        
        
    elseif direction==11
        %                 Output=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(this.BinNeighbor(binNo,11),7);
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(this.BinNeighbor(binNo,11),8);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,2) ;
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,1),direction),13)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,2),direction),13)=this.BinChildren(binNo,2);
        
        
    elseif direction==12
        %                 Output=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(this.BinNeighbor(binNo,12),3);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(this.BinNeighbor(binNo,12),4) ;
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,5),direction),14)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,6),direction),14)=this.BinChildren(binNo,6);
        
        
    elseif direction==13
        %                 Output=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(this.BinNeighbor(binNo,13),1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(this.BinNeighbor(binNo,13),2) ;
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,7),direction),11)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,8),direction),11)=this.BinChildren(binNo,8);
        
    elseif direction==14
        %                 Output=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(this.BinNeighbor(binNo,14),5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(this.BinNeighbor(binNo,14),6);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,4) ;
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,3),direction),12)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,4),direction),12)=this.BinChildren(binNo,4);
        
        
    elseif direction==15
        %                 Output=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(this.BinNeighbor(binNo,15),3);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(this.BinNeighbor(binNo,15),7);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,6) ;
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,2),direction),9)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,6),direction),9)=this.BinChildren(binNo,6);
        
        
    elseif direction==16
        %                 Output=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(this.BinNeighbor(binNo,16),1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(this.BinNeighbor(binNo,16),3);
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,6),direction),10)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,8),direction),10)=this.BinChildren(binNo,8);
        
        
    elseif direction==17
        %                 Output=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(this.BinNeighbor(binNo,17),1);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(this.BinNeighbor(binNo,17),5);
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,4),direction),7)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,8),direction),7)=this.BinChildren(binNo,8);
        
    elseif direction==18
        %                 Output=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(this.BinNeighbor(binNo,18),5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(this.BinNeighbor(binNo,18),7);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,4);
        
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,2),direction),8)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinNeighbor(this.BinChildren(binNo,4),direction),8)=this.BinChildren(binNo,4);
        
        
    end
    
    
    end
    
    
    
    
    function updateNeighborMatrixMinus(this,binNo,direction)
    
    if direction==1
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1) this.BinChildren(binNo,8) this.BinNeighbor(binNo,1)  ]';
        
        %                 this.BinNeighbor(this.BinChildren(binNo,2),16)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        %                 this.BinNeighbor(this.BinChildren(binNo,2),17)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        %                 this.BinNeighbor(this.BinChildren(binNo,4),15)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        %                 this.BinNeighbor(this.BinChildren(binNo,4),16)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        %                 this.BinNeighbor(this.BinChildren(binNo,6),17)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        %                 this.BinNeighbor(this.BinChildren(binNo,6),18)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        %                 this.BinNeighbor(this.BinChildren(binNo,8),15)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        %                 this.BinNeighbor(this.BinChildren(binNo,8),18)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        
    elseif direction==2
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7) ]';
        
        %                 this.BinNeighbor(this.BinChildren(binNo,1),8)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        %                 this.BinNeighbor(this.BinChildren(binNo,1),9)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        %                 this.BinNeighbor(this.BinChildren(binNo,3),7)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        %                 this.BinNeighbor(this.BinChildren(binNo,3),8)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        %                 this.BinNeighbor(this.BinChildren(binNo,5),9)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        %                 this.BinNeighbor(this.BinChildren(binNo,5),10)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        %                 this.BinNeighbor(this.BinChildren(binNo,7),7)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        %                 this.BinNeighbor(this.BinChildren(binNo,7),10)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        
        
    elseif direction==3
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4) ]';
        
        %                 this.BinNeighbor(this.BinChildren(binNo,1),14)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        %                 this.BinNeighbor(this.BinChildren(binNo,1),18)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        %                 this.BinNeighbor(this.BinChildren(binNo,2),10)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        %                 this.BinNeighbor(this.BinChildren(binNo,2),14)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        %                 this.BinNeighbor(this.BinChildren(binNo,3),11)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        %                 this.BinNeighbor(this.BinChildren(binNo,3),18)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        %                 this.BinNeighbor(this.BinChildren(binNo,4),11)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        %                 this.BinNeighbor(this.BinChildren(binNo,4),18)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        
        
    elseif direction==4
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4)]'  ;
        
        %                 this.BinNeighbor(this.BinChildren(binNo,5),13)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        %                 this.BinNeighbor(this.BinChildren(binNo,5),16)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        %                 this.BinNeighbor(this.BinChildren(binNo,6),8)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        %                 this.BinNeighbor(this.BinChildren(binNo,6),13)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        %                 this.BinNeighbor(this.BinChildren(binNo,7),12)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        %                 this.BinNeighbor(this.BinChildren(binNo,7),16)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        %                 this.BinNeighbor(this.BinChildren(binNo,8),8)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        %                 this.BinNeighbor(this.BinChildren(binNo,8),12)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        
        
    elseif direction==5
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
        
        %                 this.BinNeighbor(this.BinChildren(binNo,3),13)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        %                 this.BinNeighbor(this.BinChildren(binNo,3),17)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        %                 this.BinNeighbor(this.BinChildren(binNo,4),9)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        %                 this.BinNeighbor(this.BinChildren(binNo,4),13)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        %                 this.BinNeighbor(this.BinChildren(binNo,7),14)=this.BinChildren(this.BinNeighbor(binNo,direction),1);
        %                 this.BinNeighbor(this.BinChildren(binNo,7),17)=this.BinChildren(this.BinNeighbor(binNo,direction),6);
        %                 this.BinNeighbor(this.BinChildren(binNo,8),9)=this.BinChildren(this.BinNeighbor(binNo,direction),5);
        %                 this.BinNeighbor(this.BinChildren(binNo,8),14)=this.BinChildren(this.BinNeighbor(binNo,direction),2);
        
    elseif direction==6
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)   ]';
        
        %                 this.BinNeighbor(this.BinChildren(binNo,1),12)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        %                 this.BinNeighbor(this.BinChildren(binNo,1),15)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        %                 this.BinNeighbor(this.BinChildren(binNo,2),7)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        %                 this.BinNeighbor(this.BinChildren(binNo,2),12)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        %                 this.BinNeighbor(this.BinChildren(binNo,5),11)=this.BinChildren(this.BinNeighbor(binNo,direction),3);
        %                 this.BinNeighbor(this.BinChildren(binNo,5),15)=this.BinChildren(this.BinNeighbor(binNo,direction),8);
        %                 this.BinNeighbor(this.BinChildren(binNo,6),7)=this.BinChildren(this.BinNeighbor(binNo,direction),7);
        %                 this.BinNeighbor(this.BinChildren(binNo,6),11)=this.BinChildren(this.BinNeighbor(binNo,direction),4);
        
        
    elseif direction==7
        %                 Output=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinNeighbor(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinNeighbor(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,5);
        
        
    elseif direction==8
        %                 Output=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinNeighbor(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinNeighbor(binNo,8);
        
        
    elseif direction==9
        %                 Output=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinNeighbor(binNo,9);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinNeighbor(binNo,9);
        
    elseif direction==10
        %                 Output=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinNeighbor(binNo,10);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinNeighbor(binNo,10);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,3) ;
        
        
    elseif direction==11
        %                 Output=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinNeighbor(binNo,11);
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinNeighbor(binNo,11);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,2) ;
        
        
    elseif direction==12
        %                 Output=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinNeighbor(binNo,12);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinNeighbor(binNo,12) ;
        
        
    elseif direction==13
        %                 Output=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinNeighbor(binNo,13);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinNeighbor(binNo,13) ;
        
        
    elseif direction==14
        %                 Output=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
        
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinNeighbor(binNo,14);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinNeighbor(binNo,14);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,4) ;
        
        
    elseif direction==15
        %                 Output=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinNeighbor(binNo,15);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinNeighbor(binNo,15);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,6) ;
        
        
    elseif direction==16
        %                 Output=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinNeighbor(binNo,16);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinNeighbor(binNo,16);
        
        
    elseif direction==17
        %                 Output=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
        
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinNeighbor(binNo,17);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinNeighbor(binNo,17);
        
        
    elseif direction==18
        %                 Output=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinNeighbor(binNo,18);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinNeighbor(binNo,18);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,4);
        
        
    end
    
    
    end
    
    
    
    
    
    function Child=updateNeighborMatrix(this,ChildNo,direction,binNo)
    
    if ChildNo==1
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)   ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
        
    elseif ChildNo==2
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)   ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinChildren(this.BinNeighbor(binNo,10),5) -1 this.BinChildren(this.BinNeighbor(binNo,10),5) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
        
        
    elseif ChildNo==3
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)   ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
        
    elseif ChildNo==4
        
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)   ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1 this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
        
    elseif ChildNo==5
        
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
    elseif ChildNo==6
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16) ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
        
    elseif ChildNo==7
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
        
        
        
    elseif ChildNo==8
        if direction==1
            Child=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1)  this.BinChildren(binNo,8) this.BinNeighbor(binNo,1) ]';
            
        elseif direction==2
            Child=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7)  ]';
            
        elseif direction==3
            Child=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4)]';
            
        elseif direction==4
            Child=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) ]';
            
        elseif direction==5
            Child=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
            
        elseif direction==6
            Child=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)  ]';
            
        elseif direction==7
            Child=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
            
        elseif direction==8
            Child=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
            
        elseif direction==9
            Child=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
            
        elseif direction==10
            Child=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
            
        elseif direction==11
            Child=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
            
        elseif direction==12
            Child=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
            
        elseif direction==13
            Child=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
            
        elseif direction==14
            Child=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
            
        elseif direction==15
            Child=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
            
        elseif direction==16
            Child=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
            
        elseif direction==17
            Child=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
            
        elseif direction==18
            Child=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        end
    end
    
    end
    
    function Output=updateNeighborMatrixRegular(this,binNo,direction)
    
    if direction==1
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1) this.BinChildren(binNo,8) this.BinNeighbor(binNo,1)  ]';
        
    elseif direction==2
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7) ]';
        
    elseif direction==3
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4) ]';
        
    elseif direction==4
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4)]'  ;
        
    elseif direction==5
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
        
    elseif direction==6
        this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)   ]';
        
    elseif direction==7
        %                 Output=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinNeighbor(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinNeighbor(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,5);
    elseif direction==8
        %                 Output=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinNeighbor(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinNeighbor(binNo,8);
        
    elseif direction==9
        %                 Output=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinNeighbor(binNo,9);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinNeighbor(binNo,9);
        
    elseif direction==10
        %                 Output=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinNeighbor(binNo,10);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinNeighbor(binNo,10);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,3) ;
        
    elseif direction==11
        %                 Output=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinNeighbor(binNo,11);
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinNeighbor(binNo,11);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,1);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinChildren(binNo,2) ;
        
    elseif direction==12
        %                 Output=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,5);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinNeighbor(binNo,12);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinNeighbor(binNo,12) ;
        
    elseif direction==13
        %                 Output=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,7);
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinNeighbor(binNo,13);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinNeighbor(binNo,13) ;
        
    elseif direction==14
        %                 Output=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinNeighbor(binNo,14);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinNeighbor(binNo,14);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,3);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinChildren(binNo,4) ;
        
    elseif direction==15
        %                 Output=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinNeighbor(binNo,15);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinNeighbor(binNo,15);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,6) ;
        
    elseif direction==16
        %                 Output=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,6);
        this.BinNeighbor(this.BinChildren(binNo,3),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,6),direction)=this.BinNeighbor(binNo,16);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinNeighbor(binNo,16);
        
    elseif direction==17
        %                 Output=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
        this.BinNeighbor(this.BinChildren(binNo,1),direction)=this.BinChildren(binNo,4);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinNeighbor(binNo,17);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,8);
        this.BinNeighbor(this.BinChildren(binNo,8),direction)=this.BinNeighbor(binNo,17);
        
    elseif direction==18
        %                 Output=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
        this.BinNeighbor(this.BinChildren(binNo,2),direction)=this.BinNeighbor(binNo,18);
        this.BinNeighbor(this.BinChildren(binNo,4),direction)=this.BinNeighbor(binNo,18);
        this.BinNeighbor(this.BinChildren(binNo,5),direction)=this.BinChildren(binNo,2);
        this.BinNeighbor(this.BinChildren(binNo,7),direction)=this.BinChildren(binNo,4);
        
    end
    
    end
    %         function Output=updateNeighborMatrixRegular(this,binNo,direction)
    %
    %             if direction==1
    %                 Output=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,1) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinChildren(binNo,6) this.BinNeighbor(binNo,1) this.BinChildren(binNo,8) this.BinNeighbor(binNo,1)  ]';
    %
    %             elseif direction==2
    %                 Output=[this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinNeighbor(binNo,2) this.BinChildren(binNo,3) this.BinNeighbor(binNo,2) this.BinChildren(binNo,5) this.BinNeighbor(binNo,2) this.BinChildren(binNo,7) ]';
    %
    %             elseif direction==3
    %                 Output=[this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinChildren(binNo,3) this.BinChildren(binNo,4) ]';
    %
    %             elseif direction==4
    %                 Output=[this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4) this.BinNeighbor(binNo,4)]'  ;
    %
    %             elseif direction==5
    %                 Output=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5) this.BinChildren(binNo,7) this.BinChildren(binNo,8) this.BinNeighbor(binNo,5) this.BinNeighbor(binNo,5)  ]';
    %
    %             elseif direction==6
    %                 Output=[this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,1) this.BinChildren(binNo,2) this.BinNeighbor(binNo,6) this.BinNeighbor(binNo,6) this.BinChildren(binNo,5) this.BinChildren(binNo,6)   ]';
    %
    %             elseif direction==7
    %                 Output=[this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,1) this.BinNeighbor(binNo,7) -1 -1 this.BinChildren(binNo,5)   ]';
    %
    %             elseif direction==8
    %                 Output=[-1 this.BinChildren(binNo,5) -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,8) -1 this.BinNeighbor(binNo,8) -1   ]';
    %
    %             elseif direction==9
    %                 Output=[-1 this.BinChildren(binNo,3) this.BinNeighbor(binNo,9) -1 -1 this.BinChildren(binNo,7) this.BinNeighbor(binNo,9) -1   ]';
    %
    %             elseif direction==10
    %                 Output=[this.BinNeighbor(binNo,10) -1 this.BinNeighbor(binNo,10) -1 -1 this.BinChildren(binNo,1)  -1 this.BinChildren(binNo,3)   ]';
    %
    %             elseif direction==11
    %                 Output=[this.BinNeighbor(binNo,11) this.BinNeighbor(binNo,11) -1 -1 -1 -1  this.BinChildren(binNo,1) this.BinChildren(binNo,2)   ]';
    %
    %             elseif direction==12
    %                 Output=[-1 -1 this.BinChildren(binNo,5) this.BinChildren(binNo,6) this.BinNeighbor(binNo,12) this.BinNeighbor(binNo,12)  -1 -1   ]';
    %
    %             elseif direction==13
    %                 Output=[this.BinChildren(binNo,7) this.BinChildren(binNo,8) -1 -1 -1 -1  this.BinNeighbor(binNo,13) this.BinNeighbor(binNo,13)   ]';
    %
    %             elseif direction==14
    %                 Output=[-1 -1 this.BinNeighbor(binNo,14) this.BinNeighbor(binNo,14) this.BinChildren(binNo,3) this.BinChildren(binNo,4)  -1 -1   ]';
    %
    %             elseif direction==15
    %                 Output=[-1 this.BinNeighbor(binNo,15) this.BinChildren(binNo,2) -1 -1 this.BinNeighbor(binNo,15)  this.BinChildren(binNo,6) -1  ]';
    %
    %             elseif direction==16
    %                 Output=[this.BinChildren(binNo,6) -1 this.BinChildren(binNo,8) -1 -1 this.BinNeighbor(binNo,16)  -1 this.BinNeighbor(binNo,16)  ]';
    %
    %             elseif direction==17
    %                 Output=[this.BinChildren(binNo,4) -1 -1 this.BinNeighbor(binNo,17) this.BinChildren(binNo,8) -1  -1 this.BinNeighbor(binNo,17)  ]';
    %
    %             elseif direction==18
    %                 Output=[-1 this.BinNeighbor(binNo,18) -1 this.BinNeighbor(binNo,18) this.BinChildren(binNo,2) -1  this.BinChildren(binNo,4) -1  ]';
    %
    %             end
    %
    %         end
    
    
    
    function shrink(this)
    % Shrink all bins to bound only the points they contain
    % WARNING: this operation creates gaps in the final space not
    % covered by a bin. Only shrink OcTree structures when you only
    % intend to use the points used to create the tree to query the
    % tree space.
    binChildren = arrayfun(@(i)find(this.BinParents==i),1:this.BinCount,'Un',0)';
    binIsLeaf = cellfun(@isempty, binChildren);
    for i = find(binIsLeaf(:))'
        binShrink_recurse(i, true)
    end
    
        function binShrink_recurse(binNo, isLeafBin)
            % Build a list of all points that fall within one of the
            % bins to be checked, and the list of which point falls in
            % which bin.
            oldBoundaryMin = this.BinBoundaries(binNo,1:3);
            oldBoundaryMax = this.BinBoundaries(binNo,4:6);
            if isLeafBin
                % Shrink bin based on child POINTS
                ptsMask = this.PointBins==binNo;
                if ~any(ptsMask)
                    % No points, shrink the bin to infinitely small
                    proposedBoundaries = [oldBoundaryMin oldBoundaryMin];
                else
                    pts = this.Points(ptsMask,:);
                    proposedBoundaries = [...
                        max([oldBoundaryMin; min(pts,[],1)]) ...
                        min([oldBoundaryMax; max(pts,[],1)])];
                end
            else
                % Shrink bin based on child BINS
                childBoundaries = this.BinBoundaries(binChildren{binNo},:);
                proposedBoundaries = [min(childBoundaries(:,1:3),[],1) max(childBoundaries(:,4:6),[],1)];
            end
            
            if ~isequal(proposedBoundaries, [oldBoundaryMin oldBoundaryMax])
                % We just shrunk the boundary. Make it official and
                % check the parent
                this.BinBoundaries(binNo,:) = proposedBoundaries;
                parentBin = this.BinParents(binNo);
                if parentBin>0
                    binShrink_recurse(parentBin, false)
                end
            end
        end
    end
    
    function binNos = query(this, newPts, queryDepth)
    % Get the OcTree bins that new query points belong to.
    %
    % BINS = OT.query(NEWPTS) searches the OcTree object OT and
    % returns an N-by-1 vector of BINS giving the bin index in
    % which each of the N points in NEWPTS is contained. For any
    % query points outside all bins in OT, the index -1 is
    % returned.
    %
    % BINS = OT.query(NEWPTS,DEPTH) restricts the search to DEPTH
    % levels in the OT bin tree. Note that the first bin
    % (containing all other bins in OT) has DEPTH = 1.
    
    if nargin<3
        queryDepth = max(this.BinDepths);
    end
    
    numPts = size(newPts,1);
    newPts = permute(newPts,[3 2 1]);
    binNos = ones(numPts,1)*-1;
    
    binChildren = arrayfun(@(i)find(this.BinParents==i),1:this.BinCount,'Un',0)';
    binIsLeaf = cellfun(@isempty, binChildren);
    ptQuery_recurse(1:numPts, this.BinParents==0, 0)
    
        function ptQuery_recurse(newIndsToCheck_, binsToCheck, depth)
            % Build a list of all points that fall within one of the
            % bins to be checked, and the list of which point falls in
            % which bin.
            boundsToCheck = this.BinBoundaries(binsToCheck,:);
            [ptInBounds, subbinNo] = max(all(...
                bsxfun(@ge, newPts(:,:,newIndsToCheck_), boundsToCheck(:,1:3)) & ...
                bsxfun(@le, newPts(:,:,newIndsToCheck_), boundsToCheck(:,4:6))...
                ,2),[],1);
            
            if ~all(ptInBounds)
                % Special case usually when depth=0, where a point may
                % fall outside the bins entirely. This should only
                % happen once so let's fix it once and let subsequent
                % code rely on all points being in bounds
                binNos(newIndsToCheck_(~ptInBounds)) = -1;
                newIndsToCheck_(~ptInBounds) = [];
                subbinNo(~ptInBounds) = [];
            end
            binNosToAssign = binsToCheck(subbinNo);
            newIndsToAssign = newIndsToCheck_;
            binNos(newIndsToAssign) = binNosToAssign;
            
            % Allow a free exit when we reach a certain depth
            if depth>=queryDepth
                return;
            end
            
            % Otherwise, for all of the points we just placed into
            % bins, check which of the children of those bins those
            % same points fall into
            [unqbinNos, ~, unqGrpNos] = unique(binNosToAssign);
            for i = 1:length(unqbinNos)
                thisPtMask = unqGrpNos==i;
                if ~binIsLeaf(unqbinNos(i))
                    ptQuery_recurse(newIndsToCheck_(thisPtMask), binChildren{unqbinNos(i)}, depth+1)
                end
            end
            
        end
    end
    
    function h = plot(this,varargin)
    % OcTree.plot plots bin bounding boxes of an OcTree object
    %
    % H = OT.plot('name',value,...) allows you to specify any
    % properties of the bounding box lines that you would normally
    % supply to a plot(...,'name',value) command, and returns plot
    % object handles (one per bin) to H.
    hold on;
    h = zeros(this.BinCount,1);
    for i = 1:this.BinCount
        binMinMax = this.BinBoundaries(i,:);
        pts = cat(1, binMinMax([...
            1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
            1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
            nan(1,3), binMinMax([4 2 3; 4 2 6]),...
            nan(1,3), binMinMax([4 5 3; 4 5 6]),...
            nan(1,3), binMinMax([1 5 3; 1 5 6]));
        h(i) = plot3(pts(:,1),pts(:,2),pts(:,3),varargin{:});
    end
    end
    function h = plot3(this,varargin)
    % OcTree.plot plots bin bounding boxes of an OcTree
    %
    % See also OcTree.plot
    h = this.plot(varargin{:});
    end
    %         function UniformFlag(this,binNo,alpha1,alpha4,alpha5)
    
    function UniformFlag(this,binNo,alpha1,alpha2,radius)
    % alpha1=1 air , alpha2=2.25 conductor
    
    d=radius;
    
    % Left Part of almond
    a=.193333;
    b=0.064444;
    
    % Right Part of almond
    ee=4.833450;
    ff=1.61115;
    
    xmin=this.BinBoundaries(binNo,1); xmax=this.BinBoundaries(binNo,4);
    ymin=this.BinBoundaries(binNo,2); ymax=this.BinBoundaries(binNo,5);
    zmin=this.BinBoundaries(binNo,3); zmax=this.BinBoundaries(binNo,6);
    xmid=.5*(xmin+xmax); ymid=.5*(ymin+ymax); zmid=.5*(zmin+zmax);
    
    
    A=[xmin ymin zmin;
        xmax ymin zmin;
        xmin ymax zmin;
        xmax ymax zmin;
        xmin ymin zmax;
        xmax ymin zmax;
        xmin ymax zmax;
        xmax ymax zmax;
        xmid ymid zmin;
        xmid ymid zmax;
        xmid ymin zmid;
        xmid ymax zmid;
        xmin ymid zmid;
        xmax ymid zmid;
        
        ];
    
    x=A(:,1); y=A(:,2); z=A(:,3);
    dist1=zeros(14,1);
    
    for i=1:14
        
        if x(i)<0 && x(i)>-.4167*d
            dist1(i)=y(i).^2+((a/b)^2)*z(i).^2-(a*d)^2*(1-((x(i)/d)/.416667).^2);
        elseif x(i)>=0 && x(i)<0.5833*d
            dist1(i)=(y(i)/ee)^2+(z(i)/ff)^2-d^2*(sqrt(1-(x(i)/d/2.08335)^2)-0.96)^2;
        else
            dist1(i)=100 ; % nose is outside of the almond
        end
    end
    
    if any(imag(dist1)~=0)
        this.alpha(binNo)=alpha1;
    elseif all(dist1>0)
        this.alpha(binNo)=alpha1;
    elseif all(dist1<0)
        this.alpha(binNo)=alpha2;
        
    else
        this.alpha(binNo)=0;
    end
    
    %
    
    end
    
end


end