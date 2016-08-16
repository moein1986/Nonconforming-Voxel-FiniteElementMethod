function alpha=MaterialFinder_CenterPoint(BinBoundaries,binNo,alpha,alpha1,alpha2,radius)
% alpha1=1 air , alpha2=2.25 conductor
 
d=radius;

% Left Part of almond
a=.193333;
b=0.064444;

% Right Part of almond
ee=4.833450;
ff=1.61115;


A=0.5*[BinBoundaries(:,1)+BinBoundaries(:,4)  BinBoundaries(:,2)+BinBoundaries(:,5) BinBoundaries(:,3)+BinBoundaries(:,6) ] ;

             
         
  x=A(1,1); y=A(1,2); z=A(1,3);              
 
    
     if x<0 && x>-.4167*d
         dist1=y.^2+((a/b)^2)*z.^2-(a*d)^2*(1-((x/d)/.416667).^2);
     elseif x>=0 && x<0.5833*d
         dist1=(y/ee)^2+(z/ff)^2-d^2*(sqrt(1-(x/d/2.08335)^2)-0.96)^2;
     else
         dist1=100 ; % nose is outside of the almond
     end

            if imag(dist1)~=0
                alpha(binNo)=alpha1;
            elseif dist1>0
                alpha(binNo)=alpha1;
            elseif dist1<0  
               alpha(binNo)=alpha2;
         
            end

%                        

        end