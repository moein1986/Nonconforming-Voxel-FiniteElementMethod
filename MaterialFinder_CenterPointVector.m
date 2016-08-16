function alpha=MaterialFinder_CenterPointVector(BinBoundaries,alpha1,alpha2,radius)
% alpha1=1 air , alpha2=2.25 conductor
alpha=zeros(1,64);
d=radius;

% Left Part of almond
a=.193333;
b=0.064444;

% Right Part of almond
ee=4.833450;
ff=1.61115;


A=0.5*[BinBoundaries(:,1)+BinBoundaries(:,4)  BinBoundaries(:,2)+BinBoundaries(:,5) BinBoundaries(:,3)+BinBoundaries(:,6) ] ;

             
         
  x=A(:,1); y=A(:,2); z=A(:,3);              
 
  for ii=1:length(x)
      
      if x(ii)<0 && x(ii)>-.4167*d
          dist1=y(ii).^2+((a/b)^2)*z(ii).^2-(a*d)^2*(1-((x(ii)/d)/.416667).^2);
      elseif x(ii)>0 && x(ii)<0.5833*d
          dist1=(y(ii)/ee)^2+(z(ii)/ff)^2-d^2*(sqrt(1-(x(ii)/d/2.08335)^2)-0.96)^2;
      else
          dist1=100 ; % nose is outside of the almond
      end
      
      if imag(dist1)~=0
          alpha(ii)=alpha1;
      elseif dist1>0
          alpha(ii)=alpha1;
      elseif dist1<0
          alpha(ii)=alpha2;
          
      end
      
  end

        end