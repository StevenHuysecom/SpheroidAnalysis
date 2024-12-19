function [Geo] = cartesianToGeodetic( X, axis )
    %
    % Converter Cartesian coordinates to Geodetic coordinates 
    %on Triaxial Ellipsoid or Biaxial Ellipsoid or Sphere
    %
    %   (x/a)^2+(y/b)^2+(z/c)^2=1   Triaxial Ellipsoid Equation
    %    Cartesian To Geodetic  x y z ==> B L h
    
    % Parameters:
    % * X, [x y z]     - Cartesian coordinatesdata, n x 3 matrix or three n x 1 vectors
    % * axis,[a; b; c] - ellipsoid radii  [a; b; c],its axes % along [x y z] axes
    %  
    %                  For Triaxial ellipsoid ,it must be a > b > c
    %
    %                  For Biaxial ellipsoid ,it must be a = b > c
    %
    %                  For Sphere ,it must be a = b = c
    %
    % Output:
    % * Geo,[B,L,h]  -  Geodetic coordinates [latitude(deg);longitude(deg);ellipsoidal height(m)]
    % 
    %
    % Author:
    % Sebahattin Bektas, 19 Mayis University, Samsun
    % sbektas@omu.edu.tr
    
    format long
    ro=180/pi; % converter Degree to radian
    eps=0.0005; % three sholder
    a=axis(1);b=axis(2);c=axis(3);
    for NumOfPoints = 1:size(X, 1)
        x=X(1);y=X(2);z=X(3);
           
        ex2=(a^2-c^2)/a^2; ee2=(a^2-b^2)/a^2;
        
        E=1/a^2;F=1/b^2;G=1/c^2;
        
        xo=a*x/sqrt(x^2+y^2+z^2);
        yo=b*y/sqrt(x^2+y^2+z^2);
        zo=c*z/sqrt(x^2+y^2+z^2);
        
        for i=1:20
        j11=F*yo-(yo-y)*E;
        j12=(xo-x)*F-E*xo;
        
        j21=G*zo-(zo-z)*E;
        j23=(xo-x)*G-E*xo;
        
        A=[ j11   j12   0 
            j21   0   j23
            2*E*xo    2*F*yo  2*G*zo  ];
        
        sa=(xo-x)*F*yo-(yo-y)*E*xo;
        sb=(xo-x)*G*zo-(zo-z)*E*xo;
        se=E*xo^2+F*yo^2+G*zo^2-1;
        Ab=[ sa  sb  se]';
        bil=-A\Ab;
        xo=xo+bil(1);
        yo=yo+bil(2);
        zo=zo+bil(3);
        
        if max(abs(bil))<eps
            break
        end
        end
        
         fi=ro*atan(zo*(1-ee2)/(1-ex2)/sqrt((1-ee2)^2*xo^2+yo^2));
        
         l=ro*atan(1/(1-ee2)*yo/xo);
        
        h=sign(z-zo)*sign(zo)*sqrt((x-xo)^2+(y-yo)^2+(z-zo)^2);
        
        Geo(i,:)=[fi l h];
    end
end