%   Ellipsoid:
%   
%   [x, y, z] = geog2geoc(ellipsoid, B, L, latType, h)
%   [B, L, h] = geoc2geog(ellipsoid, x, y, z)
%   [A1, A2, S] = JTP2(ellipsoid, B1, B2, L1, L2)
%   [B2, L2, A2] = JTP1(ellipsoid, B1, L1, A1, S)
%   G = lat2eqdist(ellipsoid, B)
%   B = dist2eqlat(ellipsoid, G)
%   Sp = longt2dist(ellipsoid, L1, L2, B)
%   F = ellipsoidArea(ellipsoid, B1, B2, L1, L2)
%   [N, M, R, varargout] = radiusCurveture(ellipsoid, B, A)

% References:
%       Comparison of Different Algorithms between Geocentric and Geodetic Coordinates
%       Harita Dergisi Temmuz 2011 Sayı 146
%
%       Doç.Dr. FARUK YILDIRIM | AVESİS. https://avesis.ktu.edu.tr/yfaruk/dokumanlar. Erişim 28 Kasım 2020.
%
%
% for the updated version : github.com/solounextracto
% @author: 
% @date: 20202811
%
% Version: 9.9.0.1467703 (R2020b)
%
%   TODO:  * function signatures
%               * Error handling
%               * Test
classdef Ellipsoid < NivoEllipsoid
    properties (Dependent = true, Access = public)
        SemimajorAxis double
        SemiminorAxis double
        InverseFlattening double
        Eccentricity double
    end
    properties (GetAccess = public, SetAccess = private)
        SecondEccentricity = 0;
        Flattening = 0;
        ThirdFlattening = 0;
    end
    properties (Hidden = true, Access = protected)
        a = 1;
        b = 1;
        invf = inf;
        ecc = 0;
        ecc2 = 0;
    end
    properties (GetAccess = protected, Constant = true)
        AngularVelocity = 7292115e-11;
        EarthGravity = 3986004418e+8;
    end

    % GET METHODS
    methods
        function a = get.SemimajorAxis(this)
            a = this.a;
        end
        function b = get.SemiminorAxis(this)
            b = this.b;
        end
        function invf = get.InverseFlattening(this)
            invf = this.invf;
        end
        function ecc = get.Eccentricity(this)
            ecc = this.ecc;
        end
    end
    
    % SET METHODS
    methods
        function this = set.SemimajorAxis(this, a)
            this.a = a;
            this.b = (1 - this.Flattening) * a;
        end
        function this = set.SemiminorAxis(this, b)
            a_ = this.a;
            this.b = b;
            this.ecc = sqrt(a_*a_ - b*b) / a_;
            this.ecc2 = sqrt(a_*a_ - b*b) / b;
            f = (a_ - b) / a_;
            this.invf = 1 / f;
            this.Flattening = f;
            this.ThirdFlattening = f / (2 - f);
            this.SecondEccentricity = this.ecc2;
        end
        function this = set.InverseFlattening(this, invf)
            f = 1 / invf;
            this.ecc = sqrt((2 - f) * f);
            this.ecc2 = this.ecc / (1 - f);
            this.invf = invf;
            this.Flattening = f;
            this.ThirdFlattening = f / (2 - f);
            this.b = (1 - this.Flattening) * this.a;
            this.SecondEccentricity = this.ecc2;
        end
        function this = set.Eccentricity(this, ecc)
            this.ecc = ecc;
            e2 = ecc * ecc;
            this.b = (1 - f) * this.a;
            this.ecc2 = e2 - (1 - e2);
            this.invf = 1 / f;
            this.Flattening = f;
            this.ThirdFlattening = f / (2 - f);
            this.SecondEccentricity = this.ecc2;
        end
    end
    
    methods
        
        %   args:
        %       ellipsoid: ReferenceEllipsoid object
        %       B: geographic latitude
        %       L: geographic longitude
        %       latType: latitude type "geographic, geocentric, reduced"
        %       h: ellipsoid height "default: 0"
        %   returns:
        %       x: geocentric coordinate x
        %       y: geocentric coordinate y
        %       z: geocentric coordinate z
        function [x, y, z] = geog2geoc(ellipsoid, B, L, latType, h)
            if nargin < 4
                latType = 'geographic';
            end
            if nargin < 5; h = 0; end
            
            a_ = ellipsoid.SemimajorAxis;
            b_ = ellipsoid.SemiminorAxis;
            e = ellipsoid.Eccentricity; e2 = e * e;
            eu = ellipsoid.SecondEccentricity; eu2 = eu * eu;
            
            sinf = @sind;
            cosf = @cosd;
            tanf = @tand;
            atanf = @atand;
            
            switch latType
                case 'reduced'
                    B = atanf((a_/b_) * tanf(B));
                case 'geocentric'
                    B = atanf((a_*a_ / b_*b_) * tanf(B));
            end
            
            n2 = eu2 * cosf(B).^2;
            V = sqrt(1 + n2);
            W = sqrt(1 - e2 * sinf(B).^2);
            
            xv = (a_*a_ ./ (b_*V)+h) .* cosf(B) .* cosf(L);
            yv = (a_*a_ ./ (b_*V)+h) .* cosf(B) .* sinf(L);
            zv = (b_./V+h) .* sinf(B);
            
            xw = (a_./W+h) .* cosf(B) .* cosf(L);
            yw = (a_./W+h) .* cosf(B) .* sinf(L);
            zw = (b_*b_ ./ (a_*W)+h) .* sinf(B);
            
            x = round((xv + xw) / 2, 3);
            y = round((yv + yw) / 2, 3);
            z = round((zv + zw) / 2, 3);
        end
        
        %   args:
        %       "Pollard Method"
        %       ellipsoid: ReferenceEllipsoid object
        %       x: geocentric coordinate x
        %       y: geocentric coordinate y
        %       z: geocentric coordinate z
        %   returns:
        %       B: geographic latitude
        %       L: geographic longitude
        %       h: ellipsoid height
        function [B, L, h] = geoc2geog(ellipsoid, x, y, z)
            a_ = ellipsoid.SemimajorAxis;
            b_ = ellipsoid.SemiminorAxis;
            eu = ellipsoid.SecondEccentricity; eu2 = eu * eu;
            
            atanf = @atand;
            atan2f = @atan2d;
            
            x_ = x; y_ = y; z_ = z;
            z__ = (b_*z_) / sqrt(x_*x_ + y_*y_ + z_*z_);
            k_ = (a_ / b_) * (a_ / b_);
            z_old = inf;
            while abs(z__ - z_old) > 1e-7
                z_old = z__;
                k = sqrt(x_*x_ + y_*y_ + (z_ + eu2*z__)*(z_ + eu2*z__));
                l_ = x_ / k; m = y_ / k; n = (z_ + eu2*z__) / k;
                r = l_*l_ + m*m + k_*n*n;
                s = l_*x_ + m*y_ + k_*n*z_;
                t = x_*x_ + y_*y_ + k_*z_*z_ - a_*a_;
                h = (s + sqrt(s*s - r*t)) / r;
                root = sqrt(s*s - r*t);
                
                z__ = z_ - n*((s - root) / r);
                if h < 0; z__ = z_ - n*h; end
            end
            h = (s - root) / r;
            if h < 0; h = (s + root) / r; end
            
            B = atanf((z_ + eu2*z__) / sqrt(x_*x_ + y_*y_));
            L = atan2f(y_, x_);
        end
        
        %   args:
        %       "Vincenty Methods"
        %       ellipsoid: ReferenceEllipsoid object
        %       B1: latitude Pi
        %       B2: latitude Pk
        %       L1: longitude Pi
        %       L2: longitude Pk
        %   returns:
        %       A1: geodesic curve azimuth
        %       A2: geodesic curve azimuth
        %       S: edge distance
        function [A1, A2, S] = JTP2(ellipsoid, B1, B2, L1, L2)
            b_ = ellipsoid.SemiminorAxis;
            f = ellipsoid.Flattening;
            ec = ellipsoid.SecondEccentricity;
            ec2 = ec * ec;
            
            B1 = deg2rad(B1); B2 = deg2rad(B2);
            L1 = deg2rad(L1); L2 = deg2rad(L2); 
            
            u = atan((1 - f) * tan(B1));
            u_ = atan((1 - f) * tan(B2));
            w = L2 - L1 ; dL = w;
            
            dL_old = inf;
            while abs(dL - dL_old) > 1e-12
                dL_old = dL;
                snq = (sqrt((cos(u_) * sin(dL))^2 + ...
                    (cos(u) * sin(u_) - sin(u) * cos(u_) * cos(dL))^2));
                csq = sin(u) * sin(u_) + cos(u) * cos(u_) * cos(dL);
                
                q = atan2(snq, csq);
                
                Al = asin(cos(u) * cos(u_) * sin(dL) / snq);
                km = csq - (2 * sin(u) * sin(u_) / cos(Al)^2);
                c = (f / 16) * cos(Al)^2 * (4 + f * (4 - 3 * cos(Al)^2));
                dL = w + (1 - c) * f * sin(Al) * (q + c * snq * (km + c * csq * (-1 + 2 * km^2)));
            end
            
            u2_ = ec2 * cos(Al)^2;
            a__ = 1 + (u2_ / 16384) * (4096 + u2_ * (-768 + u2_ * (320 - 175 * u2_)));
            b__ = (u2_ / 1024) * (256 + u2_ * (-128 + u2_ * (74 - 47*u2_)));
            dq = b__ * snq * (km + (b__ / 4) * (csq * (-1 + 2 * km^2) - ...
                (b__ / 6) * km * (-3 + 4 * snq^2) * (-3 + 4 * km^2)));
            
            S = b_ * a__ *abs(q - dq);
            A1 = atan2(cos(u_) * sin(dL), (cos(u) * sin(u_) - sin(u) * cos(u_) * cos(dL)));
            A2 = atan2(cos(u) * sin(dL), (-sin(u) * cos(u_) + cos(u) * sin(u_) * cos(dL)));
            if A1 < pi
                A2 = rad2deg(A2 + pi);
            else
                A2 = rad2deg(A2 - pi);
            end
            A1 = rad2deg(A1);
        end
        
        %   args:
        %       "Vincenty Methods"
        %       ellipsoid: ReferenceEllipsoid object
        %       B1: latitude Pi
        %       L1: longitude Pi
        %       A1: geodesic curve azimuth
        %       S: edge distance
        %   returns:
        %       B2: latitude Pk
        %       L2: longitude Pk
        %       A2: geodesic curve azimuth
        function [B2, L2, A2] = JTP1(ellipsoid, B1, L1, A1, S)
            b_ = ellipsoid.SemiminorAxis;
            f = ellipsoid.Flattening;
            ec = ellipsoid.SecondEccentricity;
            ec2 = ec * ec;
            
            B1 = deg2rad(B1); L1 = deg2rad(L1); A1 = deg2rad(A1);
            
            u_ = atan((1 - f) * tan(B1));
            p_ = atan(tan(u_) / cos(A1));
            alpha = asin(cos(u_) * sin(A1));
            u2_ = cos(alpha) * cos(alpha) * ec2;
            
            a__ = 1 + (u2_ / 16384) * (4096 + u2_ * (-768 + u2_ * (320 - 175 * u2_)));
            b__ = (u2_ / 1024) * (256 + u2_ * (-128 + u2_ * (74 - 47 * u2_)));
            
            p = S / (b_ * a__);
            
            p_old = inf;
            while abs(p - p_old) > 1e-13
                p_old = p;
                qm2 = 2 * p_ + p;
                km = cos(qm2);
                dq = b__ * sin(p) * (km + (b__ / 4) * (cos(p) * (-1 + 2 * km * km) - ...
                    (b__ / 6) * km * (-3 + 4 * sin(p) * sin(p)) * (-3 + 4 * km * km)));
                p = S / (b_ * a__) + dq;
            end
            
            B2 = atan((sin(u_) * cos(p) + cos(u_) * sin(p) * cos(A1)) / ...
                ((1 - f) * (sin(alpha) * sin(alpha) + (sin(u_) * sin(p) - cos(u_) * cos(p) * cos(A1))^2) ^ (1 / 2)));
            L_ = atan((sin(p) * sin(A1)) / (cos(u_) * cos(p) - sin(u_) * sin(p) * cos(A1)));
            km2 = 2 * p_ + p;
            c = (f / 16) * cos(alpha) * cos(alpha) * (4 + f * (4 - 3 * cos(alpha) * cos(alpha)));
            dL = L_ - (1 - c) * f * sin(alpha) * ...
                (p + c * sin(p) * (cos(km2) + c * cos(p) * (-1 + 2 * cos(km2) * cos(km2))));
            L2 = L1 + dL;
            A2 = atan2(sin(alpha), (-sin(u_) * sin(p) + cos(u_) * cos(p) * cos(A1)));
            if A1 < pi
                A2 = rad2deg(A2 + pi);
            else
                A2 = rad2deg(A2 - pi);
            end
            B2 = rad2deg(B2);
            L2 = rad2deg(L2);
        end
        
        %   args:
        %       ellipsoid: ReferenceEllipsoid object
        %       B : latitude
        %   returns:
        %       G: meridian arc part
        function G = lat2eqdist(ellipsoid, B)
            a_  = ellipsoid.SemimajorAxis ;
            b_  = ellipsoid.SemiminorAxis ;
            c  = a_^2 / b_ ;
            eu2 = ellipsoid.SecondEccentricity^2 ;
            invro = pi / 180 ;

            sinf = @sind ;

            A_ = c * (1 - 3/4*eu2 + 45/64*eu2^2 - 175/256*eu2^3 + 11025/16384*eu2^4 - 43659/65536*eu2^5) * invro;
            B_ = c * (-3/8*eu2 + 15/32*eu2^2 - 525/1024*eu2^3 + 2205/4096*eu2^4 - 72765/131072*eu2^5) ;
            C_ = c * (15/256*eu2^2 - 105/1024*eu2^3 + 2205/16384*eu2^4 - 10395/65536*eu2^5) ;
            D_ = c * (-35/3072*eu2^3 + 315/12288*eu2^4 - 31185/786432*eu2^5) ;
            E_ = c * (315/131072*eu2^4 - 3465/524288*eu2^5 ) ;
            F_ = c * (-693/1310720*eu2^5) ;

            G = A_*B + B_*sinf(2*B) + C_*sinf(4*B) + D_*sinf(6*B) + E_*sinf(8*B) + F_*sinf(10*B);
        end
        
        %   args:
        %       ellipsoid: ReferenceEllipsoid object   
        %       G: meridian arc part
        %   returns:
        %       B: latitude
        function B = dist2eqlat(ellipsoid, G)
            a_  = ellipsoid.SemimajorAxis ;
            b_  = ellipsoid.SemiminorAxis ;
            c  = a_^2 / b_ ;
            eu2 = ellipsoid.SecondEccentricity^2 ;
            ro = 180 / pi ;

            sinf = @sind ;

            A_ = c * (1 - 3/4*eu2 + 45/64*eu2^2 - 175/256*eu2^3 + 11025/16384*eu2^4 - 43659/65536*eu2^5) / ro; %... * invro
            B__= (3/8*eu2 - 3/16*eu2^2 + 213/2048*eu2^3 - 255/4096*eu2^4 + 20861/524288*eu2^5) * ro ;
            C__= (21/256*eu2^2 - 21/256*eu2^3 + 533/8192*eu2^4 - 197/4096*eu2^5) * ro ;
            D__= (151/6144*eu2^3 - 453/12288*eu2^4 + 5019/131072*eu2^5) * ro ;
            E__= (1097/131072*eu2^4 - 1097/65536*eu2^5) * ro ;
            F__= (8011/2621440*eu2^5) * ro ;

            sigma = G / A_ ;

            B = sigma + B__*sinf(2*sigma) + C__*sinf(4*sigma) + D__*sinf(6*sigma) + ...
                E__*sinf(8*sigma) + F__*sinf(10*sigma) ;
        end
        
        %   args:
        %       ellipsoid: ReferenceEllipsoid object   
        %       L1: longitude Pi
        %       L2: longitude Pk
        %       B: latitude
        %   returns:
        %       Sp: parallel circle arc
        function Sp = longt2dist(ellipsoid, L1, L2, B)
            a_  = ellipsoid.SemimajorAxis ;
            b_  = ellipsoid.SemiminorAxis ;
            c   = a_^2 / b_ ;
            e2  = ellipsoid.Eccentricity^2 ;
            eu2 = ellipsoid.SecondEccentricity^2 ;
            
            sinf = @sind ;
            cosf = @cosd ;
            
            n2 = eu2 * cosf(B).^2 ;
            V = sqrt(1 + n2) ;
            W = sqrt(1 - e2*sinf(B).^2) ;
            
            N_ = c ./ V ;
            N__= a_ ./ W ;
            N = (N_ + N__) ./ 2 ;
            
            l = L2 - L1 ;
            
            Sp = N * l * cosf(B) * pi / 180; % pi / 180 -> radians
        end

        %   args:
        %       "Area of the ellipsoid bounded by latitude and longitude"
        %       ellipsoid: ReferenceEllipsoid object   
        %       B1: latitude Pi
        %       B2: latitude Pk
        %       L1: longitude Li
        %       L2: longitude Lk
        %   returns:
        %       F: Area
        function F = ellipsoidArea(ellipsoid, B1, B2, L1, L2)

            b_  = ellipsoid.SemiminorAxis ;
            e2  = ellipsoid.Eccentricity^2 ;
            
            Fa_ = 1 + 1/2*e2 + 3/8*e2^2 + 5/16*e2^3 + 35/128*e2^4 ;
            Fb_ = 1/6*e2 + 3/16*e2^2 + 3/16*e2^3 + 35/192*e2^4 ;
            Fc_ = 3/80*e2^2 + 1/16*e2^3 + 5/64*e2^4 ;
            Fd_ = 1/112*e2^3 + 5/256*e2^4 ;
            Fe_ = 5/2304*e2^4 ;
            
            dB = B2 - B1 ;
            dL = L2 - L1 ;
            Bm = (B1 + B2) / 2 ; % mean B
            
            sinf = @sind ;
            cosf = @cosd ;
            
            % Z = 2*pi*b^2 INTEGRAL(B1 -> B2) cos(B) / (1 - e^2*sin^2*B)^2
            % --> dB
            Z_ = 4 * pi * b_^2 * ( ...
                Fa_*sinf(dB/2)*cosf(Bm) - ...
                Fb_*sinf(3*dB/2)*cosf(3*Bm) + ...
                Fc_*sinf(5*dB/2)*cosf(5*Bm) - ...
                Fd_*sinf(7*dB/2)*cosf(7*Bm) + ...
                Fe_*sinf(9*dB/2)*cosf(9*Bm) ) ; % m^2
            
            Z = Z_ / 1e6 ; % km^2
            F = dL * Z / (2*pi) * pi / 180; % km^2
        end
        
        %   args:
        %       ellipsoid: ReferenceEllipsoid object   
        %       B: latitude
        %       A: is the radius of curvature of the azimuth normal section curve
        %   returns:
        %       N: Major Curvature radii
        %       M: Major Curvature radii
        %       R: gauss radius of curvature (average curvature)
        %       Ra: Radius of curvature of the normal section curve with azimuth
        function [N, M, R, varargout] = radiusCurveture(ellipsoid, B, A)
            a_ = ellipsoid.SemimajorAxis;
            eu = ellipsoid.SecondEccentricity; eu2 = eu * eu;
            invf_ = ellipsoid.InverseFlattening;
            c = a_ / (1 - 1 / invf_); 
            
            B = deg2rad(B);

            n2 = eu2 * cos(B)^2;
            V = sqrt(1 + n2);
 
            N = c / V;
            M = c / V^3;
            R = c / V^2;
   
            if nargin > 2
                A = deg2rad(A);
                Ra = M*N / (M*sin(A)^2 + N * cos(A)^2);
                varargout{1} = Ra;
            end
            if M > R || R > N || M > N
                warning('incorrect radius of curvature')
            end
        end
    end
end