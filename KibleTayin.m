classdef KibleTayin
    properties (Dependent = true, Access = public)
        Ellipsoid
        MiddleMeridian
    end
    properties (Hidden = true, Access = protected)
        ellipsoid ReferenceEllipsoid
        dom double
    end
    properties (Constant = true, Hidden = true)
        KABE_PHI = dms2degrees([21, 25, 21.3400]); % B
        KABE_LAMBDA = dms2degrees([39, 49, 34.1400]); % L
        R = 6373924.115; % meter
    end
    
    methods
        function this = KibleTayin(ellipsoid, dom)
            this.Ellipsoid = ellipsoid;
            this.MiddleMeridian = dom;
        end
    end
    
    % GET METHODS
    methods
        % ellipsoid parameters
        function ep = get.Ellipsoid(this)
            ep = this.ellipsoid;
        end
        % slice middle meridian.
        function smm = get.MiddleMeridian(this)
            smm = this.dom;
        end
    end
    
    % SET METHODS
    methods
        function this = set.Ellipsoid(this, ep)
            this.ellipsoid = ep;
        end
        function this = set.MiddleMeridian(this, smm)
            this.dom = smm;
        end
    end
    
    methods
        %
        function [beta, S]= kibleTayini(this, xp, yp, x, y)
            [phi, lambda] = UTM3degrees2Geographic(this, xp, yp);
            l = this.MiddleMeridian;
            R_ = this.R;
            phi_kabe = deg2rad(this.KABE_PHI);
            lambda_kabe = deg2rad(this.KABE_LAMBDA);
            
            phi = deg2rad(phi);
            lambda = deg2rad(lambda);

            deltaLambda = lambda_kabe - lambda; % kabe - p
            A = pi-atan((cos(phi)*tan(phi_kabe)-sin(phi)*cos(deltaLambda))/ ...
                (sin(deltaLambda)+1e-13))-pi/2*marking(sin(deltaLambda)+1e-13); % p - kabe
            Kappa = atan(sin(phi)*tan(lambda-l*pi/180));
            S = R_*acos(sin(phi_kabe)*sin(phi)+cos(phi_kabe)*cos(phi)*cos(deltaLambda));
            
            alpha= A - Kappa; %p-pkabe
            alpha_ = pi-atan((xp-x)/((yp-y)+1e-13))-pi/2*marking((yp-y)+1e-13); % p1-p
            
            if (alpha - alpha_) < pi
                beta = alpha - alpha_  + pi;
            else
                beta = alpha - alpha_  - pi;
            end
            beta = beta * 200 / pi; % radian to grad.
        end
        
        function [B, L, varargout] = UTM3degrees2Geographic(this, x, y)
            y = y - 5e+5;
            l = this.dom;
            f = 1 / this.Ellipsoid.InverseFlattening;
            a = this.Ellipsoid.SemimajorAxis;
            c = a / (1 - f);
            e2 = (2 * f - f^2) / (1 - f)^2;
            
            A = c*(1-3/4*e2+45/64*e2^2-175/256*e2^3+11025/...
                16384*e2^4-43659/65536*e2^5)*pi/180; % m/degrees
            B_ = (3/8*e2-3/16*e2^2+213/2048*e2^3-255/4096 *e2^4+20861/524288*e2^5)*180/pi;
            C_ = (21/ 256*e2^2-21/256*e2^3+533/8192*e2^4-197/4096*e2^5)*180/pi;
            D_ = (151/6144*e2^3-453/12288 * e2^4+5019/131072*e2^5)*180/pi;
            
            sigma = deg2rad(x / A); % degrees to radians
            Bf = deg2rad(rad2deg(sigma) + (B_ * sin(2 * sigma)) + ...
                (C_ * sin(4 * sigma)) + (D_ * sin(6 * sigma))); % degrees to radians
            
            t = tan(Bf);
            nf2 = e2 * (cos(Bf))^2;
            V = sqrt(1 + nf2);
            Nf = c / V;
            
            B1 = (180/(4*atan(1)))/(Nf*cos(Bf));
            B2 = (180/(4*atan(1)))*t*(-1-nf2)/(2*Nf^2);
            B3 = (180/(4*atan(1)))*(-1-2*t^2-nf2)/(6*Nf^3*cos(Bf));
            B4 = (180/(4*atan(1)))*t*(5+3*t^2+6*nf2-6*t^2*nf2)/(24*Nf^4);
            B5 = (180/(4*atan(1)))*(5+28*t^2+24*t^4)/(120*Nf^5*cos(Bf));
            
            c_ = (t*(180/(4*atan(1)))/Nf)*y+ ...
                (t*(180/(4*atan(1)))/(3*Nf^3))*(-1-t^2+nf2+2*nf2^2)*y^3 + ...
                (t*(180/(4*atan(1)))/(15*Nf^5))*(2+5*t^2+3*t^4)*y^5;
            
            varargout{1} = c_;
            
            B = rad2deg(Bf) + (B2 * y^2) + (B4 * y^4);
            L = l + (B1 * y) + (B3 * y^3) + (B5 * y^5); % l + l_; [Lo + l]
        end

        %
        %
        %
        %
        %
        function [A1_2, A2_1, S, varargout] = geographicJTP2(this ,phi, lambda, phi_, lambda_)
            R_ = this.R;
            phi = deg2rad(phi);
            lambda = deg2rad(lambda);
            phi_ = deg2rad(phi_);
            lambda_ = deg2rad(lambda_);

            dLambda = lambda_ - lambda;

            A1_2 = pi - atan((cos(phi) * tan(phi_) - sin(phi) * cos(dLambda)) / ...
                (sin(dLambda) + 1e-13)) - pi/2 * marking(sin(dLambda) + 1e-13);

            isA2_1 = (pi - atan((cos(dLambda) * sin(phi_) - tan(phi) * cos(phi_))/ (sin(dLambda) + 1e-13)) - ...
            pi/2 * marking(sin(dLambda) + 1e-13)) < pi;

            if isA2_1
                A2_1 = (pi - atan((cos(dLambda) * sin(phi_) - tan(phi) * cos(phi_)) / (sin(dLambda) + 1e-13)) - ...
                pi/2 * marking(sin(dLambda) + 1e-13)) + pi;
            else
                A2_1 = (pi - atan((cos(dLambda) * sin(phi_) - tan(phi) * cos(phi_)) / (sin(dLambda) + 1e-13)) - ...
                pi/2 * marking(sin(dLambda) + 1e-13)) - pi;
            end

            S = R_ * acos(sin(phi_) * sin(phi) + cos(phi_) * cos(phi) * cos(dLambda));

            A1_2 = rad2deg(A1_2);
            A2_1 = rad2deg(A2_1);
            S = round(S, 3);

            varargout{1} = rad2deg(dLambda);
        end

        %
        %
        %
        %
        %
        function [phi_, lambda_, A2_1] = geographicJTP1(this, phi, lambda, azimuth, S)
            R_ = this.R;

            phi = deg2rad(phi);
            lambda = deg2rad(lambda);
            azimuth = deg2rad(azimuth);

            phi_ = asin(sin(phi) * cos(S / R_) + cos(phi) * sin(S / R_) * cos(azimuth));
            lambda_ = atan(sin(azimuth)/(cos(phi) / tan(S / R_) - sin(phi) * cos(azimuth))) + lambda;

            if azimuth < pi
                A2_1 = (pi - atan((cos(azimuth) * cos(S / R_) - tan(phi) * sin(S / R_)) / (sin(azimuth) + 1e-13)) - ...
                pi/2 * marking(sin(azimuth) + 1e-13)) + pi;
            else
                A2_1 = (pi - atan((cos(azimuth) * cos(S / R_) - tan(phi) * sin(S / R_)) / (sin(azimuth) + 1e-13)) - ...
                pi/2 * marking(sin(azimuth) + 1e-13)) - pi;
            end

            phi_ = rad2deg(phi_);
            lambda_ = rad2deg(lambda_);
            A2_1 = rad2deg(A2_1);
        end

    end
    
end

% -- HELPER FUNCTION --
function out = marking(in)
    if in > 0
        out = 1;
    elseif in < 0
        out = -1;
    else
        out = 0;
    end
end