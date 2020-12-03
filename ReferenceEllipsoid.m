% ReferenceEllipsoid :
%
% Example:
%
% ellipsoid = ReferenceEllipsoid('grs80')
% ----
% B = 39.5 ; L = 37;
% [x, y, z] = geog2geoc(ellipsoid, B, L)
% ----
% x = 3719318.657; y = 3034764.344; z = 4185622.239;
% [B, L, h] = geoc2geog(ellipsoid, x, y, z)
% ----
% B1 = 39.5; B2 = 40; L1 = 37; L2 = 37.5;
% [A1, A2, S] = JTP2(ellipsoid, B1, B2, L1, L2)
% ----
% B1 = 39.5; L1 = 37; A1 = 141.6904861; S = 69876.5752;
% [B2, L2, A2] = JTP1(ellipsoid, B1, L1, A1, S)
% ----
% B = 39.5;
% G = lat2eqdist(ellipsoid, B)
% ----
% G = 4314212.496;
% B = dist2eqlat(ellipsoid, G)
% ----
% L1 = 37; L2 = 37.5; B = 39.5;
% Sp = longt2dist(ellipsoid, L1, L2, B)
% ----
% B1 = 39.5; B2 = 40; L1 = 37; L2 = 37.5;
% F = ellipsoidArea(ellipsoid, B1, B2, L1, L2)
% ----
% B = 39.5; A = 42.33333333;
% [N, M, R, Ra] = radiusCurveture(ellipsoid, B, A)
%
%

% References:
%       Comparison of Different Algorithms between Geocentric and Geodetic Coordinates
%       Harita Dergisi Temmuz 2011 Sayı 146
%
%       Doç.Dr. FARUK YILDIRIM | AVESİS. https://avesis.ktu.edu.tr/yfaruk/dokumanlar. Erişim 28 Kasım 2020.
%
%
% Version: 9.9.0.1467703 (R2020b)
%
% for the updated version : github.com/solounextracto
% @author: 
% @date: 20202811
% 
%               
classdef ReferenceEllipsoid < Ellipsoid
    properties
        Name string
    end
    properties (GetAccess = private, Constant = true)
        Ellipsoid = struct( ...
            "GRS80", struct( ...
                    "Name", "Geodetic Reference System 1980", ...
                    "SemimajorAxis", 6378137, ...
                    "SemiminorAxis", 6356752.314140347, ...
                    "InverseFlattening", 298.2572221008827), ...
            "HAYFORD", struct( ...
                    "Name", "Hayford 1909", ...
                    "SemimajorAxis", 6378388, ...
                    "SemiminorAxis", 6356911.94613, ...
                    "InverseFlattening", 297) ...
        );
    end
    methods
        function this = ReferenceEllipsoid(nm)
            format longG;
            name = upper(nm);
            reference = this.Ellipsoid.(name);
            this.Name = reference.Name;
            this.SemimajorAxis = reference.SemimajorAxis;
            this.SemiminorAxis = reference.SemiminorAxis;
            this.InverseFlattening = reference.InverseFlattening;
        end
    end
end