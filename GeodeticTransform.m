classdef GeodeticTransform < handle
    properties
        B
        L
        Zone
        Coordinate
        Selection
    end
    properties (Hidden = true)
        UTM
        Transform
    end
    
    methods
        function self = GeodeticTransform(varargin)
            if mod(length(varargin), 2)
                varargin(end) = [] ;
            end
            for count = 1 : 2 : length(varargin)
                self.(varargin{count}) = varargin{count + 1} ;
            end
        end
        
        function BL2UTM(self, degrees, origin)
            if nargin > 1
                self.Coordinate.degrees = degrees ;
            else
                self.Coordinate.degrees = '6' ;
            end
            for count = 1 : length(self.B)
            self.Zone = utmzone([self.B(count), self.L(count)]);
%             ellipsoid = referenceEllipsoid('international','m');
            ellipsoid = referenceEllipsoid('wgs84','m');
            self.UTM = defaultm('utm'); 
            self.UTM.zone = self.Zone;
            self.UTM.geoid = ellipsoid;
            self.UTM = defaultm(self.UTM); 
            if nargin > 2
                if self.Coordinate.degrees == '3'
                    self.UTM.scalefactor = 1 ;
                end
                self.UTM.origin = [0, origin, 0] ;
            end
            [y, x]= mfwdtran(self.UTM, self.B(count), self.L(count));
            gaussY = (y - 500000) / self.UTM.scalefactor ;
            gaussX = x / self.UTM.scalefactor;
            self.Coordinate.Gauss(count, :) = [gaussX, gaussY];
            if self.Coordinate.degrees == '6'
                y = str2double([self.Zone(1: end - 1), num2str(y)]);
            end
            self.Coordinate.UTM(count,:) = [x; y]';
            end
        end
        
        function geocentric(self, h)
            Bi=deg2rad(self.B);
            Li=deg2rad(self.L);
            ellipsoid = referenceEllipsoid('international', 'm');
            for count = 1 : length(self.B)
                [X, Y, Z] = geodetic2ecef(Bi(count),Li(count),h(count), ellipsoid); %Elipsoid yüksekliðinden jeosentrik koordinatlar
                self.Coordinate.Geocentric(count, :)=[X; Y; Z]';
            end
        end
        
        function transform(self, inverse, foward, selectiontype)
            self.Selection = selectiontype ;
            self.selectionTransform();
            t = fitgeotrans(inverse, foward, self.Transform);% tform structure ile parametreler hesaplandý...
            [xnew, ynew]= transformPointsForward(t, self.Coordinate.UTM(:,1), self.Coordinate.UTM(:,2));% Yeni noktalara koordinat verme...
            self.Coordinate.TransformCoordinate = [xnew,  ynew];
        end
        
    end
    methods (Access = private)
        function selectionTransform(self)
            if strcmpi(self.Selection, 'Helmert')
                self.Transform = 'NonreflectiveSimilarity' ;
                % ...
            end
        end
    end
end