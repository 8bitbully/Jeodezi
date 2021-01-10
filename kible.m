clearvars, clc,

% -- cami veri
B = [40, 59, 44.03];
L = [39, 40, 58.23];

B = dms2degrees(B);
L = dms2degrees(L);

% 6 dereceklik Google earth dönüşüm.
% Easting = 557432.00;
% Northing = 4538489.00;

% Easting = 4540384.68182062;
% Northing = 557458.131094429;

% 3 dereceklik matlab dönüşüm
% Easting = 557455.05421614;
% Northing = 4540305.24471597;

% Zone = '37T'; % 37T

% -- excel veri
% B = [40 59 52.192];
% L = [39 39 49.571];
% 
% B = dms2degrees(B);
% L = dms2degrees(L);

gt = GeodeticTransform('B', B, 'L', L);

BL2UTM(gt, '3', 39);

[xp_yp] = gt.Coordinate.UTM;

L0 = 39;
% xp = 4540544.663;
% yp = 555848.389;
xp = xp_yp(1);
yp = xp_yp(2);

x = 4540516.332;
y = 555820.152;

e = ReferenceEllipsoid('grs80');
k = KibleTayin(e, 39);

% [B, L] = UTM3degrees2Geographic(k, x, y);

[beta, S] = kibleTayini(k, xp, yp, x, y);