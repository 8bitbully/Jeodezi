clc;clearvars;

% Reference Ellipsoid

ellipsoid = ReferenceEllipsoid('hayford'); % grs80

%%
% geocentric

% B = 41.27737178;
% L = 39.21259856;

% reduced 

% beta = 39.08885772;
% L = 32;

% [x, y, z] = geog2geoc(ellipsoid, B, L);

%%
% 
% x = 3794779.641;
% y = 3072951.963;
% z = 4089718.725;
% 
% [b, l, h] = geoc2geog(ellipsoid, x, y, z);

%%

% JTP2:

% B1 = 40.19043611;
% L1 = 33.20710278;
% B2 = 41.74048611;
% L2 = 33.72431111;
% 
% [A1, A2, S] = JTP2(ellipsoid, B1, B2, L1, L2);
%%

% JTP1:

% B1 = 39.50765833;
% L1 = 39.35715278;
% A1 = 141.6904861;
% S = 69876.5752;
% 
% [B2, L2, A2] = JTP1(ellipsoid, B1, L1, A1, S);

%%
% G to B;

% G = 4314212.496;
% B = dist2eqlat(ellipsoid, G);

%%
% B to G;

% B = 38.0000;
% G = lat2eqdist(ellipsoid, B);

%%
% Sp

% L1 = 26;
% L2 = 45;
% B = 36;
% 
% Sp = longt2dist(ellipsoid, L1, L2, B);

%% 
% Area

% B1 = 39; B2 = 40;
% L1 = 28; L2 = 29;
% 
% F = ellipsoidArea(ellipsoid, B1, B2, L1, L2);

%%
% radius of curvature

% B = 41.27737178;
% A = 42.33333333;

% [M, N, R] = radiusCurveture(ellipsoid, B);
% [M, N, R, Ra] = radiusCurveture(ellipsoid, B, A);
