% Fonksiyon cagirici
% Tablo donduren fonksiyonlar.
% Odev icin hazırlanmıstır.
%
% [corner, edge, area, center, geo] = getFunctions()


% Version: 9.9.0.1467703 (R2020b)
% @Author:
% @date: 20201202
function [corner, edge, area, center, geo] = getFunctions()
    corner = @findCorner;
    edge = @findEdge;
    area = @findArea;
    center = @findCenterGravity;
    geo = @findGeocentric;
end

% Pafta köse koordinatlari:
% PAFTA_KOSE = [double(cidx)', B(cidx)', L(cidx)'];
function PAFTA_KOSE = findCorner(cidx, B, L, type)
    if nargin < 4
        type = 'degrees';
    end
    Latitude = B(cidx)';
    Longitude = L(cidx)';
    if strcmp(type, 'dms')
        for i = 1:length(cidx)
            dmsB = degrees2dms(Latitude(i));
            dmsL = degrees2dms(Longitude(i));
            out.Latitude{i} = [num2str(dmsB(1)), char(0176), num2str(dmsB(2)), char(39), num2str(dmsB(3)), char(34)];
            out.Longitude{i} = [num2str(dmsL(1)), char(0176), num2str(dmsL(2)), char(39), num2str(dmsL(3)), char(34)];
        end
        Latitude = out.Latitude';
        Longitude = out.Longitude';
    end
    Kose_No = cidx';

    PAFTA_KOSE = table(Kose_No, Latitude, Longitude);
end

% Kenar uzunluklari:
% KENAR = round([double(edge) ,KU.G, KU.Sp, KU.S, KU.diff], 5);
function KENAR = findEdge(e, edge, B, L)
    for i = 1 : size(edge, 1)
        KU.G(i,1) = abs(lat2eqdist(e, B(edge(i, 1))) - lat2eqdist(e, B(edge(i, 2))));
        if B(edge(i, 1)) == B(edge(i, 2))
            KU.Sp(i, 1) = longt2dist(e, L(edge(i, 1)), L(edge(i, 2)), B(edge(i, 1)));
            [~, ~, KU.S(i, 1)] = JTP2(e, B(edge(i, 1)), B(edge(i, 2)), L(edge(i, 1)), L(edge(i, 2)));
            KU.diff(i, 1) = KU.Sp(i, 1) - KU.S(i, 1);
        end
    end
    
    G = round(KU.G, 3);
    Sp = round(KU.Sp, 3);
    S = round(KU.S, 3);
    Fark = round(KU.diff, 5);
    Kenar = combineIndex(edge);

    KENAR = table(Kenar, G, Sp, S, Fark);
end

% Pafta Alani
function ALAN = findArea(e, area, B, L)
    for i = 1 : size(area,1)
        PA.F(i, 1) = ellipsoidArea(e, B(area(i, 1)), B(area(i, 3)), L(area(i, 1)), L(area(i, 3)));
        [PA.m(i, 1), PA.hektar(i, 1), PA.donum(i, 1)] = areaConvert(PA.F(i, 1));
    end

    % ALAN = round([double(area), PA.F, PA.km, PA.hektar, PA.donum], 5);
    km2 = round(PA.F , 2);
    m2 = round(PA.m , 2);
    hektar = round(PA.hektar , 2);
    donum = round(PA.donum , 2);
    Area = combineIndex(area, 'F');

    ALAN = table(Area, m2, km2, hektar, donum);
end

% Agirlik merkezi yaricaplari
% - N - R - M -
function AGIRLIK = findCenterGravity(e, area, B)
    for i = 1: size(area, 1)
        MA.Bort(i, 1) = (B(area(i, 1)) + B(area(i, 3))) / 2;
        [MA.N(i, 1), MA.M(i, 1), MA.R(i, 1)] = radiusCurveture(e, MA.Bort(i, 1));
    end
    N = round(MA.N , 3);
    R = round(MA.R , 3);
    M = round(MA.M , 3);
    Nokta_No = combineIndex(area, 'O');
    % AGIRLIK = round([double(area), MA.N, MA.M, MA.R], 3);
    AGIRLIK = table(Nokta_No, N, R, M);
end

% Jeosantrik koordinatlar
% GEOCENTRIC = [double(cidx)', x', y', z'];
function GEOCENTRIC = findGeocentric(e, cidx, B, L)
    [x, y, z] = geog2geoc(e, B(cidx), L(cidx));
    X = x'; Y = y'; Z = z';
    Kose_No = cidx';
    GEOCENTRIC = table(Kose_No, X, Y, Z);
end

% -- HELPER FUNCTIONS --

%content = [prefix, num2str([idx(i,:)])];
function str = combineIndex(idx, prefix)
    if nargin < 2
        prefix = '';
    end
    if isstring(prefix)
        prefix = char(prefix);
    end
    for i = 1 : size(idx, 1)
        content = [prefix, num2str(idx(i,:))];
        content(isspace(content)) = [];
        out.IDX{i} = content;
    end
    str = out.IDX';
end

% Alanları cevirme
function [m, hektar, donum] = areaConvert(km)
    m = km * 1e6;
    hektar = km * 1e2;
    donum =  km * 1e3;
end