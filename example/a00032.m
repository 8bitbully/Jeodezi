clearvars; clc;
% Version: 9.9.0.1467703 (R2020b)

main = localfunctions;
% 100000'lik paftanin sol alt B, L cografi koordinatlari girilir, grs80 ve
% hayford ellipsoidlerinde 100000, 50000, 25000'lik paftaların
% hesaplamalarını excel'e yazdırır.

vertexFunction(39.5, 37);

% cografi koordinatlari vertexFunction' a yazdıktan sonra, a00032.m
% scriptini calistirin. (Run or F5)
% cikan sonuclar output dosyasina gelir.

function vertexFunction(lat, long)
    ELLIPSOID_NAME = ["grs80", "hayford"];
    PAFTA_LIST = [Pafta.PAFTA_100, Pafta.PAFTA_50, Pafta.PAFTA_25];
    PAFTA_NAME = ["-100", "-50", "-25"];
    for i = 1 : length(ELLIPSOID_NAME)
        for j = 1 : length(PAFTA_NAME)
            [kose, kenar, alan, agirlik, jeosentrik] = getFunctions;
            % elipsoid secimi
            e = ReferenceEllipsoid(ELLIPSOID_NAME(i));

            % pafta secimi
            pf = PAFTA_LIST(j);

            % pafta indexleri
            [cidx, area, edge] = paftaEdge(pf);
            [B, L] = paftaParser(pf, lat, long);

            % Pafta köse koordinatlari:
            PAFTA_KOSEdegrees = kose(cidx, B, L); % degrees
            PAFTA_KOSEdms = kose(cidx, B, L, 'dms'); % degrees-minutes-seconds

            % Kenar uzunluklari:
            KENAR = kenar(e, edge, B, L);

            % Pafta Alani
            ALAN = alan(e, area, B, L);

            % Agirlik merkezi yaricaplari
            % return: - N - R - M -
            AGIRLIK = agirlik(e, area, B);

            % Jeosentrik Koordinatlar
            GEOCENTRIC = jeosentrik(e, cidx, B, L);

            writeExcel(['output/',char(ELLIPSOID_NAME(i)), char(PAFTA_NAME(j))], ...
                PAFTA_KOSEdegrees, PAFTA_KOSEdms, KENAR, ALAN, AGIRLIK, GEOCENTRIC)
        end
    end
end


% tablolari excel'e yazdirma fonksiyonu
%   args:
%       content: table // AGIRLIK, KENAR  etc.
%       filename: str // grs80-agirlik-50 etc.
%   returns:
%       excel file.
%
%   Example:
%   command window:
%   excelWriter('hayford-25', PAFTA_KOSE, KENAR, ALAN, AGIRLIK, GEOCENTRIC)
function writeExcel(filename, varargin)
    if ~isfolder('output'); mkdir('output'); end
    n = length(varargin);
%     T = table(content);
    filename = [filename, '.xlsx'];
    for i = 1 : n
%         writetable(T,filename,'Sheet','MyNewSheet','WriteVariableNames',false);
        writetable(varargin{i}, filename,'Sheet',inputname(i+1) ,'WriteVariableNames',true);
    end
end