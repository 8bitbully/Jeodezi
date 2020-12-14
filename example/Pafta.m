% Pafta Parser:
% Pafta.PAFTA_100 ; Pafta.PAFTA_50 ; Pafta.PAFTA_25 ; 
% Pafta secimi yapildiktan sonra Pafta objesini methodlarda
% kullanarak B, L koordinatlarini ve koordinat indexleri, alan indexleri
% kenar indexlerini elde edebiliriz.
% 
% """
%   Bu uygulama Jeodezi-I dersi odevi için yapilmistir ve odev cozumu icin
%   uygulanmalidir. Aksi halde yanlıs sonuclar elde edilebilir.
% """
%       2--------------3
%       |                   |
%       |   100000    |
%       |                   |
%       |                   |
%      [1]-------------4
%
%       paf = Pafta.PAFTA_100;
%
%       [B, L] = paftaParser(paf, b, l) ;
%       """
%       b: 100000 olcekli paftanin sol alt cografi koordinati
%       l:  100000 olcekli paftanin sol alt cografi koordinati
%       """
%       
%       [cidx, area, edge] = paftaEdge(paf);
%       """
%       cidx: koordinat index
%       area: alan index
%       edge: kenar index
%       """
%

%
% @author:
% @date: 20202811
%
classdef Pafta
    properties
        INDEX
    end
    enumeration
        PAFTA_100  ([
            0x02 0x03
            0x01 0x04])
        
        PAFTA_50  ([
            0x02 0x06 0x03
            0x05 0x09 0x07
            0x01 0x08 0x04])
        
        PAFTA_25  ([
            0x02 0x0C 0x06 0x0D 0x03
            0x0B 0x14 0x15 0x16 0x0E
            0x05 0x13 0x09 0x17 0x07
            0x0A 0x12 0x19 0x18 0x0F
            0x01 0x11 0x08 0x10 0x04])
    end
    methods
        function this = Pafta(pafta)
            this.INDEX = pafta;
        end
        
        function [B, L] = paftaParser(parser, b, l)
            idx = parser.INDEX;
            [m, n] = size(idx);
            b_ = repmat((b+.5:-(0.5/(m-1)):b)', [1, m]);
            l_ = repmat((l:(0.5/(m-1)):l+.5), [n, 1]);
            
            for i = 1:numel(idx)
                x = find(idx == i);
                geo.b(i, 1) = b_(x);
                geo.l(i, 1) = l_(x);
            end
            B = geo.b';
            L = geo.l';
        end
        
        function [cidx, area, edge] = paftaEdge(parser)
            idx = parser.INDEX;
            [m, ~] = size(idx);
            div = (m-1) * (m-1);
            
            if div / 4 == 1
                area = edgeParse(idx);
                edge = findEdge(area);
                reverse = area';
                cidx = reverse(:)';
            elseif div / 4 == 4
                ed1 = edgeParse(idx(3:5, 1: 3));
                ed2 = edgeParse(idx(1:3, 1: 3));
                ed3 = edgeParse(idx(1:3, 3: 5));
                ed4 = edgeParse(idx(3:5, 3: 5));
                area = [ed1; ed2; ed3; ed4];
                edge = findEdge(area);
                reverse = area';
                cidx = reverse(:)';
            elseif div == 1
                area = [1 2 3 4];
                edge = findEdge(area);
                reverse = area';
                cidx = reverse(:)';
            end
        end
    end
    
end
% HELPER FUNCTION
function out = edgeParse(pafta)
    out = flip(pafta);
    count = 1;
    for j = 1 : 2
        for i = 1 : 2
            chc.parse(count, :) = [out(i,j), out(i+1, j), out(i+1, j+1), out(i, j+1)];
            count = count + 1;
        end
    end
    change = chc.parse(4,:);
    chc.parse(4,:) = chc.parse(3, :);
    chc.parse(3, :) = change;
    out = chc.parse;
end

function edge = findEdge(in)
    out.edge = [];
    for i = 1 : size(in, 1)
        e1 = [in(i, 1) in(i, 2)];
        e2 = [in(i, 3) in(i, 4)];
        e3 = [in(i, 2) in(i, 3)];
        e4 = [in(i, 1) in(i, 4)];
        ed = [e1; e2; e3; e4];
        out.edge = [out.edge; ed];
    end
    edge = out.edge;
end