function [gdl_red] = Reduccion_CB(frec)
% Frecuencias superior de bandas de tercio
Banda_baja = [14.1;17.8;22.4;28.2;35.5;44.7;56.2;70.8;89.1;112;141;178;224;...
    282;355;447;562;708;891;1122;1413;1778;2239;2818;3548;4467;5623;7079;...
    8913;11220;14130;17780];
Banda_alta = [17.8;22.4;28.2;35.5;44.7;56.2;70.8;89.1;112;141;178;224;282;355;...
    447;562;708;891;1122;1413;1778;2239;2818;3548;4467;5623;7079;...
    8913;11220;14130;17780;22390];
% Encuentrar cuantos modos pertenecen a cada banda
counter = zeros(length(Banda_baja),1);
for i = 1: length(counter)
    counter(i) = length(frec(frec>Banda_baja(i) & frec<Banda_alta(i))); 
end
% Encontrar cuando hay mas de 5 modos
pos = find(counter>5);
% Nos quedamos con solo esas posiciones en las matrices
gdl_red = length((find(frec>=Banda_baja(pos(1)))));

end