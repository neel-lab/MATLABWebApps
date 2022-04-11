function glysgp = sgp1to2(allsgp,PTMid)
sgp1name = {'n','h','s',...
    'f','x','g'...
    ,'u','k',...
    'v','w','l','r','t'};
sgp2name = {'(??-?)HexNAc','(??-?)Hexose','(??-?)Deoxynonulosonate',...
    '(??-?)Deoxyhexose','(??-?)Pentose','(??-?)Nonulosonate',...
    '(??-?)Hexuronate','(??-?)Nonulosonate',...
    '(a2-3)Deoxynonulosonate','(a2-6)Deoxynonulosonate','(??-?)Sialolactone','(??-?)NeuNAz','(??-?)HexNAz'};
glysgp = '';
if PTMid == 1
    for i = length(allsgp):-1:1
        if ~ismember(allsgp(i),sgp1name)
            glysgp = [allsgp(i),glysgp];
        else
            for j = 1:length(sgp1name)
                if strcmpi(sgp1name{j},allsgp(i))
                    glysgp = [sgp2name{j},glysgp];
                    break
                end
            end
        end
    end
else
    if strcmpi(allsgp(1),'<') && strcmpi(allsgp(end),'>')
        glysgp = ['{(??-?)',allsgp,'}'];
    end
end