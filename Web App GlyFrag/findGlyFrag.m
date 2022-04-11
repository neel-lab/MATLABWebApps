function [Yion,Bion,Iion,Yunitind,Bunitind,Iunitind]=findGlyFrag(SmallGlyPep,breakPt,unitindex)
%GLYCANFRAG: Break input glycan/glycopeptide at break points and return fragments
%
% Syntax:
%    [NewGlyPep,Bion]=glycanFrag(SmallGlyPep,breakPt)
%
% Input: SmallGlyPep and breakPt array, denoting sites of glycan fragmentation
%
% Output: New glycopeptide that is reduced in size and cell array with Bions
%   Note: If the input is purely a glycan, then NewGlyPep will correspond to the y-ion
%   formed following either single or multiple fragmentation
%
% Example 1: for glycomics analysis
% >> SmallGlyPep='{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}';
%    breakPt=[5,26]
%    [Yion,Bion,Iion]=findGlyFrag(SmallGlyPep,breakPt)
% Answer:
%       breakPt =    5    26
%       Yion ='{n{n}}'
%       Bion ='{h}'
%       Iion='{h{h{h{h}}}{h{h{h}}{h}}}'
%
%See also glyMZCALC, FragmentGly, UQFRAGION, JOINGLYPEP, COMPILEFRAGS,GLYCANFRAG, BREAKGLYPEP,
%MULTISGPFRAG.

% Author: Sriram Neelamegham
% Date Lastly Updated: 8/11/14

nFrag=length(breakPt);
NewGlyPep=SmallGlyPep;
Bion=[];
Iion=[];
monosacind = unitindex{2}(:,1);
monosacstart = unitindex{2}(:,2);
monosacend = unitindex{2}(:,3);
Bunitind = {};
Iunitind = {};
Yunitind = zeros(0,1);
tempmonosacind = monosacind;
% Find terminal '}'
term=[regexp(SmallGlyPep,'}{'),length(SmallGlyPep)];
for i=1:length(term)
    ct=term(i);
    ctt=1;
    cluster(i,ctt)=term(i);
    while (SmallGlyPep(ct-1)=='}')
        ct=ct-1;
        ctt=ctt+1;
        cluster(i,ctt)=ct;
    end
end
for j=nFrag:-1:1   % start with outermost bond
    openBrac=1;
    closeBrac=0;
    glyBion='{';
    countTerm=breakPt(j);
    NewGlyPep(breakPt(j))='';
    bracketcount=0;
    while(openBrac~=closeBrac)
        if (NewGlyPep(breakPt(j))=='{')
            openBrac=openBrac+1;
            bracketcount = bracketcount + 1;
        end
        if (NewGlyPep(breakPt(j))=='}')
            closeBrac=closeBrac+1;
        end
        glyBion=[glyBion,NewGlyPep(breakPt(j))];
        countTerm=countTerm+1;
        NewGlyPep(breakPt(j))='';
    end
    tempBionind = find(monosacstart > breakPt(j),1,'first');
    tempsgpind = tempBionind:tempBionind + bracketcount;
    tempBionind = tempmonosacind(tempsgpind);
    tempmonosacind(tempsgpind) = [];
    
    if isempty(find(cluster == countTerm, 1))
        Iion=[Iion,cellstr(glyBion)];
        Iunitind = [Iunitind;unitindex{2}(tempBionind,:)];
    else
        Bion=[Bion,cellstr(glyBion)];
        rowDel=any((cluster==countTerm),2);
        cluster(rowDel,:)=[];
        if any(tempBionind)
            Bunitind = [Bunitind;unitindex{2}(tempBionind,:)];
        end
    end
    Yion=cellstr(NewGlyPep);
end
tempBIionmsind = [];
for i = 1:length(Bunitind)
    tempBIionmsind = [tempBIionmsind;Bunitind{i}(:,1)];
end
for i = 1:length(Iunitind)
    tempBIionmsind = [tempBIionmsind;Iunitind{i}(:,1)];
end
[~,diffind] = setdiff(monosacind,tempBIionmsind);
Yunitind = unitindex{2}(diffind,:);
end