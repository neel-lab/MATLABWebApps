function AllIons=FragmentGlyGUI(varargin)
% MULTISGPFRAG: Calculate all fragment ions formed theoretically from a
% glycan
%
% Syntax:
% Prod=FragmentGly(SmallGlyPep,form,adduct,ngFrag,z)
%
% Example:
% SmallGlyPep='B{n{h{s}}}';
% form='Me';
% ngFrag=1;
% z=1;
% adduct='Na';
%
% Input:
%   SmallGlyPep; glycan in SmallGlyPep format
%   ngFrag: # of glycan fragmentations
%   z: charge state.
%   Note:
%   This program will take a 'glycan' alone as input. glyBion and glyBion are
%   returned.
%
% Output:
%   AllFragIons: Structure containing all fragment ion data, i.e. product structure,
%   fragmentation type, ion type and m/z at charge state z etc. Note even
%   internal ions are generated in the case of multiple peptide fragmentation.
%
%  Example:
% adduct(1).name='Na';
% adduct(1).count=2;
% glyStruct(1).name='{n{n{h{h{h{h}}}{h{h{h}}{h}}}}}';
% glyStruct(1).z=2;
% glyStruct(1).adduct=adduct;
% glyStruct(1).ano='Me';
% glyStruct(1).form='Me';
% glyStruct(1).mz=glyMZCalc(glyStruct);
% ngFrag=2;
% Prod=FragmentGly(glyStruct,ngFrag);
%
% * Children function: glyMZCalc,findGlyFrag
%  joinGlyPep, glypepMW, glycanFrag, CompileFrags
%

% Author: Sriram Neelamegham

if (nargin>0)
    glyStruct = varargin{1};
end
if (nargin>1)
    nFrag=varargin{2};
    if (nargin>2)
       GlycanNew = varargin{3};
       AnomericNew=varargin{4};
       AdductNew=varargin{5};
    else
        GlycanNew = Glycan;
        AnomericNew = Anomeric;
        AdductNew = Adduct;
    end
else
    nFrag=2;
    GlycanNew = Glycan;
end

gly=glyStruct.name;
z=glyStruct(1).z;
adduct=glyStruct.adduct;
ano=glyStruct.ano;
form=glyStruct.form;
ion=glyStruct.ion;
mz=glyStruct.mz;
[st,ed] = regexp(gly,'(?<={)[^{}]+[{}]*?','start','end');
eind = zeros(0,1);  % standard empty index

% Unfragmented unit
Prod.name=gly;
Prod.z=z;
Prod.adduct=adduct;
Prod.ano=ano;
Prod.form=form;
Prod.ion = ion;
Prod.mz=mz;
Prod.nFrag=0;
Prod.type='original';
Prod.unitindex = {eind,[1:length(st);st;ed]',eind};

% collect the different types of adducts
addCount=0;
pos=[];
pos=regexp(adduct.name,',');
if (isempty(pos))
    for j=1:adduct.count
        addCount=addCount+1;
        adductType(addCount).name=adduct.name;
        adductType(addCount).count=j;
    end
else
    for j=1:adduct.count
        addCount=addCount+1;
        adductType(addCount).name=adduct.name(1:pos-1);
        adductType(addCount).count=j;
        addCount=addCount+1;
        adductType(addCount).name=adduct.name(pos+1:end);
        adductType(addCount).count=j;
    end
end

% Do glycan fragmentation
glyBion=[];                         % stores glycan B-ions that result from release of terminal or internal ions
glyYion=[];
%if ((~isempty(glyMat))&&(nFrag>0)) % fragmentation only necessary for peptides with glycan modifications
Yion=[]; Bion=[]; Iion=[];
finalYion=[]; finalBion=[]; finalIion=[];
bonds = strfind(glyStruct.name,'{');  % find location of glycosidic bonds
allind = Prod(1).unitindex;
for j=1:nFrag
    if (length(bonds)>=nFrag)         % only necessary if number of glycosidic bonds exceed nFrag
        modComb=combnk(bonds,j);
        for k=1:size(modComb,1)
            [gYion,gBion,gIion,indYion,indBion,indIion]=findGlyFrag(Prod(1).name,modComb(k,:),allind);
            for zz=1:z
                for m=1:length(adductType)
                    if (zz==sum(adductType(m).count))
                        % First, handle Y-ions
                        Yion.name=char(gYion);
                        Yion.z=zz;
                        Yion.adduct=adductType(m);
                        Yion.nFrag=j;
                        Yion.ano=ano;
                        if strcmpi(form,'OH')
                            Yion.form='Y';
                        elseif strcmpi(form,'Me')
                            Yion.form='MY';
                        else strcmpi(form,'Ac')
                            Yion.form='AY';
                        end
                        Yion.ion = ion;
                        Yion.mz=(glyMZCalc(Yion,GlycanNew,AnomericNew,AdductNew));            % This is m/z of Y-ion
                        % Second, B-ions
                        concat1='';
                        concat2='';
                        for n=1:length(gBion)
                            Bion.name=char(gBion(n));
                            Bion.z=zz;
                            Bion.adduct=adductType(m);
                            Bion.nFrag=j;
                            Bion.ano=ano;
                            if strcmpi(form,'OH')
                                Bion.form='B';
                            elseif strcmpi(form,'Me')
                                Bion.form='MB';
                            else strcmpi(form,'Ac')
                                Bion.form='AB';
                            end
                            Bion.ion = ion;
                            Bion.mz=(glyMZCalc(Bion,GlycanNew,AnomericNew,AdductNew));            % This is m/z of Y-ion
                            concat1=[concat1,'-b ',Bion.name];
                            Bion.type=['-b ',Bion.name];
                            Bion.unitindex = {eind,indBion{n},eind};
                        end
                        % I-ion
                        for n=1:length(gIion)
                            Iion.name=char(gIion(n));
                            Iion.z=zz;
                            Iion.adduct=adductType(m);
                            Iion.nFrag=j;
                            Iion.ano=ano;
                            if strcmpi(form,'OH')
                                Iion.form='I';
                            elseif strcmpi(form,'Me')
                                Iion.form='MI';
                            else strcmpi(form,'Ac')
                                Iion.form='AI';
                            end
                            Iion.ion = ion;
                            Iion.mz=(glyMZCalc(Iion,GlycanNew,AnomericNew,AdductNew));            % This is m/z of Y-ion
                            concat2=[concat2,'-i ',Iion.name];
                            Iion.type=['-i ',Iion.name];
                            Iion.unitindex = {eind,indIion{n},eind};
                        end
                        Yion.type=['-y',concat1,concat2];
                        Yion.unitindex = {eind,indYion,eind};
                        finalYion=[finalYion,Yion];
                        finalBion=[finalBion,Bion];
                        finalIion=[finalIion,Iion];
                    end
                end
            end
        end
    end
end
masslist=0;
for i=length(finalYion):-1:1
    if any(masslist==finalYion(i).mz)
        finalYion(i)=[];
    else
        masslist=[masslist,finalYion(i).mz];
    end
end
masslist=0;
for i=length(finalBion):-1:1
    if any(masslist==finalBion(i).mz)
        finalBion(i)=[];
    else
        masslist=[masslist,finalBion(i).mz];
    end
end
masslist=0;
for i=length(finalIion):-1:1
    if any(masslist==finalIion(i).mz)
        finalIion(i)=[];
    else
        masslist=[masslist,finalIion(i).mz];
    end
end

AllIons=[Prod,finalYion,finalBion,finalIion];
end