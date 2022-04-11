function [core, elements]=bracket(varargin)
% Analyze a set of glycans within glycanList cell element, and suggest
% core glycan structure in sgp format and also elements to place outside brackets
%
% Uses other functions: findGlyStr, multiSGPFrag, removeCellCol, commonGly, eraseGly

if nargin==0        % read glycan list, with options to specify nFrag and charge also if necessary
    gly{1}='{n{n{h{h{h}{h}}{h{h}{n{h}}}}}}';
    gly{2}='{n{n{h{h{h}{h}}{n}{h{h}{h}}}}}';
    z=1;
elseif nargin==1
    gly=varargin{1};
    z=1;
elseif nargin == 2
    SGP = varargin{1};
    select = varargin{2};
    gly = cell(length(select),1);
    for ii = 1:length(select)
        gly{ii}=SGP{select(ii)};
    end
    z = 1;
end
core={};
ngFrag=1;
while isempty(core)
    ngFrag=ngFrag+1;
    nMono=findGlyStr(gly,'{');       % find maxium number of monosaccharides among all the glycans
    AllFrag_B=cell(length(gly),nMono*ngFrag);  % Use this to initialize the AllFrag_B and AllFrag_Y list
    AllFrag_Y=cell(length(gly),nMono*ngFrag);
    
    for i=1:length(gly)                 % determine the list of B- and Y-ions for each of the glycans
        [AllFrag]=multiSGPFrag(gly{i},0,0,ngFrag,z);
        tempB={};
        tempY={};
        for j=1:length(AllFrag)
            if strfind(AllFrag(j).type,'-b')
                tempB=[tempB,AllFrag(j).sgp];
            end
            if strfind(AllFrag(j).type,'-y')
                tempY=[tempY,AllFrag(j).sgp];
            end
        end
        for j=1:length(tempB)
            AllFrag_B(i,j)=tempB(j);
        end
        for j=1:length(tempY)
            AllFrag_Y(i,j)=tempY(j);
        end
    end
    AllFrag_B=removeCellCol({AllFrag_B});
    AllFrag_Y=removeCellCol({AllFrag_Y});
    
    % find largest common Y-ion and call it 'core'
    list=commonGly(AllFrag_Y);
    maxGly=0;
    for i=1:length(list)
        glyLength=length(strfind(list{i},'{'));
        if maxGly<glyLength
            maxGly=glyLength;
            core=list{i};
        end
    end
end
% remove core from the orignial glycan for each element
subgly={};
for i=1:length(gly)
    remainder=eraseGly(gly{i},core,ngFrag);
    subgly=[subgly;remainder];
end

% Find out what bits and pieces remain that are common. These are the elements
elements={};
for i=1:size(subgly,2)
    master=subgly(1,i);
    boo=1;
    for j=2:size(subgly,1)
        if ~cellfun(@(x) strcmp(master,x),subgly(j,:))
            boo=0;
        end
    end
    if boo==1
        elements=[elements,master];
    else
        comp=glyComp(master);
        elements=[elements,comp];
    end
end
end

function [comp]=glyComp(gly)
list=['n','h','s','f','u'];   % change composition if there are additional monosaccharides
comp='';
for i=1:length(list)
    ct=length(strfind(gly,list(i)));
    comp=[comp,list(i),num2str(ct)];
end
end