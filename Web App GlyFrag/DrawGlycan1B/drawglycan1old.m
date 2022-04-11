function output = drawglycan1old(allsgp,varargin)
% DRAWGLYCAN: draw glycan or glycopeptide. This is the main program
%
% Syntax:
% drawglycan(allsgp,options)
% drawglycan(allsgp,optionname1,optionvalue1,...)
%
% Input:
% allsgp: string, glycan and glycopeptide string in IUPAC format.
% options: struct, user defined options listed below in Table format.
% optionname and optionvalue: user defined option name and value pairs.
%
% Output:
% N/A
%
% Available options (case sensitive):
% ________________________________________________________________________________________
% |  option name              purpose                   option value                     |
% |--------------------------------------------------------------------------------------|
% |  orientation           Orientation of the        String:'up','down','left','right'.  |
% |                        glycan structure          Default is 'right'.                 |
% |  monosacsize           Monosac. size             Number in range: [0,1].             |
%                                                    Default=0.5.                        |
% |  msperiwidth           Monosac. border           Any integer. Default=1.             |
% |                        thickness                                                     |
% |  bondwidth             Molecular bond            Any integer. Default=2.             |
% |                        thickness                                                     |
% |  perpendicularmonosac  List of perpendicular     1 x n cell array of strings.        |
% |                        monosacs.                 Default= {'Fuc','Xyl'}.             |
% |  showlink              Diplay linkage info.      String, 'yes' or 'no'.              |
% |                                                  Default = 'yes'.                    |
% |  fontsize              Linkage data font size.   Any integer. Default = 12.          |
% |  linkinfodist          Linkage text position.    Any number. Default = 0.7.          |
% |  linkinfotheta         Text orientation with     Any number in range = [-90,90].     |
% |                        respect to linkage.       Default = 30 (degrees).             |
% |  bondbreaksiglength    Length of glycosidic      Any number. Default = 0.5.          |
% |                        bond fragmentation line.                                      |
% |  structspacing         Spacing between glycans   Any number. Default = 1.            |
% |                        listed in cell array.                                         |
% |  aaspacing             Spacing between amino     Any number. Default = 0.75.         |
% |                        acids.                                                        |
% |  fileout               File saving specs.        String, output figure file          |
% |                                                  name. Default = empty.              |
% |  visible               Display figure.           String,'on' or 'off',Default='on'.  |
% ----------------------------------------------------------------------------------------
% Example:
% drawglycan('A[GalNAc(a1-3)](-N "B2" -C  "X4")AB[Xyl(b1-3)GalNAc(a1-3)]B(-C "X2.5")C[Fuc(a1-3)]CDD')
%
% Children function:
% USG, CALCGLYPOS, ESTIFIGSIZE, PAINTGLYCAN, GETFRAGINFO, REARRANGEGLYPEP,
% PLOTPEPBREAK
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%
% v1.1 update note: draw one structure at a time, the function of drawing
% multiple structure from a cell input has been cancelled.
defoptname = {'orientation','monosacsize','msperiwidth',...
    'bondwidth','perpendicularmonosac',...
    'linkinfodist','linkinfotheta','bondbreaksiglength',...
    'showlink','fontsize','workingmode',...
    'structspacing','aaspacing','fileout',...
    'visible','inputformat','displaybrackettext',...
    'figurehandle','fragmentation','curlybracketitem',...
    'posshift','zoomfactor','intfactor','prevlocation',...
    'useratio','specialoptions','linkfontsize','sortbranches'};

defoptvalue = {'left',0.5,1,...
    2,{'Fuc','Xyl','Deoxyhexose'},...
    .7,30,.5,...
    'no',8,'',...
    1,.6,'',...
    'on','IUPAC','no',...
    [],[],[],...
    [0,0],[1,1],0,[0,0,0,0],...
    1,{},16,true}; % zoomfactor must be 2xmonosacsize

defopt = usg('gen',defoptname,defoptvalue);
if ~isempty(varargin)
    if isstruct(varargin{1})
        options = usg('mod',defopt,varargin{1});
    elseif ischar(varargin{1})
        options = usg('mod',defopt,varargin(1:2:end),varargin(2:2:end));
    end
else
    options = defopt;
end

if ismember(upper(options.orientation),{'UP','RIGHT'})
    options.linkinfotheta = abs(options.linkinfotheta);
elseif ismember(upper(options.orientation),{'DOWN','LEFT'})
    options.linkinfotheta = -abs(options.linkinfotheta);
end
[pepMat,glyMat,modMat]=breakGlyPep(allsgp);  % denoting modifications should be added, either as glycan or in pep seq.
glyseq = cell(length(glyMat),1);
glypos = zeros(length(glyMat),1);
for i = 1:length(glyMat)
    glyseq{i} = glyMat(i).struct;
    glypos(i) = glyMat(i).pos;
end
pepseq = pepMat.pep;
if isempty(pepseq) && ~isempty(glyMat)
    options.workingmode = 'G';
elseif ~isempty(pepseq) && isempty(glyMat)
    options.workingmode = 'P';
else
    options.workingmode = 'GP';
end
if isempty(options.figurehandle)
    drawglycan = figure;
    set(drawglycan,'visible',options.visible);
    mainaxes = axes(drawglycan);
    options.figurehandle = mainaxes;
end
hold(options.figurehandle,'on');
switch options.workingmode
    case 'G'
        if strcmpi(options.inputformat,'SGP2')
            glysgp = allsgp;
            %% need a module to recognize sgp based special options
            specialoptions = [];
        elseif strcmpi(options.inputformat,'SGP1')
            glysgp = sgp1to2(allsgp,1);
            specialoptions = [];
        elseif strcmpi(options.inputformat,'IUPAC')
            glyseq = strtrim(allsgp);
            [cleanseq,specialoptions] = getglydressing(glyseq);
            glysgp = IUPAC2Sgp(cleanseq,2);
        end
        if ~isempty(specialoptions)  % get option position in sgp seq, because it needs conversion
            spoptfld = fieldnames(specialoptions);  % special option fields
            numMono = length(strfind(glysgp,'{'));
            for i = 1:length(spoptfld)
                tempspopt = specialoptions.(spoptfld{i});
                for j = 1:size(tempspopt,1)
                    tempspopt{j,2} = bsxfun(@minus,numMono + 1,tempspopt{j,2});  % convert option associated monosac location from IUPAC to SGP
                end
                specialoptions.(spoptfld{i}) = tempspopt;
            end
        end
        [plotinfoout,specialoptions] = calcglypos(glysgp,options,specialoptions);  % only draw 1 glycan structure, for multiple structures, call func. multiple times
        % specialoptions are merged into plotinfoout
        if ~isempty(specialoptions)
            plotinfoout = calcattachmentpos(plotinfoout,options);
        end
        if ~isempty(options.fragmentation)
            plotinfoout(1).bonddescription = options.fragmentation{2}{1};
        end
        plotinfoout.mspos(:,1) = (plotinfoout.mspos(:,1).*options.zoomfactor(1)) + (options.posshift(1,1));
        plotinfoout.mspos(:,2) = (plotinfoout.mspos(:,2).*options.zoomfactor(2)) + ((options.posshift(1,2).*(options.useratio))+(options.intfactor.*2));
        selectprevlocation = options.prevlocation(options.prevlocation(:,3) > min(plotinfoout.mspos(:,1)-(options.monosacsize*0.5)),:);
        if ~isempty(selectprevlocation)
            go_on = 1;
            while go_on == 1
                stop = 1;
                for pp = 1:size(selectprevlocation,1)
                    if (((min(plotinfoout.mspos(:,2))-(options.monosacsize*0.5)) < selectprevlocation(pp,4)) && ...
                        ((max(plotinfoout.mspos(:,2))+(options.monosacsize*0.5)) > selectprevlocation(pp,2)))
                            plotinfoout.mspos(:,2) =  plotinfoout.mspos(:,2) + (selectprevlocation(pp,4) - (min(plotinfoout.mspos(:,2))-(options.monosacsize*0.5))) + options.monosacsize;
                        stop = 0;
                    end
                end
                if stop == 1
                    go_on = 0;
                end
            end
        end
        glytcp = paintglycan1(plotinfoout,options);
        atmtcp = paintattachment(plotinfoout,options,zeros(2));
        output.tcp_ori = [min(glytcp(1,1),atmtcp(1,1)),min(glytcp(1,2),atmtcp(1,2));...
            max(glytcp(2,1),atmtcp(2,1)),max(glytcp(2,2),atmtcp(2,2))];
        output.tcp = [(min(plotinfoout.mspos(:,1))-(options.monosacsize.*0.5)),(max(plotinfoout.mspos(:,1))+(options.monosacsize.*0.5));...
                      (min(plotinfoout.mspos(:,2))-(options.monosacsize.*0.5)),(max(plotinfoout.mspos(:,2))+(options.monosacsize.*0.5))];
        output.plotinfo = plotinfoout;
%         axis equal
    case 'GP'  % default format: IUPAC
        options.orientation = 'up';
        glysgp = cellfun(@strtrim,glyseq,'uniformoutput',false);
        if strcmpi(options.inputformat,'SGP1')
            for i = 1:length(glyseq)
                glysgp{i} = sgp1to2(glyseq{i});
            end
            specialoptions = cell(size(glysgp));
        elseif strcmpi(options.inputformat,'SGP2')
            glysgp = glyseq;
            specialoptions = cell(size(glysgp));
        end
        currentlen = length(glysgp);
        for i = 1:length(modMat)
            glysgp{currentlen + i} = ['{(??-?)',modMat(i).struct(2),'}'];
            glypos(currentlen + i) = modMat(i).pos;
            tempspopt.CHAR = {[],1};
            tempspopt.BOLD = {[],1};
            specialoptions{currentlen + i} = tempspopt;
        end
        [~,ind] = sort(glypos);
        glysgp = glysgp(ind);
        glypos = glypos(ind);
        specialoptions = specialoptions(ind);
        plotinfoout = cell(size(glysgp));
        for i = 1:length(glysgp)
            if ~isempty(specialoptions{i})  % get option position in sgp seq, because it needs conversion
                spoptfld = fieldnames(specialoptions{i});  % special option fields
                numMono = length(strfind(glysgp{i},'{'));
                for j = 1:length(spoptfld)
                    tempspopt = specialoptions{i}.(spoptfld{j});
                    for k = 1:size(tempspopt,1)
                        tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                    end
                    specialoptions{i}.(spoptfld{j}) = tempspopt;
                end
            end
            plotinfoout{i} = calcglypos(glysgp{i},options,specialoptions{i});  % only draw 1 glycan structure, for multiple structures, call func. multiple times
            % specialoptions are merged into plotinfoout
            if ~isempty(specialoptions{i})
                plotinfoout{i} = calcattachmentpos(plotinfoout{i},options);
            end
        end
        [pepbreak,cleanpepseq] = getseqoptions(pepseq,'p');
        [plotinfoout,peppos] = rearrangeglypep(cleanpepseq,glypos,plotinfoout,options);
        if ~isempty(options.fragmentation)
%         if ~isempty(options.fragmentation) && ~isempty(options.fragmentation{2})
            customfraginfo = options.fragmentation{2};
            for i = 1:length(customfraginfo)
            plotinfoout{i}.bonddescription = customfraginfo{i};
            end
        end
        for i = 1:length(cleanpepseq)
%             text(options.figurehandle,peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
%                 'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
            text(peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        tcp = zeros(2,2);
        tcp(1,1) = floor(peppos(1));
        tcp(1,2) = -2;
        tcp(2,1) = ceil(peppos(end));
        for i = 1:length(plotinfoout)
            if ~isempty(plotinfoout{i})
                plotinfoout.mspos(:,1) = (plotinfoout.mspos(:,1).*options.zoomfactor(1)) + (options.posshift(1,1));
                plotinfoout.mspos(:,2) = (plotinfoout.mspos(:,2).*options.zoomfactor(2)) + (options.posshift(1,2)+(options.intfactor.*2));
                glytcp = paintglycan1(plotinfoout{i},options);
                atmtcp = paintattachment(plotinfoout{i},options,zeros(2));
            end
        end
%         axis equal
        if ~isempty(options.fragmentation)
            pepbreak = options.fragmentation{1};
        end
        pbtcp = plotpepbreak(peppos,pepbreak,glypos,options);  % this is the last step because frag marker length is related to figure size
        
        tcp = [min([glytcp(1,1),atmtcp(1,1),tcp(1,1),pbtcp(1,1)]),min([glytcp(1,2),atmtcp(1,2),tcp(1,2),pbtcp(1,2)]);...
            max([glytcp(2,1),atmtcp(2,1),tcp(2,1),pbtcp(2,1)]),max([glytcp(2,2),atmtcp(2,2),tcp(2,2),pbtcp(2,2)])];
        output.tcp = tcp;
        output.plotinfo = plotinfoout;
        
%         estimatefigsize(plotinfoout,peppos,options)
    case 'P'
        options.orientation = 'up';
        glysgp = cell(size(modMat));
        glypos = zeros(size(modMat));
        specialoptions = cell(size(modMat));
        for i = 1:length(modMat)
            glysgp{i} = ['{(??-?)',modMat(i).struct(2),'}'];
            glypos(i) = modMat(i).pos;
            tempspopt.CHAR = {[],1};
            tempspopt.BOLD = {[],1};
            specialoptions{i} = tempspopt;
        end
        plotinfoout = cell(size(glysgp));
        for i = 1:length(glysgp)
            if ~isempty(specialoptions{i})  % get option position in sgp seq, because it needs conversion
                spoptfld = fieldnames(specialoptions{i});  % special option fields
                numMono = length(strfind(glysgp{i},'{'));
                for j = 1:length(spoptfld)
                    tempspopt = specialoptions{i}.(spoptfld{j});
                    for k = 1:size(tempspopt,1)
                        tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                    end
                    specialoptions{i}.(spoptfld{j}) = tempspopt;
                end
            end
            plotinfoout{i} = calcglypos(glysgp{i},options,specialoptions{i});  % only draw 1 glycan structure, for multiple structures, call func. multiple times
            % specialoptions are merged into plotinfoout
            if ~isempty(specialoptions{i})
                plotinfoout{i} = calcattachmentpos(plotinfoout{i},options);
            end
        end
        [pepbreak,cleanpepseq] = getseqoptions(pepseq,'p');
        [plotinfoout,peppos] = rearrangeglypep(cleanpepseq,glypos,plotinfoout,options);
        if ~isempty(options.fragmentation)
            %         if ~isempty(options.fragmentation) && ~isempty(options.fragmentation{2})
            customfraginfo = options.fragmentation{2};
            for i = 1:length(customfraginfo)
                plotinfoout{i}.bonddescription = customfraginfo{i};
            end
        end
        for i = 1:length(cleanpepseq)
            %             text(options.figurehandle,peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
            %                 'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
            text(peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        tcp = zeros(2,2);
        tcp(1,1) = floor(peppos(1));
        tcp(1,2) = -2;
        tcp(2,1) = ceil(peppos(end));
        for i = 1:length(plotinfoout)
            if ~isempty(plotinfoout{i})
                plotinfoout.mspos(:,1) = (plotinfoout.mspos(:,1).*options.zoomfactor(1)) + (options.posshift(1,1));
                plotinfoout.mspos(:,2) = (plotinfoout.mspos(:,2).*options.zoomfactor(2)) + (options.posshift(1,2)+(options.intfactor.*2));
                glytcp = paintglycan1(plotinfoout{i},options);
                atmtcp = paintattachment(plotinfoout{i},options,zeros(2));
            end
        end
%         axis equal
        if ~isempty(options.fragmentation)
            pepbreak = options.fragmentation{1};
        end
        pbtcp = plotpepbreak(peppos,pepbreak,glypos,options);  % this is the last step because frag marker length is related to figure size
        
        tcp = [min([glytcp(1,1),atmtcp(1,1),tcp(1,1),pbtcp(1,1)]),min([glytcp(1,2),atmtcp(1,2),tcp(1,2),pbtcp(1,2)]);...
            max([glytcp(2,1),atmtcp(2,1),tcp(2,1),pbtcp(2,1)]),max([glytcp(2,2),atmtcp(2,2),tcp(2,2),pbtcp(2,2)])];
        output.tcp = tcp;
        output.plotinfo = plotinfoout;
end

end