function output = drawglycan1(allsgp,varargin)
% DRAWGLYCAN-SSA version (with figure coordinate output)
defoptname = {'orientation','monosacsize','msperiwidth',...
    'bondwidth','perpendicularmonosac',...
    'linkinfodist','linkinfotheta','bondbreaksiglength',...
    'showlink','fontsize','workingmode',...
    'structspacing','aaspacing','fileout',...
    'visible','inputformat','displaybrackettext',...
    'specialoptions','linkfontsize','sortbranches',...
    'posshift','zoomfactor','intfactor','prevlocation',...
    'useratio','figurehandle'};

defoptvalue = {'left',.5,1,...
    2,{'Fuc','Xyl','Deoxyhexose','Pentose'},...
    .7,30,.5,...
    'yes',8,'',...
    1,.75,'',...
    'on','iupac','no',...
    {},16,true,...
    [0,0],[1,1],0,[0,0,0,0],...
    1,[]};

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

if ismember(lower(options.orientation),{'right','down'})
    options.linkinfotheta = -options.linkinfotheta;
end
options.isaxes = true;

[glyseq,pepseq,glypos,PTMid] = distggp(allsgp,options.inputformat);  % glyseq is cell, pepseq is char
if isempty(options.workingmode)
    if isempty(pepseq)
        options.workingmode = 'G';
    else
        if ~isempty(glyseq)
            options.workingmode = 'GP';
        else
            options.workingmode = 'P';
        end
    end
end
%% handle different input format, produce a "glysgp" and "specialoptions"
glyoptions = [DrawGlycanPara.intglymodinfo,DrawGlycanPara.intglybondinfo,...
    DrawGlycanPara.extglymodinfo,DrawGlycanPara.glyidentityinfo,'CURLYBRACKET'];

if strcmpi(options.inputformat,'IUPAC')
    switch options.workingmode
        case 'G'
            [cleanseq,specialoptions] = getglydressing(strtrim(glyseq));
            glysgp = IUPAC2Sgp(cleanseq,[],2);
        case 'GP'
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            [cleanseq,specialoptions] = cellfun(@getglydressing,glyseq,'uniformoutput',false);
            glysgp = cellfun(@IUPAC2Sgp,cleanseq,num2cell(repmat(2,numel(glyseq),1)),'uniformoutput',false);
        case 'P'
            specialoptions = options.specialoptions;
    end
    if ~isempty(specialoptions)
        specialoptions = specialoptions{1};
    end
    if ~isempty(specialoptions)  % get option position in sgp seq, because it needs conversion
        if iscell(specialoptions)
            for i = 1:length(specialoptions)
                if ~isempty(specialoptions)
                    spoptfld = fieldnames(specialoptions{i});  % special option fields
                    numMono = length(strfind(glysgp{i},'{'));
                    for j = 1:length(spoptfld)
                        if ismember(upper(spoptfld{j}),glyoptions)
                            tempspopt = specialoptions{i}.(spoptfld{j});
                            for k = 1:size(tempspopt,1)
                                tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                            end
                            specialoptions{i}.(spoptfld{j}) = tempspopt;
                        end
                    end
                end
            end
        elseif isstruct(specialoptions)
            spoptfld = fieldnames(specialoptions);  % special option fields
            numMono = length(strfind(glysgp,'{'));
            for j = 1:length(spoptfld)
                if ismember(upper(spoptfld),glyoptions)
                    tempspopt = specialoptions.(spoptfld{j});
                    for k = 1:size(tempspopt,1)
                        tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                    end
                    specialoptions.(spoptfld{j}) = tempspopt;
                end
            end
        end
    end
elseif strcmpi(options.inputformat,'GLYCAM')
    if ~any(strfind(allsgp,'('))
        allsgp = glycam2iupac(allsgp);  % not tested
    end
    options.inputformat = 'IUPAC';
    [cleanseq,specialoptions] = getglydressing(strtrim(allsgp));
    glysgp = IUPAC2Sgp(cleanseq,2);
    if ~isempty(specialoptions)
        specialoptions = specialoptions{1};
    end
elseif strcmpi(options.inputformat,'SGP1')
    switch options.workingmode
        case 'G'
            glysgp = sgp1to2(allsgp,PTMid);
            specialoptions = options.specialoptions;
        case 'GP'
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            glysgp = cellfun(@sgp1to2a,glyseq,num2cell(PTMid),'UniformOutput',false);
            specialoptions = options.specialoptions;
            if isempty(specialoptions)
                specialoptions = cell(size(glysgp));
            end
        case 'P'
            specialoptions = options.specialoptions;
    end
    if ~isempty(specialoptions)
        specialoptions = specialoptions{1};
    end
elseif strcmpi(options.inputformat,'SGP1to2')
    switch options.workingmode
        case 'G'
            glysgp = sgp1to2a(allsgp);
            specialoptions = options.specialoptions;
        case 'GP'
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            glysgp = cellfun(@sgp1to2a,glyseq,'UniformOutput',false);
            specialoptions = options.specialoptions;
            if isempty(specialoptions)
                specialoptions = cell(size(glysgp));
            end
        case 'P'
            specialoptions = options.specialoptions;
    end
    if ~isempty(specialoptions)
        specialoptions = specialoptions{1};
    end
elseif strcmpi(options.inputformat,'SGP2')
    switch options.workingmode
        case 'G'
%             glysgp = getglydressing(strtrim(glyseq));
            glysgp = allsgp;
            specialoptions = options.specialoptions;
        case 'GP'
            glysgp = cellfun(@strtrim,glyseq,'UniformOutput',false);
            specialoptions = options.specialoptions;
            if isempty(specialoptions)
                specialoptions = cell(size(glysgp));
            end
    end
    if ~isempty(specialoptions)
        specialoptions = specialoptions{1};
    end
elseif strcmpi(options.inputformat,'LINUCS')
    
elseif strcmpi(options.inputformat,'BCSDB')
    
elseif strcmpi(options.inputformat,'LINEARCODE')
    
elseif strcmpi(options.inputformat,'GLYCOCT')
    %% sequence itself is a connection table, reformat and put it in specialoptions.
elseif strcmpi(options.inputformat,'WURCS')
    
elseif strcmpi(options.inputformat,'KCF')
    
elseif strcmpi(options.inputformat,'CABOSML')
    
elseif strcmpi(options.inputformat,'GLYDE')
    
else
    return
end
switch upper(options.workingmode)
    case 'G'
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
        if isempty(options.figurehandle)
            parentfigure = figure;
            options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1]);
            hold on;
            set(options.figurehandle,'visible',options.visible)
            options.isaxes = false;
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
        fig = options.figurehandle;
        glytcp = paintglycan(plotinfoout,options,[]);
        atmtcp = paintattachment(plotinfoout,options,zeros(2));
        estimatefigsize(plotinfoout,'',options);
        if all(all(atmtcp==0))
            output.tcp = [min(glytcp(1,1)),max(glytcp(2,1));...
                min(glytcp(1,2)),max(glytcp(2,2))];
        else
            output.tcp = [min(glytcp(1,1),atmtcp(1,1)),max(glytcp(2,1),atmtcp(2,1));...
                min(glytcp(1,2),atmtcp(1,2)),max(glytcp(2,2),atmtcp(2,2))];
        end
    case 'GP'  % default format: IUPAC
        options.orientation = 'up';
        plotinfoout = cell(size(glysgp));
        for i = 1:length(glysgp)
            [plotinfoout{i},specialoptions{i}] = calcglypos(glysgp{i},options,specialoptions{i});  % only draw 1 glycan structure, for multiple structures, call func. multiple times
            % specialoptions are merged into plotinfoout
            if ~isempty(specialoptions)
                plotinfoout{i} = calcattachmentpos(plotinfoout{i},options);
            end
        end
        [pepoptions,cleanpepseq] = getseqoptions(pepseq,'p');
        if ~isempty(specialoptions) && ismember('pepC',fieldnames(specialoptions))
            if isempty(pepoptions) || ~isfield(pepoptions,'C')
                pepoptions.C = specialoptions.pepC;
            else
                pepoptions.C = [pepoptions.C;specialoptions.pepC];
            end
        end
        if ~isempty(specialoptions) && ismember('pepN',fieldnames(specialoptions))
            if isempty(pepoptions) || ~isfield(pepoptions,'N')
                pepoptions.N = specialoptions.pepN;
            else
                pepoptions.N = [pepoptions.N;specialoptions.pepN];
            end
        end
        
        [plotinfoout,peppos] = rearrangeglypep(cleanpepseq,glypos,plotinfoout,options);
        if isempty(options.figurehandle)
            parentfigure = figure;
            options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1]);
            hold on;
            set(options.figurehandle,'visible',options.visible,'Units','normalized')
            options.isaxes = false;
        end
        fig = options.figurehandle;
        for i = 1:length(cleanpepseq)
            text(fig,peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        tcp = zeros(2,2);
        tcp(1,1) = peppos(1) - options.aaspacing;
        tcp(1,2) = -2;
        tcp(2,1) = peppos(end) + options.aaspacing;
        glytcp = tcp;
        atmtcp = tcp;
        for i = 1:length(plotinfoout)
            if ~isempty(plotinfoout{i})
                glytcp = paintglycan(plotinfoout{i},options,glytcp);
                atmtcp = paintattachment(plotinfoout{i},options,atmtcp);
            end
        end
    case 'P'
        [pepoptions,cleanpepseq] = getseqoptions(pepseq,'p');
        peppos = options.aaspacing * [0:length(cleanpepseq)-1];
        tcp = zeros(2);
        tcp(1,1) = peppos(1) - options.aaspacing;
        tcp(1,2) = -2;
        tcp(2,1) = peppos(end) + options.aaspacing;
        glytcp = zeros(2);
        atmtcp = zeros(2);
        if isempty(options.figurehandle)
            parentfigure = figure;
            options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1]);
            hold on;
            set(options.figurehandle,'visible',options.visible,'Units','normalized')
            options.isaxes = false;
        end
        fig = options.figurehandle;
        if ~isempty(specialoptions) && ismember('pepC',fieldnames(specialoptions))
            pepoptions.C = specialoptions.pepC;
        end
        if  ~isempty(specialoptions) && ismember('pepN',fieldnames(specialoptions))
            pepoptions.N = specialoptions.pepN;
        end
        for i = 1:length(cleanpepseq)
            text(fig,peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        glytcp = zeros(2);
        atmtcp = zeros(2);
        tcp = zeros(2);
end
if ~options.isaxes
    set(fig,'position',[0,0,1,1],'visible','off')
end
% set(fig,'visible','off')
axis equal
if strcmpi(options.workingmode,'gp') || strcmpi(options.workingmode,'p')
    pbtcp = plotpepbreak(peppos,pepoptions,glypos,options);  % this is the last step because frag marker length is related to figure size
    if isfield(pepoptions,'NGMOD')
        ngmods = pepoptions.NGMOD;
        for i = 1:size(ngmods,1)
            text(fig,peppos(ngmods{i,2}),-2,ngmods{i,1},'FontName','Helvetica','VerticalAlignment','bottom',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
    end
    tcp = [min([glytcp(1,1),atmtcp(1,1),tcp(1,1),pbtcp(1,1)]),min([glytcp(1,2),atmtcp(1,2),tcp(1,2),pbtcp(1,2)]);...
        max([glytcp(2,1),atmtcp(2,1),tcp(2,1),pbtcp(2,1)]),max([glytcp(2,2),atmtcp(2,2),tcp(2,2),pbtcp(2,2)])];
    output.tcp = tcp;
    if strcmpi(options.workingmode,'gp')
        set(fig,'XLim',[tcp(1,1),tcp(2,1)],'YLim',[tcp(1,2),tcp(2,2)+options.monosacsize]);
    elseif strcmpi(options.workingmode,'p')
        set(fig,'XLim',[peppos(1) - options.aaspacing,peppos(i) + options.aaspacing],'YLim',[-2.5,1]);
    end
end

if ~isempty(options.fileout)
    set(gcf, 'InvertHardcopy', 'off')
    saveas(gcf,options.fileout)
end
if strcmpi(options.visible,'off')
    p=get(options.figurehandle,'Parent');
    close(p);
end
end