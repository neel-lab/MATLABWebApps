function tcp = paintattachment(plotinfoout,options,prevtcp)
%% adduct, anomer curlybracket,repeat
mspos = plotinfoout.mspos;
bondmap = plotinfoout.bondmap;
tr = DrawGlycanPara.tipradius;
fig = options.figurehandle;
tcp = prevtcp;

%% 1.cbitem
%% analysis of curly bracet add-on structure
if isfield(plotinfoout,'CBITEM')
    cbitem = plotinfoout.CBITEM;
    waypoints = cbitem.bracketwaypoints;
    plot(fig,waypoints(2:3,1),waypoints(2:3,2),'k')
    plot(fig,waypoints(6:7,1),waypoints(6:7,2),'k')
    mainmspos = plotinfoout.mspos;
    switch lower(options.orientation)
        case 'up'
            drawcurve(fig,pi/2,pi,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(fig,pi*3/2,pi*2,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(fig,pi,pi*3/2,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(fig,0,pi/2,waypoints(8,1),waypoints(8,2),tr)  % IV
        case 'down'
            drawcurve(fig,pi,pi*3/2,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(fig,0,pi/2,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(fig,pi/2,pi,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(fig,pi*3/2,pi*2,waypoints(8,1),waypoints(8,2),tr)  % IV
        case 'left'
            drawcurve(fig,pi,pi*3/2,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(fig,0,pi/2,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(fig,pi*3/2,pi*2,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(fig,pi/2,pi,waypoints(8,1),waypoints(8,2),tr)  % IV
        case 'right'
            drawcurve(fig,pi*3/2,pi*2,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(fig,pi/2,pi,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(fig,pi,pi*3/2,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(fig,0,pi/2,waypoints(8,1),waypoints(8,2),tr)  % IV
    end
    cboption = options;
    cboption.workingmode = 'g';
    cb_tcp = paintglycan(cbitem.cbplotinfosto,cboption,tcp);
    if all(all(cb_tcp==0))
        tcp = [min(tcp(1,1)),min(tcp(1,2));...
            max(tcp(2,1)),max(tcp(2,2))];  % update 2 corner pos.
    else
        tcp = [min([cb_tcp(:,1);tcp(1,1)]),min([cb_tcp(:,2);tcp(1,2)]);...
            max([cb_tcp(:,1);tcp(2,1)]),max([cb_tcp(:,2);tcp(2,2)])];  % update 2 corner pos.

    end
end

%% 2. repeat
if isfield(plotinfoout,'RS')
    repeatstart = plotinfoout.RS;
    repeatend = plotinfoout.RE;
    if strcmpi(options.orientation,'UP')
        textpos = 'end';
    elseif strcmpi(options.orientation,'DOWN')
        textpos = 'start';
        inter = repeatend;
        repeatend = repeatstart;
        repeatstart = inter;
    elseif strcmpi(options.orientation,'LEFT')
        textpos = 'start';
    elseif strcmpi(options.orientation,'RIGHT')
        textpos = 'end';
        inter = repeatend;
        repeatend = repeatstart;
        repeatstart = inter;
    end
    for i = 1:size(repeatend,1)
        thismonosacpos = mspos(repeatend{i,2},:);
        prevmonosacind = find(bondmap(:,repeatend{i,2}),1,'first');
        if ~isempty(prevmonosacind)
            prevmonosacpos = mspos(prevmonosacind,:);
        else  % only at the beginning of the structure, it's the nature of tree structure
            switch lower(options.orientation)
                case 'up'
                    prevmonosacpos = thismonosacpos - [0,1];
                case 'down'
                    prevmonosacpos = thismonosacpos + [0,1];
                case 'left'
                    prevmonosacpos = thismonosacpos + [1,0];
                case 'right'
                    prevmonosacpos = thismonosacpos - [1,0];
            end
        end
        if thismonosacpos(1)-prevmonosacpos(1) ~= 0
            alpha = atan((thismonosacpos(2)-prevmonosacpos(2))/(thismonosacpos(1)-prevmonosacpos(1)));
        else
            alpha = pi/2;
        end
        xdif = thismonosacpos(1)-prevmonosacpos(1);
        ydif = thismonosacpos(2)-prevmonosacpos(2);
        msdistance = sqrt(xdif^2+ydif^2);
        if xdif > 0 && ydif <= 0
            alpha = alpha + 2*pi;
        elseif xdif <= 0 && ydif < 0
            alpha = alpha + pi;
        elseif xdif < 0 && ydif >= 0
            alpha = alpha + pi;
        end  % fix alpha value
        %         centx = prevmonosacpos(1)*.3 + thismonosacpos(1)*.7;
        %         centy = prevmonosacpos(2)*.3 + thismonosacpos(2)*.7;
        centx = prevmonosacpos(1)*.45 + thismonosacpos(1)*.55;
        centy = prevmonosacpos(2)*.45 + thismonosacpos(2)*.55;
        linestart = [centx + options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy + options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        lineend = [centx - options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy - options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        plot(fig,[linestart(1),lineend(1)],[linestart(2),lineend(2)],'Color',[.5 .5 .5],'linewidth',options.bondwidth*.75);
        bracketlinefix = [options.bondbreaksiglength/ 4 * cos( 0 /180*pi+alpha),options.bondbreaksiglength/ 4 * sin( 0 /180*pi+alpha)];
        bracketlineend1 = linestart + bracketlinefix;
        bracketlineend2 = lineend + bracketlinefix;
        plot(fig,[bracketlineend1(1),linestart(1)],[bracketlineend1(2),linestart(2)],'Color',[.5 .5 .5],'linewidth',options.bondwidth*.75);
        plot(fig,[bracketlineend2(1),lineend(1)],[bracketlineend2(2),lineend(2)],'Color',[.5 .5 .5],'linewidth',options.bondwidth*.75);
        if ~isempty(repeatend{i,1})
            switch textpos
                case 'start'
                    brackettextpos = [linestart(1) + (linestart(1) - thismonosacpos(1))* 0.4,linestart(2) + (linestart(2) - thismonosacpos(2))* 0.4];
                    text(fig,brackettextpos(1),brackettextpos(2),repeatend{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
                case 'end'
                    brackettextpos = [lineend(1) + (lineend(1) - thismonosacpos(1))* 0.4,lineend(2) + (lineend(2) - thismonosacpos(2))* 0.4];
                    text(fig,brackettextpos(1),brackettextpos(2),repeatend{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
            end
        end
    end
    
    
    for i = 1:size(repeatstart,1)
        thismonosacpos = mspos(repeatstart{i,2},:);
        nextmonosacind = find(bondmap(repeatstart{i,2},:),1,'first');
        if ~isempty(nextmonosacind)
            nextmonosacpos = mspos(nextmonosacind,:);
        else  % end of antenna
            prevmonosacind = find(bondmap(:,repeatstart{i,2}),1,'first');
            if ~isempty(prevmonosacind)
                prevmonosacpos = mspos(prevmonosacind,:);
            else  % only at the beginning of the structure, it's the nature of tree structure
                switch lower(options.orientation)
                    case 'up'
                        prevmonosacpos = thismonosacpos - [0,1];
                    case 'down'
                        prevmonosacpos = thismonosacpos + [0,1];
                    case 'left'
                        prevmonosacpos = thismonosacpos + [1,0];
                    case 'right'
                        prevmonosacpos = thismonosacpos - [1,0];
                end
            end
            nextmonosacpos = 2*thismonosacpos - prevmonosacpos;
        end
        if thismonosacpos(1)-nextmonosacpos(1) ~= 0
            alpha = atan((nextmonosacpos(2)-thismonosacpos(2))/(nextmonosacpos(1)-thismonosacpos(1)));
        else
            alpha = pi/2;
        end
        xdif = nextmonosacpos(1) - thismonosacpos(1);
        ydif = nextmonosacpos(2) - thismonosacpos(2);
        msdistance = sqrt(xdif^2+ydif^2);
        if xdif > 0 && ydif <= 0
            alpha = alpha + 2*pi;
        elseif xdif <= 0 && ydif < 0
            alpha = alpha + pi;
        elseif xdif < 0 && ydif >= 0
            alpha = alpha + pi;
        end  % fix alpha value
        centx = nextmonosacpos(1)*.45 + thismonosacpos(1)*.55;
        centy = nextmonosacpos(2)*.45 + thismonosacpos(2)*.55;
        linestart = [centx + options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy + options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        lineend = [centx - options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy - options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        plot(fig,[linestart(1),lineend(1)],[linestart(2),lineend(2)],'Color',[.5 .5 .5],'linewidth',options.bondwidth*.75);
        bracketlinefix = [options.bondbreaksiglength/ 4 * cos( 0 /180*pi+alpha),options.bondbreaksiglength/ 4 * sin( 0 /180*pi+alpha)];
        brackettextfix = [options.bondbreaksiglength/ 4 * cos( 90 /180*pi+alpha),options.bondbreaksiglength/ 4 * sin( 90 /180*pi+alpha)];
        bracketlineend1 = linestart - bracketlinefix;
        bracketlineend2 = lineend - bracketlinefix;
        plot(fig,[bracketlineend1(1),linestart(1)],[bracketlineend1(2),linestart(2)],'Color',[.5 .5 .5],'linewidth',options.bondwidth*.75);
        plot(fig,[bracketlineend2(1),lineend(1)],[bracketlineend2(2),lineend(2)],'Color',[.5 .5 .5],'linewidth',options.bondwidth*.75);
        if ~isempty(repeatstart{i,1})
            switch textpos
                case 'start'
                    brackettextpos = [linestart(1) + (linestart(1) - thismonosacpos(1))* 0.4,linestart(2) + (linestart(2) - thismonosacpos(2))* 0.4];
                    text(fig,brackettextpos(1),brackettextpos(2),repeatstart{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
                case 'end'
                    brackettextpos = [lineend(1) + (lineend(1) - thismonosacpos(1))* 0.4,lineend(2) + (lineend(2) - thismonosacpos(2))* 0.4];
                    %                     text(fig,lineend(1)-1.0*brackettextfix(1),lineend(2)-1.0*bracketlinefix(2),repeatend{i,1},'HorizontalAlignment','center');
                    text(fig,brackettextpos(1),brackettextpos(2),repeatstart{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
            end
        end
    end
    
end

%% 3.bracket
if isfield(plotinfoout,'bracketspos')
    bracketspos = plotinfoout.bracketspos;
    adduct = plotinfoout.ADDUCT;
    adducttextpos = plotinfoout.adducttextpos;
    plot(fig,bracketspos(1:4,1),bracketspos(1:4,2),'k');
    plot(fig,bracketspos(5:8,1),bracketspos(5:8,2),'k');
    adducttext = '';
    for i = 1:size(adduct,1)
        adducttext = [adducttext,adduct{i,1},' '];
    end
    adducttext = adducttext(1:end - 1);
    text(fig,adducttextpos(1),adducttextpos(2),adducttext,'fontsize',options.fontsize * 1.5);
    currentax = gca;
    set(currentax,'Units','pixels')
    width = currentax.Position(3);
    pixelperunit = width/(tcp(2,1) - tcp(1,1));
    tcp(1,1) = min(tcp(1,1),min(bracketspos(:,1)));
    tcp(2,1) = max(tcp(2,1),adducttextpos(1) + (length(adducttext)*options.fontsize * 1.5 * DrawGlycanPara.letterwidthpxl)*.41/pixelperunit);
    tcp(2,2) = max(tcp(2,2),adducttextpos(2));
    set(currentax,'Units','normalized')
end
end

function drawcurve(figurehandle,a,b,h,k,r)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians,
% b is end of arc in radians,
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Author:  Matt Fig
t = linspace(a,b);
x = r*cos(t) + h;
y = r*sin(t) + k;
plot(figurehandle,x,y,'k');
end