function drawbondmap(bondmap)
h=figure;

distance = zeros(size(bondmap,1),1);
firstbranch = 0;
for i = 1:size(bondmap,1)
    distance(logical(bondmap(i,:))) = distance(i) + 1;
    if firstbranch == 0
        if sum(bondmap(i,:)) > 1
            firstbranch = i;
        end
    end
end
mspos = drawrawtree(bondmap,distance,[]);
ax = axes(h);
hold on;
options.figurehandle = ax;
options.monosacsize = 0.5;
options.bondwidth = 1;
mspos = branchequalizer(mspos,bondmap);
fix = mspos(firstbranch,2);
mspos(firstbranch:end,2) = mspos(firstbranch:end,2) - fix;
tmp = mspos(:,1);
mspos(:,1) = -mspos(:,2);
mspos(:,2) = tmp;
mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
mspos(:,1) = -mspos(:,1);
paintlink(mspos,bondmap,'',options)
end