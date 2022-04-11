function [glyMZ,comp] = glyMZCalc(varargin)
% GLYMW: Returns m/z for glycan 
%
% Syntax:
%      glyMZ = glyMZCalc(glyStruct)
%
%   Input:
%     glyStruct: A structure describing the glycan
%
%   Output:
%     glyMZ: monoisotopic glycan M/Z (including charge state consideration)
%
% Examples:
%   adduct.name='Na';
%   adduct.count=2;
%   glyStruct.name='{n{n{h{h{h{h}}}{h{h{h}}{h}}}}}';
%   glyStruct.z=2;
%   glyStruct.adduct=adduct;
%   glyStruct.ano='Me';
%   glyStruct.form='Me';
%   glyMZCalc(glyStruct)
%
% Date Last updated : 01/08/17  by Sriram Neelamegham

% Test data
% adduct(1).name='Na';
% adduct(1).count=2;
% glyStruct(1).name='{n{n{h{h{h{h}}}{h{h{h}}{h}}}}}';
% glyStruct(1).z=2;
% glyStruct(1).adduct=adduct;
% glyStruct(1).ano='Me';
% glyStruct(1).form='Me';

if (nargin>0)
    glyStruct=varargin{1};
end

if (nargin>1)
    GlycanNew=varargin{2};   % This part of the code is not written but it will yield glycan molecular composition when done
    AnomericNew=varargin{3};
    AdductNew=varargin{4};
else
    GlycanNew = Glycan;
    AnomericNew = Anomeric;
    AdductNew = Adduct;
end

format longg;
glyMW=0;

gly=glyStruct.name;
ion=glyStruct.ion;
z=glyStruct.z;
ano=glyStruct.ano;
adduct=glyStruct.adduct;
form=glyStruct.form;

switch ion
    case 'Positive'
        % mass of adduct
        pos=regexp(adduct.name,',', 'once');
        if isempty(pos)
            adductMass=AdductNew.mass(adduct.name)*adduct.count(1);
        else
            alladduct = strsplit(adduct.name,',');
            if mod(length(alladduct),2)
                error('Adduct information error')
            end
            adductMass = 0;
            for i = 1:length(alladduct)
                adductMass = adductMass + AdductNew.mass(alladduct{i})*adduct.count(i);
            end
        end
        % read monosaccharide mass from glycan class depending on the form of the monosaccharide
        if any(strcmpi(form,{'Me','MY','MB','MI'}))         % for all methylated ions
            typesofmono = GlycanNew.glycanMSMeMap.keys;
        elseif any(strcmpi(form,{'Ac','AY','AB','AI'}))     % for all acetylated ions
            typesofmono = GlycanNew.glycanMSAcetylMap.keys;
        elseif any(strcmpi(form,{'OH','Y','B','I'}))      % for non-derivatized ions
            typesofmono = GlycanNew.glycanMSMap.keys;
        end
        numtypesofmono = length(typesofmono);
        brac = regexp(gly,"[{}]");
        monoresidues = cell((length(brac)/2),1);
        count = 1;
        for ii = 1:(length(brac)-1)
            if brac(ii+1)-brac(ii) > 1
                monoresidues{count} = gly((brac(ii)+1):(brac(ii+1)-1));
                count = count + 1;
            end
        end
        monoresidues = monoresidues(~cellfun('isempty',monoresidues));
        for i = 1 : length(monoresidues)                  % add all internal masses of monosaccharides
            nummonoresidues = 1; %length(strfind(gly,typesofmono{i}));
            if any(strcmpi(form,{'Me','MY','MB','MI'}))
                glyMW = glyMW+nummonoresidues*GlycanNew.glycanMSMeMap(monoresidues{i});
            elseif any(strcmpi(form,{'Ac','AY','AB','AI'}))
                glyMW = glyMW+nummonoresidues*GlycanNew.glycanMSAcetylMap(monoresidues{i});
            else % non-derivatized monosaccharide
                glyMW = glyMW+nummonoresidues*GlycanNew.glycanMSMap(monoresidues{i});
            end
        end
        % set the reducing end mass
        if (strcmpi(ano,'P')&&any(strcmpi(form,{'Me','MY','MB','MI'})))
            ano='PperM';
        end
        redEnd=AnomericNew.mass(ano);
        % non-reducing end addition
        if any(strcmpi(form,{'Me','MY','Ac','AY','OH','Y'}))
            if any(strcmpi(form,{'MY','AY','Y','OH'}))  % Methyl at terminal
                TotEnd=redEnd+1.0078246;    % Add H at the non-reducing end
            elseif any(strcmpi(form,{'Me'}))
                TotEnd=redEnd+15.0234738;    % Add +CH3 at the non-reducing end
            elseif any(strcmpi(form,{'Ac'}))
                TotEnd=redEnd+43.0183879;     % ?+COCH3 at non-reducing end
            end
        end
        % non-reducing end
        if any(strcmpi(form,{'B','MB','AB'}))
            if any(strcmpi(form,{'MB'}))  % Methyl
                TotEnd=15.0234738;        % +CH3at the non-reducing end
            elseif any(strcmpi(form,{'AB'}))  % Acetyl
                TotEnd=43.0183879;             % ?+COCH3 at the non-reducing end
            else                            % plain B-ion
                TotEnd=1.0078246;           % +H at the non-reducing end
            end
            TotEnd=TotEnd-1.0078246;    % Add H at reducing end
        end
        if any(strcmpi(form,{'MI','AI','I'}))
            TotEnd=0;               % loss of H at red-end and gain of H at other end compensate
        end
        glyFull=glyMW+TotEnd;
        glyMZ=(glyFull+adductMass)/z;
        comp='';
        
    case 'Negative'
        z = -z;
        % mass of adduct
        AdductNew.charge('H') = -1;
        pos=regexp(adduct.name,',', 'once');
        if isempty(pos)
            adductMass=AdductNew.mass(adduct.name)*adduct.count(1)*AdductNew.charge(adduct.name);
        else
            alladduct = strsplit(adduct.name,',');
            if mod(length(alladduct),2)
                error('Adduct information error')
            end
            adductMass = 0;
            for i = 1:length(alladduct)
                adductMass = adductMass + (AdductNew.mass(alladduct{i})*adduct.count(i)*AdductNew.charge(alladduct{i}));
            end
        end
        % account for positive and negative adducts
        % read monosaccharide mass from glycan class depending on the form of the monosaccharide
        if any(strcmpi(form,{'Me','MY','MB','MI'}))         % for all methylated ions
            typesofmono = GlycanNew.glycanMSMeMap.keys;
        elseif any(strcmpi(form,{'Ac','AY','AB','AI'}))     % for all acetylated ions
            typesofmono = GlycanNew.glycanMSAcetylMap.keys;
        elseif any(strcmpi(form,{'OH','Y','B','I'}))      % for non-derivatized ions
            typesofmono = GlycanNew.glycanMSMap.keys;
        end
        numtypesofmono = length(typesofmono);
        brac = regexp(gly,"[{}]");
        monoresidues = cell((length(brac)/2),1);
        count = 1;
        for ii = 1:(length(brac)-1)
            if brac(ii+1)-brac(ii) > 1
                monoresidues{count} = gly((brac(ii)+1):(brac(ii+1)-1));
                count = count + 1;
            end
        end
        monoresidues = monoresidues(~cellfun('isempty',monoresidues));
        for i = 1 : length(monoresidues)                  % add all internal masses of monosaccharides
            nummonoresidues = 1; %length(strfind(gly,typesofmono{i}));
            if any(strcmpi(form,{'Me','MY','MB','MI'}))
                glyMW = glyMW+nummonoresidues*GlycanNew.glycanMSMeMap(monoresidues{i});
            elseif any(strcmpi(form,{'Ac','AY','AB','AI'}))
                glyMW = glyMW+nummonoresidues*GlycanNew.glycanMSAcetylMap(monoresidues{i});
            else % non-derivatized monosaccharide
                glyMW = glyMW+nummonoresidues*GlycanNew.glycanMSMap(monoresidues{i});
            end
        end
        % set the reducing end mass
        if (strcmpi(ano,'P')&&any(strcmpi(form,{'Me','MY','MB','MI'})))
            ano='PperM';
        end
        redEnd=AnomericNew.mass(ano);
        % non-reducing end addition
        if any(strcmpi(form,{'Me','MY','Ac','AY','OH','Y'}))
            if any(strcmpi(form,{'MY','AY','Y','OH'}))  % Methyl at terminal
                TotEnd=redEnd+1.0078246;    % Add H at the non-reducing end
            elseif any(strcmpi(form,{'Me'}))
                TotEnd=redEnd+15.0234738;    % Add +CH3 at the non-reducing end
            elseif any(strcmpi(form,{'Ac'}))
                TotEnd=redEnd+43.0183879;     % ?+COCH3 at non-reducing end
            end
        end
        % non-reducing end
        if any(strcmpi(form,{'B','MB','AB'}))
            if any(strcmpi(form,{'MB'}))  % Methyl
                TotEnd=15.0234738;        % +CH3at the non-reducing end
            elseif any(strcmpi(form,{'AB'}))  % Acetyl
                TotEnd=43.0183879;             % ?+COCH3 at the non-reducing end
            else                            % plain B-ion
                TotEnd=1.0078246;           % +H at the non-reducing end
            end
            TotEnd=TotEnd-1.0078246;    % Add H at reducing end
        end
        if any(strcmpi(form,{'MI','AI','I'}))
            TotEnd=0;               % loss of H at red-end and gain of H at other end compensate
        end
        % Calculate final mz
        glyFull=glyMW+TotEnd;
        glyMZ=(glyFull+adductMass)/(-z);
        comp='';
end
end

% % if any(strcmpi(form,{'Me','MY','Ac','AY','full','Y'}))
% %     if strcmpi(ano,'Bn')
% %         redEnd=redEnd+107.0496863;    % +C7H7O at reducing end (one H less than Benzyl alc.)
% %     elseif strcmpi(ano,'SNAP')
% %         redEnd=redEnd+173.0424941;    % C11H9O
% %     elseif strcmpi(ano,'ONAP')
% %         redEnd=redEnd+157.0653355;    % C11H9S
% %     elseif strcmpi(ano,'2AB')
% %         redEnd=redEnd+137.0714831;    % C7H9ON2    
% %     else                                      % ano='-'
% %         if any(strcmpi(form,{'Me','MY'}))  % Methyl at anomeric position
% %             redEnd=31.0183879;             % +OCH3 at anomer
% %         elseif any(strcmpi(form,{'Ac','AY'}))  % Acetyl at anomeric position
% %             redEnd=59.013302;             % ? +OCOCH3 at anomer
% %         else                               % it is either full or Y
% %             redEnd=17.0027387;         % +OH at reducing end
% %         end
% %     end
% %     if any(strcmpi(form,{'MY','AY','Y','full'}))  % Methyl at terminal
% %     TotEnd=redEnd+1.0078246;    % Add H at the non-reducing end
% %    elseif any(strcmpi(form,{'Me'}))
% %     TotEnd=redEnd+15.0234738;    % Add +CH3 at the non-reducing end
% %         else any(strcmpi(form,{'Ac'}))
% %     TotEnd=redEnd+43.0183879;     % ?+COCH3 at non-reducing end
% %         end
% % end


% if strcmpi(gly(1),'B')
%     ano='Bn';
% elseif strcmpi(gly(1),'S')
%     ano='SNAP';
% elseif strcmpi(gly(1),'O')
%     ano='ONAP';
% elseif strcmpi(gly(1),'2')
%     ano='2AB';
% else
%     ano='-';   % free end
% end
% 
% if (nargin>1)
%     form = varargin{1};
% elseif isempty(form)
%     form='Me';                                  % default is permetylated glycan mass
% end

