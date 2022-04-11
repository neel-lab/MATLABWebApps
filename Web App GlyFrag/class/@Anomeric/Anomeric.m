classdef Anomeric
    % Anomeric: An object storing the properties of the anomeric entity
    %
    %   Syntax:
    %       anoMass = Anomeric.mass(anoText)
    %       monoFormula = Glycan.anoformulaMap(mono1let)
    %
    %   Input:
    %       anoText: Text describing the anomeric entity
    %
    %   Output;
    %       anoMass:  Anomeric entitity mass
    %       anoFormula: Anomeric end formula
    %
    %   Examples:
    %       Mass    = Anomeric.mass('free')
    %       Formula = Anomeric.anoformulaMap('Bn')
    %
    %   Note:                 Anomeric end table:
    %                          __________________
    %                            code | type
    %                          __________________
    %                            free | free anomeric end
    %                             red | Reduced
    %                              Me | Methylated
    %                              Ac | Acetylated
    %                           MeRed | perMethylated Reduced
    %                              Bn | Benzyl gp
    %                             SBn | S-Benzyl gp.
    %                            ONAP | ONAP
    %                            SNAP | SNAP
    %                              AB | 2AB (2-aminobenzamide)
    %                          ABperM | permethylated 2AB
    %                            N2AB | 2-amino-N,N-dimethylbenzamide
    %                        N2ABperM | perMethylated 2-amino-N,N-dimethylbenzamide
    %                               P | Procainamide
    %                           PperM | permethylated Procainamide
    %                             pNP | para-Nitrophenol
    %                           Allyl | Allyl group
    %                          _________________
    %
    %See also Glycan, Aminoacid.
    
    %Author: Sriram Neelamegham
    %Date Lastly Updated: 4/23/2021
    
    properties (Constant)
        free      = struct('C',0,'H',1,'O',1,'N',0,'S',0,'P',0);
        ol        = struct('C',0,'H',2,'O',2,'N',0,'S',0,'P',0);
        Me        = struct('C',1,'H',3,'O',1,'N',0,'S',0,'P',0);
        Ac        = struct('C',2,'H',3,'O',2,'N',0,'S',0,'P',0);
        perMol    = struct('C',1,'H',3,'O',2,'N',0,'S',0,'P',0);
        Bn        = struct('C',7,'H',7,'O',1,'N',0,'S',0,'P',0);
        SBn       = struct('C',7,'H',7,'O',0,'N',0,'S',1,'P',0);
        ONAP      = struct('C',11,'H',9,'O',1,'N',0,'S',0,'P',0);
        SNAP      = struct('C',11,'H',9,'O',0,'N',0,'S',1,'P',0);
        AA       = struct('C',7,'H',8,'O',2,'N',1,'S',0,'P',0);
        AB       = struct('C',7,'H',9,'O',1,'N',2,'S',0,'P',0);
        ABperM   = struct('C',11,'H',17,'O',1,'N',2,'S',0,'P',0);
        N2AB = struct('C',9,'H',13,'O',1,'N',2,'S',0,'P',0);
        N2ABperM = struct('C',11,'H',17,'O',1,'N',2,'S',0,'P',0);
        P = struct('C',13,'H',22,'O',1,'N',3,'S',0,'P',0);
        PperM = struct('C',17,'H',31,'O',1,'N',3,'S',0,'P',0);
        pNP = struct('C',6,'H',5,'O',3,'N',1,'S',0,'P',0);
        allyl = struct('C',3,'H',5,'O',0,'N',0,'S',0,'P',0);
        anofullname = {'free','red','Me','Ac','MeRed','Bn','SBn','ONAP','SNAP','AA','AB','ABperM','N2AB','N2ABperM','P','PperM','pNP','allyl'};
        ano1let     = {'-','ol','M','A','pol','B','T','O','S','N','2','2perM','N2','N2perM','P','PperM','pNP','allyl'};  % not used in program-- for future development
        anoformula =  {struct('H',1,'O',1)...
            struct('C',1,'H',3,'O',1)...
            struct('C',1,'H',3,'O',1)...
            struct('C',2,'H',3,'O',2)...
            struct('C',1,'H',3,'O',2)...
            struct('C',7,'H',7,'O',1)...
            struct('C',7,'H',7,'S',1)...
            struct('C',11,'H',9,'O',1)...
            struct('C',11,'H',9,'S',1)...
            struct('C',7,'H',8,'O',2,'N',1)...
            struct('C',7,'H',9,'O',1,'N',2)...                % This is 2AB + 1H
            struct('C',11,'H',17,'O',1,'N',2)...              % VdS: C10H14O1N2? Reply: C11H17O1N2 is correct. Because ring is opened, we need to consider the methylation of that OH. 
            struct('C',9,'H',13,'O',1,'N',2,'S',0,'P',0)...
            struct('C',11,'H',17,'O',1,'N',2,'S',0,'P',0)...
            struct('C',13,'H',22,'O',1,'N',3,'S',0,'P',0)...  %
            struct('C',17,'H',31,'O',1,'N',3,'S',0,'P',0)...
            struct('C',6,'H',5,'O',3,'N',1,'S',0,'P',0)...
            struct('C',3,'H',5,'O',0,'N',0,'S',0,'P',0)}; 
        
        anoformulaMap = containers.Map({'free','red','Me','Ac','MeRed','Bn','SBn','ONAP','SNAP','2AA','2AB','2ABperM','N2AB','N2ABperM','P','PperM','pNP','allyl'},...
            {struct('H',1,'O',1)...
            struct('H',1,'O',1)...
            struct('C',1,'H',3,'O',1)...
            struct('C',2,'H',3,'O',2)...
            struct('C',1,'H',3,'O',2)...
            struct('C',7,'H',7,'O',1)...
            struct('C',7,'H',7,'S',1)...
            struct('C',11,'H',9,'O',1)...
            struct('C',11,'H',9,'S',1)...
            struct('C',7,'H',8,'O',2,'N',1)...
            struct('C',7,'H',9,'O',1,'N',2)...
            struct('C',11,'H',17,'O',1,'N',2)...
            struct('C',9,'H',13,'O',1,'N',2,'S',0,'P',0)...
            struct('C',11,'H',17,'O',1,'N',2,'S',0,'P',0)...
            struct('C',13,'H',22,'O',1,'N',3,'S',0,'P',0)...
            struct('C',17,'H',31,'O',1,'N',3,'S',0,'P',0)...
            struct('C',6,'H',5,'O',3,'N',1,'S',0,'P',0)...
            struct('C',3,'H',5,'O',0,'N',0,'S',0,'P',0)});  
        
        
        mass = containers.Map({'free','red','Me','Ac','MeRed','Bn','SBn',...
            'ONAP','SNAP','2AA','2AB','2ABperM','N2AB',...
            'N2ABperM','P','PperM','pNP','allyl'},...
            {17.0027387,19.0186187,31.0183879,59.013302,47.0608479,107.0496863,123.0268449,...
            157.0653355,173.0424941,138.0554988,137.0714831,193.1340799,165.10279,...
            193.1340799,236.17629,292.23889,139.0269391,41.0391230}); 
        charge = containers.Map({'free','red','Me','Ac','MeRed','Bn','SBn',...
            'ONAP','SNAP','2AA','2AB','2ABperM','N2AB',...
            'N2ABperM','P','PperM','pNP','allyl'},...
            {0,0,0,0,0,0,0,...
            0,0,0,0,0,0,...
            0,0,1,0,0});
        
    end
end