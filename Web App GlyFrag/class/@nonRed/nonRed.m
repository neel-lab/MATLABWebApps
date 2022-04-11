classdef nonRed
    % Anomeric: An object storing the properties of the non-reducing end
    %    
    %   Syntax:
    %       nonRedMass = nonRed.mass(nonRedText) 
    %       nonRedFormula = nonRed.formula(nonRedText) 
    %
    %   Input:
    %       nonRedText: Text describing the nonRed entity
    %       
    %   Output;
    %       nonRedMass:  nonRed entitity mass
    %       nonRedFormula: nonRed formula
    %    
    %   Examples:
    %       Mass    = nonRed.mass('free')
    %       Formula = nonRed.formula('free')
    %                      
    %   Note:                    nonRed end table:
    %                          __________________
    %                            code | type
    %                          __________________
    %                           free  | free anomeric end
    %                              Me | Methylated
    %                              Ac | Acetylated
    %                          _________________
    %
    %See also Glycan, Anomeric, Aminoacid.
    
    %Author: Sriram Neelamegham
    %Date Lastly Updated: 1/5/2017
    
      properties (Constant)
             free      = struct('C',0,'H',1,'O',0,'N',0,'S',0,'P',0);
             Me        = struct('C',1,'H',3,'O',0,'N',0,'S',0,'P',0);
             Ac        = struct('C',2,'H',3,'O',1,'N',0,'S',0,'P',0);
            nonRedfullname = {'free','Me','Ac'}
            nonRed1let     = {'-','M','A'};
            nonRedformula = {struct('H',1)...
                           struct('C',1,'H',3)...
                           struct('C',2,'H',3,'O',1)};      
            formula = containers.Map({'free','Me','Ac'},...
                           {struct('H',1)...
                           struct('C',1,'H',3)...
                           struct('C',2,'H',3,'O',1)}); 
           mass = containers.Map({'free','Me','Ac'},...
           {1.0078246,15.0234738,43.0183879}); 
      end
end