function [w, M]=getPropensities(params, y)

    G1=y(1);
    R12=y(2);

    %Total number of plasmids
    m=G1+R12;
    
    % Compute propensities of each reaction
    if m>params.pcn || m==0
        w(1)=0;
        w(2)=0;
    else
        w(1)=  G1/m; % (G1/m)*(1-G1/params.mu);
        w(2)=  R12/m; %(R12/m)*(1-R12/params.mu);        
    end
    
    
    %Matrix of stoichiometries (Nreactions x Nspecies).
    %Each row gives the stoichiometry of a reaction.
    %     G1  R12
    M = [  1   0 ;  % Plasmid G1 replicates
           0   1 ];  % Plasmid R12 replicates
    
    
    
    
    
    