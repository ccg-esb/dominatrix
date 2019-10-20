function params=defineParams


    %%%%%%%%%%%%%%%
    
    params.nodes=[1:2];
    params.name_nodes={'WT','MUT'};

    
    %%%%%%%%%%
    params.pcn=20; %Mean plasmid copy number
    params.rho=.05; %Plasmid replication rate
   
    params.iniPlasmids = [params.pcn; 0];  %Initial plasmid copy number
    
    params.mut_rate=1e-7; %Mutation rate

    