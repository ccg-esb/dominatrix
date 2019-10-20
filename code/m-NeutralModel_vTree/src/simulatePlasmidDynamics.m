function [time, pG1, pR12]=simulatePlasmidDynamics(params)


pG1=[];
pR12=[];
time=[];

%Initial conditions
pG1_0=params.iniPlasmids(1);
pR12_0=params.iniPlasmids(2);
for n=1:params.numDivisions
    
    %Plasmid replication events
    params.iniPlasmids = [pG1_0; pR12_0];
    [rep_time, rep_y]=simulateReplication(params, params.iniPlasmids, params.T);
    
    %Division event (randomly segregate plasmids)
    p1=rep_y(end,1);
    p2=rep_y(end,2);
    
    pG1_0=binornd(floor(p1),0.5,1);
    pR12_0=binornd(floor(p2),0.5,1);
    
    %disp(['Division event: (',num2str(p1),', ',num2str(p2),') -> (',num2str(pG1_0),', ',num2str(pR12_0),'): freq=',num2str(freq)]);
    
    
    %Save variables
    pG1=[pG1 interp1(rep_time, rep_y(:,1), rep_time(2:end))];
    pR12=[pR12 interp1(rep_time, rep_y(:,2), rep_time(2:end))];
    
    time=[time (n-1)*params.T+rep_time(2:end)];
end