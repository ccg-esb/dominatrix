function [t, A]=simulateReplication(params, A0, T)

    if nargin<2
        T=1;
    end
    
    numNodes = size(A0,1);

    numSteps = 1000;
    A = zeros(numNodes,numSteps); 
    t = zeros(1,numSteps);

    % set the initial values
    m = 1;
    A(:,1) = A0;
    t(1) = 0;

    % Gillespie SSA
    
    while (t(m) < T)
      %******* (a)

      r = rand(1,2); % 1-by-2 vector of two random numbers
      %******* (b)

      [w, M]=getPropensities(params, A(:,m));
      w0 = sum(w);
      w0_1 = w0*(1- sum(A(:,m))/params.pcn);
      
      %******* (c)
      if ~isinf(w0_1) && w0_1>0
        tau = r(1) / w0_1; %uniform distribution
      else
        tau=Inf; 
      end

      %******(d)
      if r(2)<w(1) 
          i=1;
      else
          i=2;
      end

      % i-th reaction occurs
      %fprintf('\n%d: A(m)=%d t(m)=%f i=%d\n',m,A(1,m),t(m),i);
      %disp(['Plasmid ', params.name_nodes{i},' Replicates ']);
      Omega=1;  %Plasmid copy number increment 
      A(:,m+1) = A(:,m) + (1/Omega)*M(i,:)';
      t(m+1) = t(m) + params.rho*tau;
      %fprintf('-> t(m+1)=%f A(m+1)=%d  w=[%g,%g]  tau=%g\n',t(m+1),A(1,m+1),w(1),w(2), tau);

      m = m + 1;

      % memory management
      if (m >= numSteps)
        %disp(['Doubling memory for data. (Reallocation.)  t(',num2str(m),')=',num2str(t(m))]);
        %tic
        Aaux = zeros(numNodes,2*numSteps);
        taux = zeros(1,2*numSteps);
        Aaux(:,1:numSteps) = A;
        taux(1:numSteps) = t;
        A = Aaux;
        t = taux;
        clear Aaux taux;
        numSteps = 2*numSteps;
        %fprintf('  done. [%fsec]\n', toc);
      end
      
    end % while

    % cutting the zeros at the end of arrays
    A = A(:,1:m);
    t = t(1:m);

    % postprocessing
    if (t(m) > T)
      t(m) = T;
    end
    for j=1:numNodes
      if (A(j,m) < 0)
        A(j,m) = 0;
      end
    end
    A=A';



