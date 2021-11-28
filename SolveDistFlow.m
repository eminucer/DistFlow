
function [y] = SolveDistFlow(Grid, PlotGraph, MaxIter, tol)
% This function solves the distribution power flow equations for a given
% radial distribution network based on DistFlow equations introduced in the Baran&Wu paper(1989).
%
% Inputs
%   ->Grid      : a structure contaning grid parameters
%   ->PlotGraph : Binary variable to display informative figures
%   ->MaxIter   : Maximum number of iterations
%   ->tol       : error tolerance (0.0001 recommended)
% 
% Outputs
%   ->y         : output contaning the solved line powers, voltages and phase angles.

function [x_next] = F(x_prev,rk,xk,pk,qk)
    P_prev = x_prev(1);
    Q_prev = x_prev(2);
    V_prev = x_prev(3);

    P_next = P_prev-pk-rk*(P_prev^2 + Q_prev^2)/V_prev^2;
    Q_next = Q_prev-qk-xk*(P_prev^2 + Q_prev^2)/V_prev^2;
    V_next = V_prev^2-2*rk*P_prev-2*xk*Q_prev+(rk^2+xk^2)*(P_prev^2 + Q_prev^2)/V_prev^2;

    x_next(1) = P_next;
    x_next(2) = Q_next;
    x_next(3) = sqrt(V_next);
end

function J = Jacobian(x_prev,rk,xk)
P_prev = x_prev(1);
Q_prev = x_prev(2);
V_prev = x_prev(3);

zk = rk^2 + xk^2;

J = [1-2*rk*P_prev/(V_prev^2) -2*rk*Q_prev/(V_prev^2) rk*(P_prev^2 + Q_prev^2)/(V_prev^4); ...
     -2*xk*P_prev/(V_prev^2) 1-2*xk*Q_prev/(V_prev^2) xk*(P_prev^2 + Q_prev^2)/(V_prev^4); ...
     -2*(rk-zk*P_prev/(V_prev^2)) -2*(xk-zk*Q_prev/(V_prev^2)) 1-zk*(P_prev^2 + Q_prev^2)/(V_prev^4)];
end

    N = Grid.N;
    L = Grid.L;
    pN = Grid.pN;
    qN = Grid.qN;
    rk = Grid.rk;
    xk = Grid.xk;
    Vo = Grid.Vo;
    
    Phase = cell(N+1,1);
    Phase{1} = 0;

    Ainv = -L;
    A = inv(Ainv);

    FromNode = [];
    ToNode = [];
    FromNode(1) = 0;
    ToNode(1) = 1;

    for i=1:N-1
        a = find(A(:,i) == 1);
        FromNode = [FromNode repmat(i,1,length(a))];
        ToNode = [ToNode a'];
    end

    FromNode = FromNode +1;
    ToNode = ToNode + 1;
    GR = digraph(FromNode,ToNode);
    eLabels = cell(1,N);
    nLabels = cell(1,N+1);
    for i=1:N
        eLabels{i} = strcat('l_{',num2str(FromNode(i)-1),'-',num2str(ToNode(i)-1),'}');
    end
    for i=1:N+1
        nLabels{i} = strcat(num2str(i-1));
    end

    if PlotGraph
        figure('Renderer', 'painters', 'Position', [250 250 900 350])  
        p=plot(GR,'NodeLabel',nLabels,'LineW',2);
%        p.NodeFontSize = 14;
%        p.EdgeFontSize = 14;
        p.Marker = 's';
        p.NodeColor = 'r';
    end

    FinalNodes = find(sum(abs(A),1) == 1); % Nodes with no further branches/children (terminal nodes)
    Paths = cell(length(FinalNodes),1);
    for i=1:length(FinalNodes)
        Paths{i} = [0, find(Ainv(FinalNodes(i),:) == -1)];
    end

    NumOfPaths = length(Paths);
    len = zeros(NumOfPaths,2);

    for i=1:NumOfPaths
        len(i,:) = [length(Paths{i}),i];
    end

    [~,idx] = sort(len,'descend');
    
    pL = zeros(N,1);
    qL = zeros(N,1);

    LinDistFlow_solution_V = Vo^2*ones(N,1) - 2*Ainv*diag(rk)*(Ainv)'*pN - 2*Ainv*diag(xk)*(Ainv)'*qN;
    LinDistFlow_solution_V = [Vo^2 ; LinDistFlow_solution_V];
    LinDistFlow_solution_V = sqrt(LinDistFlow_solution_V);

    P_terminal = zeros(NumOfPaths,1);

    for iter = 1:MaxIter
    UpdatedVoltages = nan(N+1,1);    
    UpdatedVoltages(1) = Vo;
    
    counter = 0;
    
        for i= idx(:,1)' 
            counter = counter + 1;
            
            UpdatedNodes = Paths{i}(find(~isnan(UpdatedVoltages(Paths{i}(2:end)+1)))+1);
            if isempty(UpdatedNodes)
                LastNode = Paths{i}(1);
            else
                LastNode = UpdatedNodes(end);
            end
            RemainingNodes = Paths{i}(find(Paths{i} == LastNode)+1:end);
            % If LastNode == RemainingNodes => We already solved this branch

            NN = length(RemainingNodes);

            FirstNode = RemainingNodes(1);

            B = abs(Ainv);
            B = B + -diag(ones(N,1));

            % Initialize
            ConnectedNodes = cell(length(Paths{i})-1,1);
            AggregatedNodes = cell(length(Paths{i})-1,1);

            for j=1:length(Paths{i})-1
                ConnectedNodes{j} = find(B(:,Paths{i}(j+1))==1);
            end

            for j=1:length(Paths{i})-1
                if Paths{i}(j+1) == Paths{i}(end)
                    AggregatedNodes{j} = [];
                else
                    AggregatedNodes{j} = setdiff(ConnectedNodes{j},[ConnectedNodes{j+1}; Paths{i}(j+2)]);
                end
                AggregatedNodes{j,2} = Paths{i}(j+1);
            end

            StartFrom = find(cell2mat(AggregatedNodes(:,2)) == FirstNode);
            %We assume FirstNode's voltage is known/updated.
            pNN = zeros(NN,1);
            qNN = zeros(NN,1);
            
            NumOfNodesToSolve = StartFrom:StartFrom+NN-1;
            for j=1:length(NumOfNodesToSolve)
                ContributingNodes = AggregatedNodes{NumOfNodesToSolve(j),1};
                if isempty(ContributingNodes)
                        pNN(j) =  pN(AggregatedNodes{NumOfNodesToSolve(j),2});
                        qNN(j) =  qN(AggregatedNodes{NumOfNodesToSolve(j),2});
                else    
                    if pL(ContributingNodes(1)) == 0 %if not updated (include qL as well)
                        for k=1:length(ContributingNodes)
                            pNN(j) = pNN(j) + pN(ContributingNodes(k));
                            qNN(j) = qNN(j) + qN(ContributingNodes(k));
                        end
                        %add itself
                        pNN(j) = pNN(j) + pN(AggregatedNodes{NumOfNodesToSolve(j),2});
                        qNN(j) = qNN(j) + qN(AggregatedNodes{NumOfNodesToSolve(j),2});
                    else
                        pNN(j) = pNN(j) + pL(ContributingNodes(1)) + pN(AggregatedNodes{NumOfNodesToSolve(j),2});
                        qNN(j) = qNN(j) + qL(ContributingNodes(1)) + qN(AggregatedNodes{NumOfNodesToSolve(j),2});
                    end
                end
            end
                       
            xx = cell(NN+1,1); % State variable (N+1)x1
            J = cell(NN,1);   % Jacobian matrix NxN

            if pL(FirstNode) == 0
                xx{1} = [sum(pNN); ...
                         sum(qNN); ...
                         UpdatedVoltages(LastNode+1)]; % LastNode = FirstNode in un-updated path
            else
                xx{1} = [pL(FirstNode);
                         qL(FirstNode);
                         UpdatedVoltages(LastNode+1)];
            end

            % Forward sweep to find PLs and QLs
            rkk = rk(RemainingNodes);
            xkk = xk(RemainingNodes);

            for ii=1:NN
                xx{ii+1} = (F(xx{ii},rkk(ii),xkk(ii),pNN(ii),qNN(ii)))';
                pL(RemainingNodes(ii)) = xx{ii}(1);
                qL(RemainingNodes(ii)) = xx{ii}(2);
                UpdatedVoltages(RemainingNodes(ii)+1) = xx{ii+1}(3);
                J{ii} = Jacobian(xx{ii},rkk(ii),xkk(ii));
            end
            P_terminal(counter) = xx{ii+1}(1);

            JJ = J{NN}(1:2,1:3); % J(N) -> 2x3 matrix

            for ii=NN-1:-1:1
               JJ = JJ * J{ii};
            end

            %x_next = xx{1}-JJ\xx{NN+1}; 
            x_next = xx{1}(1:3)-JJ\xx{NN+1}(1:2); 
            
            pL(FirstNode) = x_next(1);
            qL(FirstNode) = x_next(2);

        end   
        
        error = norm(abs(P_terminal),1);
        
        fprintf('Iteration: %d\nTerminal node error: %f', iter, error)
        if error > tol
            fprintf(' > tolerance = %f -> N\n', tol);
        else
            fprintf('< tolerance = %f -> Y\n\nTotal Num of Iterations: %d\n', tol, iter);
            if 1
                fprintf('\nNode voltages and phase angles:\n');
                for i=1:N+1
                    if i<= N
                        Phase{i+1} = Phase{i}-angle(UpdatedVoltages(i)^2 - (conj(rk(i)+1j*xk(i)))*(pL(i)+1j*qL(i))); 
                    end
                    fprintf('V[%d] = %.4f /_ %.2f\n',i-1, (abs(UpdatedVoltages(i))/Vo), rad2deg(Phase{i}));
                end
            end
           break;
        end
        fprintf('\n');

    end

    figure('Renderer', 'painters', 'Position', [250 250 900 350]) 
    hold on; 
    grid on;
    plot(0:N,UpdatedVoltages/Vo,'-o','LineW',2)
    plot(0:N,LinDistFlow_solution_V/Vo, 'LineW',2)
    xticks([0:N])
    xlabel('Nodes')
    ylabel('Voltage')
    set(gca,'FontSize',13)
    legend('DistFlow','LinDistFlow')
    xlim([0,N]);
    
    y = cell(4,1);
    y{1} = pL;
    y{2} = qL;
    y{3} = UpdatedVoltages;
    y{4} = cell2mat(Phase);
    y{5} = P_terminal;

end
