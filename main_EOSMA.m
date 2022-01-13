clear;clc
close all
tic

Algorithm_name='EOSMA';
optimum=[100, 1100, 700, 1900, 1700, 1600, 2100, 2200, 2400, 2500];

runs=30;

for num=1:10
    
    for i=1:runs
        
        lb=-100;  ub=100;  dim=20;
        SearchAgents_no=100; % Number of search agents
        Max_Nfes=1000000;    % Maximum number of function evaluations
        [Destination_fitness,bestPositions,Convergence_curve]=EOSMA(SearchAgents_no,Max_Nfes,lb,ub,dim,num);

        fitness(i,:)=Destination_fitness;
        Positions(i,:)=bestPositions;
        curve(i,:)=Convergence_curve;

    end
    fitness=fitness-optimum(num);
    [~,Index]=sort(fitness);
    bestPosition=Positions(Index(1),:);
    curve=curve-optimum(num);
    curve(curve<1e-8)=0;
    Convergence_curve=sum(curve,1)/runs;
    
    Convergence_curve=Convergence_curve(1,1:10:end);
    
    disp([Algorithm_name,':function ',num2str(num)]);
    disp(['Mean Fitness ',num2str(Convergence_curve(end))]);

end

Time=toc;

%Developer: Shihong Yin