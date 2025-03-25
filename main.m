clear all
clc

SearchAgents_no=30; % Number of search agents
Max_iteration=500; % Maximum numbef of iterations

figure('Position',[100 200 500 350])
Function_name=sprintf('%s%d', 'F', 1);
% Load details of the selected benchmark function
fobj = @Chung_Reynolds;
lb = -100;
ub = 100;
dim = 30;
[Best_score,Best_pos,ESHO_cg_curve]=ESHO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

%Draw objective space
semilogy(ESHO_cg_curve,'Color','r','LineWidth',1.5)

title(Function_name)
xlabel('Iteration count');
ylabel('Fitness value');
axis tight 
grid off  
box on 



