% main training
n_problems = 500;
n_params = 5;
todo = "training";
tot_time = 0.0;

range_beta = [1, 5, 10, 15, 20];
range_force = [1:0.1:10];

for i = 1:n_problems
    tstart = tic;
    fprintf('\n..........Problem number %d........\n', i);
    rng(i,'twister');
    v = ones([n_params,1]);
    v(2) = range_beta(randi(length(range_beta),1));
    v(5) = range_force(randi(length(range_force),1));
    higamod_call(v, todo);
    tstop = toc(tstart);
    tstop = datevec(tstop./(60*60*24));
    tot_time = tot_time + tstop;
    fprintf('\n..........time needed for problem number %d is %f seconds..........\n', i, tstop(6));    
end

for i = length(tot_time):-1:2
    tot_time(i-1) = floor(tot_time(i)/60);
    tot_time(i) = tot_time(i) - tot_time(i-1)*60;
end

fprintf('\n..........total time needed is %d hours %d minutes %f seconds..........\n', tot_time(4), tot_time(5), tot_time(6));

%% main testing

% n_problems = 10;
% n_params = 5;
% todo = "testing";
% tot_time = 0.0;
% 
% a = 1;
% b = 10;
% 
% for i = 1:n_problems
%     tstart = tic;
%     fprintf('\n..........Problem number %d........\n', i);
%     rng(i,'twister');
%     v = (b-a).*rand(n_params,1) + a;
%     higamod_call(v, todo);
%     tstop = toc(tstart);
%     tstop = datevec(tstop./(60*60*24));
%     tot_time = tot_time + tstop;
%     fprintf('\n..........time needed for problem number %d is %f seconds..........\n', i, tstop(6));    
% end
% 
% for i = length(tot_time):-1:2
%     tot_time(i-1) = floor(tot_time(i)/60);
%     tot_time(i) = tot_time(i) - tot_time(i-1)*60;
% end
% 
% fprintf('\n..........total time needed is %d hours %d minutes %f seconds..........\n', tot_time(4), tot_time(5), tot_time(6));

%% main testing for 1 problem
todo = "testing";
n_params = 5;

tot_time = 0.0;

fprintf('\n.......... Problem ........\n');
tstart = tic;
v = ones([n_params,1]);
v(2) = 18;
v(5) = 3.4;
[plotStruct, obj_solverIGA, numbVerNodes, cmin, cmax] = higamod_call(v, todo);
tstop = toc(tstart);
tstop = datevec(tstop./(60*60*24));
tot_time = tot_time + tstop;
fprintf('\n..........time needed for problem is %f seconds..........\n', tstop(6)); 

for i = length(tot_time):-1:2
    tot_time(i-1) = floor(tot_time(i)/60);
    tot_time(i) = tot_time(i) - tot_time(i-1)*60;
end

fprintf('\n..........total time needed is %d hours %d minutes %f seconds..........\n', tot_time(4), tot_time(5), tot_time(6));

%% plotting the solution obtained through the NN

tstart = tic;

plot_solution_ML(plotStruct.dimModalBasis,...
    plotStruct.liftCoeffA,plotStruct.liftCoeffB,...
    obj_solverIGA.domainLimit_inX,plotStruct.stepMeshX, ...
    plotStruct.u,plotStruct.label_upBoundDomain,plotStruct.label_downBoundDomain, ...
    plotStruct.coefficientForm,plotStruct.simulationCase,...
    plotStruct.degreePolySplineBasis,plotStruct.continuityParameter,...
    plotStruct.space,plotStruct.refDomain1D,plotStruct.map, ...
    numbVerNodes, cmin, cmax);

% plot_solution_ML2(plotStruct.dimModalBasis,...
%     plotStruct.liftCoeffA, plotStruct.liftCoeffB,...
%     obj_solverIGA.domainLimit_inX, plotStruct.stepMeshX, ...
%     plotStruct.label_upBoundDomain, ...
%     plotStruct.label_downBoundDomain, plotStruct.coefficientForm,...
%     plotStruct.simulationCase, plotStruct.degreePolySplineBasis,...
%     plotStruct.continuityParameter, plotStruct.space,...
%     plotStruct.refDomain1D, plotStruct.map, ...
%     numbVerNodes, cmin, ...
%     cmax);

tstop = toc(tstart);
tstop = datevec(tstop./(60*60*24));
fprintf('\n..........time needed for NN to solve the problem is %f seconds..........\n', tstop(6)); 

