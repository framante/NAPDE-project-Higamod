% main training
n_problems = 200;
n_params = 5;
todo = "training";
tot_time = 0.0;
for i = 1:n_problems
    tstart = tic;
    fprintf('\n..........Problem number %d........\n', i);
    v = i*normrnd(0, 10, [n_params,1]);
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
n_problems = 10;
n_params = 5;
todo = "testing";
tot_time = 0.0;
for i = 1:n_problems
    tstart = tic;
    fprintf('\n..........Problem number %d........\n', i);
    v = i*normrnd(0, 10, [n_params,1]);
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

%% main testing for 1 problem
todo = "testing";
n_params = 5;

tot_time = 0.0;

fprintf('\n.......... Problem ........\n');
tstart = tic;
%v = normrnd(0, 10, [n_params,1]);
v = ones(n_params,1);
[plotStruct, obj_solverIGA, numbVerNodes] = higamod_call(v, todo);
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

plot_solution_ML(plotStruct.dimModalBasis,...
    plotStruct.liftCoeffA,plotStruct.liftCoeffB,...
    obj_solverIGA.domainLimit_inX,plotStruct.stepMeshX, ...
    plotStruct.u,plotStruct.label_upBoundDomain,plotStruct.label_downBoundDomain, ...
    plotStruct.coefficientForm,plotStruct.simulationCase,...
    plotStruct.degreePolySplineBasis,plotStruct.continuityParameter,...
    plotStruct.space,plotStruct.refDomain1D,plotStruct.map, ...
    numbVerNodes);
