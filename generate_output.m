function [] = generate_output(solution, todo)
    
    if todo == "training"
        fid = fopen('SolutionOutput.txt', 'a+');
    elseif todo == "testing"
       fid = fopen('TestOutput.txt', 'w');
    % 'w' to overwrite, 'a+' not to overwrite
    end
    [m,n]=size(solution);
    
    %% prints the solution of the problem evaluated over a grid (n FEM nodes x n QUAD nodes)
    
    for i=1:m
        for j=1:n
            fprintf(fid,' %12.8f ',solution(i,j));
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);
%     sz = size(solution);
%     %fprintf('\n\nthe solution is given by %d rows -> FEM nodes and %d columns -> QUADRATURE nodes\n\n', sz(1), sz(2));

    %% prints the solution fixing a FEM node for all the QUAD nodes 
    % (i.e. vector (1 x n QUAD nodes))
    
%     i = 5; %chosen FEM node
%     for j=1:n
%         fprintf(fid,' %12.8f ',solution(i,j));
%     end
%     fprintf(fid,'\n');
%         
%     fclose(fid);
   % sz = size(solution);
   % fprintf('\nsolution printed over FEM node %d for all the %d QUAD nodes\n', i , sz(2));
end

