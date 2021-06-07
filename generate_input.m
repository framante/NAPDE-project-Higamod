function [] = generate_input(xnodes, ynodes, Coeff_forma, Dati, todo)
% (xnodes, ynodes) are two row vectors
    if todo == "training"
       fid = fopen('SolutionInput.txt', 'a+');
    elseif todo == "testing"
       fid = fopen('TestInput.txt', 'w'); 
    end
    
    sz1 = size(xnodes);
    sz2 = size(ynodes);
    
    %disp(sz1);
    %disp(sz2);

    params = struct('mu', Coeff_forma.mu, 'beta1', Coeff_forma.beta1, 'beta2', Coeff_forma.beta2, ...
       'sigma', Coeff_forma.sigma, 'force', Dati.force);
    numb_params = numel(fieldnames(params));
    cell_params = { params.mu, params.beta1, params.beta2, params.sigma, params.force};

    %% this way you get a tensor (n FEM nodes x n QUAD nodes) x n params
    % where each matrix (n FEM node x n QUAD nodes) is referred to a single param

%     for k = 1: numb_params
%         for i=1:sz1(2)
%             for j=1:sz2(2)
%                fprintf(fid,'%12.8f',cell_params{k}(i,j));
%             end
%             fprintf(fid,'\n');
%         end
%     end

    %% this way you get a tensor (n params x n QUAD nodes) x n FEM nodes
    % where each matrix (n params x n QUAD nodes) is referred to a single FEM node
    % THIS IS THE VERSION TO BE USED IN ORDER TO COLLECT DATA FOR THE CNN

    for i=1:sz1(2)
        for k = 1: numb_params
            for j=1:sz2(2)
               fprintf(fid,' %12.8f ',cell_params{k}(i,j));
            end
            fprintf(fid,'\n');
        end
    end

    fclose(fid);


    %% this way you get a matrix ( n params x n QUAD nodes) for a fixed FEM node
%     i = 5; % chosen FEM node
% 
%     for k = 1: numb_params
%         for j=1:sz2(2)
%             fprintf(fid,' %12.8f ',cell_params{k}(i,j));
%         end
%         fprintf(fid,'\n');
%     end
% 
%     fclose(fid);
%     %fprintf('\nparams printed over FEM node %d for all the %d QUAD nodes\n', i , sz2(2));

end