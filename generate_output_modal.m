function [] = generate_output_modal(u, todo, size_mb, cutx, hx)
    
    if todo == "training"
        fid = fopen('SolutionOutput.txt', 'a+');
    elseif todo == "testing"
       fid = fopen('TestOutput.txt', 'w');
    % 'w' to overwrite, 'a+' not to overwrite
    end
    
    %% SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne*(p - k) + k + 1;           % NUMBER of Control Point on the X Direction
                                        % for the IGA Mesh
    
    %% prints a matrix corresponding to the modal coefficients for each FEM node (n modes x n FEM nodes)
    
    for m=1:size_mb
        for j=1:nnx
            fprintf(fid,' %12.8f ', u((m-1) * nnx + j));
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);
% remember that here u = [u1, ...., um]
% where ui is a vector of nnx, one value for each FEM node


end