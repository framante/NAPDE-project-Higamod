function [] = print_solution_IGA_scatter( size_mb,a_ril,b_ril,cutx,...
                                          hx, u,bc_up,bc_down,Coeff_forma,...
                                          caso,p,k,space,geometry,map, Dati, ...
                                          numbVerNodes, todo)
    % p = degreePolySplineBasis
    % k = continuityParameter

    % IMPORT CLASSES
            
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
                
    %% SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne*(p - k) + k + 1;           % NUMBER of Control Points on the X Direction
                                        % for the IGA Mesh
    %disp(nnx);

   %% SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = numbVerNodes;

   %% CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

%     augVerNodes = linspace(0,1,M);

    % Vertical direction
            
    obj_gaussLegendre_2 = IntegrateHandler();
    obj_gaussLegendre_2.numbQuadNodes = M;
    [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2); 
    
    objVertQuadRule = IntegrateHandler();

    objVertQuadRule.leftBoundInterval = 0;
    objVertQuadRule.rightBoundInterval = 1;
    objVertQuadRule.inputNodes = verGLNodes;
    objVertQuadRule.inputWeights = verWeights;

    [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule); % nodes and weights rescaled in the desired interval
    %disp(augVerNodes);

    %% INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(nnx,M);
    sol_x = zeros(nnx,M);
    sol_y = zeros(nnx,M);

    import Core.BasisHandler

    obj_newModalBasis = BasisHandler();

    obj_newModalBasis.dimModalBasis         = size_mb;
    obj_newModalBasis.evalNodesY            = verGLNodes;
    %obj_newModalBasis.evalWeightsY          = augVerWeights;
    obj_newModalBasis.labelUpBoundCond      = bc_up{1};
    obj_newModalBasis.labelDownBoundCond    = bc_down{1};
    obj_newModalBasis.coeffForm             = Coeff_forma;

    [coeffModalBase,coeffModalBaseDer] = newModalBasis(obj_newModalBasis);

    %---------------------------------------------------------------------%
    % Note:
    % The coefficients 'coeffModalBase' and 'coeffModalBaseDer' are
    % vectors corresponding respectively to the coefficients of the
    % modal bases and derivative of the modal basis on the points
    % assigned in the Y direction.
    %---------------------------------------------------------------------%

    %% SETTING OF THE ISOGEOMETRIC MESH IN THE X DIRECTION
    
    xx3 = {linspace(cutx(1),cutx(2),nnx)};
    evalNodesX = linspace(cutx(1),cutx(2),nnx);
    
    solVect = zeros(nnx * M,1);
    solMat = zeros(M,nnx);
    sol = zeros(nnx, M);

    for h = 1:nnx
        for k = 1:M
            for imb = 1:size_mb
                sol(h,k)   = sol(h,k)   + u(h+(imb-1)*nnx)*coeffModalBase(k,imb);
            end
        sol(h,k) = sol(h,k)+ a_ril(1) * augVerNodes(k) + b_ril(1);
        end
    end
    
    finalSol = sol;
    
    for ii = 1:M
        solMat(ii,:) = finalSol((ii-1)*nnx + 1 : ii * nnx);
    end
    
  %%   PRINT
    xnodes = evalNodesX;
    ynodes = augVerNodes';
    generate_output(solMat', todo);
    generate_input(xnodes, ynodes, Coeff_forma, Dati, todo);
    
  %%   PLOT
  
%     finalSol = sol;
%     for ii = 1:M
%         solMat(ii,:) = finalSol((ii-1)*nnx + 1 : ii * nnx);
%     end
%     
% 
%     X = mapOut(evalNodesX,augVerNodes,map,1);
%     Y = mapOut(evalNodesX,augVerNodes,map,2);
%     
%     Nx = 103;
%     Ny = 103;
%     
%     x_eval = linspace(0,1,Nx);
%     y_eval = linspace(0,1,Ny);
%     Xeval  = mapOut(x_eval,y_eval,map,1);
%     Yeval  = mapOut(x_eval,y_eval,map,2);
%     
%     higaSol = griddata(X,Y,solMat,Xeval,Yeval);
%     [numbX,numbY] = size(higaSol);
%     
%     for ii = 1:numbX
%         for jj = 1:numbY
%             if(isnan(higaSol(ii,jj)))
%                 higaSol(ii,jj) = 0;
%             end
%         end
%     end

  %%   %--------------------------------------------------%
    % CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
    %--------------------------------------------------%
    
%     minX = min(min(Xeval));
%     minY = min(min(Yeval));
%     maxX = max(max(Xeval));
%     maxY = max(max(Yeval));
%     
%     scale = 0.01;
    
  %%   CONTOUR PLOT FOR THE SOLUTION
    
%     figure;
   
    %[~,~] = contourf(Xeval,Yeval,higaSol,20);    

%     [~,~] = contourf(Xeval,Yeval,higaSol,20,'edgecolor','none');    
%     colormap(jet);
%     cmin = min(min(higaSol));
%     cmax = max(max(higaSol));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     daspect([1 1 scale])
%     pbaspect([1 1 scale])
%     set(gca, 'FontSize', 14)
%     hold on;
    
  %%   SURF PLOT OF THE SOLUTION
%     
%     figure;
%     mesh(Xeval,Yeval,higaSol);
%     colormap(jet);
%     cmin = min(min(higaSol));
%     cmax = max(max(higaSol));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY cmin cmax]);
%     axis equal
%     daspect([1 1 scale])
%     pbaspect([1 1 scale])
%     set(gca, 'FontSize', 14)
%     

  %%   SURF PLOT OF THE CUTS IN THE SOLUTION
%     
%     xCut1 = (Nx-1)/2 + 1;   % Position of the centerline
%     xCut2 = (Nx-3)/4 + 1;   % Position of the first quarter line
%     xCut3 = (Nx-3)*3/4 + 2; % Position of the third quarter line
%     
%     solX1 = higaSol(xCut1,:);   % Solution at the centerline
%     solX2 = higaSol(xCut2,:);   % Solution at the first quarter line
%     solX3 = higaSol(xCut3,:);   % Solution at the third quarter line
%     
%     figure;
%     plot(Xeval(xCut1,:),solX1,'b','LineWidth',2); hold on;
%     plot(Xeval(xCut2,:),solX2,'r','LineWidth',2); 
%     plot(Xeval(xCut3,:),solX3,'g','LineWidth',2);
%     cmin = 0.02; % min(min(higaSol));
%     cmax = 0.02; % max(max(higaSol));
%     axis([minX maxX cmin cmax]);
%     axis equal
%     set(gca, 'FontSize', 14)
%     grid on
%     daspect([1 scale scale])
%     pbaspect([1 scale scale])
%     
    errL2 = 0;
    errH1 = 0;

end
