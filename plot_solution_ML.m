function [] = print_solution_ML( size_mb, a_ril, b_ril, cutx, ...
                                 hx, u, bc_up, bc_down, Coeff_forma, ...
                                 caso, p, k, space, geometry, map, ...
                                 numbVerNodes)
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
    
    solMat_transpose = readmatrix('output_ML.txt');
    solMat = solMat_transpose';
    
    approx= zeros(M, nnx);
    
    for m = 1:M
        approx(m, 1:nnx)  = sp_eval(solMat(m,(1:nnx)), space, geometry, xx3);
    end
    % COMPUTATION OF THE NURBS BASE
    solMat  = approx;    
    
  %%  PLOT 
    X = mapOut(evalNodesX,augVerNodes,map,1);
    Y = mapOut(evalNodesX,augVerNodes,map,2);
    
    Nx = 103;
    Ny = 103;
    
    x_eval = linspace(0,1,Nx);
    y_eval = linspace(0,1,Ny);
    Xeval  = mapOut(x_eval,y_eval,map,1);
    Yeval  = mapOut(x_eval,y_eval,map,2);
    
    higaSol = griddata(X,Y,solMat,Xeval,Yeval);
    [numbX,numbY] = size(higaSol);
    
    for ii = 1:numbX
        for jj = 1:numbY
            if(isnan(higaSol(ii,jj)))
                higaSol(ii,jj) = 0;
            end
        end
    end

  %%   %--------------------------------------------------%
    % CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
    %--------------------------------------------------%
    
    minX = min(min(Xeval));
    minY = min(min(Yeval));
    maxX = max(max(Xeval));
    maxY = max(max(Yeval));
    
    scale = 0.01;
    
  %%   CONTOUR PLOT FOR THE SOLUTION
    
    figure;
   
    %[~,~] = contourf(Xeval,Yeval,higaSol,20);    

    [~,~] = contourf(Xeval,Yeval,higaSol,20,'edgecolor','none');    
    colormap(jet);
    cmin = min(min(higaSol));
    cmax = max(max(higaSol));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    daspect([1 1 scale])
    pbaspect([1 1 scale])
    set(gca, 'FontSize', 14)
    hold on;
    
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
