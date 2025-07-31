clear all ; close all ; clc ;

%% Define geomtry and material
N = 8; h = 1; lyr = 4; phi = -pi/8;
% N-gon prism; height of each layer; number of inter-layer planes; twisting
% angle of the prism in each layer
MaxIcr = 500 ; blam = 0.032; 
% Maximum increment number & initial load factor
Kf = 1e-3; Kb = Kf; E0 = 5e3; Abar = 0.1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 45; limrht = 315;
% Left and right limits for the linear range of rotational stiffness
[Node, Panel] = GetDiSym(N,h,lyr,phi);
BarMater = @(Ex)Ogden(Ex, E0); % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);
% Define rotational spring constitutive

% GetDiSym
% Ogden
% EnhancedLinear


%% Set boundary condition
indsupp = find(Node(:,3)<0.01);
nsupp = numel(indsupp);
Supp = [          indsupp(1), 1, 1, 1;
                  indsupp(2), 1, 1, 1;
        indsupp(3:end), zeros(nsupp-2,1)+1, zeros(nsupp-2,1)+1, ones(nsupp-2,1);];
m = size(Node,1);
indp = find( abs( Node(:,3) - max(Node(:,3 ) ) ) < 1e-5 ) ; %indp = indp(3) ;


%% Set boundary condition (around line 25)
% Add displacement control parameters
total_displacement = -0.8;  % Total compression displacement
disp_increment = total_displacement / MaxIcr;  % Per-step increment
control_dof = 3*indp(1);  % Control z-displacement of first top node

% Remove or reduce load
%Load = [];  % No external loads for displacement control
%OR keep small reference load: Load = [indp, 0*ones(npp,1), 0*ones(npp,1), -0.01*ones(npp,1)];
npp = numel(indp) ;
Load = [indp, 0*ones(npp,1), 0*ones(npp,1), -0.0001*ones(npp,1)];

%% Perform analysis (around line 35)
disp_params.control_dof = control_dof;
disp_params.disp_increment = disp_increment;
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar,disp_params);
%[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);

%%

truss.U0 = zeros(3*size(truss.Node,1),1) ;
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr) ;
U_his = real(U_his);
LF_his = real(LF_his);

%%

figure()
PlotOri(Node,angles.Panel,truss.Trigl);
axis equal; axis off;
camproj('perspective')
light
view(117,18)
rotate3d on

% PrepareData 
% PathAnalysis 

%% Visualize simulation

instdof = -indp(1)*3;

interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001,LF_his,instdof,[-inf inf -inf inf])

% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001)

% VisualFold,


%% Plot stored energy vs. pseudo time

% Red line is the total profile. Between red and cyan is the folding
% energy. Between cyan and magenta is the portion of energy for bending. 
% Below magenta is the stretching energy of bars.

STAT = PostProcess(Data,truss,angles);  
figure()
plot(1:size(U_his,2),STAT.PE.strain,'r-','linewidth',2);
grid on
hold on
plot(1:size(U_his,2),STAT.bend.PE+STAT.bar.PE,'c-');
plot(1:size(U_his,2),STAT.bar.PE,'m-');
xlabel('Increment Number (Pseudo-time)','fontsize',14);
ylabel('Stored Energy','fontsize',14);

%% Plot final configuration

Ux = U_his(:,end);

Nodew = truss.Node;
Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 
figure()
PlotOri(Nodew,angles.Panel,truss.Trigl);
axis equal; axis off;
camproj('perspective')
light
view(117,18)
rotate3d on


%% function : 

function [truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarCM,RotSpring,kpf,kpb,Abar,disp_params)

    [Bend] = findbend(Panel, Node);
    [Fold, Bdry, Trigl] = findfdbd(Panel,Bend);
    Bars = [Bend(:,1:2);Fold(:,1:2);Bdry];
    [B, L] = dirc3d(Node,Bars);

    if size(Supp,1) == 0
        rs = []; 
    else
        rs = [reshape([Supp(:,1)*3-2,Supp(:,1)*3-1,Supp(:,1)*3]',[],1),...
              reshape(Supp(:,2:4)',[],1)];
        rs(rs(:,2)==0,:)=[]; rs = rs(:,1);
    end
    
    if numel(Abar)==1
        Abar = Abar*ones(size(Bars,1),1);
    end
    
    pf0 = zeros(size(Fold,1),1); 
    for i = 1:size(Fold,1), pf0(i) = FoldKe(Node,Fold(i,:),kpf,0); end ;
    
    pb0 = zeros(size(Bend,1),1); 
    for i = 1:size(Bend,1), pb0(i) = FoldKe(Node,Bend(i,:),kpb,0); end ;
    
    m = size(Node,1);
    F = zeros(3*m,1);
    indp = Load(:,1);
    F(3*indp-2) = Load(:,2); 
    F(3*indp-1) = Load(:,3); 
    F(3*indp) = Load(:,4);
    
    truss.CM = BarCM ; 
    truss.Node = Node ;
    truss.Bars = Bars ;
    truss.Trigl = Trigl ;
    truss.B = B ; 
    truss.L = L ;
    truss.FixedDofs = unique(rs) ;
    truss.A = Abar ; 
    angles.CM = RotSpring ;
    angles.fold = Fold ;
    angles.bend = Bend ;
    angles.kpf = kpf*ones(1,size(Fold,1)) ;
    angles.kpb = kpb*ones(1,size(Bend,1)) ; 
    angles.pf0 = pf0' ;
    angles.pb0 = pb0'*0+pi ;
    angles.Panel = Panel ;

    if nargin > 9 && ~isempty(disp_params)
        truss.control_dof = disp_params.control_dof;
        truss.disp_increment = disp_params.disp_increment;
        truss.displacement_control = true;
    else
        truss.displacement_control = false;
    end
    

end 

function [Uhis,load_his,Data] = PathAnalysis(truss,angles,F,b_lambda,MaxIcr)

    tol = 1e-6; MaxIter = 50;
    Node = truss.Node;
    AllDofs = [1:3*size(Node,1)];
    U = truss.U0;
    
    Uhis = zeros(3*size(Node,1),MaxIcr);
    Data.Exbar = zeros(size(truss.Bars,1),MaxIcr);
    Data.FdAngle = zeros(size(angles.fold,1),MaxIcr); Data.LFdAngle = Data.FdAngle;
    Data.BdAngle = zeros(size(angles.bend,1),MaxIcr); Data.LBdAngle = Data.BdAngle;
    
    % Modified DOF handling for displacement control
    if truss.displacement_control
        FreeDofs = setdiff(AllDofs, [truss.FixedDofs; truss.control_dof]);
    else
        FreeDofs = setdiff(AllDofs, truss.FixedDofs);
    end
    
    lmd = 0; icrm = 0; MUL = [U,U];
    load_his = zeros(MaxIcr,1);
    F = F(:,1);
    
    % Add convergence failure counter
    convergence_failures = 0;
    max_failures = 5;
    
    while icrm < MaxIcr
        icrm = icrm + 1;
        iter = 0; err = 1;
        
        % For displacement control, prescribe displacement directly
        if truss.displacement_control
            target_disp = icrm * truss.disp_increment;
            U(truss.control_dof) = target_disp;
            fprintf('icrm = %d, prescribed_disp = %6.4f\n', icrm, target_disp);
        else
            fprintf('icrm = %d, lambda = %6.4f\n', icrm, lmd);
        end
        
        while err > tol && iter < MaxIter
            iter = iter + 1;
            [IF,K] = GlobalK_edu_ver(U,Node,truss,angles);
            
            if truss.displacement_control
                % For displacement control: solve for reaction at controlled DOF
                R_free = -IF(FreeDofs);  % Residual forces at free DOFs only
                dU_free = K(FreeDofs,FreeDofs) \ R_free;
                
                % Update only free DOFs
                dUt = zeros(size(U));
                dUt(FreeDofs) = dU_free;
                U(FreeDofs) = U(FreeDofs) + dUt(FreeDofs);
                
                % Keep prescribed displacement fixed
                U(truss.control_dof) = target_disp;
                
                % Calculate reaction force at controlled DOF
                reaction_force = IF(truss.control_dof);
                
                err = norm(dUt(FreeDofs));
                fprintf('    iter = %d, err = %6.4f, reaction = %6.4f\n', iter, err, reaction_force);
                
            else
                % Original force control method
                R = lmd*F - IF;
                MRS = [F,R];
                MUL(FreeDofs,:) = K(FreeDofs,FreeDofs) \ MRS(FreeDofs,:);
                dUp = MUL(:,1); dUr = MUL(:,2);
                
                if iter == 1, dUr = 0*dUr; end
                dlmd = nlsmgd(icrm,iter,dUp,dUr,b_lambda);
                
                dUt = dlmd*dUp + dUr;
                U = U + dUt;
                err = norm(dUt(FreeDofs));
                lmd = lmd + dlmd;
                fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n', iter, err, dlmd);
            end
            
            % Enhanced divergence detection
            if err > 1e8 || any(isnan(U)) || any(isinf(U))
                disp('Divergence detected!'); 
                break; 
            end
        end
        
        % Enhanced convergence failure handling
        if iter >= MaxIter
            convergence_failures = convergence_failures + 1;
            fprintf('Max iterations reached (failure %d/%d) - ', convergence_failures, max_failures);
            
            if convergence_failures >= max_failures
                disp('Too many convergence failures - stopping analysis');
                break;
            else
                disp('reducing step size and retrying');
                if truss.displacement_control
                    % For displacement control, reduce displacement increment
                    truss.disp_increment = truss.disp_increment * 0.5;
                    fprintf('New displacement increment: %6.4f\n', truss.disp_increment);
                else
                    % For force control, reduce constraint radius
                    b_lambda = b_lambda / 2;
                end
                icrm = icrm - 1;  % Retry this increment
                if icrm > 0
                    U = Uhis(:,icrm);
                    if ~truss.displacement_control
                        lmd = load_his(icrm);
                    end
                end
                continue;
            end
        end
        
        % Original adaptive step size logic (for force control only)
        if ~truss.displacement_control
            if iter > 15
                b_lambda = b_lambda/2;
                disp('Reduce constraint radius!')
                icrm = icrm-1;
                U = Uhis(:,max(icrm,1));
                lmd = load_his(max(icrm,1));
                continue;
            elseif iter < 3
                disp('Increase constraint radius!')
                b_lambda = b_lambda*1.5;
            end
        end
        
        % Store results
        Uhis(:,icrm) = U;
        if truss.displacement_control
            load_his(icrm) = IF(truss.control_dof);  % Store reaction force
        else
            load_his(icrm) = lmd;  % Store load factor
        end
        
        [Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei] = GetData(U,Node,truss,angles);
        Data.Exbar(:,icrm) = Exbari;
        Data.FdAngle(:,icrm) = FdAnglei; Data.LFdAngle(:,icrm) = LFdAnglei;
        Data.BdAngle(:,icrm) = BdAnglei; Data.LBdAngle(:,icrm) = LBdAnglei;
        
        % Reset convergence failure counter on successful increment
        convergence_failures = 0;
    end
    
    % Clean up unused array elements
    icrm = icrm + 1;
    Uhis(:,icrm:end) = [];
    load_his(icrm:end,:) = [];
    Data.Exbar(:,icrm:end) = [];
    Data.FdAngle(:,icrm:end) = []; Data.LFdAngle(:,icrm:end) = [];
    Data.BdAngle(:,icrm:end) = []; Data.LBdAngle(:,icrm:end) = [];

end

function [IF,K] = GlobalK_edu_ver(Ui,Node,truss,angles)

    Nn = size(Node,1); % Nn = numbers of nodes 
    IFb = zeros(3*Nn,1); IFp = IFb; % what is IFb and IFp ? Why is it 3 times the number of nodes ?
    indi = zeros(36*size(truss.Bars,1),1); indj = indi; kentry = indi; % where does the 36 comes from? What are indi and indj?

    Nodenw(:,1) = Node(:,1)+Ui(1:3:end); % Compute new node position by adding displacement components 
    Nodenw(:,2) = Node(:,2)+Ui(2:3:end);
    Nodenw(:,3) = Node(:,3)+Ui(3:3:end);

    % Next sections we have a clear 
    
    for bel = 1:size(truss.Bars,1) 
        eDof = [(-2:0)+(truss.Bars(bel,1)*3),(-2:0)+(truss.Bars(bel,2)*3)]';
        [~,Rbe,Kbe] = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel));
        IFb(eDof) = IFb(eDof)+Rbe;
        I=repmat(eDof,1,6); J=I';
        indi(36*(bel-1)+1:36*bel) = I(:);
        indj(36*(bel-1)+1:36*bel) = J(:); 
        kentry(36*(bel-1)+1:36*bel) = Kbe(:);
    end

    Kb = sparse(indi,indj,kentry,3*Nn,3*Nn);
    
    indi = zeros(144*size(angles.bend,1),1); indj = indi; kentry = indi;
    Lbend = truss.L(1:size(angles.bend,1));
    for del = 1:size(angles.bend,1)
        eDof = reshape([3*angles.bend(del,:)-2;...
                        3*angles.bend(del,:)-1;...
                        3*angles.bend(del,:)],12,1);
        bend = angles.bend(del,:);
        [~,Rpe,Kpe] = FoldKe(Nodenw,bend,angles.kpb,angles.pb0(del),Lbend(del),angles.CM);
        IFp(eDof) = IFp(eDof)+Rpe;
        I=repmat(eDof,1,12); J=I';
        indi(144*(del-1)+1:144*del) = I(:);
        indj(144*(del-1)+1:144*del) = J(:); 
        kentry(144*(del-1)+1:144*del) = Kpe(:);
    end;
    Kbd = sparse(indi,indj,kentry,3*Nn,3*Nn);
    if isempty(Kbd), Kbd = zeros(3*Nn); end
    
    indi = zeros(144*size(angles.fold,1),1); indj = indi; kentry = indi;
    Lfold = truss.L(size(angles.bend,1)+1:size(angles.bend,1)+size(angles.fold,1));
    for fel = 1:size(angles.fold,1)
        eDof = reshape([3*angles.fold(fel,:)-2;...
                        3*angles.fold(fel,:)-1;...
                        3*angles.fold(fel,:)],12,1);
        fold = angles.fold(fel,:);
        [~,Rpe,Kpe] = FoldKe(Nodenw,fold,angles.kpf,angles.pf0(fel),Lfold(fel),angles.CM);
        IFp(eDof) = IFp(eDof)+Rpe;
        I=repmat(eDof,1,12); J=I';
        indi(144*(fel-1)+1:144*fel) = I(:);
        indj(144*(fel-1)+1:144*fel) = J(:); 
        kentry(144*(fel-1)+1:144*fel) = Kpe(:);
    end;
    Kfd = sparse(indi,indj,kentry,3*Nn,3*Nn);
    
    IF = IFb+IFp;
    K = Kb+Kbd+Kfd;
    K = (K+K')/2;

end 



function dl = nlsmgd_disp(step, ite, dup, dur, cmp, control_dof, target_disp_incr)
    global dupp1 sinal dupc1
    
    if step==1 && ite==1
        sinal = 1;
        dupp1 = dup;
        dupc1 = dup;
    end
    
    if ite==1
        sinal = sinal * sign(dupp1(control_dof));
        dupp1 = dup;
        dupc1 = dup;
        % Direct displacement increment control
        if abs(dup(control_dof)) > 1e-12  % Avoid division by zero
            dl = target_disp_incr / dup(control_dof);
        else
            dl = 0;
        end
    else
        % Newton correction to maintain displacement constraint
        if abs(dup(control_dof)) > 1e-12
            dl = -dur(control_dof) / dup(control_dof);
        else
            dl = 0;
        end
    end
end


function [] = PlotOri(Node,Panel,Trigl,lstyle,alfa,color)
    if nargin < 4
        lstl = '-'; al = 1; cc = 'g';
    else
        lstl = lstyle; al = alfa; cc = color;
    end
    
    if ~isempty(Trigl)
        patch('faces', Trigl, 'vertices', Node, 'facecolor', cc, ...
              'linestyle', 'none', 'facelighting', 'flat', 'edgecolor', (1-al)*[1 1 1]);
    end
    
    hold on;
    Panelsize = cellfun(@numel,Panel);
    Ptri = cell(sum(Panelsize==3),1);
    Pquad = cell(sum(Panelsize==4),1);
    flg3 = find(Panelsize==3);
    flg4 = find(Panelsize==4);
    for i = 1:numel(flg3), Ptri{i} = Panel{flg3(i)}; end
    for j = 1:numel(flg4), Pquad{j} = Panel{flg4(j)}; end
    
    if ~isempty(Trigl)
        patch('faces', cell2mat(Ptri), 'vertices', Node, 'facecolor', 'none', ...
              'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
        patch('faces', cell2mat(Pquad), 'vertices', Node, 'facecolor', 'none', ...
              'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
    else
        patch('faces', cell2mat(Ptri), 'vertices', Node, 'facecolor', cc, ...
              'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
        patch('faces', cell2mat(Pquad), 'vertices', Node, 'facecolor', cc, ...
              'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
    end
end 