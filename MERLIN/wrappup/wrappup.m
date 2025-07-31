
clear all ; close all ; clc ;

%% Define geomtry and material

N = 6; h = 1; lyr = 4; phi = 2*pi/8;
% N-gon prism; height of each layer; number of inter-layer planes; twisting
% angle of the prism in each layer
MaxIcr = 180; blam = 0.032; 
% Maximum increment number & initial load factor
Kf = 1e-3; Kb = Kf; E0 = 5e3; Abar = 0.1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 90; limrht = 310;
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

npp = numel(indp) ;

Lo = 0.001 ; %0.5 
Load = [ indp , 0*ones(npp,1) , 0*ones(npp,1) , -Lo*ones(npp,1) ; ] ;

%% Prepare data 

[truss, angles, F] = PrepareData( Node , Panel , Supp , Load , BarMater , RotSpring , Kf , Kb , Abar ) ;
%[truss, angles, F] = PrepareData( Node , Panel , Supp , Load , BarMater , RotSpring , Kf , Abar ) ;

%% Perform analysis

truss.U0 = zeros(3*size(truss.Node,1),1) ;
MaxSteps = 100 ;
max_disp = 0.1

 
%[U_his,LF_his,Data] = PathAnalysis_Displacement(truss, angles, indp, max_disp, MaxIcr) ; 
%%
truss.U0 = zeros(3*size(truss.Node,1),1) ;

[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr) ;

U_his = real(U_his);
LF_his = real(LF_his);

% %%
% 
% figure()
% PlotOri(Node,angles.Panel,truss.Trigl);
% axis equal; axis off;
% camproj('perspective')
% light
% view(117,18)
% rotate3d on
% 


%% Plot final configuration

Ux = U_his(:,end);

Nodew = truss.Node;
Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 

%%
figure()
PlotOri(Nodew,angles.Panel,truss.Trigl);
axis equal; axis off;
camproj('perspective')
light
view(117,18)
rotate3d on

%%


%% Visualize simulation

instdof = -indp(1)*3;
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001,LF_his,instdof,[-inf inf -inf inf])

%% function : 

function [truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarCM,RotSpring,kpf,kpb,Abar) %PrepareData(Node,Panel,Supp,Load,BarCM,RotSpring,kpf, Abar)

    %[Bend] = findbend(Panel, Node);
    %[Fold, Bdry, Trigl] = findfdbd(Panel,Bend);
    [Fold, Bdry, Trigl] = findfdbd(Panel);
    %Bars = [Bend(:,1:2);Fold(:,1:2);Bdry];
    Bars = [Fold(:,1:2);Bdry];
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
    
    % pb0 = zeros(size(Bend,1),1); 
    % for i = 1:size(Bend,1), pb0(i) = FoldKe(Node,Bend(i,:),kpb,0); end ;
    % 
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
    %angles.bend = Bend ;
    angles.kpf = kpf*ones(1,size(Fold,1)) ;
    % angles.kpb = kpb*ones(1,size(Bend,1)) ; 
    angles.pf0 = pf0' ;
%    angles.pb0 = pb0'*0+pi ;
    angles.Panel = Panel ;
end 

function [Uhis, R_his, Data] = PathAnalysis_Displacement(truss, angles, indp, max_disp, MaxIcr)

    Node = truss.Node;
    Ndof = 3 * size(Node,1);
    U = truss.U0;

    % Prescribed Z DOFs
    dof_z = 3 * indp;

    % Other DOFs
    all_dofs = 1:Ndof;
    free_dofs = setdiff(all_dofs, dof_z);

    % Displacement path
    disp_step = max_disp / MaxIcr ;
    disp_vals = - disp_step * ( 1:MaxIcr ) ;  % negative Z

    % Init logs
    Uhis = zeros(Ndof, MaxIcr) ;
    R_his = zeros(MaxIcr, 1);  % reaction force at Z DOFs
    Data.Exbar = zeros(size(truss.Bars,1), MaxIcr) ;
    Data.FdAngle = zeros(size(angles.fold,1), MaxIcr) ;
    Data.LFdAngle = Data.FdAngle ; 

    for icrm = 1:MaxIcr
        fprintf("Step %d / %d: displacement = %.4f\n", icrm, MaxIcr, disp_vals(icrm));

        % Impose current displacement on Z DOFs
        U(dof_z) = disp_vals(icrm);

        % Compute internal force and stiffness
        [IF, K] = GlobalK_edu_ver(U, Node, truss, angles);

        % Solve for the free DOFs to satisfy equilibrium
        R = -IF;  % residual (no external force)
        U(free_dofs) = K(free_dofs, free_dofs) \ R(free_dofs);

        % Log full displacement
        Uhis(:, icrm) = U;

        % Log reaction force at Z DOFs
        R_his(icrm) = sum(IF(dof_z));

        % Log energy and deformation data
        [Exbari, FdAnglei, LFdAnglei] = GetData(U, Node, truss, angles);
        Data.Exbar(:,icrm) = Exbari;
        Data.FdAngle(:,icrm) = FdAnglei;
        Data.LFdAngle(:,icrm) = LFdAnglei;
    end
end



function [Uhis,load_his,Data] = PathAnalysis(truss,angles,F,b_lambda,MaxIcr)

    tol = 1e-6; MaxIter = 50 ; 
    Node = truss.Node ;

    AllDofs = [ 1:3*size(Node,1) ] ;

    %%%%%%%%%%%%%%%%%%%%%%%%
    max_disp = -0.6 ; 
    max_disp_step = max_disp / MaxIcr ; 
    current_disp = max_disp_step ; 
    U_add = truss.U0 ; 
    indp = find( abs( Node(:,3) - max(Node(:,3 ) ) ) < 1e-5 ) ; %indp = indp(3) ;
    npp = numel(indp) ;
    disp = [ indp , 0*ones(npp,1) , 0*ones(npp,1) , -current_disp*ones(npp,1) ; ] ;
    % 
    % indp = disp(:,1);
    % U_add(3*indp-2) = disp(:,2); 
    % U_add(3*indp-1) = disp(:,3); 
    % U_add(3*indp) = disp(:,4);

    %%%%%%%%%%%%%%%%%%%%%%%%
    U = truss.U0;

    Uhis = zeros(3*size(Node,1),MaxIcr) ;
    Data.Exbar = zeros(size(truss.Bars,1),MaxIcr) ; 
    Data.FdAngle = zeros(size(angles.fold,1),MaxIcr) ; Data.LFdAngle = Data.FdAngle ;
    %Data.BdAngle = zeros(size(angles.bend,1),MaxIcr) ; Data.LBdAngle = Data.BdAngle ;
    FreeDofs = setdiff(AllDofs,truss.FixedDofs) ;

    lmd = 0 ; icrm = 0 ; MUL = [ U , U ] ; dU = U ; 
    load_his = zeros( MaxIcr , 1 ) ;

    F = F(:,1) ;

    while icrm<MaxIcr

        icrm = icrm+1 ;
        iter = 0; err = 1 ;
        fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd) ;

        while err>tol && iter<MaxIter
            iter = iter+1 ; 
            % There are two versions of the function for assembling global
            % stiffness matrix:
            % 'GlobalK_fast_ver': speed-optimized version;
            % 'GlobalK_edu_ver': easy-to-read educational version.
            % [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);

            %[IF,K] = GlobalK_edu_ver(U,Node,truss,angles) ;

            % indp = disp(:,1);
            % U(3*indp-2) = disp(:,2); 
            % U(3*indp-1) = disp(:,3); 
            U(3*indp) = U(3*indp) + max_disp_step ; 

            [IF,K] = GlobalK_edu_ver(U,Node,truss, angles) ;

            R = -IF;  % residual (no external force)
            dU(FreeDofs) = K(FreeDofs, FreeDofs) \ R(FreeDofs);
            U = U + dU ; 

            R = lmd*F-IF;   MRS = [F,R] ;
            % %R = F-IF;   MRS = [F,R] ;
            MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:) ;
            dUp = MUL(:,1); dUr = MUL(:,2) ;
            if iter==1, dUr = 0*dUr; end
            % 
            dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda) ;
            dUt = dlmd*dUp+dUr ; 
            % %dUt=dUp+dUr
            % 
            U = U+dUt ;
            err = norm(dUt(FreeDofs)) ; 

            %err = norm(dU(FreeDofs)) ; 
            lmd = lmd + dlmd ;
            fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n', iter , err , dlmd ) ;
            %fprintf('    iter = %d, err = %6.4f', iter , err ) ;
            if err > 1e8, disp('Divergence!') ; break ; end
        end
    

        if iter>15
            b_lambda = b_lambda/2;
            %disp('Reduce constraint radius!')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));
            %lmd = load_his(max(icrm,1));

        elseif iter<3
            %disp('Increase constraint radius!')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            load_his(icrm) = lmd; 
            %[Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei] = GetData(U,Node,truss,angles);
            [Exbari,FdAnglei,LFdAnglei] = GetData(U,Node,truss,angles);
            %[Exbari] = GetData(U, Node, truss) ; 
            Data.Exbar(:,icrm) = Exbari ; 
            Data.FdAngle(:,icrm) = FdAnglei; Data.LFdAngle(:,icrm) = LFdAnglei;
            %Data.BdAngle(:,icrm) = BdAnglei; Data.LBdAngle(:,icrm) = LBdAnglei;
        else
            Uhis(:,icrm) = U;
            load_his(icrm) = lmd; 
            %[Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei] = GetData(U,Node,truss,angles);
            %[Exbari,FdAnglei,LFdAnglei] = GetData(U,Node,truss,angles);
            [Exbari,FdAnglei,LFdAnglei] = GetData(U,Node,truss,angles);
            %[Exbari] = GetData(U, Node, truss) ; 
            Data.Exbar(:,icrm) = Exbari; 
            Data.FdAngle(:,icrm) = FdAnglei; Data.LFdAngle(:,icrm) = LFdAnglei;
            %Data.BdAngle(:,icrm) = BdAnglei; Data.LBdAngle(:,icrm) = LBdAnglei;
        end
        % 

        % Uhis(:,icrm) = U;
        % load_his(icrm) = lmd; 
        % %[Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei] = GetData(U,Node,truss,angles);
        % %[Exbari,FdAnglei,LFdAnglei] = GetData(U,Node,truss,angles);
        % [Exbari,FdAnglei,LFdAnglei] = GetData(U,Node,truss,angles);
        % %[Exbari] = GetData(U, Node, truss) ; 
        % Data.Exbar(:,icrm) = Exbari; 
        % Data.FdAngle(:,icrm) = FdAnglei; Data.LFdAngle(:,icrm) = LFdAnglei;
        % %Data.BdAngle(:,icrm) = BdAnglei; Data.LBdAngle(:,icrm) = LBdAnglei;

        current_disp = current_disp + max_disp_step ;
        disp = [ indp , 0*ones(npp,1) , 0*ones(npp,1) , -current_disp*ones(npp,1) ; ] ;
        
        %U = U + U_add ; 

    end
    
    icrm = icrm+1;
    Uhis(:,icrm:end) = [];

    load_his(icrm:end,:) = [];
    Data.Exbar(:,icrm:end) = []; 

    Data.FdAngle(:,icrm:end) = []; Data.LFdAngle(:,icrm:end) = [];
    %Data.BdAngle(:,icrm:end) = []; Data.LBdAngle(:,icrm:end) = [];


end 

function [IF,K] = GlobalK_edu_ver(Ui,Node,truss, angles) %[IF,K] = GlobalK_edu_ver(Ui,Node,truss,angles)

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
    
    % indi = zeros(144*size(angles.bend,1),1); indj = indi; kentry = indi;
    % Lbend = truss.L(1:size(angles.bend,1));
    % for del = 1:size(angles.bend,1)
    %     eDof = reshape([3*angles.bend(del,:)-2;...
    %                     3*angles.bend(del,:)-1;...
    %                     3*angles.bend(del,:)],12,1);
    %     bend = angles.bend(del,:);
    %     [~,Rpe,Kpe] = FoldKe(Nodenw,bend,angles.kpb,angles.pb0(del),Lbend(del),angles.CM);
    %     IFp(eDof) = IFp(eDof)+Rpe;
    %     I=repmat(eDof,1,12); J=I';
    %     indi(144*(del-1)+1:144*del) = I(:);
    %     indj(144*(del-1)+1:144*del) = J(:); 
    %     kentry(144*(del-1)+1:144*del) = Kpe(:);
    % end;
    % 
    % Kbd = sparse(indi,indj,kentry,3*Nn,3*Nn);
    % if isempty(Kbd), Kbd = zeros(3*Nn); end
    
    indi = zeros(144*size(angles.fold,1),1); indj = indi; kentry = indi;
    % Lfold = truss.L(size(angles.bend,1)+1:size(angles.bend,1)+size(angles.fold,1));
    Lfold = truss.L(1:size(angles.fold,1));
    for fel = 1:size(angles.fold,1)
        eDof = reshape([3*angles.fold(fel,:)-2;...
                        3*angles.fold(fel,:)-1;...
                        3*angles.fold(fel,:)],12,1) ;
        fold = angles.fold(fel,:) ;
                 
                     % BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel)); 
        [~,Rpe,Kpe] = FoldKe(Nodenw,fold,angles.kpf,angles.pf0(fel),Lfold(fel),angles.CM) ;

        IFp(eDof) = IFp(eDof)+Rpe ;
        fprintf('    eDof = %d', eDof)
        disp(Rpe)
        I=repmat(eDof,1,12); J=I';
        indi(144*(fel-1)+1:144*fel) = I(:) ;
        indj(144*(fel-1)+1:144*fel) = J(:) ; 
        kentry(144*(fel-1)+1:144*fel) = Kpe(:) ;

    end ;

    Kfd = sparse(indi,indj,kentry,3*Nn,3*Nn);
    
    
    IF = IFb +IFp ;
    %K = Kb+Kbd+Kfd;
    %K = Kb ; %+ Kfd ;
    K = Kb + Kfd ; 
    K = (K+K')/2;

end 



function dl=nlsmgd(step,ite,dup,dur,cmp)
    % Modified generalized displacement control method.
    
    % Define the incremental load factor signal
    global dupp1 sinal
    
    if step==1 && ite==1
        sinal=1;
        dupp1=dup;
    end
    
    if ite==1
        sinal=sinal*sign(dot(dupp1,dup));
        dupp1=dup;
    end
    
    % Calculate the incremental load factor
    global dupc1 numgsp
    
    if ite==1
        if step==1
            dl=cmp;
            numgsp=dot(dup,dup);
            dupc1=dup;
        else
            gsp=numgsp/dot(dup,dup);
            dl=sinal*cmp*sqrt(gsp);
            dupc1=dup;
        end
    else
        dl=-dot(dupc1,dur)/dot(dupc1,dup);
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

function [Node, Panel] = GetDiSym(N,h,lyr,phi)

    if numel(phi) == 1
        rotangle = (0:lyr-1)*phi;
    else
        rotangle = phi;
    end
    rdl = zeros(lyr,N);
    for i = 1:lyr
        rdl(i,:) = linspace(rotangle(i),2*pi/N*(N-1)+rotangle(i),N);
    end
    Xcood = cos(reshape(rdl',[],1));
    Ycood = sin(reshape(rdl',[],1));
    Zcood = h*reshape(repmat(0:lyr-1,N,1),[],1);
    Node = [Xcood,Ycood,Zcood];
    
    PMat = zeros(2*N*(lyr-1),3);
    for i = 1:lyr-1
        PMat((1:N)+2*(i-1)*N,:) = [1:N;(1:N)+N;mod(1:N,N)+N+1]'+(i-1)*N;
        PMat((N+1:2*N)+2*(i-1)*N,:) = [1:N;mod(1:N,N)+1;mod(1:N,N)+N+1]'+(i-1)*N;
    end
    Panel = mat2cell(PMat,ones(size(PMat,1),1),3);

end 

function [fold, bdry, Trigl] = findfdbd(Panel)%findfdbd(Panel,bend)

    Nn = max(cellfun(@max,Panel)); % Panel : constains list of nodes making one panel
    % triangularization
    Panelsize = cellfun(@numel,Panel);
    Ptri = cell(sum(Panelsize==3),1);
    flg = find(Panelsize==3);

    for i = 1:sum(Panelsize==3), Ptri{i} = Panel{flg(i)}; end

    %Triglraw = [bend(:,[1,2,3]);bend(:,[1,2,4]);cell2mat(Ptri)];
    Triglraw = [cell2mat(Ptri)];
    Triglraw = sort(Triglraw,2);
    Trigl = unique(Triglraw ,'rows');

    % formulate connectivity matrix
    Comm = sparse(Nn,size(Trigl,1));

    for i=1:size(Trigl,1), Comm(Trigl(i,:),i) = true; end;
    % search for fold lines
    Ge = Comm'*Comm;
    [mf, me] = find(triu(Ge==2)); % triangular meshes that share two common nodes
    fold = zeros(length(mf),4);

    for i=1:length(mf)
        [link,ia,ib] = intersect(Trigl(mf(i),:),Trigl(me(i),:));
        oftpa = setdiff(1:3,ia);
        oftpb = setdiff(1:3,ib);
        fold(i,:) = [link,Trigl(mf(i),oftpa),Trigl(me(i),oftpb)];
    end

    % fdandbd = sort(fold(:,1:2),2);
    % onlybd = sort(bend(:,1:2),2);
    % [~,ibd] = intersect(fdandbd,onlybd,'rows');
    % fold(ibd,:) = [];
    
    % search for boundaries
    Edge = sort( [ Trigl(:,1) Trigl(:,2) ; Trigl(:,2) Trigl(:,3); Trigl(:,3) Trigl(:,1)],2) ;
    [u,~,n] = unique( Edge ,'rows' ) ;
    counts = accumarray( n(:) , 1 ) ;
    bdry = u( counts==1 , : ) ;
end

function [bend] = findbend(Panel, Node)

    bend = zeros(length(Panel),4) ;
    
    for i = 1:length(Panel)
        if numel(Panel{i}) == 4
            L1 = norm( ( Node(Panel{i}(1),:) - Node(Panel{i}(3),:)) ) ;
            L2 = norm( ( Node(Panel{i}(4),:) - Node(Panel{i}(2),:)) ) ;
            if L1>L2, lclbend = [2,4,1,3] ;
            else lclbend = [1,3,2,4] ; end ;
            bend(i,:) = Panel{i}(lclbend) ;
        end
    end ;

    bend(sum(bend,2)==0,:)=[] ;  
end

function [B, L] = dirc3d(Node,Ele)
    Ne = size(Ele,1); Nn = size(Node,1);
    D = [Node(Ele(:,2),1)-Node(Ele(:,1),1), Node(Ele(:,2),2)-Node(Ele(:,1),2),...
        Node(Ele(:,2),3)-Node(Ele(:,1),3)];
    L = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
    D = [D(:,1)./L D(:,2)./L D(:,3)./L];
    B = sparse(repmat((1:Ne)',1,6),[3*Ele(:,1)-2 3*Ele(:,1)-1 3*Ele(:,1),...
               3*Ele(:,2)-2 3*Ele(:,2)-1 3*Ele(:,2)],[D -D],Ne,3*Nn);
    B = -B;
end

function [he,Rhe,Khe] = FoldKe(Cood, List, kpi, h0, L0, CM)

    rkj = [Cood(List(2),:)-Cood(List(1),:)]'; 
    rij = [Cood(List(3),:)-Cood(List(1),:)]'; 
    rkl = [Cood(List(2),:)-Cood(List(4),:)]'; 
    rmj = icross(rij,rkj); rnk = icross(rkj,rkl);

    sgn = ((abs(rnk'*rij)>1e-8)*sign(rnk'*rij)+(abs(rnk'*rij)<=1e-8)*1);
    he = real(acos(rmj'*rnk/(norm(rmj)*norm(rnk)))); 

    he = real(sgn*he);
    if he<0 
        he = 2*pi+he; 
    end;
    
    if nargin > 4     % RotSpring(he, h0, kpi, L0, limlft, limrht)
        [Rspr, Kspr] = CM(he,h0,kpi,L0);
    
        if nargout>1
            di = norm(rkj)/(rmj'*rmj)*rmj;
            dl = -norm(rkj)/(rnk'*rnk)*rnk;
            dj = (rij'*rkj/(rkj'*rkj)-1)*di-rkl'*rkj/(rkj'*rkj)*dl;
            dk = -rij'*rkj/(rkj'*rkj)*di+(rkl'*rkj/(rkj'*rkj)-1)*dl;
        
            Jhe = [dj;dk;di;dl];
            Rhe = Rspr*Jhe;
        end
        
        if nargout > 2
            dii = -norm(rkj)/(rmj'*rmj)^2*((rmj*cross(rkj,rmj)')+(rmj*cross(rkj,rmj)')');
            
            dtempij = -norm(rkj)/(rmj'*rmj)^2*(rmj*(cross(rij-rkj,rmj))'+(cross(rij-rkj,rmj))*rmj');
            dij = -rmj*rkj'/(rmj'*rmj*norm(rkj))+dtempij;
            
            dtempik = norm(rkj)/(rmj'*rmj)^2*(rmj*(cross(rij,rmj))'+(cross(rij,rmj))*rmj');
            dik = rmj*rkj'/(rmj'*rmj*norm(rkj))+dtempik;
            
            dil = zeros(3);
            
            dll = norm(rkj)/(rnk'*rnk)^2*(rnk*cross(rkj,rnk)'+(rnk*cross(rkj,rnk)')');
            
            dtemplk = norm(rkj)/(rnk'*rnk)^2*(rnk*(cross(rkl-rkj,rnk))'+(cross(rkl-rkj,rnk))*rnk');
            dlk = -rnk*rkj'/(rnk'*rnk*norm(rkj))+dtemplk;
            
            dtemplj = norm(rkj)/(rnk'*rnk)^2*(rnk*(cross(rnk,rkl))'+(rnk*(cross(rnk,rkl))')');
            dlj = rnk*rkj'/(rnk'*rnk*norm(rkj))+dtemplj;
            
            dT1jj = 1/(rkj'*rkj)*((-1+2*rij'*rkj/(rkj'*rkj))*rkj-rij);
            dT2jj = 1/(rkj'*rkj)*(2*rkl'*rkj/(rkj'*rkj)*rkj-rkl);
            djj = di*dT1jj'+(rij'*rkj/(rkj'*rkj)-1)*dij-(dl*dT2jj'+rkl'*rkj/(rkj'*rkj)*dlj);
            
            dT1jk = 1/(rkj'*rkj)*(-2*rij'*rkj/(rkj'*rkj)*rkj+rij);
            dT2jk = 1/(rkj'*rkj)*((1-2*rkl'*rkj/(rkj'*rkj))*rkj+rkl);
            djk = di*dT1jk'+(rij'*rkj/(rkj'*rkj)-1)*dik-(dl*dT2jk'+rkl'*rkj/(rkj'*rkj)*dlk);
            
            dT1kk = dT2jk;
            dT2kk = dT1jk;
            dkk = dl*dT1kk'+(rkl'*rkj/(rkj'*rkj)-1)*dlk-(di*dT2kk'+rij'*rkj/(rkj'*rkj)*dik);
            
            Hp = [ djj , djk , dij', dlj';
                   djk', dkk , dik', dlk';
                   dij , dik , dii , dil ;
                   dlj , dlk , dil', dll];
                  
                                          
            Khe = (Kspr*(Jhe*Jhe')+Rspr*Hp);
        end
    end
end


function c = icross(a,b)
    c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
         a(3,:).*b(1,:)-a(1,:).*b(3,:)
         a(1,:).*b(2,:)-a(2,:).*b(1,:)];

end

% [~,Rbe,Kbe] = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel));
function [Ex,Rbe,Kbe] = BarKe(u,B,L,CM,A)

    du = u(1:3)-u(4:6);
    Du = [du;-du];
    Ex = (B*u/L+0.5*(du'*du)/L^2);
    [Sx,Et] = CM(Ex);  
    Fx = Sx*A;
    if nargout>1
        Rbe = Fx*(B'+Du/L); 
    end

    if nargout>2
        Kel = B'*B;
        Kg = Fx/L*[eye(3),-eye(3);-eye(3),eye(3)];
        K1 = ((Du*B)+(Du*B)')/L;
        K2 = (Du*Du')/L^2;
        Kbe = Et*A/L*(Kel+K1+K2)+Kg;
    %     Kbe = HA/L*(Kel)+Kg;
    end
end


function [Sx, Et, Wb] = Ogden(Ex, C0)
    % Ogden hyperelastic constitutive model for bar elements
    alfa = [5,1]; % Specify parameteres
    pstr = real(sqrt(2*Ex+1));
    C0 = (pstr<1)*1*C0+(~(pstr<1))*1*C0;
    Et = C0/(alfa(1)-alfa(2)).*((alfa(1)-2)*pstr.^(alfa(1)-4)-(alfa(2)-2)*pstr.^(alfa(2)-4));
    Sx = C0/(alfa(1)-alfa(2)).*(pstr.^(alfa(1)-2)-pstr.^(alfa(2)-2));
    if nargout>2
        Wb = C0/(alfa(1)-alfa(2)).*((pstr.^alfa(1)-1)/alfa(1)-(pstr.^alfa(2)-1)/alfa(2));
    end
end

function [Rspr, Kspr, Espr] = EnhancedLinear(he,h0,kpi,L0,limlft,limrht)
    limlft = limlft/180*pi; partl = pi/limlft;
    limrht = limrht/180*pi; partr = pi/(2*pi-limrht);
    % limlft: theta_1: left partition point
    % limrht: theta_2: right partition point
    if numel(kpi)==1, kpi = kpi*ones(size(he)); end;
    Rspr = zeros(size(he)); Kspr = Rspr; 
    
    Lind = he<limlft; Rind = he>limrht; Mind = ~(Lind|Rind);
    Rspr(Lind) = kpi(Lind).*real(limlft-h0(Lind))+kpi(Lind).*tan(partl/2*(he(Lind)-limlft))/(partl/2);
    Kspr(Lind) = kpi(Lind).*sec(partl/2*(he(Lind)-limlft)).^2;
    Rspr(Rind) = kpi(Rind).*real(limrht-h0(Rind))+kpi(Rind).*tan(partr/2*(he(Rind)-limrht))/(partr/2);
    Kspr(Rind) = kpi(Rind).*sec(partr/2*(he(Rind)-limrht)).^2;
    Rspr(Mind) = kpi(Mind).*real(he(Mind)-h0(Mind));
    Kspr(Mind) = kpi(Mind);
    Rspr = L0.*Rspr; Kspr = L0.*Kspr;
    
    if nargout>2
        Espr = zeros(size(he));
        Espr(Lind) = 0.5*kpi(Lind).*real(h0(Lind)-limlft).^2+kpi(Lind).*real(h0(Lind)-limlft).*(limlft-he(Lind))-4*kpi(Lind)/partl^2.*log(abs(cos(partl/2*(limlft-he(Lind)))));
        Espr(Rind) = 0.5*kpi(Rind).*real(limrht-h0(Rind)).^2+kpi(Rind).*real(limrht-h0(Rind)).*(he(Rind)-limrht)-4*kpi(Rind)/partr^2.*log(abs(cos(partr/2*(he(Rind)-limrht))));
        Espr(Mind) = 0.5*kpi(Mind).*real(he(Mind)-h0(Mind)).^2;
        Espr = L0.*Espr;
    end
end

function [Exbar,FdAngle,LFd] = GetData(Ui,Node,truss,angles) %[Exbar,FdAngle,BdAngle,LFd,LBd] = GetData(Ui,Node,truss,angles) % [Exbar] = GetData(Ui,Node,truss, angles) 
    Exbar = zeros(size(truss.Bars,1),1) ; 
    FdAngle = zeros(size(angles.fold,1),1) ; LFd = FdAngle ;
    % BdAngle = zeros(size(angles.bend,1),1) ; LBd = BdAngle ;
    Nodenw = Node ;
    Nodenw(:,1) = Node(:,1)+Ui(1:3:end) ;
    Nodenw(:,2) = Node(:,2)+Ui(2:3:end) ;
    Nodenw(:,3) = Node(:,3)+Ui(3:3:end) ;

    for bel = 1:size(truss.Bars,1) 
        eDof = [(-2:0)+(truss.Bars(bel,1)*3),(-2:0)+(truss.Bars(bel,2)*3)]' ;
        Exbar(bel) = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel)) ;
    end
    
    % for del = 1:size(angles.bend,1)
    %     bend = angles.bend(del,:) ;
    %     BdAngle(del) = FoldKe(Nodenw,bend,angles.kpb,angles.pb0(del)) ;
    %     LBd(del) = norm(Nodenw(bend(2),:)-Nodenw(bend(1),:)) ;
    % end ;
    
    for fel = 1:size(angles.fold,1)
        fold = angles.fold(fel,:);
        FdAngle(fel) = FoldKe(Nodenw,fold,angles.kpf,angles.pf0(fel));
        LFd(fel) = norm(Nodenw(fold(2),:)-Nodenw(fold(1),:));
    end
end

%% Plotting

function VisualFold(U_his,truss,angles,recordtype,filename,pausetime,LF_his,instdof,axislim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record the simulation to file if needed:                                %
    % recordtype = 'none': do not save the simulation                         %
    % recordtype = 'video': save simulation in MP4 format;                    %
    % recordtype = 'imggif': save simulatrion in GIF format.                  %
    % pausement: pause time between each frame (in seconds).                  %
    %            If recordtype = 'none': use a small number such as 0.0001;   %
    %            Otherwise, pausetime = 1/fps;                                %
    % If input data does not include 'load_his', 'instdof', 'axislim', the    %
    % function does not plot load vs. displacement diagram.                   %
    % axislim: Axis limits (bounding box) for load vs. displacement diagram.  %
    %          format: [xmin,xmax,ymin,ymax].                                 %
    % instdof: specify the DOF of interest for displacement measure.          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Node = truss.Node; Trigl = truss.Trigl;
    Panel = angles.Panel;
    if nargin<7
        if strcmp(recordtype,'video')
            opengl('hardware')
            writerObj = VideoWriter(filename,'MPEG-4');
            writerObj.FrameRate = 1/pausetime;
            open(writerObj);
        elseif strcmp(recordtype,'imggif')
            filename = [filename '.gif'];
        else
            disp('Not recording');
        end;
        f1 = figure('units','pixels');
        f1.Color = 'w';
        hold on
        for i = 1:size(U_his,2)
            U = U_his(:,i);
            clf
    %         view(10,30)
            view(35,30)
            Nodew = Node;
            Nodew(:,1) = Node(:,1)+U(1:3:end);
            Nodew(:,2) = Node(:,2)+U(2:3:end);
            Nodew(:,3) = Node(:,3)+U(3:3:end);
    
            PlotOri(Node,Panel,Trigl,'-',0.3,'none');
            PlotOri(Nodew,Panel,Trigl);
            axis equal; axis off; 
            light
            pause(pausetime);
            if strcmp(recordtype,'imggif')
                frame = getframe(f1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);                
                if i == 1;
                    imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', pausetime);
                end  
            elseif strcmp(recordtype,'video')
                frame = getframe(f1);
                writeVideo(writerObj,frame);
            end
        end
        hold off;
        if strcmp(recordtype,'video'), close(writerObj); end
        
    else
        if size(LF_his,2)>1, LF_his = sum(LF_his,2); end;
        if strcmp(recordtype,'video')
            opengl('hardware')
            writerObj = VideoWriter(filename,'MPEG-4');
            writerObj.FrameRate = 1/pausetime;
            open(writerObj);
        elseif strcmp(recordtype,'imggif')
            filename = [filename '.gif'];
        else
            disp('Not recording');
        end;
        f1 = figure('units','pixels');
        f1.Color = 'w';
        hold on
        for i = 1:size(U_his,2)
            U = U_his(:,i);
            clf
            view(10,30)
            Nodew = Node;
            Nodew(:,1) = Node(:,1)+U(1:3:end);
            Nodew(:,2) = Node(:,2)+U(2:3:end);
            Nodew(:,3) = Node(:,3)+U(3:3:end);
    
            PlotOri(Node,Panel,Trigl,'-',0.3,'none');
            PlotOri(Nodew,Panel,Trigl);
            axis equal; 
            pause(pausetime);
            if strcmp(recordtype,'imggif')
                frame = getframe(f1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);                
                if i == 1;
                    imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', pausetime);
                end  
            elseif strcmp(recordtype,'video')
                frame = getframe(f1);
                writeVideo(writerObj,frame);
            end
        end
        hold off
        if strcmp(recordtype,'video'), close(writerObj); end
        
        if strcmp(recordtype,'video')
            opengl('hardware')
            writerObj = VideoWriter([filename '_dispvslambda'],'MPEG-4');
            writerObj.FrameRate = 1/pausetime;
            open(writerObj);
        elseif strcmp(recordtype,'imggif')
            filename = [filename 'dispvslambda' '.gif'];
        end;
        f2 = figure('units','pixels');
        f2.Color = 'w';
        dsp = sign(instdof)*U_his(abs(instdof),:);
        for i = 1:numel(LF_his)
            clf;
            plot(dsp(1:i),LF_his(1:i),'b-','linewidth',2);
            hold on;
            plot(dsp(i),LF_his(i),'ro','linewidth',2);
            axis(axislim);
            xlabel('displacement','fontsize',14)
            ylabel('load factor','fontsize',14)
            pause(pausetime)
            if strcmp(recordtype,'imggif')
                frame = getframe(f2);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);                
                if i == 1;
                    imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', pausetime);
                end  
            elseif strcmp(recordtype,'video')
                frame = getframe(f2);
                writeVideo(writerObj,frame);
            end
        end
        hold off
        if strcmp(recordtype,'video'), close(writerObj); end
    end
end
