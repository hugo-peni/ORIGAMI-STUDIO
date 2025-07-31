
clear all ; close all ; clc ;

%% Define geomtry and material
N = 8; h = 0.5; lyr = 10; phi = -pi/8;
% N-gon prism; height of each layer; number of inter-layer planes; twisting
% angle of the prism in each layer
MaxIcr = 180; blam = 0.032; 
% Maximum increment number & initial load factor
Kf = 1e-7; Kb = Kf; E0 = 5e4; Abar = 0.1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 5; limrht = 355;
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
        indsupp(3:end), ones(nsupp-2,1), ones(nsupp-2,1), ones(nsupp-2,1);];

Total_Displacement = -4 ; 
ForcedDisplacement = Total_Displacement / MaxIcr ; 

if( ForcedDisplacement ~= 0)
    use_displacement_control = 0 ; 
    LoadMagnitude = -0.000 ;
    indPrescribed = find( abs( Node(:,3) - max(Node(:,3 ) ) ) < 1e-5 ) ;
    nPrescribed = numel(indPrescribed);
    PrescribedDisplacement1 = [          indPrescribed(1), 1e-10, 1e-10, ForcedDisplacement;
                      indPrescribed(2), 1e-10, 1e-10, ForcedDisplacement;
            indPrescribed(3:end), ones(nPrescribed-2,1)*1e-10, ones(nPrescribed-2,1)*1e-10, ForcedDisplacement*ones(nPrescribed-2,1);];

    indPrescribed = find( abs( Node(:,3) - 0.5  ) < 1e-5 ) ;
    nPrescribed = numel(indPrescribed);
    PrescribedDisplacement2 = [          indPrescribed(1), 1e-10, 1e-10, 0;
                      indPrescribed(2), 1e-10, 1e-10, 0;
            indPrescribed(3:end), ones(nPrescribed-2,1)*1e-10, ones(nPrescribed-2,1)*1e-10, 0*ones(nPrescribed-2,1);];
    
    PrescribedDisplacement5 = [PrescribedDisplacement1; PrescribedDisplacement2] ;
    for i = 1:0.5:3.5
        indPrescribed = find( abs( Node(:,3) - i  ) < 1e-5 ) ;
        nPrescribed = numel(indPrescribed);
        PrescribedDisplacement3 = [          indPrescribed(1), 1e-10, 1e-10, 0;
                      indPrescribed(2), 1e-10, 1e-10, 0;
            indPrescribed(3:end), ones(nPrescribed-2,1)*1e-10, ones(nPrescribed-2,1)*1e-10, 0*ones(nPrescribed-2,1);];

        PrescribedDisplacement = [PrescribedDisplacement5; PrescribedDisplacement3] ; 
        PrescribedDisplacement5 = PrescribedDisplacement ; 
    end
else
    use_displacement_control = 1 ; 
    PrescribedDisplacement = [] ; 
    LoadMagnitude = -0.4885 ; %-0.0001 ;
end 



m = size(Node,1);
indp = find( abs( Node(:,3) - max(Node(:,3 ) ) ) < 1e-5 ) ; %indp = indp(3) ;
npp = numel(indp) ;
Load = [ indp , 0*ones(npp,1) , 0*ones(npp,1) , LoadMagnitude*ones(npp,1) ; ] ;

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,PrescribedDisplacement, Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1) ;
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr,use_displacement_control) ;
U_his = real(U_his);
LF_his = real(LF_his);

% figure()
% PlotOri(Node,angles.Panel,truss.Trigl);
% axis equal; axis off;
% camproj('perspective')
% light
% view(117,18)
% rotate3d on

% PrepareData 
% PathAnalysis 

%% Visualize simulation

instdof = -indp(1)*3;

interv = 1; endicrm = size(U_his,2);
%VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001,LF_his,instdof,[-inf inf -inf inf])
VisualFold_3(U_his(:,1:interv:endicrm),truss,angles,0.0001, LF_his)
% 
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001)

% VisualFold,


%% Plot stored energy vs. pseudo time

% Red line is the total profile. Between red and cyan is the folding
% energy. Between cyan and magenta is the portion of energy for bending. 
% Below magenta is the stretching energy of bars.

% STAT = PostProcess(Data,truss,angles);  
% figure()
% plot(1:size(U_his,2),STAT.PE.strain,'r-','linewidth',2);
% grid on
% hold on
% plot(1:size(U_his,2),STAT.bend.PE+STAT.bar.PE,'c-');
% plot(1:size(U_his,2),STAT.bar.PE,'m-');
% xlabel('Increment Number (Pseudo-time)','fontsize',14);
% ylabel('Stored Energy','fontsize',14);

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

function [truss, angles, F] = PrepareData(Node,Panel, PrescribedDisplacement , Supp , Load,BarCM,RotSpring,kpf,kpb,Abar)

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

    if size(PrescribedDisplacement,1) == 0
        Prescribed_rs = []; 
    else
        Prescribed_rs = [reshape([PrescribedDisplacement(:,1)*3-2,PrescribedDisplacement(:,1)*3-1,PrescribedDisplacement(:,1)*3]',[],1),...
              reshape(PrescribedDisplacement(:,2:4)',[],1)];
        Prescribed_rs(Prescribed_rs(:,2)==0,:)=[]; Prescribed_rs = Prescribed_rs(:,1);
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

    Displacement_field = zeros(3*m,1) ;  
   
    if( size(PrescribedDisplacement, 1 ) > 0 )
        indPrescribedDisplacement = PrescribedDisplacement(:,1);
        Displacement_field(3*indPrescribedDisplacement-2) = PrescribedDisplacement(:,2);
        Displacement_field(3*indPrescribedDisplacement-1) = PrescribedDisplacement(:,3);
        Displacement_field(3*indPrescribedDisplacement) = PrescribedDisplacement(:,4);
    end

    truss.CM = BarCM ; 
    truss.Node = Node ;
    truss.Bars = Bars ;
    truss.Trigl = Trigl ;
    truss.B = B ; 
    truss.L = L ;
    truss.FixedDofs = unique(rs) ;
    truss.PrescribedDofs = unique(Prescribed_rs) ;
    truss.A = Abar ; 
    truss.PrescribedDisplacement = Displacement_field ; 


    angles.CM = RotSpring ;
    angles.fold = Fold ;
    angles.bend = Bend ;
    angles.kpf = kpf*ones(1,size(Fold,1)) ;
    angles.kpb = kpb*ones(1,size(Bend,1)) ; 
    angles.pf0 = pf0' ;
    angles.pb0 = pb0'*0+pi ;
    angles.Panel = Panel ;


end 

function [Uhis,load_his,Data] = PathAnalysis(truss,angles,F, b_lambda,MaxIcr, use_displacement_control)

    tol = 1e-6; MaxIter = 100 ; 
    Node = truss.Node ;

    AllDofs = [ 1:3*size(Node,1) ] ;

    U = truss.U0;

    Uhis = zeros(3*size(Node,1),MaxIcr) ;
    Data.Exbar = zeros(size(truss.Bars,1),MaxIcr) ; 
    Data.FdAngle = zeros(size(angles.fold,1),MaxIcr) ; 
    Data.LFdAngle = Data.FdAngle ;
    Data.BdAngle = zeros(size(angles.bend,1),MaxIcr) ; 
    Data.LBdAngle = Data.BdAngle ;
    Data.BarForces = zeros(size(truss.Bars,1), MaxIcr);      % Bar internal forces
    Data.TotalReactionForce = zeros(3, MaxIcr);

    FixedDofs = [ truss.FixedDofs ; truss.PrescribedDofs ] ; 
    FreeDofs = setdiff(AllDofs,FixedDofs) ;

    lmd = 0; icrm = 0; MUL = [U,U] ;
    dUt = U ; 
    load_his = zeros(MaxIcr,1) ;

    F = F(:,1) ;
    
    
    while icrm<MaxIcr

        icrm = icrm+1 ;
        iter = 0; err = 1 ;
        fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd) ;

        U = U + truss.PrescribedDisplacement ; 

        while err>tol && iter<MaxIter
            iter = iter+1 ; 
            % There are two versions of the function for assembling global
            % stiffness matrix:
            % 'GlobalK_fast_ver': speed-optimized version;
            % 'GlobalK_edu_ver': easy-to-read educational version.
            % [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
            [IF,K] = GlobalK_edu_ver(U,Node,truss,angles) ;

            if use_displacement_control == 1 
                R = lmd*F-IF;   MRS = [F,R] ;
                if( icrm == MaxIcr-1 )
                    disp('hello world')
                end
                MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:) ;
                dUp = MUL(:,1); dUr = MUL(:,2) ;
                if iter==1, dUr = 0*dUr; end
                dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda) ;
                dUt = dlmd*dUp+dUr ; 
            else
                R = -IF ; 
                dUt(FreeDofs,:) = K(FreeDofs, FreeDofs) \ R(FreeDofs, :) ;
            end

            U = U+dUt ;

            err = norm(dUt(FreeDofs)) ; 

            if use_displacement_control == 1 
                lmd = lmd + dlmd ;
                fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n', iter , err , dlmd ) ;
            else 
                fprintf('    iter = %d, err = %6.4f\n', iter , err ) ;
            end 

            if err > 1e8, disp('Divergence!'); break; end
        end


        if iter>15
            if use_displacement_control == 1
                disp('Reduce constraint radius!')
                icrm = icrm-1;
                U = Uhis(:,max(icrm,1));
                b_lambda = b_lambda/2;
                lmd = load_his(max(icrm,1));
            end
        elseif iter<3
            disp('Increase constraint radius!')
            Uhis(:,icrm) = U;
            if use_displacement_control == 1
                b_lambda = b_lambda*1.5;
                load_his(icrm) = lmd;
            end
            [Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei, BarForcesi, TotalReactionForcei, ReactionForcesi] = GetData(U,Node,truss,angles);
            Data.Exbar(:,icrm) = Exbari; 
            Data.FdAngle(:,icrm) = FdAnglei; 
            Data.LFdAngle(:,icrm) = LFdAnglei;
            Data.BdAngle(:,icrm) = BdAnglei;
            Data.LBdAngle(:,icrm) = LBdAnglei;
            Data.BarForces(:,icrm) = BarForcesi ; 
            Data.TotalReactionForce(:,icrm) = TotalReactionForcei ; 
        else
            Uhis(:,icrm) = U;
            if use_displacement_control == 1
                load_his(icrm) = lmd;
            end
            [Exbari,FdAnglei,BdAnglei,LFdAnglei,LBdAnglei, BarForcesi, TotalReactionForcei, ReactionForcesi] = GetData(U,Node,truss,angles);
            Data.Exbar(:,icrm) = Exbari; 
            Data.FdAngle(:,icrm) = FdAnglei; 
            Data.LFdAngle(:,icrm) = LFdAnglei;
            Data.BdAngle(:,icrm) = BdAnglei;
            Data.LBdAngle(:,icrm) = LBdAnglei;
            Data.BarForces(:,icrm) = BarForcesi ; 
            Data.TotalReactionForce(:,icrm) = TotalReactionForcei ; 
        end

    end
   
    icrm = icrm+1;
    Uhis(:,icrm:end) = [];

    load_his(icrm:end,:) = [];
    Data.Exbar(:,icrm:end) = []; 

    Data.FdAngle(:,icrm:end) = []; Data.LFdAngle(:,icrm:end) = [];
    Data.BdAngle(:,icrm:end) = []; Data.LBdAngle(:,icrm:end) = [];
    Data.BarForces(:,icrm:end) =[] ; 
    Data.TotalReactionForce(:,icrm:end) = [] ; 

    if icrm == 100 
        disp('hugo')
    end

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

function [Exbar,FdAngle,BdAngle,LFd,LBd, BarForces, TotalReactionForce, ReactionForces] = GetData(Ui,Node,truss,angles)

    Exbar = zeros(size(truss.Bars,1),1); 
    FdAngle = zeros(size(angles.fold,1),1); 
    LFd = FdAngle;
    BdAngle = zeros(size(angles.bend,1),1); 
    LBd = BdAngle;

    Nodenw = Node;
    Nodenw(:,1) = Node(:,1)+Ui(1:3:end);
    Nodenw(:,2) = Node(:,2)+Ui(2:3:end);
    Nodenw(:,3) = Node(:,3)+Ui(3:3:end);

    for bel = 1:size(truss.Bars,1) 
        eDof = [(-2:0)+(truss.Bars(bel,1)*3),(-2:0)+(truss.Bars(bel,2)*3)]';
        Exbar(bel) = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel));
    end
    
    for del = 1:size(angles.bend,1)
        bend = angles.bend(del,:);
        BdAngle(del) = FoldKe(Nodenw,bend,angles.kpb,angles.pb0(del));
        LBd(del) = norm(Nodenw(bend(2),:)-Nodenw(bend(1),:));
    end
    
    for fel = 1:size(angles.fold,1)
        fold = angles.fold(fel,:);
        FdAngle(fel) = FoldKe(Nodenw,fold,angles.kpf,angles.pf0(fel));
        LFd(fel) = norm(Nodenw(fold(2),:)-Nodenw(fold(1),:));
    end

    % NEW: Calculate bar internal forces
    BarForces = zeros(size(truss.Bars,1),1);
    for i = 1:size(truss.Bars,1)
        eDof = [(-2:0)+(truss.Bars(i,1)*3), (-2:0)+(truss.Bars(i,2)*3)]';
        [~, Rbe, ~] = BarKe(Ui(eDof), truss.B(i,eDof), truss.L(i), truss.CM, truss.A(i));
        
        % Calculate axial force magnitude (force along bar direction)
        bar_direction = (Nodenw(truss.Bars(i,2),:) - Nodenw(truss.Bars(i,1),:))';
        bar_direction = bar_direction / norm(bar_direction);
        
        % Extract nodal forces and project onto bar direction
        F_node1 = Rbe(1:3);
        axial_force = dot(F_node1, bar_direction);
        BarForces(i) = axial_force;
    end

    % Calculate reaction forces
    [IF, ~] = GlobalK_edu_ver(Ui, Node, truss, angles);
    ReactionForces = IF(truss.PrescribedDofs);
    
    % Calculate total reaction force components
    total_fx = 0; total_fy = 0; total_fz = 0;
    
    for i = 1:length(truss.PrescribedDofs)
        dof = truss.PrescribedDofs(i);
        direction = mod(dof-1, 3) + 1; % 1=x, 2=y, 3=z
        
        if direction == 1
            total_fx = total_fx + ReactionForces(i);
        elseif direction == 2
            total_fy = total_fy + ReactionForces(i);
        elseif direction == 3
            total_fz = total_fz + ReactionForces(i);
        end
    end
    
    TotalReactionForce = [total_fx, total_fy, total_fz];

end

function [TotalReactionForce, ReactionForces] = GetReactionForces(U, Node, truss, angles)
    % Calculate the global internal force vector
    [IF, ~] = GlobalK_edu_ver(U, Node, truss, angles);
    
    % Extract reaction forces at prescribed displacement DOFs
    ReactionForces = IF(truss.PrescribedDofs);
    
    % Calculate total force in each direction
    % Assuming prescribed displacements are in Z-direction (vertical compression)
    total_fx = 0; total_fy = 0; total_fz = 0;
    
    for i = 1:length(truss.PrescribedDofs)
        dof = truss.PrescribedDofs(i);
        direction = mod(dof-1, 3) + 1; % 1=x, 2=y, 3=z
        
        if direction == 1
            total_fx = total_fx + ReactionForces(i);
        elseif direction == 2
            total_fy = total_fy + ReactionForces(i);
        elseif direction == 3
            total_fz = total_fz + ReactionForces(i);
        end
    end
    
    TotalReactionForce = [total_fx, total_fy, total_fz];
end

function [Ex,Rbe,Kbe] = BarKe(u,B,L,CM,A)

    du = u(1:3)-u(4:6);
    Du = [du;-du];
    Ex = (B*u/L+0.5*(du'*du)/L^2);
    [Sx,Et] = CM(Ex);  
    Fx = Sx*A;
    if nargout>1, Rbe = Fx*(B'+Du/L); end;
    if nargout>2
        Kel = B'*B;
        Kg = Fx/L*[eye(3),-eye(3);-eye(3),eye(3)];
        K1 = ((Du*B)+(Du*B)')/L;
        K2 = (Du*Du')/L^2;
        Kbe = Et*A/L*(Kel+K1+K2)+Kg;
    %     Kbe = HA/L*(Kel)+Kg;
    end
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
    end ;
    
    if nargin > 4
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
        dii = -norm(rkj)/(rmj'*rmj)^2*( ( rmj*cross(rkj,rmj)' ) + ( rmj*cross(rkj,rmj)' )' );
        
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
