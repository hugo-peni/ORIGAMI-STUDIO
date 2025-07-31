clear all ; close all ; clc ;

%% Define geomtry and material
N = 8; h = 1; lyr = 4; phi = -pi/8;
% N-gon prism; height of each layer; number of inter-layer planes; twisting
% angle of the prism in each layer
MaxIcr = 180; blam = 0.032; 
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
npp = numel(indp) ;
Load = [ indp , 0*ones(npp,1) , 0*ones(npp,1) , -1*ones(npp,1) ; ] ;

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1) ;
%[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr) ;
%U_his = real(U_his);
%LF_his = real(LF_his);

figure()
PlotOri(Node,angles.Panel,truss.Trigl);
axis equal; axis off;
camproj('perspective')
light
view(117,18)
rotate3d on

% PrepareData 
% PathAnalysis 


%% function : 

function [truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarCM,RotSpring,kpf,kpb,Abar)

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


function [fold, bdry, Trigl] = findfdbd(Panel,bend)

    Nn = max(cellfun(@max,Panel)); % Panel : constains list of nodes making one panel
    % triangularization
    Panelsize = cellfun(@numel,Panel);
    Ptri = cell(sum(Panelsize==3),1);
    flg = find(Panelsize==3);

    for i = 1:sum(Panelsize==3), Ptri{i} = Panel{flg(i)}; end

    Triglraw = [bend(:,[1,2,3]);bend(:,[1,2,4]);cell2mat(Ptri)];
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

    fdandbd = sort(fold(:,1:2),2);
    onlybd = sort(bend(:,1:2),2);
    [~,ibd] = intersect(fdandbd,onlybd,'rows');
    fold(ibd,:) = [];
    
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

