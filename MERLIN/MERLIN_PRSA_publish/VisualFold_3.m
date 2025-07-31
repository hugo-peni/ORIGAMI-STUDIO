function VisualFold_3(U_his,truss,angles,pausetime,LF_his)

    Node = truss.Node; Trigl = truss.Trigl;
    Panel = angles.Panel;

    if size(LF_his,2)>1, LF_his = sum(LF_his,2); end;
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
        end