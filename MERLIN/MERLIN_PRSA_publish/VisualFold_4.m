function VisualFold_4(U_his,truss,angles, pausetime,LF_his,instdof,axislim)

    Node = truss.Node; Trigl = truss.Trigl;
    Panel = angles.Panel;

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
    end
    hold off