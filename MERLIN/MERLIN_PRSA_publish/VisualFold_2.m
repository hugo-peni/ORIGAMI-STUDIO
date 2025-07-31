function VisualFold_2(U_his,truss,angles,recordtype,filename,pausetime,LF_his,instdof,axislim)

    Node = truss.Node; Trigl = truss.Trigl;
    Panel = angles.Panel;

    if size(LF_his,2)>1, LF_his = sum(LF_his,2); end;
        % if strcmp(recordtype,'video')
        %     opengl('hardware')
        %     writerObj = VideoWriter(filename,'MPEG-4');
        %     writerObj.FrameRate = 1/pausetime;
        %     open(writerObj);
        % elseif strcmp(recordtype,'imggif')
        %     filename = [filename '.gif'];
        % else
        %     disp('Not recording');
        % end;
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
            % if strcmp(recordtype,'imggif')
            %     frame = getframe(f1);
            %     im = frame2im(frame);
            %     [imind,cm] = rgb2ind(im,256);                
            %     if i == 1;
            %         imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
            %     else
            %         imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', pausetime);
            %     end  
            % elseif strcmp(recordtype,'video')
            %     frame = getframe(f1);
            %     writeVideo(writerObj,frame);
            % end
        end
        % hold off
        % if strcmp(recordtype,'video'), close(writerObj); end
        % 
        % if strcmp(recordtype,'video')
        %     opengl('hardware')
        %     writerObj = VideoWriter([filename '_dispvslambda'],'MPEG-4');
        %     writerObj.FrameRate = 1/pausetime;
        %     open(writerObj);
        % elseif strcmp(recordtype,'imggif')
        %     filename = [filename 'dispvslambda' '.gif'];
        % end;

        % f2 = figure('units','pixels');
        % f2.Color = 'w';
        % dsp = sign(instdof)*U_his(abs(instdof),:);
        % for i = 1:numel(LF_his)
        %     clf;
        %     plot(dsp(1:i),LF_his(1:i),'b-','linewidth',2);
        %     hold on;
        %     plot(dsp(i),LF_his(i),'ro','linewidth',2);
        %     axis(axislim);
        %     xlabel('displacement','fontsize',14)
        %     ylabel('load factor','fontsize',14)
        %     pause(pausetime)
        %     if strcmp(recordtype,'imggif')
        %         frame = getframe(f2);
        %         im = frame2im(frame);
        %         [imind,cm] = rgb2ind(im,256);                
        %         if i == 1;
        %             imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
        %         else
        %             imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', pausetime);
        %         end  
        %     elseif strcmp(recordtype,'video')
        %         frame = getframe(f2);
        %         writeVideo(writerObj,frame);
        %     end
        % end
        % hold off
        %if strcmp(recordtype,'video'), close(writerObj); end