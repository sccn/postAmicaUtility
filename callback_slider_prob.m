function callback_slider_prob(sliderHandle,~,sRate,epoched,epoch_length,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2,Linewidth)

figh = gcf;
currentFrame = sliderHandle.Value;
data = get(figh,'userdata');
windowlength = data{2};
data(2) = [];
data = data{1};
if gca ~= ax1
    axes(ax1);
end


obj = findobj('parent',figh,'Tag','ToggleEvent');
showEvent = get(obj,'Value');
if ~showEvent
    events2plot = {};
end

%childrenAx1 = get(ax1,'children');

if ~epoched
    
    T = ceil(windowlength);
    if sRate*(currentFrame + T)+1 > size(data,1)
%         for i = 1:size(data,2)-2
%             set(childrenAx1(end-size(data,2)+2+i),'ydata',data((sRate*currentFrame+1:end),i))
%             set(childrenAx1(end-size(data,2)+2+i),'xdata',data((sRate*currentFrame+1:end),end))
%         end
        
        plot(ax1,data((sRate*currentFrame+1:end),end),data((sRate*currentFrame+1:end),1:end-2),'linewidth',Linewidth)
        %[currentFrame ceil(size(data,1)/sRate)]
        set(ax1,'xlim',[currentFrame ceil(size(data,1)/sRate)]);
        xlim = get(ax1,'xlim');
        set(ax1,'xtick',linspace(xlim(1),xlim(2),11))
        
    else
%         for i = 1:size(data,2)-2
%             set(childrenAx1(end-size(data,2)+2+i),'ydata',data((sRate*currentFrame:sRate*(currentFrame + T))+1,i))
%             set(childrenAx1(end-size(data,2)+2+i),'xdata',data((sRate*currentFrame:sRate*(currentFrame + T))+1,end))
%         end
        plot(ax1,data((sRate*currentFrame:sRate*(currentFrame + T))+1,end),data((sRate*currentFrame:sRate*(currentFrame + T))+1,1:end-2),'linewidth',Linewidth)
        set(ax1,'xlim',[currentFrame currentFrame+T]);
        set(ax1,'xtick',linspace(currentFrame,currentFrame+T,11));
    end
    
    xlim = get(ax1,'xlim');
    set(ax1,'ylim',[-0.1 1.1]);
    
    ylabel(ax1,'Probability of Model Being Active')
    xlabel(ax1,fastif(epoched,'Trials','Time(sec)'))
    hold(ax1,'on');
    %draw events
    x_coord = [];
    if ~isempty(events2plot)
        count = 1;
        
        for i = 1:length(latency)
            if latency(i)>sRate*xlim(2)
                break;
            else
                if latency(i)>=currentFrame*sRate & sum(strcmp({events(i).type},events2plot))>0
                    x_coord(count) = latency(i)/sRate;
                    plot(ax1,[x_coord(count) x_coord(count)],[-0.1 1.1],'linestyle','- -')
                    if ~ischar(events(i).type)
                        text(x_coord(count),1.12,num2str(events(i).type),'rotation',90);
                    else
                        text(x_coord(count),1.12,events(i).type,'rotation',90);
                    end
                    count = count + 1;
                end
            end
        end
    end
    hold(ax1,'off')
    %childrenAxLikelihood = get(axLikelihood,'children');
    if sRate*(currentFrame + T)+1 > size(data,1)
%         set(childrenAxLikelihood(i),'ydata',data((sRate*currentFrame+1:end),end-1))
%         set(childrenAxLikelihood(i),'xdata',data((sRate*currentFrame+1:end),end))
        plot(axLikelihood,data((sRate*currentFrame+1:end),end),data((sRate*currentFrame+1:end),end-1),'linewidth',Linewidth)
        set(axLikelihood,'xlim',[currentFrame ceil(size(data,1)/sRate)]);
        
        set(axLikelihood,'xtick',linspace(xlim(1),xlim(2),11))
    else
%         set(childrenAxLikelihood,'ydata',data((sRate*currentFrame:sRate*(currentFrame + T))+1,end-1))
%         set(childrenAxLikelihood,'xdata',data((sRate*currentFrame:sRate*(currentFrame + T))+1,end))
        plot(axLikelihood,data((sRate*currentFrame:sRate*(currentFrame + T))+1,end),data((sRate*currentFrame:sRate*(currentFrame + T))+1,end-1),'linewidth',Linewidth);
        set(axLikelihood,'xlim',[currentFrame currentFrame+T]);
        set(axLikelihood,'xtick',linspace(currentFrame,currentFrame+T,11));
    end
    
    set(axLikelihood,'ylim',[ylim1 ylim2]);
    
    
    xlabel(axLikelihood,fastif(epoched,'Trials','Time(sec)'))
    hold(axLikelihood,'on');
    
    %draw events for second plot
    for i = 1:length(x_coord)
        plot(axLikelihood,[x_coord(i) x_coord(i)],[ylim1 ylim2],'linestyle','- -')
    end
    
    hold(axLikelihood,'off')
    if T ~= windowlength
        set(findobj('parent',gcf,'tag','WindowSizeTag'),'string',num2str(T))
    end
    
    sliderHandle.UnitIncrement = ceil(T/10);
    
    sliderHandle.Maximum = ceil(size(data,1)/sRate-T + ceil(T/10));
    sliderHandle.VisibleAmount = sliderHandle.UnitIncrement;
    
    
else
    
    T = ceil(windowlength);
    
    if epoch_length*(currentFrame + T) > size(data,1)
        
        currentFrame = size(data,1)/epoch_length - T;
        sliderHandle.Value = currentFrame;
    end
    
    plot(ax1,data((currentFrame*epoch_length+1:epoch_length*(currentFrame + T)),end),data((currentFrame*epoch_length+1:epoch_length*(currentFrame + T)),1:end-2),'linewidth',Linewidth)
    hold(ax1,'on');
    set(ax1,'ylim',[-0.1 1.1]);
    set(ax1,'xlim',[currentFrame+0.5 currentFrame+T+0.5]);
    lim = get(ax1,'xlim');
    set(ax1,'xtick',floor(lim(1):lim(2)))
    ylabel(ax1,'Probability of Model Being Active')
    xlabel(ax1,fastif(epoched,'Trials','Time(sec)'))
    
    %draw epoch seperating lines
    for i = currentFrame+1:currentFrame+T
        plot(ax1,[i+0.5 i+0.5],[-0.1 1.1],'color',[1 0.5 0.25],'linestyle','- -');
    end
    %-------------------------------
    windowInMs = 1000*windowInS;
    x_coord = [];
    if ~isempty(events2plot)
        % draw event lines
        
        C = 1/(windowInMs(2)-windowInMs(1));
        count = 1;
        
        for i = currentFrame+1:currentFrame+T
            lengthevents = length(events(i).event);
            
            for j = 1:lengthevents
                
                if lengthevents == 1
                    if iscell(events(i).eventtype)
                        str = num2str(events(i).eventtype{1});
                        eventlatency = events(i).eventlatency{1};
                    else
                        str = num2str(events(i).eventtype);
                        eventlatency = events(i).eventlatency;
                    end
                else
                    str = num2str(events(i).eventtype{j});
                    eventlatency = [events(i).eventlatency{j}];
                end
                
                if sum(strcmp(str,events2plot))>0
                    x_coord(count) = C*(eventlatency-windowInMs(1)) + (i - 0.5);
                    plot(ax1,[x_coord(count) x_coord(count)],[-0.1 1.1],'linestyle','- -')
                    text(x_coord(count),1.12,str,'rotation',90);
                    count = count + 1;
                end
            end
        end
        
    end
    hold(ax1,'off');
    
    plot(axLikelihood,data((currentFrame*epoch_length+1:epoch_length*(currentFrame + T)),end),data((currentFrame*epoch_length+1:epoch_length*(currentFrame + T)),end-1),'linewidth',Linewidth)
    set(axLikelihood,'ylim',[ylim1 ylim2]);
    set(axLikelihood,'xlim',[currentFrame+0.5 currentFrame+T+0.5]);
    set(axLikelihood,'xtick',floor(lim(1):lim(2)))
    xlabel(axLikelihood,fastif(epoched,'Trials','Time(sec)'))
    hold(axLikelihood,'on');
    
    for i = currentFrame+1:currentFrame+T
        plot(axLikelihood,[i+0.5 i+0.5],[ylim1 ylim2],'color',[1 0.5 0.25],'linestyle','- -');
    end
    
    % draw event lines
    
    for i = 1:length(x_coord)
        plot(axLikelihood,[x_coord(i) x_coord(i)],[ylim1 ylim2],'linestyle','- -');
    end
    
    hold(axLikelihood,'off');
    if T ~= windowlength
        set(findobj('parent',gcf,'tag','WindowSizeTag'),'string',num2str(T))
    end
    
    sliderHandle.UnitIncrement = ceil(T/5);
    sliderHandle.Maximum = size(data,1)/epoch_length-T + ceil(T/5);
    sliderHandle.VisibleAmount = sliderHandle.UnitIncrement;
    
end
ylabel(axLikelihood,'Log-likelihood of data under most probable model')






end