function scrollplot_fast(y,L,models2plot,windowlength,sRate,events, windowInS, events2plot,smoothlength)
if nargin<9
    smoothlength = [];
end
if nargin<6
    events = {};
    events2plot = {};
end

epoched = (size(y,3) > 1);

if nargin<4 && epoched
    sRate = [];
end

if size(y,1) ~= length(models2plot)
    y = y(models2plot,:,:);
end

figure;
figh = gcf;
set(gcf,'position',[634 94 898 1002]);
set(gcf,'Toolbar','auto')
set(gcf,'name','AMICA model probabilities')
try, icadefs;
catch
    BACKCOLOR = [.93 .96 1];
end;
set(gcf,'color',BACKCOLOR);
Linewidth = 2;

ax = gca;
pos = get(gca,'Position');
set(gca,'Position',[pos(1)-0.05 pos(2)+0.33 pos(3)-0.1 pos(4)-0.35]);

pos = get(gca,'Position');
sliderPosition = [pos(1) 0.03 pos(3) 0.02];
%sliderPosition = [pos(1) pos(2)-0.12 pos(3) 0.03];

axLikelihood = axes('Position',[sliderPosition(1) sliderPosition(2)+sliderPosition(4)+0.05 pos(3) pos(2)-(sliderPosition(2)+sliderPosition(4)+0.05)-0.06]);
axes(ax);
ylim1  = min(min(L))-1;
ylim2 = max(max(L))+1;
latency = [];
if ~epoched
    
    T = windowlength; %number of seconds to be displayed in each frame.
    t = 1:size(y,2);
    t = (t - 1)/sRate;
    set(figh,'userdata',{[y' L' t'] T})
    %set(figh,'userdata',{[y;L;t] T})
    plot(t(1:sRate*T),y(:,1:sRate*T),'linewidth',Linewidth);
    
    set(gca,'xlim',[0 T]);
    set(gca,'xtick',0:2:T)
    
    
    sliderHandle = uicontrol('style','slider',...
        'units','normalized','position',sliderPosition,...
        'min',0,'max',size(y,2));
    
    
    
    jScrollBar = findjobj(sliderHandle);
    jScrollBar.Minimum = 0;
    jScrollBar.UnitIncrement = T/10;
    
    jScrollBar.Maximum = ceil(size(y,2)/sRate-T + ceil(T/10));
    jScrollBar.VisibleAmount = jScrollBar.UnitIncrement;
    
    size2 = 0;
    
    %draw event lines
    currentFrame = 0;
    x_coord = [];
    count = 1;
    if ~isempty(events2plot)
        latency = [events.latency];
        for i = 1:length(latency)
            if latency(i)>sRate*(currentFrame + T)
                break;
            else
                if latency(i)>=currentFrame*sRate & sum(strcmp({events(i).type},events2plot))>0
                    x_coord(count) = latency(i)/sRate;
                    hl = line([x_coord(count) x_coord(count)],[-0.1 1.1]);
                    set(hl,'LineStyle','- -');
                    ht = text(x_coord(count),1.12,num2str(events(i).type));
                    set(ht,'rotation',90);
                    count = count + 1;
                end
            end
        end
    end
    
    %----------------
    
    plot(axLikelihood,t(1:sRate*T),L(:,1:sRate*T),'linewidth',Linewidth);
    set(axLikelihood,'ylim',[ylim1 ylim2]);
    set(axLikelihood,'xlim',[0 T]);
    set(axLikelihood,'xtick',0:2:T);
    xlabel(axLikelihood,fastif(epoched,'Trials','Time(sec)'))
    hold(axLikelihood,'on');
    
    %draw events
    
    for i = 1:length(x_coord)
        plot(axLikelihood,[x_coord(i) x_coord(i)],[ylim1 ylim2],'linestyle','- -')
    end
    hold(axLikelihood,'off')
    
else
    
    
    T = windowlength; %number of epochs to be shown
    size2 = size(y,2);
    size3 = size(y,3);
    t = 1:size2*size3;
    t = (t + size2-ceil(size2/2))/size2;
    y = reshape(y,size(y,1),size2*size3);
    L = reshape(L,size(L,1),size2*size3);
    set(figh,'userdata',{[y' L' t'] T})
    plot(t(1:T*size2),y(:,1:T*size2),'linewidth',Linewidth);
    hold on;
    set(gca,'xlim',[0.5 T+0.5]);
    lim = get(gca,'xlim');
    set(gca,'xtick',floor(lim(1):lim(2)))
    
    
    sliderHandle = uicontrol('style','slider',...
        'units','normalized','position',sliderPosition,...
        'min',0,'max',size(y,2));
    
    
    jScrollBar = findjobj(sliderHandle);
    jScrollBar.Minimum = 0;
    jScrollBar.UnitIncrement = 1;
    jScrollBar.Maximum = size3-T + 1;
    jScrollBar.VisibleAmount = jScrollBar.UnitIncrement;
    
    
    %draw epoch seperating lines
    for i = 1:T-1
        plot([i+0.5 i+0.5],[-0.1 1.1],'color',[1 0.5 0.25],'linestyle','- -')
        
    end
    
    x_coord = [];
    windowInMs = 1000*windowInS;
    if ~isempty(events2plot)
        % draw event lines
        currentFrame = 0;
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
                    plot(ax,[x_coord(count) x_coord(count)],[-0.1 1.1],'linestyle','- -')
                    ht = text(x_coord(count),1.12,str);
                    set(ht,'rotation',90);
                    count = count + 1;
                end
            end
        end
        hold off;
    end
    
    plot(axLikelihood,t(1:T*size2),L(:,1:T*size2),'linewidth',Linewidth);
    hold(axLikelihood,'on');
    set(axLikelihood,'ylim',[ylim1 ylim2]);
    set(axLikelihood,'xlim',[0.5 T+0.5]);
    lim = get(axLikelihood,'xlim');
    set(axLikelihood,'xtick',floor(lim(1):lim(2)))
    xlabel(axLikelihood,fastif(epoched,'Trials','Time(sec)'))
    
    % draw epoch seperating lines
    for i = 1:T-1
        plot(axLikelihood,[i+0.5 i+0.5],[ylim1 ylim2],'color',[1 0.5 0.25],'linestyle','- -')
    end
    
    
    % draw event lines
    for i = 1:length(x_coord)
        plot(axLikelihood,[x_coord(i) x_coord(i)],[ylim1 ylim2],'linestyle','- -');
    end
    hold(axLikelihood,'off');
    
end
ylabel(axLikelihood,'Log-likelihood of data under most probable model')
axes(ax);
set(gca,'ylim',[-0.1 1.1]);
ax1 = gca;
ylabel('Probability of Model Being Active')
xlabel(fastif(epoched,'Trials','Time(sec)'))
str = {};
for i = 1:length(models2plot)
    newstr = ['Model ' num2str(models2plot(i))];
    str{i} = newstr;
end

axlegend = axes('position',[pos(1)+pos(3)+0.05 pos(2) 0.06 pos(4)]);
set(axlegend,'Tag','legendaxes');
plot(axlegend,1,models2plot');
set(axlegend,'visible','off')
legend(str,'location','NorthOutside')
%legend(str,'position',[pos(1)+pos(3)+0.07 pos(2)+0.614 0.1 0.1]);

% Smoothing Window Length Text
if ~isempty(smoothlength)
    smoothLengthText = ['Smoothing Window Length: ' num2str(smoothlength) ' sec.'];
    yl = get(gca,'ylim');
    text(-1,(3*yl(1)+yl(2))/4,smoothLengthText)
    
    if epoched
        
    end
end
% -----------------------------

% Window Size Adjustment
posWindowSize = [0.835 0.3 0.07 0.045];
posToggleEvent = [0.835 0.25 0.15 0.03];
posWindowSizeLabel = [0.795 0.35 0.17 0.02];
uicontrol('Parent',figh, ...
    'Units', 'normalized', ...
    'BackgroundColor',[1 1 1], ...
    'Position', posWindowSize, ...
    'Style','edit', ...
    'Tag','WindowSizeTag',...
    'string',num2str(T),...
    'Callback', {@callback_windowSize,jScrollBar,sRate,epoched,size2,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2});

uicontrol('Parent',figh, ...
    'Units', 'normalized', ...
    'BackgroundColor',BACKCOLOR, ...
    'Position', posWindowSizeLabel, ...
    'Style','text', ...
    'Tag','WindowSizeLabelTag',...
    'string','Adjust Window Size')


toggleEventHandle = uicontrol('Parent',figh, ...
    'Units', 'normalized', ...
    'BackgroundColor',BACKCOLOR, ...
    'Position', posToggleEvent, ...
    'Style','checkbox', ...
    'Tag','ToggleEvent',...
    'string','Toggle Events',...
    'Callback', {@callback_toggleEvent,jScrollBar,sRate,epoched,size2,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2,Linewidth});

set(toggleEventHandle,'value',1);

%-----------------------


jScrollBar.AdjustmentValueChangedCallback = {@callback_slider_prob,sRate,epoched,size2,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2,Linewidth};

end

function callback_windowSize(~,~,jScrollBar,sRate,epoched,epoch_length,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2)
obj = findobj('parent',gcbf,'Tag','WindowSizeTag'); T = get(obj,'string'); T = str2num(T);
data = get(gcf,'userdata');
data{2} = T;
set(gcf,'userdata',data);
if ~epoched
    jScrollBar.UnitIncrement = ceil(T/10);
    jScrollBar.Maximum = size(data{1},1)/sRate-T + ceil(T/10);
    jScrollBar.VisibleAmount = jScrollBar.UnitIncrement;
else
    jScrollBar.UnitIncrement = ceil(T/5);
    jScrollBar.Maximum = size(data{1},1)/epoch_length-T + ceil(T/5);
    jScrollBar.VisibleAmount = jScrollBar.UnitIncrement; 
    
end
end

function callback_toggleEvent(~,~,jScrollBar,sRate,epoched,size2,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2,Linewidth)
callback_slider_prob(jScrollBar,1,sRate,epoched,size2,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2,Linewidth);
end




