% pop_modprobplot() - postAmicaUtility plugin. Plot multiple model
%                     probability as time series in scrolling pages.

% Author
% Ozgur Baklan.

% Hisotry
% 12/26/2016 Makoto. Help editied.

% Copyright (C) Ozgur Baklan, SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG com] = pop_modprobplot(EEG,models2plot,smoothed,smoothlength,store,varargin)
if nargin<1
    help pop_modprobplot
    return
end

if ~isfield(EEG.etc,'amica') || isempty(EEG.etc.amica)
    error('No AMICA solution found. You should first load AMICA components');
end
com = '';
if nargin<6
events2plot = {};
end
if nargin<2
    plotopts = 'Smoothed probabilities|Unsmoothed Probabilities';
    guititle = 'Plot Model Probabilities -- pop_modprobplot()';
    if isfield(EEG.etc.amica,'smooth_length')
        smoothlength = num2str(EEG.etc.amica.smooth_length);
    else
        smoothlength = '2';
    end
    comsmooth = ['if get(gcbo,''value'') == 2, set(findobj(''Tag'',''smoothlengthTag''),''enable'',''off'');'...
        'set(findobj(''Tag'',''saveTag''),''enable'',''off'');else, set(findobj(''Tag'',''smoothlengthTag''),''enable'',''on'');' ...
        'set(findobj(''Tag'',''saveTag''),''enable'',''on'');end;'];
    
    comtype = [ 'if ~isfield(EEG.event, ''type'')' ...
        '   errordlg2(''No type field'');' ...
        'else' ...
        '   if isnumeric(EEG.event(1).type),' ...
        '        [tmps,tmpstr] = pop_chansel(unique([ EEG.event.type ]));' ...
        '   else,' ...
        '        [tmps,tmpstr] = pop_chansel(unique({ EEG.event.type }));' ...
        '   end;' ...
        '   if ~isempty(tmps)' ...
        '       set(findobj(''parent'', gcbf, ''tag'', ''eventtype''), ''string'', tmpstr);' ...
        '   end;' ...
        'end;' ...
        'clear tmps tmpstr;' ];
    
    
    allmodels = ['1:' num2str(EEG.etc.amica.num_models)];
    
    
    uilist = {{'style' 'text' 'string' 'Probabilities of which models'} ...
        {'style' 'edit' 'string' allmodels } ...
        {'style' 'text' 'string' 'Plot Option: '} ...
        {'style' 'popupmenu' 'string' plotopts 'value' 1 'callback' comsmooth} ...
        {'style' 'text' 'string' 'Smoothing Hanning Window Length (sec)'} ...
        {'style' 'edit' 'tag' 'smoothlengthTag' 'string' smoothlength 'enable' 'on'} ...
        {'style' 'checkbox' 'string' 'Save new smoothed probabilites under EEG.etc.amica' 'value' 0 'tag' 'saveTag'} ...
        {'style' 'pushbutton' 'string' 'Mark events of type..' 'callback' comtype} ...
        {'style' 'edit' 'string' '' 'tag' 'eventtype'}};
    uigeom = {[1 1] [1 1] [1 1] [1] [1 1]};
    
    result = inputgui(uigeom, uilist, 'pophelp(''pop_modprobplot'')', guititle, [], 'normal');
    if isempty(result)
        return
    end
    
    models2plot = eval(['[' result{1} ']']);
    if min(models2plot)<1 && max(models2plot) > EEG.etc.amica.num_models
        error('Models chosen are not in the proper range');
    end
    
    smoothed = (result{2} == 1);
    store = result{4};
    events = result{5};
    
    if smoothed
        smoothlength = str2num(result{3});
    end
    
    if ~isempty(events)
        events2plot = parsetxt(events);
        
        if isnumeric(EEG.event(1).type)
            eeg_events = unique([EEG.event.type]);
            for i = 1:length(events2plot)
                if isempty(find(eeg_events == str2num(events2plot{i}), 1))
                    disp(['No event called ''' events2plot{i} '''' 'is found, ignoring that event'])
                    events2plot(i) = [];
                end
            end
        else
            eeg_events = unique({EEG.event.type});
            events2plot(ismember(events2plot,eeg_events) == 0) = [];
            
        end
    else
        events2plot = {};
    end
else
    events2plot = varargin;
end

if ~smoothed
    smoothlength = 0;
end
if ~isfield(EEG.etc.amica,'LLt')
    if isfield(EEG.etc.amica,'Lht')
        EEG.etc.amica.LLt = EEG.etc.amica.Lht;
    end
end

if size(EEG.data,2) ~= size(EEG.etc.amica.LLt,2) || size(EEG.data,3) ~= size(EEG.etc.amica.LLt,3)
    recompute_prob = questdlg('Data size and AMICA probabilities do not match. Would you like to recompute the probabilities?', ...
        'Recompute probabilities?','Yes','No','No');
    if strcmpi(recompute_prob,'Yes')
        EEG.etc.amica.LLt = findLLt(EEG,EEG.etc.amica);
        EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
        [~,EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlength);
       
    end
end

 EEG = eeg_checkamica(EEG);
[v2plot, llt2plot] = modprobplot(EEG,models2plot,smoothlength,events2plot);


if store && smoothed
    EEG.etc.amica.v_smooth = v2plot;
    EEG.etc.amica.LLt_smooth = llt2plot;
    EEG.etc.amica.smooth_length = smoothlength;
end

if ~isempty(events2plot)
    com = sprintf('%s = pop_modprobplot(%s,%s,%d,%d,%d,%s)',inputname(1),inputname(1),vararg2str(models2plot),smoothed,smoothlength,store,vararg2str(events2plot));
else
    com = sprintf('%s = pop_modprobplot(%s,%s,%d,%d,%d)',inputname(1),inputname(1),vararg2str(models2plot),smoothed,smoothlength,store);
end











