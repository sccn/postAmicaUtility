% pop_ozeegplot() - postAmicaUtility plugin. eegplot() is modified to show 
%                   alteration of dominant multiple models represented probabilistically. 

% Author
% Ozgur Baklan.
% Arnaud Delorme for the original pop_eegplot().

% Hisotry
% 12/26/2016 Makoto. Help editied.

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [com] = pop_ozeegplot( EEG, icacomp, superpose, reject, topcommand, varargin);

com = '';
if nargin < 1
    help pop_ozeegplot;
    return;
end;
if nargin < 2
    icacomp = 1;
end;
if nargin < 3
    superpose = 0;
end;
if nargin < 4
    reject = 1;
end;
if icacomp == 0
    if isempty( EEG.icasphere )
        disp('Error: you must run ICA first'); return;
    end;
end;

if nargin < 3 & EEG.trials > 1
    
    % which set to save
    % -----------------
    uilist       = { { 'style' 'text' 'string' 'Add to previously marked rejections? (checked=yes)'} , ...
        { 'style' 'checkbox' 'string' '' 'value' 1 } , ...
        { 'style' 'text' 'string' 'Reject marked trials? (checked=yes)'} , ...
        { 'style' 'checkbox' 'string' '' 'value' 0 } };
    result = inputgui( { [ 2 0.2] [ 2 0.2]} , uilist, 'pophelp(''pop_ozeegplot'');', ...
        fastif(icacomp==0, 'Manual component rejection -- pop_ozeegplot()', ...
        'Reject epochs by visual inspection -- pop_ozeegplot()'));
    size_result  = size( result );
    if size_result(1) == 0 return; end;
    
    if result{1}, superpose=1; end;
    if ~result{2}, reject=0; end;
    
end;

if EEG.trials > 1
    if icacomp == 1 macrorej  = 'EEG.reject.rejmanual';
        macrorejE = 'EEG.reject.rejmanualE';
    else			macrorej  = 'EEG.reject.icarejmanual';
        macrorejE = 'EEG.reject.icarejmanualE';
    end;
    elecrange = [1:EEG.nbchan];
    colrej = EEG.reject.rejmanualcol;
    rej  = eval(macrorej);
    rejE = eval(macrorejE);
    
    eeg_rejmacro; % script macro for generating command and old rejection arrays
    
else % case of a single trial (continuous data)
    %if icacomp,
    %    	command = ['if isempty(EEG.event) EEG.event = [eegplot2event(TMPREJ, -1)];' ...
    %         'else EEG.event = [EEG.event(find(EEG.event(:,1) ~= -1),:); eegplot2event(TMPREJ, -1, [], [0.8 1 0.8])];' ...
    %         'end;'];
    %else, command = ['if isempty(EEG.event) EEG.event = [eegplot2event(TMPREJ, -1)];' ...
    %         'else EEG.event = [EEG.event(find(EEG.event(:,1) ~= -2),:); eegplot2event(TMPREJ, -1, [], [0.8 1 0.8])];' ...
    %         'end;'];
    %end;
    %if reject
    %   command = ...
    %   [  command ...
    %      '[EEG.data EEG.xmax] = eegrej(EEG.data, EEG.event(find(EEG.event(:,1) < 0),3:end), EEG.xmax-EEG.xmin);' ...
    %      'EEG.xmax = EEG.xmax+EEG.xmin;' ...
    %   	'EEG.event = EEG.event(find(EEG.event(:,1) >= 0),:);' ...
    %      'EEG.icaact = [];' ...
    %      'EEG = eeg_checkset(EEG);' ];
    eeglab_options; % changed from eeglaboptions 3/30/02 -sm
    if reject == 0, command = [];
    else,
        command = ...
            [  '[EEGTMP LASTCOM] = eeg_eegrej(EEG,eegplot2event(TMPREJ, -1));' ...
            'if ~isempty(LASTCOM),' ...
            '  [ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
            '  if ~isempty(tmpcom),' ...
            '     EEG = eegh(LASTCOM, EEG);' ...
            '     eegh(tmpcom);' ...
            '     eeglab(''redraw'');' ...
            '  end;' ...
            'end;' ...
            'clear EEGTMP tmpcom;' ];
        if nargin < 4
            res = questdlg2( strvcat('Mark stretches of continuous data for rejection', ...
                'by dragging the left mouse button. Click on marked', ...
                'stretches to unmark. When done,press "REJECT" to', ...
                'excise marked stretches (Note: Leaves rejection', ...
                'boundary markers in the event table).'), 'Warning', 'Cancel', 'Continue', 'Continue');
            if strcmpi(res, 'Cancel'), return; end;
        end;
    end;
    eegplotoptions = { 'winlength', 5, 'events', EEG.event };
    if ~isempty(EEG.chanlocs) & icacomp
        eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs };
    end;
end;

if EEG.nbchan > 100
    disp('pop_ozeegplot() note: Baseline subtraction disabled to speed up display');
    eegplotoptions = { eegplotoptions{:} 'submean' 'off' };
end;
if isfield(EEG.etc,'amica') 
    if isfield(EEG.etc.amica,'v_smooth') && size(EEG.etc.amica.v_smooth,2) == size(EEG.data,2) && size(EEG.etc.amica.v_smooth,3) == size(EEG.data,3)
        disp('Using EEG.etc.amica.v_smooth to determine most probable models')
        v_smooth = EEG.etc.amica.v_smooth;
    else
        if size(EEG.etc.amica.LLt,2) == size(EEG.data,2) && size(EEG.etc.amica.LLt,3) == size(EEG.data,3)
            uilist = {{'style' 'text' 'string' 'Length of probability smoothing window (in sec.)' } ...
                      {'style' 'edit' 'string' '2'}};
             result = inputgui({[1 1]},uilist);
                 if isempty(result)
                     return;
                 else
                     smoothlength = str2num(result{1})
                 end
                 if isfield(EEG.etc.amica,'smooth_length') && EEG.etc.amica.smooth_length == smoothlength
                     
                     EEG = eeg_checkamica(EEG);
                     v_smooth = EEG.etc.amica.v_smooth;
                 else
                    v_smooth = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlength);
                    
                 end
        else
            error('Size of probabilities does not match the size of data')
        end
        
    end
    [~, winners] = max(v_smooth,[],1);
    if icacomp == 1
        %for index = 1:2:length(eegplotoptions)
        %    if strcmpi(eegplotoptions{index}, 'winlength') eegplotoptions{index+1} = 20; end;
        %end;
        %eegplotoptions = { eegplotoptions{:} 'spacing', 100 };
        
        ozeegplot( EEG.data,'modelactive',winners, 'srate', EEG.srate, 'title', 'Scroll channel activities with active AMICA models-- ozeegplot()', ...
            'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}, varargin{:});
        
        %eeg_multieegplot( EEG.data, [], [], oldrej, oldrejE, 'title', 'Scroll channel activities -- eegplot()', 'srate', ...
        %	      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command);
    else
        tmpdata = eeg_getdatact(EEG, 'component', [1:size(EEG.icaweights,1)]);
        for i=1:length(EEG.icaweights(:,1));
            chans(i).labels=sprintf('%s%s','comp',num2str(i));
            chans(i).badchan=EEG.reject.gcompreject(i);
        end
        
        ozeegplot( tmpdata, 'modelactive',winners,'srate', EEG.srate, 'title', 'Scroll component activities with active AMICA models-- ozeegplot()', ...
            'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}, varargin{:});
        %eeg_multieegplot( tmpdata, [], [], oldrej, oldrejE, 'title', 'Scroll component activities -- eegplot()', 'srate', ...
        %	      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command',
        %	      command);
    end;
else
    error('No AMICA solution found');
end
com = [ com sprintf('pop_ozeegplot( %s, %d, %d, %d);', inputname(1), icacomp, superpose, reject) ];
return;
