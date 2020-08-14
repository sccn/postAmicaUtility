% pop_moderp() - postAmicaUtility plugin. pop_erpimage() is modified
%                to show event-related alteration of dominant models in AMICA results. 

% Author
% Ozgur Baklan.
% Arnaud Delorme for the original pop_erpimage().

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



function varargout = pop_moderp( EEG, typeplot, channel, projchan, titleplot, smooth, decimate, sortingtype, sortingwin, sortingeventfield, varargin)
 
varargout{1} = '';
if nargin < 1
   help pop_erpimage;
   return;
end;
typeplot = 1;

if typeplot == 0 & isempty(EEG.icasphere)
   error('no ICA data for this set, first run ICA');
end;   
if EEG.trials == 1
   error('erpimage of one trial cannot be plotted');
end;   

if nargin < 2	
	typeplot = 1; %1=signal; 0=component
end;
lastcom = [];
if nargin < 3
	popup = 1;
else
	popup = isstr(channel) | isempty(channel);
	if isstr(channel)
		lastcom = channel;
	end;
end;
mod = 1;
if popup
	% get contextual help
	% -------------------
    clear functions;
    erpimagefile = which('erpimage.m');
	[txt2 vars2] = gethelpvar(erpimagefile);
	txt  = { txt2{:}};
	vars = { vars2{:}};
	% [txt vars] = gethelpvar('erpimopt.m');
	% txt  = { txt{:} txt2{:}};
	% vars = { vars{:} vars2{:}};
    
    opt.topo         = getkeyval(lastcom, 'topo', 'present', 1);
    opt.fieldname    = getkeyval(lastcom, 10);
    opt.type         = getkeyval(lastcom, 8);
    opt.renorm       = getkeyval(lastcom, 'renorm','', 'no');
    opt.erp          = fastif(getkeyval(lastcom, '''erp''', 'present', 1), 'on', 'off');
    opt.cbar         = fastif(getkeyval(lastcom, 'cbar', 'present', 1), 'on', 'off');
    opt.nosort       = fastif(getkeyval(lastcom, 'nosort', 'present', 0), 'on', 'off');
    opt.noplot       = fastif(getkeyval(lastcom, 'noplot', 'present', 0), 'on', 'off');
    opt.plotamps     = fastif(getkeyval(lastcom, 'plotamps', 'present', 0), 'on', 'off');
    opt.index        = str2num(getkeyval(lastcom,3,[],'1'));
    opt.smoothing    = str2num(getkeyval(lastcom, 6, [], int2str(min(max(EEG.trials-5,0), 0))));
    opt.downsampling = str2num(getkeyval(lastcom, 7, [], '1'));
    opt.caxis        = [0.5 EEG.etc.amica.num_models+0.5]; %str2num(getkeyval(lastcom, 'caxis'));
    opt.eventrange   = str2num(getkeyval(lastcom, 9));
    opt.align        = str2num(getkeyval(lastcom, 'align'));
    opt.phasesort    = str2num(getkeyval(lastcom, 'phasesort'));
    opt.coher        = str2num(getkeyval(lastcom, 'coher'));
    opt.spec         = str2num(getkeyval(lastcom, 'spec'));
    opt.vert         = str2num(getkeyval(lastcom, 'vert'));
    opt.limits       = str2num(getkeyval(lastcom, 'limits'));
    opt.limits       = [ opt.limits NaN NaN NaN NaN NaN NaN NaN NaN NaN ]; opt.limits = opt.limits(1:9);
    opt.coher        = [ opt.coher NaN NaN NaN NaN NaN NaN NaN NaN NaN ];  opt.coher  = opt.coher(1:3);
    opt.phasesort    = [ opt.phasesort NaN NaN NaN NaN NaN NaN NaN NaN NaN ]; opt.phasesort = opt.phasesort(1:4);
    if isnan(opt.limits(1)), opt.limits(1:2) = 1000*[EEG.xmin EEG.xmax]; end;
    
	commandphase = [ 'if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''phase''),''string'')),' ...
					 '   if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''coher''),''string'')), ' ...
					 '      set(findobj(''parent'', gcbf, ''tag'', ''coher''), ''string'', get(findobj(''parent'', gcbf, ''tag'', ''phase''),''string''));' ...
					 'end; end;' ];
	commandcoher = [ 'if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''coher''),''string'')),' ...
					 '   if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''phase''),''string'')), ' ...
					 '      set(findobj(''parent'', gcbf, ''tag'', ''phase''), ''string'', get(findobj(''parent'', gcbf, ''tag'', ''coher''),''string''));' ...
					 'end; end;' ];
	
	commandfield = ['if isempty(EEG.event)' ...
				   '   errordlg2(''No events'');' ...
				   'else' ...
				   '   tmpfieldnames = fieldnames(EEG.event);' ...
				   '   [tmps,tmpv] = listdlg2(''PromptString'', ''Select fields'', ''SelectionMode'',''single'',''ListString'', tmpfieldnames);' ...
				   '   if tmpv' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''field''), ''string'', tmpfieldnames{tmps});' ...
				   '   end;' ...
				   'end;' ...
                   'if isempty(get(findobj(''parent'', gcbf, ''tag'', ''type''), ''string'')),' ...
                   '   warndlg2(''Do not forget to select an event type in the next edit box'');' ...
                   'end;' ...
				   'clear tmps tmpv tmpfieldnames;' ];
	commandtype = [ 'if ~isfield(EEG.event, ''type'')' ...
				   '   errordlg2(''No type field'');' ...
				   'else' ...
                   '   if isnumeric(EEG.event(1).type),' ...
				   '        [tmps,tmpstr] = pop_chansel(unique([ EEG.event.type ]));' ...
				   '   else,' ...
                   '        [tmps,tmpstr] = pop_chansel(unique({ EEG.event.type }));' ...
                   '   end;' ...
				   '   if ~isempty(tmps)' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''type''), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
                   'if isempty(get(findobj(''parent'', gcbf, ''tag'', ''field''), ''string'')),' ...
                   '   warndlg2(''Do not forget to select an event type in the previous edit box'');' ...
                   'end;' ...
				   'clear tmps tmpv tmpstr tmpfieldnames;' ];
	
	geometry = { [1 1 0.1 0.8 2.1] [1 1 1 1 1] [1 1 1 1 1] [1 1 1 1 1] [1] [1] [1 1 1 0.8 0.8 1.2] [1 1 1 0.8 0.8 1.2] [1] [1] ...
				 [1.6 1.7 1.2 1 .5] [1.6 1.7 1.2 1 .5] [1] [1] [1.5 1 1 1 1] [1.5 1 1 1 1] [1] [1] [1.5 1 1 2.2] [1.5 1 1 2.2]};
    uilist = { { 'Style', 'text', 'visible', fastif(mod,'off','on'), 'string', fastif(typeplot, 'Channel', 'Component(s)'), 'fontweight', 'bold'} ...
			   { 'Style', 'edit', 'string', num2str(opt.index), 'tag', 'chan', 'visible', fastif(mod,'off','on')} { } ...
			   { 'Style', 'text', 'string', 'Figure title', 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', '', 'tag', 'title'  } ...
			   ...
			   { 'Style', 'text', 'string', 'Smoothing', 'fontweight', 'bold', 'tooltipstring', context('avewidth',vars,txt) } ...
			   { 'Style', 'edit', 'string', num2str(opt.smoothing), 'tag', 'smooth' } ...
			   { 'Style', 'checkbox', 'string', 'Plot scalp map', 'tooltipstring', 'plot a 2-d head map (vector) at upper left', ...
				 'value', fastif(mod,0,opt.topo), 'tag', 'plotmap', 'visible', fastif(mod,'off','on')} { } { } ...
			   { 'Style', 'text', 'string', 'Downsampling', 'fontweight', 'bold', 'tooltipstring', context('decimate',vars,txt) } ...
			   { 'Style', 'edit', 'string', num2str(opt.downsampling), 'tag', 'decimate' } ...
			   { 'Style', 'checkbox', 'string', 'Plot ERP', 'tooltipstring', context('erp',vars,txt), 'value', fastif(mod,0,fastif(strcmpi(opt.erp, 'on'), 1,0)), 'tag', 'erp', 'visible', fastif(mod,'off','on')} ...
			   { 'Style', 'text', 'string', fastif(typeplot, 'ERP limits (uV)','ERP limits'), 'tooltipstring', [ 'Plotting limits for ERP trace [min_uV max_uV]' 10 '{Default: ERP data limits}'], 'visible', fastif(mod,'off','on')} ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.limits(3)), num2str(opt.limits(3:4)), ''), 'tag', 'limerp', 'visible', fastif(mod,'off','on')} ...
			   { 'Style', 'text', 'string', 'Time limits (ms)', 'fontweight', 'bold', 'tooltipstring',  'Select time subset in ms' } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.limits(1)), num2str(opt.limits(1:2)), ''), 'tag', 'limtime' } ...
			   { 'Style', 'checkbox', 'string', 'Plot colorbar','tooltipstring', context('caxis',vars,txt), 'value', fastif(strcmpi(opt.cbar, 'on'), 1,0), 'tag', 'cbar'} ...
			   { 'Style', 'text', 'string', 'Color limits (see Help)','tooltipstring', context('caxis',vars,txt) } ...
			   { 'Style', 'edit', 'string', num2str(opt.caxis), 'tag', 'caxis' } ...
			   {} ...
			   { 'Style', 'text', 'string', 'Sort/align trials by epoch event values', 'fontweight', 'bold'} ...
			   { 'Style', 'pushbutton', 'string', 'Epoch-sorting field', 'callback', commandfield, ...
				 'tooltipstring', 'Epoch-sorting event field name (Ex: latency; default: no sorting):' } ...
			   { 'Style', 'pushbutton', 'string', 'Event type(s)', 'callback', commandtype, 'tooltipstring', ['Event type(s) subset (default: all):' 10 ...
                '(See ''/Edit/Edit event values'' for event types)']  } ...
			   { 'Style', 'text', 'string', 'Event time range', 'tooltipstring', [ 'Sorting event window [start, end] in milliseconds (default: whole epoch):' 10 ...
												  'events are only selected within this time window (can be usefull if several' 10 ...
												  'events of the same type are in the same epoch, or for selecting trials with given response time)']} ...
			   { 'Style', 'text', 'string', 'Rescale', 'tooltipstring', 'Rescale sorting variable to plot window (yes|no|a*x+b)(Ex:3*x+2):' } ...
			   { 'Style', 'text', 'string', 'Align', 'tooltipstring',  context('align',vars,txt) } ...
			   { 'Style', 'checkbox', 'string', 'Don''t sort by value', 'tooltipstring', context('nosort',vars,txt), 'value', fastif(strcmpi(opt.nosort, 'on'), 1,0), 'tag', 'nosort' } ...
			   { 'Style', 'edit', 'string', opt.fieldname, 'tag', 'field' } ...
			   { 'Style', 'edit', 'string', opt.type, 'tag', 'type' } ...
			   { 'Style', 'edit', 'string', num2str(opt.eventrange), 'tag', 'eventrange' } ...
			   { 'Style', 'edit', 'string', opt.renorm, 'tag', 'renorm' } ...
			   { 'Style', 'edit', 'string', num2str(opt.align),  'tag', 'align' } ...
			   { 'Style', 'checkbox', 'string', 'Don''t plot values', 'tooltipstring', context('noplot',vars,txt), 'value', fastif(strcmpi(opt.noplot, 'on'), 1,0), 'tag', 'noplot' } ...
			   {} ...
			   { 'Style', 'text', 'string', 'Sort trials by phase', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Frequency (Hz | minHz maxHz)', 'tooltipstring', ['sort by phase at maximum-power frequency' 10 ...
               'in the data within the range [minHz,maxHz]' 10 '(overrides frequency specified in ''coher'' flag)']  } ...
			   { 'Style', 'text', 'string', 'Percent low-amp. trials to ignore', 'tooltipstring', ['percent of trials to reject for low' ...
			    'amplitude. Else,' 10 'if prct is in the range [-100,0] -> percent to reject for high amplitude'] } ...
			   { 'Style', 'text', 'string', 'Window center (ms)', 'tooltipstring', 'Center time of the n-cycle window'  } ...
			   { 'Style', 'text', 'string', 'Wavelet cycles', 'tooltipstring', 'cycles per wavelet window'  } {}...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.phasesort(3)),num2str(opt.phasesort(3:4)),'') 'tag', 'phase', 'callback', commandphase } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.phasesort(2)),num2str(opt.phasesort(2)),''), 'tag', 'phase2' } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.phasesort(1)),num2str(opt.phasesort(1)),''), 'tag', 'phase3' } ...
			   { 'Style', 'text', 'string', '        3'  } {}...
			   {} ...
			   { 'Style', 'text', 'string', 'Inter-trial coherence options', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Frequency (Hz | minHz maxHz)', 'tooltipstring', [ '[freq] -> plot erp plus amp & coher at freq (Hz)' 10 ...
               '[minHz maxHz] -> find max in frequency range' 10 '(or at phase freq above, if specified)']} ...
			   { 'Style', 'text', 'string', 'Signif. level (<0.20)', 'tooltipstring', 'add coher. signif. level line at alpha (alpha range: (0,0.1])' } ...
			   { 'Style', 'text', 'string', 'Amplitude limits (dB)'  } ...
			   { 'Style', 'text', 'string', 'Coher limits (<=1)'  } ...
			   { 'Style', 'checkbox', 'string', 'Image amps', 'tooltipstring',  context('plotamps',vars,txt), 'value', fastif(strcmpi(opt.plotamps, 'on'), 1,0), 'tag', 'plotamps' } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.coher(1)), num2str(opt.coher(1:2)), ''), 'tag', 'coher', 'callback', commandcoher } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.coher(3)), num2str(opt.coher(3)),   ''), 'tag', 'coher2'  } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.limits(5)), num2str(opt.limits(5:6)), ''), 'tag', 'limamp' } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.limits(7)), num2str(opt.limits(7:8)), ''), 'tag', 'limcoher' } ...
               {'style', 'text', 'string', '   (Requires signif.)' } ... 
			   {} ...
			   { 'Style', 'text', 'string', 'Other options', 'fontweight', 'bold'} ...
			   { 'Style', 'text', 'string', 'Plot spectrum (minHz maxHz)','tooltipstring',  context('spec',vars,txt)} ...
			   { 'Style', 'text', 'string', 'Baseline ampl. (dB)', 'tooltipstring', 'Use it to fix baseline amplitude' } ...
			   { 'Style', 'text', 'string', 'Mark times (ms)','tooltipstring',  context('vert',vars,txt)} ...
			   { 'Style', 'text', 'string', 'More options (see >> help erpimage)' } ...
			   { 'Style', 'edit', 'string', num2str(opt.spec), 'tag', 'spec' } ...
			   { 'Style', 'edit', 'string', fastif(~isnan(opt.limits(9)), num2str(opt.limits(9)), ''), 'tag', 'limbaseamp' } ...
			   { 'Style', 'edit', 'string', num2str(opt.vert), 'tag', 'vert' } ...
			   { 'Style', 'edit', 'string', '', 'tag', 'others' } ...
			};
	if typeplot == 0 % add extra param for components
		geometry = { [1 1 0.1 0.8 2.1] geometry{:} };
		uilist   = { { } { } { } { } { } uilist{:}};
		uilist{1} = uilist{6};
		uilist{2} = uilist{7};
		uilist{6} = { 'Style', 'text', 'string', 'Project to channel #', 'fontweight', 'bold','tooltipstring', ['Project component(s) to data channel' 10 ...
												  'This allows plotting projected component activity at one channel in microvolts'] };
		uilist{7} = { 'Style', 'edit', 'string', getkeyval(lastcom, 4), 'tag', 'projchan' };
	end;
    [oldres a b res] = inputgui( geometry, uilist, 'pophelp(''pop_erpimage'');', ...
							 fastif(mod,'AMICA: Event Related Model Activity -- pop_erpimage()',fastif( typeplot, 'Channel ERP image -- pop_erpimage()', 'Component ERP image -- pop_erpimage()')));
	if isempty(oldres), return; end;

	% first rows
	% ---------
	channel   	 = eval( [ '[' res.chan ']' ]);
	titleplot    = res.title;
	if isfield(res, 'projchan'), projchan = str2num(res.projchan); else, projchan = []; end;
    opt = [];
	if ~isempty(res.others)
        try,
            tmpcell = eval( [ '{' res.others '}' ] );
            opt = struct( tmpcell{:} );
        catch, error('Additional options ("More options") requires ''key'', ''val'' arguments');
        end;
	end;
	if ~typeplot & isempty(projchan)
        opt.yerplabel = '';
    else
        opt.yerplabel = '\muV' ;
	end;
	if isempty(titleplot)
        if typeplot==1 % if channel plot
            if ~isempty(EEG.chanlocs) % if channel information exist
                  titleplot = [ EEG.chanlocs(channel).labels ];
            else, titleplot = [ int2str(channel) ];
            end
        else
            titleplot = [ 'Comp. ' int2str(channel) ];
            if ~isempty(projchan),
                if ~isempty(EEG.chanlocs) % if channel plot
                      titleplot = [ titleplot ' -> ' EEG.chanlocs(projchan).labels ];
                else, titleplot = [ titleplot ' -> Chan. ' int2str(projchan) ];
                end
            end;
        end
    end;
	smooth       = eval(res.smooth);
    if res.plotmap
		if isfield(EEG.chanlocs, 'theta')
            if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end;
			if typeplot == 0
				     opt.topo = [ ' { mean(EEG.icawinv(:,[' int2str(channel) ']),2) EEG.chanlocs EEG.chaninfo } '];
			else     opt.topo = [ ' { [' int2str(channel) '] EEG.chanlocs EEG.chaninfo } '];
			end;	
		end;
	end;
	
	decimate     = eval( res.decimate );
	if res.erp
		opt.erp = 'on';
	end;
	
	% finding limits
	% --------------
	limits(1:8)  = NaN;
	if ~isempty(res.limerp)
		limits(3:4) = eval( [ '[' res.limerp ']' ]); 
	end;
	if ~isempty(res.limtime) % time limits
		if ~strcmp(res.limtime, num2str(1000*[EEG.xmin EEG.xmax]))
			limits(1:2) = eval( [ '[' res.limtime ']' ]);
		end;
	end;
	if ~isempty(res.limamp)
		limits(5:6) = eval( [ '[' res.limamp ']' ]);
	end;
	if ~isempty(res.limcoher)
		limits(7:8) = eval( [ '[' res.limcoher ']' ]);
	end;
	if ~isempty(res.limbaseamp)
		limits(9) = eval( res.limbaseamp ); %bamp
	end;
	if ~all(isnan(limits))
        opt.limits = limits;
	end;
	
	% color limits
	% --------------
	if res.cbar
		opt.cbar = 'on';
	end;
	if res.caxis
		opt.caxis = str2num(res.caxis);
	end;
	
	% event rows
	% ----------
	if res.nosort
		opt.nosort = 'on';
	end;
	try, sortingeventfield = eval( res.field ); catch, sortingeventfield = res.field; end;
    if ~isempty(res.type)
       if strcmpi(res.type(1),'''')
            sortingtype = eval( [ '{' res.type '}' ] );
       else sortingtype = parsetxt( res.type );
       end;
    end
	sortingwin   = eval( [ '[' res.eventrange ']' ] );
	if ~isempty(res.field) & ~strcmp(res.renorm, 'no')
		opt.renorm = res.renorm;
	end;
	if ~isempty(res.align)
		opt.align = str2num(res.align);
	end;
	if res.noplot
		opt.noplot = 'on';
	end;

	% phase rows
	% ----------
	tmpphase = [];
	if ~isempty(res.phase)
		tmpphase = eval( [ '[ 0 0 ' res.phase ']' ]);
	end;
	if ~isempty(res.phase2)
		tmpphase(2) = eval( res.phase2 );
	end;
	if ~isempty(res.phase3)
		tmpphase(1) = eval( res.phase3 );
	end;
	if ~isempty(tmpphase)
		opt.phasesort = tmpphase;
	end;
	
	% coher row
	% ----------
	tmpcoher = [];
	if res.plotamps
		opt.plotamps = 'on';
	end;
	if ~isempty(res.coher)
		tmpcoher = eval( [ '[' res.coher ']' ]);
	end;
	if ~isempty(res.coher2)
		if length(tmpcoher) == 1
			tmpcoher(2) = tmpcoher(1);
		end;
		tmpcoher(3) = eval( res.coher2 );
	end;
	if ~isempty(tmpcoher)
		opt.coher = tmpcoher;
	end;

	% options row
	% ------------
	if ~isempty(res.spec)
		opt.spec = eval( [ '[' res.spec ']' ]);
	end;
	if ~isempty(res.vert)
		opt.vert = eval( [ '[' res.vert ']' ]);
	end;
	figure;
    options = '';
else
	options = '';
	if nargin < 4
		projchan = [];
	end;
	if nargin < 5
		titleplot = ' ';
	end;
	if nargin < 6
		smooth = 5;
	end;
	if nargin < 7
		decimate = 0;
	end;
	if nargin < 8
		sortingtype = [];
	end;
	if nargin < 9
		sortingwin = [];
	end;
	if nargin < 10
		sortingeventfield = [];
	end;
    %options = vararg2str(varargin); % NO BECAUSE OF THE CHANNEL LOCATION
    %                                  PROBLEM BELOW
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else  
		  if ~iscell( varargin{ i } )
		      options = [ options ',' vararg2str({varargin{i}}) ];
		  else
		      options = [ options ', { [' num2str(varargin{ i }{1}') ']'' EEG.chanlocs EEG.chaninfo }' ];
		  end;    
		end;
	end;	
end;
try, icadefs; set(gcf, 'color', BACKCOLOR,'Name',' erpimage()'); catch, end;

% testing inputs
% --------------
if typeplot == 0 & length(channel) > 1 & isempty(projchan)1
	error('A channel must be selected to plot (the sum of) several component projections');
end;

% find sorting latencies
% ---------------------
typetxt = '';
if ~isempty(sortingeventfield)
    %events = eeg_getepochevent( EEG, sortingtype, sortingwin, sortingeventfield);
	events = sprintf('eeg_getepochevent( EEG, %s)', vararg2str({sortingtype, sortingwin, sortingeventfield}));
	
    % generate text for the command
    % -----------------------------
    for index = 1:length(sortingtype)
        if isstr(sortingtype{index})
            typetxt = [typetxt ' ''' sortingtype{index} '''' ];
        else
            typetxt = [typetxt ' ' num2str(sortingtype{index}) ];
        end;
    end;    
% $$$ 	% renormalize latencies if necessary
% $$$ 	% ----------------------------------
% $$$ 	switch lower(renorm)
% $$$ 	    case 'yes',
% $$$ 	         disp('Pop_erpimage warning: *** sorting variable renormalized ***');
% $$$ 	         events = (events-min(events)) / (max(events) - min(events)) * ...
% $$$ 	                        0.5 * (EEG.xmax*1000 - EEG.xmin*1000) + EEG.xmin*1000 + 0.4*(EEG.xmax*1000 - EEG.xmin*1000);
% $$$ 	    case 'no',;
% $$$ 	    otherwise,
% $$$ 	        locx = findstr('x', lower(renorm))
% $$$ 	        if length(locx) ~= 1, error('Pop_erpimage error: unrecognize renormalazing formula'); end;
% $$$ 	        eval( [ 'events =' renorm(1:locx-1) 'events' renorm(locx+1:end) ';'] );
% $$$ 	end;
else
	events = 'ones(1, EEG.trials)*EEG.xmax*1000';
    %events = ones(1, EEG.trials)*EEG.xmax*1000;
    sortingeventfield = '';
end;           

if typeplot == 1
	tmpsig = ['mean(EEG.data([' int2str(channel) '], :),1)'];
else
    % test if ICA was computed or if one has to compute on line
    % ---------------------------------------------------------
    tmpsig = [ 'eeg_getdatact(EEG, ''component'', [' int2str(channel) '])' ];
end;


    tmpsig = ['find_winners(EEG.etc.amica.v_smooth)'];


% outputs
% -------
outstr = '';
if ~popup
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the data and generate output command
% --------------------------------------------
if isempty( options )
    if isfield(opt, 'topo')
        tmptopo = opt.topo;
        opt = rmfield(opt, 'topo');
    else
        tmptopo = '';
    end;
    fields = fieldnames(opt);
    values = struct2cell(opt);
    params = { fields{:}; values{:} };
    options = [ ',' vararg2str( { params{:} } ) ];
    tmpind = find( options == '\' ); options(tmpind(1:2:end)) = [];
    if ~isempty(tmptopo), options = [ options ',''topo'',' tmptopo ]; end;
end;


EEG = eeg_checkamica(EEG);
answer1 = questdlg2('Would you like to use smooth model probabilities to assign most probable models?','Smooth Probabilities?','Yes','No','');
newSmoothing = false;

if isempty(answer1)
   return;
end

if strcmpi(answer1,'No')
    tmpsig = 'find_winners(EEG.etc.amica.v)';
else
    if isfield(EEG.etc.amica,'v_smooth')
        answer2 = questdlg2('Would you like to use already calculated smooth probabilities under AMICA structure?','','Yes','No','');
        if isempty(answer2)
            return;
        end
        
        if strcmpi(answer2,'Yes')
            tmpsig = 'find_winners(EEG.etc.amica.v_smooth)';
        else
            tmpsig = 'find_winners(EEG.etc.amica.v_smooth)';
            newSmoothing = true;
        end
        
    else
        newSmoothing = true;
        
    end
end

if newSmoothing
    prompt = {'Enter the smoothing Hanning window size: '};
    name = 'Input for probability smoothing for AMICA';
    numlines = 1;
    defaultanswer = {'2'};
    smoothlength = inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(smoothlength) 
        return   
    end 
    if isempty(str2num(smoothlength{1}))
        disp('Smoothing window size is not numeric!'); return;
    end
    disp('Smoothing probabilities... (To save smooth probability results, Use AMICA -> Smooth Model Probabilities)')
    [EEG.etc.amica.v_smooth,EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,str2num(smoothlength{1}));
    EEG.etc.amica.smooth_length = str2num(smoothlength{1});
    
end

% varargout{1} = sprintf('figure; pop_erpimage(%s,%d,%d,''%s'',%d,%d,{%s},[%s],''%s'',''%s''%s);', inputname(1), typeplot, channel, titleplot, smooth, decimate, typetxt, int2str(sortingwin), sortingeventfield, renorm, options);
popcom = sprintf('figure; pop_moderp(%s,%d, [%s],[%s],''%s'',%d,%d,{%s},[%s],''%s'' %s);', inputname(1), typeplot, int2str(channel), int2str(projchan), titleplot, smooth, decimate, typetxt, int2str(sortingwin), sortingeventfield, options);

    com = sprintf('%s erpimage( %s, %s, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), ''%s'', %d, %d %s);', outstr, tmpsig, events, titleplot, smooth, decimate, options);

    %com = sprintf('%s erpimage( %s, %s, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), ''%s'', %d, %d %s);', outstr, tmpsig, events, titleplot, smooth, decimate, options);


disp('Command executed by pop_erpimage:');
disp(' '); disp(com); disp(' ');
eval(com)

if popup
	varargout{1} = popcom; % [10 '% Call: ' com];
end;

% Color axis should be integers
obj = findobj('parent',gcf,'tag','cbar');
if ~isempty(obj)
    colorYlim = get(obj,'YLim');
    colorYlabel = get(obj,'YtickLabel');
    lim1 = str2num(colorYlabel(1,:));
    lim2 = str2num(colorYlabel(end,:));
    newYTickLabel = num2str((ceil(lim1):floor(lim2))');
    
    firstTick = colorYlim(1) + (ceil(lim1)-lim1)/(lim2-lim1)*(colorYlim(2)-colorYlim(1));
    lastTick = colorYlim(2) - (lim2-floor(lim2))/(lim2-lim1)*(colorYlim(2)-colorYlim(1));
    set(obj,'YTick',linspace(firstTick,lastTick,size(newYTickLabel,1)))
    set(obj,'YTickLabel',newYTickLabel)
end
%------------
return;

% get contextual help
% -------------------
function txt = context(var, allvars, alltext);
	loc = strmatch( var, allvars);
	if ~isempty(loc)
		txt= alltext{loc(1)};
	else
		disp([ 'warning: variable ''' var ''' not found']);
		txt = '';
	end;

