% eegplugin_postAmicaUtility() - postAmicaUtility plugin 
%
% Usage:
%   >> eegplugin_postAmicaUtility(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   To create a new plugin, simply create a file beginning with "eegplugin_"
%   and place it in your eeglab folder. It will then be automatically 
%   detected by eeglab. See also this source code internal comments.
%   For eeglab to return errors and add the function's results to 
%   the eeglab history, menu callback must be nested into "try" and 
%   a "catch" strings. For more information on how to create eeglab 
%   plugins, see http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Ozgur Balkan and Makoto Miyakoshi. SCCN, INC, UCSD.

% History
% 11/07/2018 Makoto. Updated.
% 12/26/2016 Makoto. Created.

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2016 Ozgur Balkan and Makoto Miyakoshi
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA



function vers = eegplugin_postAmicaUtility(fig, trystrs, catchstrs)
    
    vers = 'postAmicaUtility1.01';
    if nargin < 3
        error('eegplugin_postAmicaUtility requires 3 arguments');
    end;
    
    % Find tools menu.
    menu = findobj(fig, 'tag', 'tools'); 
        % tag can be 
        % 'import data'  -> File > import data menu
        % 'import epoch' -> File > import epoch menu
        % 'import event' -> File > import event menu
        % 'export'       -> File > export
        % 'tools'        -> tools menu
        % 'plot'         -> plot menu
    
    % Create menu callback commands.
    comLoad          = [trystrs.no_check   '[EEG, LASTCOM] = pop_loadmodout(EEG);      [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);' catchstrs.store_and_hist ];
    comSmoothProb    = [trystrs.check_data '[EEG, LASTCOM] = pop_smoothprobamica(EEG); [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);' catchstrs.store_and_hist ];
    comChangeWeights = [trystrs.no_check   '[EEG, LASTCOM] = pop_changeweights(EEG);   [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);' catchstrs.store_and_hist];
    comPlotModelProb = [trystrs.no_check   '[EEG, LASTCOM] = pop_modprobplot(EEG);     [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);' catchstrs.store_and_hist];
    comSelectModel   = [trystrs.check_data '[EEG, LASTCOM] = pop_selectmodel(EEG);     [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);' catchstrs.new_and_hist ];
    comPlotModelPmi  = [trystrs.check_data '[EEG, LASTCOM] = pop_modPMI(EEG);          [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);' catchstrs.store_and_hist ];
    comTopoHistPlot  = [trystrs.check_data 'LASTCOM = pop_topohistplot(EEG, 0);' catchstrs.add_to_hist];
    comBackgEegplot  = [trystrs.check_data 'LASTCOM = pop_ozeegplot(EEG,1,1,1);' catchstrs.add_to_hist];
    comBackgCompplot = [trystrs.check_data 'LASTCOM = pop_ozeegplot(EEG,0,1,1);' catchstrs.add_to_hist];
    comModErp        = [trystrs.check_data 'LASTCOM = pop_moderp(EEG,1);'        catchstrs.add_to_hist];
    
    % Create menus.
    submenu = uimenu( menu, 'Label', 'post AMICA utility', 'separator', 'on');
    uimenu( submenu, 'Label', 'Load AMICA components',                                 'CallBack', comLoad);
    uimenu( submenu, 'Label', 'Change current ICA weights (beta)',                     'CallBack', comChangeWeights, 'separator', 'on');
    uimenu( submenu, 'Label', 'Smooth model probabilities',                            'CallBack', comSmoothProb);
    uimenu( submenu, 'Label', 'Plot pairwise mutual information (Infomax compatible)', 'CallBack', comPlotModelPmi,  'separator', 'on');
    uimenu( submenu, 'Label', 'Plot Model Probability time series',                    'CallBack', comPlotModelProb);
    uimenu( submenu, 'Label', '[For epoched] Plot event-related dominant model image', 'CallBack', comModErp);
    uimenu( submenu, 'Label', 'Plot PDF of Gaussian mixtures for each IC',             'CallBack', comTopoHistPlot,  'separator', 'on');
    uimenu( submenu, 'Label', 'Scroll channel data with AMICA model (beta)',           'CallBack', comBackgEegplot);
    uimenu( submenu, 'Label', 'Scroll IC activation with AMICA model (beta)',          'CallBack', comBackgCompplot);
    uimenu( submenu, 'Label', '[For continuous] Select data using a dominant model',   'CallBack', comSelectModel);  
