% pop_smoothprobamica() - postAmicaUtility plugin. (Re)computes smoothed probabilities and
%                         likelihoods (v_smooth, LLt_smooth).

% Author: Ozgur Balkan. SCCN, INC, UCSD.

% History
% 07/24/2018 Makoto. eegh problem fixed (Thanks Kevin Verdiere!)
% 12/26/2016 Makoto. Modified.

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Ozgur Balkan
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

function [EEG, com] = pop_smoothprobamica(EEG, smoothlength)

com = '';

% From runamica15()
%   mods.W(:,:,h) is the unmixing matrix for model h.
%   mods.S is the sphering matrix.
%   mods is a structure containing the output components and density models.
%   mods.A(:,:,h) is the components for model h.
%   mods.varord(:,h) is the index order of the components in variance order.
%   mods.Lht is the likelihood of time points for each model (if set).
%   mods.LL is the history of the log likelihood over iterations.
%   mods.c(:,h) is the center for model h.

if ~isfield(EEG.etc,'amica') || ~isfield(EEG.etc.amica,'LLt')
    error('Please check if probabilities are under EEG.etc.amica')
end

if nargin < 2
    prompt = {'Enter the smoothing Hanning window size (in sec): '};
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
    fprintf('Smoothing probabilities...')
    [EEG.etc.amica.v_smooth,EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,str2num(smoothlength{1}));
    fprintf('Done \n')
    EEG.etc.amica.smooth_length = str2num(smoothlength{1});
    smoothlength = str2num(smoothlength{1});
    
else
    fprintf('Smoothing probabilities...')
    [EEG.etc.amica.v_smooth,EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlength);
    fprintf('Done \n')
    EEG.etc.amica.smooth_length = smoothlength;
end

com = sprintf('%s = pop_smoothprobamica(%s,%d)',inputname(1),inputname(1),smoothlength);
