% provides consistency between LLt and v, and between LLt_smooth and v
% stored under EEG.etc.amica
function EEG = eeg_checkamica(EEG)

if isfield(EEG.etc,'amica') && ~isempty(EEG.etc.amica)
    if isfield(EEG.etc.amica,'LLt') && ~isempty(EEG.etc.amica.LLt) 
        
        if size(EEG.data,2) ~= size(EEG.etc.amica.LLt,2) || size(EEG.data,3) ~= size(EEG.etc.amica.LLt,3)
            EEG.etc.amica.LLt = findLLt(EEG,EEG.etc.amica);
            EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
        else
            if size(EEG.data,2) ~= size(EEG.etc.amica.v,2) || size(EEG.data,3) ~= size(EEG.etc.amica.v,3)
                 EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
            end
        end
        
        
        
    end
    if isfield(EEG.etc.amica,'LLt_smooth') && ~isempty(EEG.etc.amica.LLt_smooth)
        if ~isfield(EEG.etc.amica,'v_smooth') | ~isequal(size(EEG.etc.amica.LLt_smooth),size(EEG.etc.amica.v_smooth))
            EEG.etc.amica.v_smooth = LLt2v(EEG.etc.amica.LLt_smooth);
        end
    end
end