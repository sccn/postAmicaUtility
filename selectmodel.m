function [EEG]= selectmodel(EEG, model, type)
startpoints = [];
endpoints = [];

if ~isequal(size(EEG.data,2),size(EEG.etc.amica.LLt,2)) || ~isequal(size(EEG.data,3),size(EEG.etc.amica.LLt,3))
        recompute_prob = questdlg('Data size and AMICA probabilities do not match. Would you like to recompute the probabilities?', ...
        'Recompute probabilities?','Yes','No','No');
     emptyStr = '';
     if strcmpi(recompute_prob,'No') || strcmpi(recompute_prob,emptyStr)
         disp('Cannot do the rejection without correct size of proabability matrix');
         return;
     else
        EEG.etc.amica.LLt = findLLt(EEG,EEG.etc.amica);
        EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
        EEG = eeg_checkamica(EEG);
     end
end

prompt = {'Enter the smoothing Hanning window size: '};
name = 'Input for probability smoothing for AMICA';
numlines = 1;
defaultanswer = {'2'};
smoothlength = inputdlg(prompt,name,numlines,defaultanswer);
[v_smooth,~] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,str2num(smoothlength{1}));


for i = 1:EEG.nbchan
    EEG.icaweights(i,:) = EEG.etc.amica.W(i,:,model);
    EEG.icasphere(i,:)  = EEG.etc.amica.S(i,:);
    EEG.icawinv(i,:)    = EEG.etc.amica.A(i,:,model);
end

disp('Calculating...')


if type == 1  % reject the data where the model is not active
    if EEG.trials > 1
        trials2keep = [];
        for i = 1:EEG.trials
            [dummy in] = max(sum(v_smooth(:,:,i),2));
            if model == in
                trials2keep = [trials2keep i];
            end
        end
        
        EEG = pop_select(EEG,'trial',trials2keep);
        
    else
        nepochs = 0; %number of epochs that will be excluded
        rejecting = 0;
        startpoints = [];
        endpoints = [];
        for i = 1:size(EEG.data,2)
            
            [y in] = max(v_smooth(:,i));
            if in ~= model
                if ~rejecting
                    nepochs = nepochs + 1;
                    startpoints(nepochs) = i;
                    rejecting = 1;
                end
            else
                if rejecting == 1
                    endpoints(nepochs) = i-1;
                    rejecting = 0;
                end
                
            end
            
        end
        if rejecting == 1
            endpoints(end + 1) = size(EEG.data,2);
        end
        
        
        if nepochs == 1 && endpoints == i && startpoints == 1
            error('Selected model is never active. Original data will be preserved.');
            startpoints = [];
            endpoints = [];
            
            return;
        end
        if nepochs == 0
            disp('Model is always active. Original data will be preserved, ICA weights for the model is loaded.');
            
            return;
        end
        
        %-------Rejecting the data where model is not active------------------
        temp_start = startpoints;
        temp_end = endpoints;
        EEG = pop_select(EEG,'nopoint',[startpoints' endpoints']);
%         for i = 1:nepochs
%             newEEG = pop_ozselect(newEEG,'nopoint',[temp_start(i) temp_end(i)]);
%             shift = temp_end(i)-temp_start(i)+1;
%             temp_start = temp_start - shift;
%             temp_end = temp_end - shift;
%         end
        %---------------------------------------------------------------------
    end
    EEG.setname = [EEG.setname ' - Model ' num2str(model)];
    
else
    if type == 2
        %weight the data w.r.t to the probability of the model chosen
        if isequal(EEG.etc.amica.v_smooth(model,:,:),zeros(size(EEG.etc.amica.v_smooth(model,:,:))))
           disp('Selected model has probability of 0 through the data, old dataset is preserved') ;
        else
            if EEG.trials > 1
                for i = 1:EEG.trials
                    for j = 1:size(EEG.data,2)
                        EEG.data(:,j,i) = EEG.data(:,j,i)*EEG.etc.amica.v_smooth(model,j,i);
                        
                    end
                end
                
            else
                for i = 1:size(EEG.data,2)
                    EEG.data(:,i) = EEG.data(:,i)*EEG.etc.amica.v_smooth(model,i);
                end
                
            end
            EEG.etc.amica.weighted = 1;
        end
        EEG.setname = [EEG.setname ' - Model ' num2str(model) ' Weighted'];
    else
        if type == 3
            
        end
        
    end
    
end









