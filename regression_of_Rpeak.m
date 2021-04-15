clc; clear;
ts=[-1000 2000];     % <--- input in ms
tecg=[-30 30];
cd('D:\eegexp27/behavioral_excFAoutliers')
pbfile=dir; pb=pbfile([pbfile.isdir]); pb(1:2)=[]; subnum=numel(pb);
save2folder='D:\eegexp27\eegdata9\filter05\reref\regress_ecg\Rpeak_late400_detmiss';
mkdir(save2folder); cd(save2folder);
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
conds={'Rpeak_hit', 'Rpeak_miss'};
condsnum=length(conds);
%%
eeglab
EEGpath='D:\eegexp27\eegdata9\filter05';
for s=1:subnum
    EEGname=['VP' pb(s).name '.set'];
    EEG = pop_loadset(EEGname, EEGpath);
    ECG = pop_select( EEG,'channel',{'E'});
    EEG = pop_subcomp(EEG, [EEG.comprej], 0); %remove bad ica components
    EEG = pop_select( EEG,'nochannel',{'VEOG', 'E'}); % remove eye and ecg electrode
    
    %rereference to common average
    EEG = eeg_checkset(EEG,'eventconsistency');
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref( EEG, []);
    
    % add ECG back to EEG structure and replace initialReference
    EEG.data(end,:)= ECG.data;
    EEG.chanlocs(1,end).labels = 'ECG';
    EEG.chanlocs(1,end).ref = 'FCZ';
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    load(['D:\eegexp27\eegdata9\Rpeak_late400_detmiss\rlatency' pb(s).name '.mat']); % Rpeak latency preceding hits and misses
    r_hit=cell2mat(rlatency.hit); r_miss=cell2mat(rlatency.miss);
    numtrial=min([length(r_hit) length(r_miss)]); r_hit=r_hit(1:numtrial); r_miss=r_miss(1:numtrial); %choose equal number of trials
    
    EEG = eeg_addnewevents(EEG, {r_hit}, {'Rpeak_hit'});
    EEG=eeg_checkset(EEG,'eventconsistency');
    EEG = eeg_addnewevents(EEG, {r_miss}, {'Rpeak_miss'});
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    EEG = pop_epoch(EEG,{'Rpeak_hit', 'Rpeak_miss'},ts/1000,'newname','1 pruned with ICA epochs','epochinfo','yes');
    EEG = pop_rmbase(EEG,[-100 0]);
    
    [~,ind1]=min(abs(EEG.times-tecg(1)));
    [~, ind2]=min(abs(EEG.times-tecg(2)));
    
    eegc=EEG.data(:,ind1:ind2,:);
    % EKG channel index
    ek = strcmpi({EEG.chanlocs.labels},'ECG');
    
    % initialize residual data
    EEG.rdata = zeros(size(EEG.data));
    
    % loop over trials
    for triali=1:EEG.trials
        
        % build the least-squares model as intercept and EKG from this trial
        Xeeg = [ ones(EEG.pnts,1) EEG.data(ek,:,triali)' ];
        X = [ ones(size(eegc,2),1) eegc(ek,:,triali)' ];
        
        % compute regression coefficients for all channels simultaneously
        b = (X'*X) \ (X'*eegc(:,:,triali)');
        
        % predicted data
        yHat = Xeeg*b;
        
        % new data are the residuals after projecting out the best EKG fit
        EEG.rdata(:,:,triali) = ( EEG.data(:,:,triali)' - yHat )';
    end
    
%     a visual aside for subject average:
%     chan2plot = 15;
%     plot(EEG.times,mean(EEG.data(chan2plot,:,:),3), EEG.times,yHat(:,chan2plot)', EEG.times, mean_res(chan2plot,:)','linew',3)
%     legend({'C4';'Fit';'resid'})
%     plot(EEG.times, mean_res(:,chan2plot),'linew',3)
%     plot(EEG.times,EEG.data(chan2plot,:,triali), EEG.times,yHat(:,chan2plot), EEG.times,EEG.data(chan2plot,:,triali)'-yHat(:,chan2plot),'linew',3)
%     legend({'C4';'Fit';'resid'})
%     plot(EEG.times,mean(EEG.rdata(chan2plot,:,:),3),'linew',3)
%     
    EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
    EEG.data = EEG.rdata;
    EEG = rmfield(EEG, 'rdata');
    
    for c=1:condsnum
        [EEG_cond,ind] = pop_selectevent(EEG,'type',conds{c});
        EEG_cond = pop_editset(EEG_cond, 'subject', pb(s).name, 'condition', conds{c});
        %pop_saveset(EEG_cond, 'filename', ['VP',pb(s).name,conds{c}], 'filepath', save2folder);
        
        %convert eeglab to fieldtrip
        data = eeglab2fieldtrip(EEG_cond, 'timelockanalysis','chan_loc');
        data.dimord = 'chan_time';
        data = ft_timelockbaseline(cfg2, data);
        data.fsample = EEG_cond.srate;
        enum = size(EEG_cond.data,3);
        grandtrials{s,c} = data;
        epochnum(s,c) = enum;
        clear EEG_cond data enum
        clear EEG_cond
    end
end
save('grandtrials','grandtrials','conds','subnum','epochnum')

