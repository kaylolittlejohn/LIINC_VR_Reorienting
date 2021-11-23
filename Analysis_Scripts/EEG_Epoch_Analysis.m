%% NiDyN Data Analysis - EEG P300 Epoch Analysis
close all
clear all

%% Parameters
subject_number = 14;
target_category = 3; % Found in Subject_Notes.txt
condition = 'eye'; % eye = Fixed condition; free = Free condition
eeg_srate = 2048;
num_eeg_chan = 89;
P300_channels = [35 45 28 27]; % Subject dependent
removed_channels = [13 22 24 55 57 59 60]; % Subject dependent
removed_ICA_components = [1 2]; % Subject dependent

%% Load Data
% Load all files of each subject/condition combination
sPath = fullfile('C:\Users\pawan\Documents\NiDyN\Summer 2019 Data\s14');
files = dir(sPath);
data = {};
for file = 1:length(files)
    fn = files(file).name;
    if contains(fn,condition)
        data{file} = load_xdf(fullfile(sPath,fn));
    end
end

%% Separate Out EEG Data
eeg_data = {};
eeg_time = {};
for i = 1:length(data)   
    for j = 1:length(data{i})
        if strcmp(data{i}{j}.info.type,'EEG')
            eeg_data{i} = data{i}{j}.time_series;
            eeg_time{i} = data{i}{j}.time_stamps;
        end
    end
end

%% Filter EEG Data from .5-50 Hz - Butterworth Fourth Order
% Design the butterworth filter
[bb, aa] = butter(2,[0.5 50]./(eeg_srate/2));
% Check the filter magnitude and phase response
% freqz ( bb, aa, 10000, eeg_srate )

eeg_filt = {};
for i = 1:length(eeg_data)
    if size(eeg_data{i},1) == num_eeg_chan
        for chan = 1:num_eeg_chan
            eeg_filt{i}(chan,:) = filtfilt(bb,aa,double(eeg_data{i}(chan,:)));
        end
    end
end

%% Separate Out Unity Data
unity_data = {};
unity_time = {};
for i = 1:length(data)   
    for j = 1:length(data{i})
        if strcmp(data{i}{j}.info.type,'object_info')
            unity_data{i} = data{i}{j}.time_series;
            unity_time{i} = data{i}{j}.time_stamps;
        end
    end
end

%% Match Unity Time to EEG Time and Add Event Markers to the Last Row of EEG Data
eeg_w_events = {};
for i = 1:length(data)
    if ~isempty(unity_data{i})
        eeg_w_events{i} = cat(1,eeg_filt{i},zeros(1,length(eeg_filt{i})));
        for unity_ind = 1:length(unity_data{i})
            if unity_data{i}(1,unity_ind) ~= 0
                [~,eeg_ind] = min(abs(unity_time{i}(unity_ind) - eeg_time{i}));
                eeg_w_events{i}(end,eeg_ind) = unity_data{i}(1,unity_ind);
            end
        end
    end
end

%% Concatenate All EEG Data (w/ Events) Together
eeg_events_all = [];
for i = 1:length(eeg_w_events)
    eeg_events_all = cat(2,eeg_events_all,eeg_w_events{i});
end

%% EEGLAB 
% Select only the EEG channels + events for EEGLAB
eeg_events = eeg_events_all(2:65,:); 
eeg_events = cat(1,eeg_events,eeg_events_all(end,:));


% Input data into EEGLAB
eeglab redraw; 
EEG = pop_importdata('dataformat','array','nbchan',0,'data','eeg_events','srate',eeg_srate,'pnts',0,'xmin',0);
EEG = pop_chanevent(EEG,65,'edge','leading','edgelen',0); 
EEG = pop_editset(EEG, 'chanlocs', 'C:\Users\pawan\Documents\MATLAB\Tanner\biosemi64.sph'); 
EEG = eeg_checkset( EEG );

% Remove the same channels that was removed in EEG_ICA_Analysis 
EEG = pop_select( EEG, 'nochannel',removed_channels);
EEG = eeg_checkset( EEG );

% Re-reference EEG data
EEG = pop_reref(EEG,[]);
EEG = eeg_checkset (EEG);

% Load the saved ICA weights
load(sprintf('s%i_ICA_Weights.mat',subject_number));
EEG.icaweights = ICA_weights;
EEG.icasphere = ICA_sphere;
EEG = eeg_checkset( EEG );

% Remove ICA components
EEG = pop_subcomp( EEG, removed_ICA_components, 0);
EEG = eeg_checkset( EEG );

% Epoch the EEG data and remove baseline
EEG = pop_epoch( EEG, {  }, [-0.5 1.0], 'epochinfo', 'yes'); 
EEG = pop_rmbase( EEG, [-200    0]); 
EEG = eeg_checkset( EEG );

% Epoch rejection
[EEG, aThresh] = pop_eegthresh(EEG,1,1:EEG.nbchan,-200,200,EEG.xmin, EEG.xmax - 1/EEG.srate,1,0);
[EEG, ~, ~, iProb] = pop_jointprob(EEG,1,1:EEG.nbchan,5,5,1,0);
[EEG, ~, ~, iKurt] = pop_rejkurt(EEG,1,1:EEG.nbchan,5,5,1,0);
EEG = eeg_rejsuperpose (EEG,1,1,1,1,1,1,1,1);
EEG = pop_rejepoch(EEG,find(EEG.reject.rejglobal),0);
EEG = eeg_checkset( EEG );

%% Separate Epoched Data into Targets and Distractors
eeg_epoch = EEG.data;
event = EEG.event;
for i=1:size(event,2)
    event_type(i) = event(i).type;   
end

eeg_epoch_filt_targ = eeg_epoch(:,:,event_type == target_category);
eeg_epoch_filt_dist = eeg_epoch(:,:,event_type ~= target_category);

%% Plot Average Data for Fz, Cz, Pz, POZ
x_axis = linspace(-500,1000,size(eeg_epoch_dist, 2));
channel_title = {'Fz','Cz','Pz','POz'};
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:4
    subplot(2,2,i)
    channel = P300_channels(i);
    Dist = shadedErrorBar(x_axis,mean(eeg_epoch_dist(channel,:,:),3),std(eeg_epoch_dist(channel,:,:),[],3),'-b',1);
    hold on
    Targ = shadedErrorBar(x_axis,mean(eeg_epoch_targ(channel,:,:),3),std(eeg_epoch_targ(channel,:,:),[],3),'-r',1);
    plot(x_axis,median(eeg_epoch_dist(channel,:,:),3),'b.')
    plot(x_axis,median(eeg_epoch_targ(channel,:,:),3),'r.')
    plot(x_axis,mean(eeg_epoch_targ(channel,:,:),3)-mean(eeg_epoch_dist(channel,:,:),3),'--k','LineWidth',2)
    title(channel_title{i})
    xlabel('Time (ms)')
    legend([Dist.mainLine, Targ.mainLine], 'distractors', 'targets', 'Location', 'SouthWest')
end
fn_save = sprintf('%s_%s_EEG_P300.png',subject_number,condition);
saveas(gcf,fn_save)





