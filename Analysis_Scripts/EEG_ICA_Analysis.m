%% NiDyN Data Analysis - EEG ICA Analysis
close all
clear all

%% Parameters
eeg_srate = 2048;
subject_number = 8;
num_eeg_chan = 89;

%% Load Data
sPath = fullfile('C:\Users\pawan\Documents\NiDyN\Summer 2019 Data\s8');
files = dir(sPath);
data = {};
for file = 1:length(files)
    fn = files(file).name;
    if contains(fn,'.xdf')
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

%% Filter EEG Data from 2-50 Hz - Butterworth Fourth Order
% Design the butterworth filter
[bb, aa] = butter(4,[2 50]./(eeg_srate/2));
eeg_filt = {};
for i = 1:length(eeg_data)
    if size(eeg_data{i},1) == num_eeg_chan
        for chan = 1:num_eeg_chan
            eeg_filt{i}(chan,:) = filtfilt(bb,aa,double(eeg_data{i}(chan,:)));
        end
    end
end

%% Concatenate All EEG Data Together
eeg_filt_all = [];
for i = 1:length(eeg_filt)
    eeg_filt_all = cat(2,eeg_filt_all,eeg_filt{i});
end

%% Downsample the EEG Data to 256 Hz
eeg_dsrate = 256; 
eeg_dsfilt_all = [];
for chan = 1:num_eeg_chan
    temp = downsample(eeg_filt_all(chan,:),eeg_srate/eeg_dsrate);
    eeg_dsfilt_all = cat(1,eeg_dsfilt_all,temp);
end

%% EEGLAB
% Select the first 64 EEG channels for EEGLAB
eeg_dsfilt_eeg = eeg_dsfilt_all(2:65,:);

% Input data into EEGLAB
eeglab redraw;
EEG = pop_importdata('dataformat','array','nbchan',64,'data','eeg_dsfilt_eeg','srate',256,'pnts',0,'xmin',0);
EEG = pop_editset(EEG, 'chanlocs', 'C:\Users\pawan\Documents\NiDyN\Summer 2019 Data\biosemi64.sph');
EEG = eeg_checkset( EEG );
EEG = pop_rejchan(EEG,'threshold',5,'norm','on','measure','kurt');
EEG = pop_rejchan (EEG,'threshold',4,'norm','on','measure','prob');
EEG = eeg_checkset( EEG );
EEG = eeg_regepochs(EEG, 0.5, [0 0.5], NaN);
EEG = eeg_checkset( EEG );
[EEG, aThresh] = pop_eegthresh(EEG,1,1:EEG.nbchan,-150,150,EEG.xmin, EEG.xmax - 1/EEG.srate,1,0);
[EEG, ~, ~, iProb] = pop_jointprob(EEG,1,1:EEG.nbchan,4,4,1,0);
[EEG, ~, ~, iKurt] = pop_rejkurt(EEG,1,1:EEG.nbchan,4,4,1,0);
EEG = eeg_rejsuperpose (EEG,1,1,1,1,1,1,1,1);
EEG = pop_rejepoch(EEG,find(EEG.reject.rejglobal),0);
EEG = eeg_epoch2continuous(EEG);
EEG = eeg_checkset (EEG);
EEG = pop_reref(EEG,[]);
EEG = eeg_checkset(EEG);
EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',size((EEG.data),1));
EEG = eeg_checkset(EEG);

% Save the ICA weights and sphere
ICA_weights = EEG.icaweights;
ICA_sphere = EEG.icasphere;
save(sprintf('s%i_ICA_Weights.mat',subject_number),'ICA_weights','ICA_sphere');

ALLEEG = EEG;

eeglab redraw
