%% NiDyN Data Analysis - Pupil Dilation Analysis
close all
clear all

%% Parameters
eeg_srate = 2048;
eye_srate = 120;
unity_srate = 75;
target_category = 4; %camera piano laptop boat
subject_number = 16;
num_eeg_chan = 89;
condition = 'free';

%% Load Data
sPath = fullfile('C:\Users\Pawan Lapborisuth\Documents\MATLAB\NiDyN\Summer 2019 Data\s16');
files = dir(sPath);
data = {};
run_counter = 0;
for file = 1:length(files)
    fn = files(file).name;
    if contains(fn,condition)
        data_file = load_xdf(fullfile(sPath,fn));
        if ~isempty(data_file)
            run_counter = run_counter + 1;
            data{run_counter} = data_file;
        end
    end
end
fprintf('Successfully loaded subject %i data %s condition with %i runs',subject_number,condition,run_counter);

%% Separate Out Pupil Data
pupil_data = {};
pupil_time = {};
for i = 1:length(data)
    for j = 1:length(data{i})
        if strcmp(data{i}{j}.info.name,'NEDE_Pupil')
            pupil_data{i} = data{i}{j}.time_series;
            pupil_time{i} = data{i}{j}.time_stamps;
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

%% Separate Out Head Movement Data
head_data = {};
head_time = {};
for i = 1:length(data)
    for j = 1:length(data{i})
        if strcmp(data{i}{j}.info.name,'NEDE_HeadPosition')
            head_data{i} = data{i}{j}.time_series;
            head_time{i} = data{i}{j}.time_stamps;
        end
    end
end

%% Separate Out Left and Right Pupil Data and Their Validity Signal
left_pupil = {};
left_valid = {};
right_pupil = {};
right_valid = {};
for i = 1:length(pupil_data)
    if ~isempty(pupil_data{i})
        left_pupil{i} = pupil_data{i}(1,:);
        left_valid{i} = pupil_data{i}(2,:);
        right_pupil{i} = pupil_data{i}(6,:);
        right_valid{i} = pupil_data{i}(7,:);
    end
end

%% Separate Out Right and Left Gaze Data
% Will use right gaze for analysis but have left gaze data just in case
right_gaze_data = {};
right_gaze_time = {};
left_gaze_data = {};
left_gaze_time = {};
for i = 1:length(data)
    for j = 1:length(data{i})
        if strcmp(data{i}{j}.info.name,'NEDE_RightGazeData')
            right_gaze_data{i} = data{i}{j}.time_series;
            right_gaze_time{i} = data{i}{j}.time_stamps;
        elseif strcmp(data{i}{j}.info.name,'NEDE_LeftGazeData')
            left_gaze_data{i} = data{i}{j}.time_series;
            left_gaze_time{i} = data{i}{j}.time_stamps;
        end
    end
end

%% Match Head Movement Data with Unity Data and Add Event Markers
% Eye condition should result in close to zero markers created
head_rotation_events = {};

% Loop through all runs to pull out head rotation data and add event
% markers to the second row
for i = 1:length(head_data)
    head_rotation_events{i} = zeros(2,length(head_data{i}));
    % Separate out head rotation data
    head_rotation_events{i}(1,:) = head_data{i}(5,:);
    % Add event markers to head rotation data
    for unity_ind = 1:length(unity_data{i})
        if unity_data{i}(1,unity_ind) ~= 0
            [~,head_ind] = min(abs(unity_time{i}(unity_ind) - head_time{i}));
            head_rotation_events{i}(end,head_ind) = unity_data{i}(1,unity_ind);
        end
    end
end


%% Match Right Horizontal Gaze Data with Unity Data and Add Event Markers
right_horizontal_gaze_events = {};
for i = 1:length(right_gaze_data)
    right_horizontal_gaze_events{i} = zeros(2,length(right_gaze_data{i}));
    % Separate out horizontal gaze data
    right_horizontal_gaze_events{i}(1,:) = right_gaze_data{i}(1,:);
    % Extrapolate through NaNs in gaze data
    right_horizontal_gaze_events{i}(1,:) = fillgaps(right_horizontal_gaze_events{i}(1,:));
    % Add event markers to horizontal gaze data
    for unity_ind = 1:length(unity_data{i})
        if unity_data{i}(1,unity_ind) ~= 0
            [~,gaze_ind] = min(abs(unity_time{i}(unity_ind) - right_gaze_time{i}));
            right_horizontal_gaze_events{i}(end,gaze_ind) = unity_data{i}(1,unity_ind);
        end
    end
end


%% Zero Out Left and Right Pupil Data for Invalid Signal
for i = 1:length(pupil_data)
    if ~isempty(pupil_data{i})
        for ind = 1:length(pupil_data{i})
            if left_valid{i}(ind) == 0
                left_pupil{i}(ind) = NaN;
            end
            if right_valid{i}(ind) == 0
                right_pupil{i}(ind) = NaN;
            end
        end
    end
end

%% Match Pupil Data with Unity Data and Add Event Markers
left_pupil_events = {};
right_pupil_events = {};
for i = 1:length(unity_data)
    if ~isempty(unity_data{i})
        left_pupil_events{i} = cat(1,left_pupil{i},zeros(1,length(left_pupil{i})));
        right_pupil_events{i} = cat(1,right_pupil{i},zeros(1,length(right_pupil{i})));
        for unity_ind = 1:length(unity_data{i})
            if unity_data{i}(1,unity_ind) ~= 0
                [~,pupil_ind] = min(abs(unity_time{i}(unity_ind) - pupil_time{i}));
                left_pupil_events{i}(end,pupil_ind) = unity_data{i}(1,unity_ind);
                right_pupil_events{i}(end,pupil_ind) = unity_data{i}(1,unity_ind);
            end
        end
    end
end



%% Convert Pupil Data into Percentage Change from Mean Value by Run
left_pupil_events_pc = {};
right_pupil_events_pc = {};
for i = 1:length(unity_data)
    if ~isempty(unity_data{i})
        left_pupil_events_pc{i} = left_pupil_events{i};
        left_pupil_events_pc{i}(1,:) = (left_pupil_events{i}(1,:)-nanmean(left_pupil_events{i}(1,:)))/nanmean(left_pupil_events{i}(1,:))*100;
        right_pupil_events_pc{i} = right_pupil_events{i};
        right_pupil_events_pc{i}(1,:) = (right_pupil_events{i}(1,:)-nanmean(right_pupil_events{i}(1,:)))/nanmean(right_pupil_events{i}(1,:))*100;
    end
end

%% Cut Left and Right Pupil Data into Epoch for Targets and Distractors
% Left eye
left_pupil_targ = {};
left_pupil_dist = {};
for i = 1:length(pupil_data)
    if ~isempty(pupil_data)
        left_pupil_targ{i} = [];
        left_pupil_dist{i} = [];
        for ind = 1:length(pupil_data{i})
            if left_pupil_events_pc{i}(2,ind) == target_category
                if ind+3*eye_srate < length(left_pupil_events_pc{i})
                    left_pupil_targ{i} = cat(1,left_pupil_targ{i},left_pupil_events_pc{i}(1,ind-0.2*eye_srate:ind+3*eye_srate));
                end
            elseif left_pupil_events_pc{i}(2,ind) ~= target_category && left_pupil_events_pc{i}(2,ind) ~= 0
                if ind+3*eye_srate < length(left_pupil_events_pc{i})
                    left_pupil_dist{i} = cat(1,left_pupil_dist{i},left_pupil_events_pc{i}(1,ind-0.2*eye_srate:ind+3*eye_srate));
                end
            end     
        end
    end
end

% Right eye - Note: Do gaze events with right eye because it is derived
% from right horizontal gaze but can be applied to both pupil data
%%%%NOTE: Pupil is cut into -0.2-3 seconds epoch and gaze into 0-1.5 seconds
right_pupil_targ = {};
right_pupil_dist = {};
right_horizontal_gaze_targ = {};
right_horizontal_gaze_dist = {};
head_rotation_targ = {};
head_rotation_dist = {};
for i = 1:length(pupil_data)
    if ~isempty(pupil_data)
        right_pupil_targ{i} = [];
        right_pupil_dist{i} = [];
        right_horizontal_gaze_targ{i} = [];
        right_horizontal_gaze_dist{i} = [];
        head_rotation_targ{i} = [];
        head_rotation_dist{i} = [];
        for ind = 1:length(pupil_data{i})
            if right_pupil_events_pc{i}(2,ind) == target_category
                if ind+3*eye_srate < length(right_pupil_events_pc{i})
                    right_pupil_targ{i} = cat(1,right_pupil_targ{i},right_pupil_events_pc{i}(1,ind-0.2*eye_srate:ind+3*eye_srate));
                    right_horizontal_gaze_targ{i} = cat(1,right_horizontal_gaze_targ{i},right_horizontal_gaze_events{i}(1,ind:ind+1.5*eye_srate-1));
                    head_rotation_targ{i} = cat(1,head_rotation_targ{i},head_rotation_events{i}(1,ind:ind+1.5*eye_srate-1));
                end
            elseif right_pupil_events_pc{i}(2,ind) ~= target_category && right_pupil_events_pc{i}(2,ind) ~= 0
                if ind+3*eye_srate < length(right_pupil_events_pc{i})
                    right_pupil_dist{i} = cat(1,right_pupil_dist{i},right_pupil_events_pc{i}(1,ind-0.2*eye_srate:ind+3*eye_srate));
                    right_horizontal_gaze_dist{i} = cat(1,right_horizontal_gaze_dist{i},right_horizontal_gaze_events{i}(1,ind:ind+1.5*eye_srate-1));
                    head_rotation_dist{i} = cat(1,head_rotation_dist{i},head_rotation_events{i}(1,ind:ind+1.5*eye_srate-1));                    
                end
            end
        end
    end
end

%% Concatenate All Target and Distractor Epochs Together From All Runs
% Pupil Data
left_pupil_targ_all = [];
left_pupil_dist_all = [];
right_pupil_targ_all = [];
right_pupil_dist_all = [];
for i = 1:length(pupil_data)
    left_pupil_targ_all = cat(1,left_pupil_targ_all,left_pupil_targ{i});
    left_pupil_dist_all = cat(1,left_pupil_dist_all,left_pupil_dist{i});
    right_pupil_targ_all = cat(1,right_pupil_targ_all,right_pupil_targ{i});
    right_pupil_dist_all = cat(1,right_pupil_dist_all,right_pupil_dist{i});
end

% Gaze Data
right_horizontal_gaze_targ_all = [];
right_horizontal_gaze_dist_all = [];
for i = 1:length(right_horizontal_gaze_events)
    right_horizontal_gaze_targ_all = cat(1,right_horizontal_gaze_targ_all,right_horizontal_gaze_targ{i});
    right_horizontal_gaze_dist_all = cat(1,right_horizontal_gaze_dist_all,right_horizontal_gaze_dist{i});
end

% Head Data
head_rotation_targ_all = [];
head_rotation_dist_all = [];
for i = 1:length(head_rotation_events)
    head_rotation_targ_all = cat(1,head_rotation_targ_all,head_rotation_targ{i});
    head_rotation_dist_all = cat(1,head_rotation_dist_all,head_rotation_dist{i});  
end

%% ONLY for Free Condition Add Head Data to Horizontal Gaze Data
if strcmp(condition,'free')
    for i = 1:size(right_horizontal_gaze_targ_all,1)
        right_horizontal_gaze_targ_all(i,:) = right_horizontal_gaze_targ_all(i,:) + head_rotation_targ_all(i,:);    
    end
    for j = 1:size(right_horizontal_gaze_dist_all,1)
        right_horizontal_gaze_dist_all(i,:) = right_horizontal_gaze_dist_all(i,:) + head_rotation_dist_all(i,:);
    end
end

%% Piece Wise Linear Modeling of Right Horizontal Gaze Movement
max_change_pts = 4;
% For plotting
how_many_trials = 20;
% Targets
saccade_fixation_index_target = zeros(size(right_horizontal_gaze_targ_all,1),4);
saccade_fixation_index_target_validity = zeros(size(right_horizontal_gaze_targ_all,1),1);
for trial = 1:size(right_horizontal_gaze_targ_all,1)
    y = right_horizontal_gaze_targ_all(trial,:);
    [ind,residuals] = findchangepts(y,'Statistic','linear','MaxNumChanges',max_change_pts);
    saccade_fixation_valid = 0;
    if length(ind) == max_change_pts
        % Marking trials where piecewise modeling does not work based on
        % saccades direction -> The two saccades must be of opposite
        % directions
        saccade_slope = zeros(1,2);
        saccade_count = 0;
        for saccade_ind = 1:length(ind)
            if saccade_ind == 1 || saccade_ind == 3
                saccade_count = saccade_count + 1;
                y_pw = y(ind(saccade_ind):ind(saccade_ind+1));
                x_pw = linspace(1,length(y_pw),length(y_pw));
                p = polyfit(x_pw,y_pw,1);
                saccade_slope(saccade_count) = p(1);
            end
        end
        % Marking trials where piecewise modeling does not work based on
        % saccades length -> Fixation length must be longer than the length
        % of the shorter saccades
        saccade_fixation_length = diff(ind);
        fixation_length = saccade_fixation_length(2);
        saccade_length = min([saccade_fixation_length(1) saccade_fixation_length(2)]);
        if saccade_length < fixation_length
            saccade_fixation_valid = 1;
        end
        if sum(sign(saccade_slope)) == 0 && saccade_fixation_valid == 1
            saccade_fixation_index_target_validity(trial) = 1;
        end
        saccade_fixation_index_target(trial,:) = ind;
    end
    % Example plots
    if trial <= how_many_trials
        figure
        findchangepts(y,'Statistic','linear','MaxNumChanges',max_change_pts)
        if saccade_fixation_index_target_validity(trial) ~= 1
            title('BAD TRIAL')
        end
    end
end

% Distractors
saccade_fixation_index_distractor = zeros(size(right_horizontal_gaze_dist_all,1),4);
saccade_fixation_index_distractor_validity = zeros(size(right_horizontal_gaze_dist_all,1),1);
for trial = 1:size(right_horizontal_gaze_dist_all,1)
    y = right_horizontal_gaze_dist_all(trial,:);
    [ind,residuals] = findchangepts(y,'Statistic','linear','MaxNumChanges',max_change_pts);
    saccade_fixation_valid = 0;
    if length(ind) == max_change_pts
        % Marking trials where piecewise modeling does not work based on
        % saccades direction -> The two saccades must be of opposite
        % directions
        saccade_slope = zeros(1,2);
        saccade_count = 0;
        for saccade_ind = 1:length(ind)
            if saccade_ind == 1 || saccade_ind == 3
                saccade_count = saccade_count + 1;
                y_pw = y(ind(saccade_ind):ind(saccade_ind+1));
                x_pw = linspace(1,length(y_pw),length(y_pw));
                p = polyfit(x_pw,y_pw,1);
                saccade_slope(saccade_count) = p(1);
            end
        end
        % Marking trials where piecewise modeling does not work based on
        % saccades length -> Fixation length must be longer than the length
        % of the shorter saccades
        saccade_fixation_length = diff(ind);
        fixation_length = saccade_fixation_length(2);
        saccade_length = min([saccade_fixation_length(1) saccade_fixation_length(2)]);
        if saccade_length < fixation_length
            saccade_fixation_valid = 1;
        end
        if sum(sign(saccade_slope)) == 0 && saccade_fixation_valid == 1
            saccade_fixation_index_distractor_validity(trial) = 1;
        end
        saccade_fixation_index_distractor(trial,:) = ind;
    end
    % Example plots
    if trial <= how_many_trials
        figure
        findchangepts(y,'Statistic','linear','MaxNumChanges',max_change_pts)
        if saccade_fixation_index_distractor_validity(trial) ~= 1
            title('BAD TRIAL')
        end
    end
end
%% Baseline Correct Each Epoch From -200 to 0 ms
% Left eye
left_pupil_targ_base = zeros(size(left_pupil_targ_all,1),size(left_pupil_targ_all,2));
for i = 1:size(left_pupil_targ_all,1)
    left_pupil_targ_base(i,:) = left_pupil_targ_all(i,:) - mean(left_pupil_targ_all(i,1:1+0.2*eye_srate)); 
end
left_pupil_dist_base = zeros(size(left_pupil_dist_all,1),size(left_pupil_dist_all,2));
for i = 1:size(left_pupil_dist_all,1)
    left_pupil_dist_base(i,:) = left_pupil_dist_all(i,:) - mean(left_pupil_dist_all(i,1:1+0.2*eye_srate)); 
end
% Right eye
right_pupil_targ_base = zeros(size(right_pupil_targ_all,1),size(right_pupil_targ_all,2));
for i = 1:size(right_pupil_targ_all,1)
    right_pupil_targ_base(i,:) = right_pupil_targ_all(i,:) - mean(right_pupil_targ_all(i,1:1+0.2*eye_srate)); 
end
right_pupil_dist_base = zeros(size(right_pupil_dist_all,1),size(right_pupil_dist_all,2));
for i = 1:size(right_pupil_dist_all,1)
    right_pupil_dist_base(i,:) = right_pupil_dist_all(i,:) - mean(right_pupil_dist_all(i,1:1+0.2*eye_srate)); 
end

%% Plot Example of Baseline Corrected Right Pupil Percentage Change by Epoch to Observe Patterns
%{
x_axis = linspace(-0.2,3,eye_srate*3.2+1);
how_many_trials = 10;

% Targets
for i = 1:how_many_trials
    figure
    plot(x_axis,right_pupil_targ_base(i,:));
    title(sprintf('Baseline Corrected Right Pupil Percentage Change for Target Trial %i \n Subject %i %s Condition',i,subject_number,condition))
    xlim([-0.2 3])
end

% Distractors
for i = 1:how_many_trials
    figure
    plot(x_axis,right_pupil_dist_base(i,:));
    title(sprintf('Baseline Corrected Right Pupil Percentage Change for Distractor Trial %i \n Subject %i %s Condition',i,subject_number,condition))
    xlim([-0.2 3])
end
%}

%% Clean Up Pupil Data by Removing 10 Samples Surrounding Peaks (>5% change)
% Left eye
left_pupil_targ_clean = left_pupil_targ_base;
for i = 1:size(left_pupil_targ_base,1)
    peak_ind = find(abs(diff(left_pupil_targ_base(i,:))) >= 5);
    for j = 1:length(peak_ind)
        if peak_ind(j)-5 > 1 && peak_ind(j)+5 < size(left_pupil_targ_base,2)
            left_pupil_targ_clean(i,peak_ind(j)-5:peak_ind(j)+5) = NaN;
        elseif peak_ind(j)+5 > size(left_pupil_targ_base,2)
            left_pupil_targ_clean(i,peak_ind(j):end) = NaN;
            left_pupil_targ_clean(i,end) = 0;
        end
    end
end

left_pupil_dist_clean = left_pupil_dist_base;
for i = 1:size(left_pupil_dist_base,1)
    peak_ind = find(abs(diff(left_pupil_dist_base(i,:))) >= 5);
    for j = 1:length(peak_ind)
        if peak_ind(j)-5 > 1 && peak_ind(j)+5 < size(left_pupil_dist_base,2)
            left_pupil_dist_clean(i,peak_ind(j)-5:peak_ind(j)+5) = NaN;
        elseif peak_ind(j)+5 > size(left_pupil_dist_base,2)
            left_pupil_dist_clean(i,peak_ind(j):end) = NaN;
            left_pupil_dist_clean(i,end) = 0; 
        end
    end
end

% Right eye
right_pupil_targ_clean = right_pupil_targ_base;
for i = 1:size(right_pupil_targ_base,1)
    peak_ind = find(abs(diff(right_pupil_targ_base(i,:))) >= 5);
    for j = 1:length(peak_ind)
        if peak_ind(j)-5 > 1 && peak_ind(j)+5 < size(right_pupil_targ_base,2)
            right_pupil_targ_clean(i,peak_ind(j)-5:peak_ind(j)+5) = NaN;
        elseif peak_ind(j)+5 > size(right_pupil_targ_base,2)
            right_pupil_targ_clean(i,peak_ind(j):end) = NaN;
            right_pupil_targ_clean(i,end) = 0; 
        end
    end
end

right_pupil_dist_clean = right_pupil_dist_base;
for i = 1:size(right_pupil_dist_base,1)
    peak_ind = find(abs(diff(right_pupil_dist_base(i,:))) >= 5);
    for j = 1:length(peak_ind)
        if peak_ind(j)-5 > 1 && peak_ind(j)+5 < size(right_pupil_dist_base,2)
            right_pupil_dist_clean(i,peak_ind(j)-5:peak_ind(j)+5) = NaN;
        elseif peak_ind(j)+5 > size(right_pupil_dist_base,2)
            right_pupil_dist_clean(i,peak_ind(j):end) = NaN;
            right_pupil_dist_clean(i,end) = 0; 
        end
    end
end

%% Plot Example of Peak-Removed Right Pupil Percentage Change by Epoch to Observe Patterns
x_axis = linspace(-0.2,3,eye_srate*3.2+1);
how_many_trials = 10;

% Targets
for i = 1:how_many_trials
    figure
    hold on
    plot(x_axis,right_pupil_targ_base(i,:),'bx','LineWidth',2);
    %plot(x_axis,right_pupil_targ_clean(i,:));
    xline((saccade_fixation_index_target(i,1))/eye_srate,'r--');
    xline((saccade_fixation_index_target(i,2))/eye_srate,'r--');
    xline((saccade_fixation_index_target(i,3))/eye_srate,'r--');
    xline((saccade_fixation_index_target(i,4))/eye_srate,'r--');
    yLimits = get(gca,'YLim');
    h1 = fill([(saccade_fixation_index_target(i,1))/eye_srate (saccade_fixation_index_target(i,1))/eye_srate (saccade_fixation_index_target(i,2))/eye_srate (saccade_fixation_index_target(i,2))/eye_srate],[yLimits(2) yLimits(1) yLimits(1) yLimits(2)],'r');
    h2 = fill([(saccade_fixation_index_target(i,3))/eye_srate (saccade_fixation_index_target(i,3))/eye_srate (saccade_fixation_index_target(i,4))/eye_srate (saccade_fixation_index_target(i,4))/eye_srate],[yLimits(2) yLimits(1) yLimits(1) yLimits(2)],'r');
    set(h1,'facealpha',.5)
    set(h2,'facealpha',.5)
    title(sprintf('Right Pupil Percentage Change for Target Trial %i with Marked Saccades\n Subject %i %s Condition',i,subject_number,condition))
    xlim([-0.2 3])
end

% Distractors
for i = 1:how_many_trials
    figure
    hold on
    plot(x_axis,right_pupil_dist_base(i,:),'bx','LineWidth',2);
    %plot(x_axis,right_pupil_dist_clean(i,:));
    xline((saccade_fixation_index_distractor(i,1))/eye_srate,'r--');
    xline((saccade_fixation_index_distractor(i,2))/eye_srate,'r--');
    xline((saccade_fixation_index_distractor(i,3))/eye_srate,'r--');
    xline((saccade_fixation_index_distractor(i,4))/eye_srate,'r--');
    yLimits = get(gca,'YLim');
    h1 = fill([(saccade_fixation_index_distractor(i,1))/eye_srate (saccade_fixation_index_distractor(i,1))/eye_srate (saccade_fixation_index_distractor(i,2))/eye_srate (saccade_fixation_index_distractor(i,2))/eye_srate],[yLimits(2) yLimits(1) yLimits(1) yLimits(2)],'r');
    h2 = fill([(saccade_fixation_index_distractor(i,3))/eye_srate (saccade_fixation_index_distractor(i,3))/eye_srate (saccade_fixation_index_distractor(i,4))/eye_srate (saccade_fixation_index_distractor(i,4))/eye_srate],[yLimits(2) yLimits(1) yLimits(1) yLimits(2)],'r');
    set(h1,'facealpha',.5)
    set(h2,'facealpha',.5)
    title(sprintf('Right Pupil Percentage Change for Distractor Trial %i with Marked Saccades\n Subject %i %s Condition',i,subject_number,condition))
    xlim([-0.2 3])
end

%% Interpolate Through Missing Values
% Left eye
left_pupil_targ_interp = left_pupil_targ_clean;
for i = 1:size(left_pupil_targ_clean,1)
    left_pupil_targ_interp(i,:) = fillmissing(left_pupil_targ_clean(i,:),'linear');
end
left_pupil_dist_interp = left_pupil_dist_clean;
for i = 1:size(left_pupil_dist_clean,1)
    left_pupil_dist_interp(i,:) = fillmissing(left_pupil_dist_clean(i,:),'linear');
end

% Right eye
right_pupil_targ_interp = right_pupil_targ_clean;
for i = 1:size(right_pupil_targ_clean,1)
    right_pupil_targ_interp(i,:) = fillmissing(right_pupil_targ_clean(i,:),'linear');
end
right_pupil_dist_interp = right_pupil_dist_clean;
for i = 1:size(right_pupil_dist_clean,1)
    right_pupil_dist_interp(i,:) = fillmissing(right_pupil_dist_clean(i,:),'linear');
end

%% Plot Example of Interpolated Right Pupil Percentage Change by Epoch to Observe Patterns
x_axis = linspace(-0.2,3,eye_srate*3.2+1);
how_many_trials = 10;

% Targets
figure
hold on
for i = 1:how_many_trials
    plot(x_axis,right_pupil_targ_interp(i,:));
end
title(sprintf('Interpolated Right Pupil Percentage Change for Target Trials \n Subject %i %s Condition',subject_number,condition))
xlim([-0.2 3])

% Distractors
figure
hold on
for i = 1:how_many_trials
    plot(x_axis,right_pupil_dist_interp(i,:));
end
title(sprintf('Interpolated Right Pupil Percentage Change for Distractor Trials \n Subject %i %s Condition',subject_number,condition))
xlim([-0.2 3])

%% Find Mean Value Across Trials and Plot Pupil Data
x_axis = linspace(-200,3000,size(left_pupil_targ_base,2));
% Left eye
figure
hold on
Dist = shadedErrorBar(x_axis,nanmean(left_pupil_dist_interp,1),nanstd(left_pupil_dist_interp,1),'-b',1);
Targ = shadedErrorBar(x_axis,nanmean(left_pupil_targ_interp,1),nanstd(left_pupil_targ_interp,1),'-r',1);
plot(x_axis,nanmean(left_pupil_dist_base,1),'b')
plot(x_axis,nanmean(left_pupil_targ_base,1),'r')
ylabel('Percentage Change')
xlabel('Time (ms)')
xlim([-200 3000])
title(sprintf('Left Eye Pupil Data Subject %i %s condition',subject_number,condition))
legend([Dist.mainLine, Targ.mainLine], 'Distractors', 'Targets', 'Location', 'SouthWest')
saveas(gcf,sprintf('s%i_pupil_left_%s.jpg',subject_number,condition))

% Right eye
figure
hold on
Dist = shadedErrorBar(x_axis,nanmean(right_pupil_dist_interp,1),nanstd(right_pupil_dist_interp,1),'-b',1);
Targ = shadedErrorBar(x_axis,nanmean(right_pupil_targ_interp,1),nanstd(right_pupil_targ_interp,1),'-r',1);
plot(x_axis,nanmean(right_pupil_dist_base,1),'b')
plot(x_axis,nanmean(right_pupil_targ_base,1),'r')
ylabel('Percentage Change')
xlabel('Time (ms)')
xlim([-200 3000])
title(sprintf('Right Eye Pupil Data Subject %i %s condition',subject_number,condition))
legend([Dist.mainLine, Targ.mainLine], 'Distractors', 'Targets', 'Location', 'SouthWest')
saveas(gcf,sprintf('s%i_pupil_right_%s.jpg',subject_number,condition))
