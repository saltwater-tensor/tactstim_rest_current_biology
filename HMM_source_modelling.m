% 
clear
global GlobalData
% This script detects the multiple stimulus instances for the Tactstim
% measurement. The stimulus event is called '16' 

database = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\data';
original_sFiles = {};
sublist = {};
subjects_name = {};
database_files = dir(database);
sFile_numbers = 1;
sublist = 1;

%%
for subjects = 6:35 %6:length(database_files)-1
    
    sub = fullfile(database,database_files(subjects).name);
    subjects_list{sublist} = sub;
    subjects_name{sublist} = database_files(subjects).name;
    sublist = sublist + 1;
    subject_dir = dir(sub);
    for sbj = 5:length(subject_dir)
        
       if ((~isempty(regexpi(subject_dir(sbj).name,'notch','once')))) &&...
        ((~isempty(regexpi(subject_dir(sbj).name,'taskblock','once'))) || ...
        (~isempty(regexpi(subject_dir(sbj).name,'task_block','once'))) )&&...
        (~isempty(regexpi(subject_dir(sbj).name,'@','once')))&&...
        (isempty(regexpi(subject_dir(sbj).name,'pre','once')))&&...
        (isempty(regexpi(subject_dir(sbj).name,'post','once')))
        
            raw_data_dir = dir(fullfile(database,database_files(subjects).name,subject_dir(sbj).name));
            for rdd = 3:length(raw_data_dir)
                
                if ((~isempty(regexpi(raw_data_dir(rdd).name,'notch','once')))) &&...
                    (~isempty(regexpi(raw_data_dir(rdd).name,'data','once')))
                    
                    file_name = fullfile(database_files(subjects).name,subject_dir(sbj).name,...
                        raw_data_dir(rdd).name);
                    file_name_brainstorm = strrep(file_name,'\','/'); 
                    original_sFiles{sFile_numbers} = file_name_brainstorm;
                    sFile_numbers = sFile_numbers + 1;
                    break
                end
            end
            
       end
        
        
    end
end

original_sFiles_raw = original_sFiles;

%% Cleaning up repetitive data blocks saved in brainstorm

% Please edit day 1 of subject 1 MANUALLY as it does not follow naming conventions
previous_block_number = [];
previous_subject = [];
previous_merge_level = '(0)';
previous_file_to_keep = [];
files_to_keep.number = [];
files_to_keep.name = {};
files_to_keep.name_full = {};
previous_block_name = {};

for osfl = 7:length(original_sFiles)
    
    current_file_to_keep = osfl;
    current_block_name = original_sFiles{osfl};
    current_block_name = strsplit(current_block_name,'/');
    current_block_name = current_block_name{end};
    expression1 = '\w*taskblock\w*';
    matchStr1 = regexp(current_block_name,expression1,'match');
    matchStr2 = regexp(current_block_name,expression1,'split');
    matchStr1_split = strsplit(matchStr1{1},'_');
    matchStr2_split = strsplit(matchStr2{end},'_');
    
    
    if length(matchStr1_split) >= 6
        current_block_number = strjoin(matchStr1_split(1,[4,5,6]),'_');
        current_subject = matchStr1_split(3);
    
%     elseif length(matchStr1_split) == 7
%         current_block_number = strjoin(matchStr1_split(1,[4,5,6]),'_');
%         current_subject = matchStr1_split(3);
        
    end
    
    if length(matchStr2_split) == 3
        current_merge_level = matchStr2_split{2};
    elseif length(matchStr2_split) == 2
        current_merge_level = matchStr2_split{1};
        if  contains(current_merge_level,'-')
            current_merge_level = '(0)';
        end
    elseif length(matchStr2_split) == 1
        current_merge_level = matchStr2_split{1};
        if  contains(current_merge_level,'.mat')
            current_merge_level = '(0)';
        end
    end
    clc
    current_subject
    current_block_number
    current_merge_level
    
    if ~isempty(previous_subject)
    str_check1_previous = [previous_subject previous_block_number previous_merge_level];
    str_check1_current = [current_subject current_block_number current_merge_level];
    if strcmp(convertCharsToStrings(strjoin(str_check1_previous(1:3),'_')),...
            convertCharsToStrings(strjoin(str_check1_current(1:3),'_')))
        % This is problematic as you cannot have two exactly similar files
        files_to_keep.number = [files_to_keep.number;previous_file_to_keep];
        str_file = [previous_subject previous_block_number previous_merge_level];
        files_to_keep.name{length(files_to_keep.number),1} =  strjoin(str_file(1:3),'_');
        files_to_keep.name_full{length(files_to_keep.number),1} = previous_block_name;
        
        files_to_keep.number = [files_to_keep.number;current_file_to_keep];
        files_to_keep.name{length(files_to_keep.number),1} = strjoin(str_check1_current(1:3),'_');
        files_to_keep.name_full{length(files_to_keep.number),1} = current_block_name;
        continue
    end 
    end
    if strcmpi(convertCharsToStrings(current_subject),previous_subject)
        
        if strcmpi(convertCharsToStrings(current_block_number),previous_block_number)
            
                if ~strcmpi(convertCharsToStrings(current_merge_level),previous_merge_level)
                    
                    if str2num(current_merge_level(2)) > str2num(previous_merge_level(2))
                        previous_merge_level = current_merge_level;
                        current_merge_level = [];
                        previous_file_to_keep = current_file_to_keep;
                        current_file_to_keep = [];
                        previous_block_name = current_block_name;
                        current_block_name = {};
                        
                        if osfl == length(original_sFiles)
                            files_to_keep.number = [files_to_keep.number;previous_file_to_keep];
                            str_file = [previous_subject previous_block_number previous_merge_level];
                            files_to_keep.name{length(files_to_keep.number),1} =  strjoin(str_file(1:3),'_');
                            files_to_keep.name_full{length(files_to_keep.number),1} = previous_block_name;
                        end
                    end
                end           
        else
            
            % A change in block means that the highest merge level must
            % have been found for the current block and the file should be
            % saved
            files_to_keep.number = [files_to_keep.number;previous_file_to_keep];
            str_file = [previous_subject previous_block_number previous_merge_level];
            files_to_keep.name{length(files_to_keep.number),1} =  strjoin(str_file(1:3),'_');
            files_to_keep.name_full{length(files_to_keep.number),1} = previous_block_name;
            
            previous_block_number = current_block_number;
            current_block_number = [];
            previous_merge_level = current_merge_level;
            current_merge_level = [];
            previous_file_to_keep = current_file_to_keep;
            current_file_to_keep = [];
            previous_block_name = current_block_name;
            current_block_name = {};
            
        end
        
    else
        if isempty(previous_subject)
            previous_subject = current_subject;
            current_subject = [];
            previous_block_number = current_block_number;
            current_block_number = [];
            previous_merge_level = current_merge_level;
            current_merge_level = [];
            previous_file_to_keep = current_file_to_keep;
            current_file_to_keep = [];
        elseif ~isempty(previous_subject)
            
            % A change in subject means that the last block which was
            % processed it's file number should be saved
            
            files_to_keep.number = [files_to_keep.number;previous_file_to_keep];
            str_file = [previous_subject previous_block_number previous_merge_level];
            files_to_keep.name{length(files_to_keep.number),1} =  strjoin(str_file(1:3),'_');
            files_to_keep.name_full{length(files_to_keep.number),1} = previous_block_name;
            
            previous_subject = current_subject;
            current_subject = [];
            previous_block_number = current_block_number;
            current_block_number = [];
            previous_merge_level = current_merge_level;
            current_merge_level = [];
            previous_file_to_keep = current_file_to_keep;
            current_file_to_keep = [];
            previous_block_name = current_block_name;
            current_block_name = {};
            
        end
    end
end

to_delete = [];
for ftk = 2:length(files_to_keep.number)
    ftk_n_1 = files_to_keep.number(ftk-1);
    ftk_n_2 = files_to_keep.number(ftk);
    
    if ftk_n_2 < ftk_n_1
        to_delete = [to_delete;ftk];
        
    end
end
files_to_keep.number(to_delete) = [];
files_to_keep.name(to_delete) = [];
files_to_keep.name_full(to_delete) = [];

% Remove
O = original_sFiles(1,files_to_keep.number);

% Final Safety check
safety_signal = NaN(length(files_to_keep.number),1);
for ftk = 1:length(files_to_keep.number) 
    
    O_file = O{ftk};
    O_file = strsplit(O_file,'/');
    O_file = O_file{3};
    
    chck_file = files_to_keep.name_full{ftk};
    
    if sum(double(chck_file == O_file)) == length(O_file)
        safety_signal(ftk) = 1;
    end
    
end
[safety_nan] = isnan(safety_signal);
original_sFiles_1 = original_sFiles_raw(1,[3,6]);
original_sFiles = original_sFiles_raw(1,files_to_keep.number);
original_sFiles = [original_sFiles_1, original_sFiles];

%% Import event based pre stim and response epochs for neural and behavioral data
% respectively

% Importing in the sequence coded below and using the eventname argument as used
% ensures that the neural and behavioral data epochs are in SYNC!!! (CRITICAL)


% IMPORT PRE STIM EPOCH FOR BOTH RESPONSE 1 AND RESPONSE 2 TRIALS
bst_report('Start', original_sFiles);
% Process: Import MEG/EEG: Events
sFiles = bst_process('CallProcess', 'process_import_data_event', original_sFiles, [], ...
    'subjectname', [], ...
    'condition',   '', ...
    'eventname',   'Stim_StimResponse1,Stim_StimResponse2', ...
    'timewindow',  [], ...
    'epochtime',   [-8, -0.1], ... %time in seconds
    'createcond',  0, ...
    'ignoreshort', 0, ... %will be handled in HMM_dataset_creator script
    'usectfcomp',  0, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% IMPORT RESPONSE 1 EPOCHS
% Process: Import MEG/EEG: Events
sFiles1 = bst_process('CallProcess', 'process_import_data_event', original_sFiles, [], ...
    'subjectname', [], ...
    'condition',   [], ...
    'eventname',   'Stim_StimResponse1', ...
    'timewindow',  [], ...
    'epochtime',   [-0.5, 0.5], ... %time in seconds
    'createcond',  0, ...
    'ignoreshort', 0, ...
    'usectfcomp',  0, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% IMPORT RESPONSE 2 EPOCHS
% % Process: Import MEG/EEG: Events
sFiles2 = bst_process('CallProcess', 'process_import_data_event', original_sFiles, [], ...
    'subjectname', [], ...
    'condition',   [], ...
    'eventname',   'Stim_StimResponse2', ...
    'timewindow',  [], ...
    'epochtime',   [-0.5, 0.5], ... %time in seconds
    'createcond',  0, ...
    'ignoreshort', 0, ...
    'usectfcomp',  0, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile,'Q:\tactstim_resting\code');

%% Source modeling

% Calculate head model for all the subjects listed
missing_noise_file = {};
missing_data_file = {};
missing_noise = 1;
missing_data = 1;
for headmodel_sub_list = 1:length(original_sFiles)
    
    current_file = original_sFiles{headmodel_sub_list};
    current_file = strsplit(current_file,'/');
    current_subject = current_file{1};
    current_file = current_file{2};
    current_file = current_file(5:end);
    date_current_file = strsplit(current_file,'_');
    date_current_file = date_current_file{3};
    stim_file_folder = fullfile(database,current_subject,current_file);
    stim_files = dir(stim_file_folder);
    
    headmodel_files = [];
    hmfiles = [];
    hmfiles_counter = 1;
    for hmfiles = 5:length(stim_files)
        
        % IMPORTANT PLEASE READ THIS!!!! DONT IGNORE!!
        % Although headmodel is independent of the data files used for
        % source modelling we will restrict ourselves to files that contain
        % the PRE STIM NEURAL DATA. This will be helpful to track the
        % variable called as sFiles which is used by Brainstorm to
        % subsequently calculate Data covariance and source model
        
        if ((~isempty(regexpi(stim_files(hmfiles).name,'data','once')))) &&...
            (isempty(regexpi(stim_files(hmfiles).name,'trial(\d*)_02','once'))) &&...
            (~isempty(regexpi(stim_files(hmfiles).name,'trial(\d*)','once')))
            
            file_name = fullfile(current_subject,current_file,stim_files(hmfiles).name) ;
            file_name_brainstorm = strrep(file_name,'\','/'); 
            headmodel_files{hmfiles_counter} = file_name_brainstorm;
            hmfiles_counter = hmfiles_counter + 1;

        end

    end
    
%% Process: Compute head model
    sFiles = bst_process('CallProcess', 'process_headmodel', headmodel_files, [], ...
    'Comment',     '', ...
    'sourcespace', 1, ...  % Cortex surface
    'meg',         3, ...  % Overlapping spheres
    'eeg',         3, ...  % OpenMEEG BEM
    'ecog',        2, ...  % OpenMEEG BEM
    'seeg',        2, ...  % OpenMEEG BEM
    'openmeeg',    struct(...
         'BemFiles',     {{}}, ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemCond',      [1, 0.0125, 1], ...
         'BemSelect',    [1, 1, 1], ...
         'isAdjoint',    0, ...
         'isAdaptative', 1, ...
         'isSplit',      0, ...
         'SplitLength',  4000), ...
    'channelfile', '');
%% Check for noise covariance file  and copy if not present

% This section will produce error if your empty room file has any ot
filenames = {GlobalData.DataBase.ProtocolStudies.Study.Name};
OutputFiles = [];
if isempty(find(ismember({stim_files.name}, {'noisecov_full.mat'})))
        
        subject_folder = fullfile(database,current_subject);
        subject_folder_dir = dir(subject_folder);
%         index = find(ismember({subject_folder_dir.name}, {['emptyroom_' date_current_file '_notch']}));
        index = -1;
        for i_s_f = 1:length(subject_folder_dir)
            
            if ((~isempty(regexpi(subject_folder_dir(i_s_f).name,['empty(\w*)_' date_current_file '_notch'],'once')))) &&...
                    ((isempty(regexpi(subject_folder_dir(i_s_f).name,'@','once')))) 
                index = i_s_f;
                break
            end
        end
        
        if index == -1
            for i_s_f = 1:length(subject_folder_dir)
            
                if ((~isempty(regexpi(subject_folder_dir(i_s_f).name,['empty(\w*)_'  'notch'],'once')))) &&...
                        ((isempty(regexpi(subject_folder_dir(i_s_f).name,'@','once')))) 
                    index = i_s_f;
                    break
                end
                
            end
        end
        
        
        
        
        emptyroom_dir = dir(fullfile(database,current_subject,subject_folder_dir(index).name));
        index_noise_file = (find(ismember({emptyroom_dir.name}, {'noisecov_full.mat'})));
        
        if ~isempty(index_noise_file)
            
            emptyroom_name = subject_folder_dir(index).name; 
            Target = (find(ismember(filenames, {current_file})));
            iSrcStudy = (find(ismember(filenames, {emptyroom_name})));
            try
                OutputFiles = db_set_noisecov(iSrcStudy, Target, 0, 1);
            catch
                continue
            end
            
        else
            emptyroom_indices = find(startsWith({subject_folder_dir.name},{'emptyroom'}));
            for emp_ind = 1:length(emptyroom_indices)
                emptyroom_dir = dir(fullfile(database,current_subject,subject_folder_dir(emptyroom_indices(emp_ind)).name));
                index_noise_file_new = (find(ismember({emptyroom_dir.name}, {'noisecov_full.mat'})));
                
                if ~isempty(index_noise_file_new)
                    
                    emptyroom_name = subject_folder_dir(emptyroom_indices(emp_ind)).name;
   
                    Target = (find(ismember(filenames, {current_file})));
                    iSrcStudy = (find(ismember(filenames, {emptyroom_name})));
                    
                    if ~isempty(Target) && ~isempty(iSrcStudy)
                        OutputFiles = db_set_noisecov(iSrcStudy, Target, 0, 1);
                    else
                       missing_data_file{missing_data} = current_file;
                       missing_data = missing_data + 1;
                    end


                end
                
            end
            
        end
        
        if isempty(OutputFiles)
            missing_noise_file{missing_noise} = current_file;
            missing_noise = missing_noise + 1;
%             error('Cannot continue with this subject noise covariance is missing, please calculate')
        end
end



%% Process: Compute covariance (noise or data)

% CRITICAL: To ensure that the correct sFiles are being used 
% Check the field sFiles.FileName
% The last part of all filenames should be of the type > e.g.:
% data_Stim_StimResponse1_trial001.mat 
% IT SHOULD NOT BE: data_Stim_StimResponse1_trial001_002.mat
% _002 are BEHAVIORAL FILES !

    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline',       [], ...
    'datatimewindow', [], ...
    'sensortypes',    'MEG', ...
    'target',         2, ...  % Data covariance      (covariance over data time window)
    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       0, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    3);  % Keep    

% Process: Compute sources [2018]
    sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
    'output',  1, ...  % Kernel only: shared
    'inverse', struct(...
         'Comment',        'MN: MEG ALL', ...
         'InverseMethod',  'minnorm', ...
         'InverseMeasure', 'amplitude', ...
         'SourceOrient',   {{'fixed'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      0.5, ...
         'WeightLimit',    10, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       0.1, ...
         'SnrMethod',      'fixed', ...
         'SnrRms',         1e-06, ...
         'SnrFixed',       3, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));
     
     
% Process: Compute sources [2018]
    sFiles = bst_process('CallProcess', 'process_inverse_2018',  sFiles, [], ...
    'output',  1, ...  % Kernel only: shared
    'inverse', struct(...
         'Comment',        'PNAI: MEG ALL', ...
         'InverseMethod',  'lcmv', ...
         'InverseMeasure', 'nai', ...
         'SourceOrient',   {{'fixed'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      0.5, ...
         'WeightLimit',    10, ...
         'NoiseMethod',    'median', ...
         'NoiseReg',       0.1, ...
         'SnrMethod',      'rms', ...
         'SnrRms',         1e-06, ...
         'SnrFixed',       3, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));     


end

%%
% HMM_dataset_creator