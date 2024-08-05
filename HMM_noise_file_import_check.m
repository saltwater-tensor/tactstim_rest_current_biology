% 
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
for subjects = [6]%6:length(database_files)-1
    
    sub = fullfile(database,database_files(subjects).name);
    subjects_list{sublist} = sub;
    subjects_name{sublist} = database_files(subjects).name;
    sublist = sublist + 1;
    subject_dir = dir(sub);
    for sbj = 5:length(subject_dir)
        
       if ((~isempty(regexpi(subject_dir(sbj).name,'notch','once')))) &&...
        (~isempty(regexpi(subject_dir(sbj).name,'empty(\w*)','once')))&&...
        (~isempty(regexpi(subject_dir(sbj).name,'@','once')))
        
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


%% 
% Process: Import MEG/EEG: Time
% Compute noise covariance

for sbj = 1:length(original_sFiles)
    
    subj_nm = strsplit(original_sFiles{sbj},'/');
    subj_nm = subj_nm{1};
    sFiles = bst_process('CallProcess', 'process_import_data_time', original_sFiles{sbj}, [], ...
        'subjectname', subj_nm, ...
        'condition',   '', ...
        'timewindow',  [], ...
        'split',       0, ...
        'ignoreshort', 0, ...
        'usectfcomp',  0, ...
        'usessp',      1, ...
        'freq',        [], ...
        'baseline',    []);
    
    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline',       [], ...
    'datatimewindow', [], ...
    'sensortypes',    'MEG', ...
    'target',         1, ...  % Noise covariance      (covariance over data time window)
    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       0, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    3);  % Keep   
    
end


