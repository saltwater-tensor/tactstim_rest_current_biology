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
        ((~isempty(regexpi(subject_dir(sbj).name,'taskblock','once'))) || ...
        (~isempty(regexpi(subject_dir(sbj).name,'task_block','once'))) )&&...
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

%% RENAMING AND RELABELLING EVENTS 
% CLASSIFY AND TAG AS RESPONSE 1 VS RESPONSE 2 EVENTS
% DOES NOT DIFFERENTIATE BETWEEN CATCH TRIALS AND ACTUAL TWO STIM TRIALS
% DOES NOT DIFFERENTIATE BASED ON SOA
% Start a new report
bst_report('Start', original_sFiles);

% Process: Detect multiple responses
sFiles = bst_process('CallProcess', 'process_evt_multiresp', original_sFiles, [], ...
    'responses', '16', ...
    'dt',        0.5, ...
    'action',    1, ...  % Keep only the first event
    'rename',    1);

% Process: Rename events
sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
    'src',   '16,Multiple', ...
    'dest',  'Stim,Stim_2');

% Process: Rename events
sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
    'src',   '4,8,512,1024', ...
    'dest',  'First,Second,Index,Middle');

% RESPONSE 1
% Process: Combine stim/response
sFiles = bst_process('CallProcess', 'process_evt_combine', sFiles, [], ...
    'combine', ['First_FirstIndex , ignore , First, Index' 10 'Second_SecondMiddle , ignore , Second, Middle' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 ''], ...
    'dt',      2);

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', 'First_FirstIndex,Second_SecondMiddle', ...
    'newname',  'Response1');

% RESPONSE 2
% Process: Combine stim/response
sFiles = bst_process('CallProcess', 'process_evt_combine', sFiles, [], ...
    'combine', ['First_FirstMiddle , ignore , First , Middle' 10 'Second_SecondIndex , ignore , Second , Index' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 ''], ...
    'dt',      2);

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', 'First_FirstMiddle,Second_SecondIndex', ...
    'newname',  'Response2');

% Process: Combine stim/response
sFiles = bst_process('CallProcess', 'process_evt_combine', sFiles, [], ...
    'combine', ['Stim_StimResponse1 , ignore , Stim , Response1' 10 'Stim_StimResponse2 , ignore , Stim , Response2' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 ''], ...
    'dt',      2);

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', 'Stim_StimResponse1,Stim_StimResponse2', ...
    'newname',  'Stim_final');

% Process: Combine stim/response
sFiles = bst_process('CallProcess', 'process_evt_combine', sFiles, [], ...
    'combine', ['Stim_StimResponse1 , ignore , Stim , Response1' 10 'Stim_StimResponse2 , ignore , Stim , Response2' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 '' 10 ''], ...
    'dt',      2);

ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile,'Q:\tactstim_resting\code');
    
%% 

