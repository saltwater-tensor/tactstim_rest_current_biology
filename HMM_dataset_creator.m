%% HMM dataset creator

source_model = 'PNAI';
trial_type = 'preStim_rest';
database = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\data';
database_anat = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\';

%% Source space projection
use_source_space = 1;
% Project to default anatomy
project_to_default_anat = 1;
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';

%% Atlas parcellation registered in brainstorm
% THESE ROIS are for the MINDBOGGLE ATLAS

ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];

if ~exist('atlas_name_check','var')
    atlas_name_check = input('Enter atlas name for the above ROIs');
end

data_subject= [];
T = [];
BEHAVIOR = [];
subject_total_trials = [];
subject_bad_trials = [];
current_subject = [];
old_subject = [];
iSubject = 0;

%%
for i_sub_sFile= 1:length(original_sFiles)
    
    current_file = original_sFiles{i_sub_sFile};
    current_file = strsplit(current_file,'/');
    current_subject = current_file{1};
    current_file = current_file{2};
    current_file = current_file(5:end);
    date_current_file = strsplit(current_file,'_');
    date_current_file = date_current_file{3};
    stim_file_folder = fullfile(database,current_subject,current_file);
    stim_files = dir(stim_file_folder);
    
%     splitname = strsplit(subjects_name{isubject},'_');
%     joined_name = strjoin(splitname(1:2),' ');
    
    
    Data=[];
    

%     stim_folder_name = fullfile(subjects_list{isubject},trial_folder_name);
%     stim_files = dir(stim_folder_name);
    ImagingKernel = [];
    
    c_sub = convertCharsToStrings(current_subject);
    o_sub = convertCharsToStrings(old_subject);
    if ~strcmpi(c_sub,o_sub)
        iSubject = iSubject + 1;
        combined_across_trials = []; 
        total_trials_count = 0;
        bad_trials_count = 0;
        good_trials_count = 0;
        current_subject_disp = strsplit(current_subject,'_');
        current_subject_disp = strjoin(current_subject_disp(1:2));
        h = waitbar(0,['Processing Subject ' current_subject_disp ]);
        old_subject = current_subject;
    end
    
%% For loading shared imaging kernel once
    
    for i_trial = 5:length(stim_files)
           
        if ((~isempty(regexpi(stim_files(i_trial).name,source_model,'once'))))
            
            Kernel = load(fullfile(stim_file_folder,stim_files(i_trial).name));
            ImagingKernel = Kernel.ImagingKernel;
            surface_file = Kernel.SurfaceFile;
            break
        end
           
    end
    
    
    
%% Data files for individual trials
    for i_trial = 5:length(stim_files)
        sRate=[];
        
       waitbar(i_trial/length(stim_files))
        
        if ((~isempty(regexpi(stim_files(i_trial).name,'data','once')))) &&...
            (isempty(regexpi(stim_files(i_trial).name,'trial(\d*)_02','once'))) &&...
            (~isempty(regexpi(stim_files(i_trial).name,'trial(\d*)','once')))
            
                total_trials_count = total_trials_count + 1;
                good_trials_count = good_trials_count + 1;
                file_name = fullfile(stim_file_folder,stim_files(i_trial).name) ;
                data = load(file_name);
                recordings = data.F(Kernel.GoodChannel(1:end),:);
                data.F = recordings;
                clear recordings
                
                % The corresponding behavioral file should be named as
                % follows
                behavior_name = strsplit(stim_files(i_trial).name,'.');
                behavior_name = behavior_name{1};
                behavior_name = [behavior_name '_02.mat'];
                

    %% Bad trial detection
            events = data.Events;
            bad_trial_flag = 0;
            for ev = 1:length(events)
                event_type = convertCharsToStrings(events(ev).label);
                if (strcmpi(event_type,"BAD"))
                    bad_trial_flag = 1;
                    bad_trials_count = bad_trials_count + 1;
                    good_trials_count = good_trials_count - 1;
                    break
                end
            end

            if bad_trial_flag
                continue
            end
    %% If it is a good trial then find behavior metrics
    
    % At each trial for every subject the behavior struct is instantiated
   Subj_behavior = [];
   if ~isempty(find(ismember({stim_files.name}, {behavior_name})))
       
       behavior_index = find(ismember({stim_files.name}, {behavior_name}));
       behavior_file = load(fullfile(stim_file_folder,stim_files(behavior_index).name));
       behavior_events = behavior_file.Events;
       
       % The event based import for the behavioral file is centered at/ is
       % 0 at the first stimulus pulse which is labelled Stim_final
       % The second stimulus pulse is labelled as Stim_2
       % If Stim_2 and Stim_final both are present i.e. the current trial
       % is not a catch trial then the subtraction will work 
       % Otherwise SOA_ms = -1 indicating a catch trial
       
       %For details please refer to HMM_behavior_extraction script
       % Section: RENAMING AND RELABELLING EVENTS
       try
           SOA_ms = 1000*(behavior_events(find(ismember({behavior_events.label}, {'Stim_2'}))).times);
       catch
           SOA_ms = -1;           
       end
       
       Response = behavior_file.Comment;
       Response = strsplit(Response,' '); 
       Trial_number = Response{2};
       Response = strsplit(Response{1},'_');
       Response = strsplit(Response{2},'R');
       Response = ['R' Response{2}];
       
       % Storing all behavior relevant information for this TRIAL
       Subj_behavior.SOA_ms = SOA_ms;
       Subj_behavior.Response = Response;
       Subj_behavior.Trial_number = Trial_number;
       Subj_behavior.original_file_location = fullfile(stim_file_folder,stim_files(behavior_index).name);
       Subj_behavior.comment = 'Remove _02 from the file location to get the location of pre stim file' ;
   else
       
       bad_trials_count = bad_trials_count + 1;
       good_trials_count = good_trials_count - 1;
       continue
   end
    %% Downsampling 
            Time = data.Time;
            downsample_rate = 250;
            if ~isempty(data)
                sRate = round(abs(1 / (data.Time(2) - data.Time(1)))); %new sampling rate This is the sampling rate that we obtain after BRAINSTORM preprocessing.
            end

            [data.F, data.Time] = process_resample('Compute', data.F,data.Time,downsample_rate,'resample-cascade');

    %% Source space projection
            if use_source_space
             try
                  Data = ImagingKernel * data.F;
                  if strcmpi...
                                (convertCharsToStrings(source_model), ...
                                convertCharsToStrings('MNE'))
                    % Rescale using minmax scaling the MNE data
                    % This applies only to MNE source space data
                    
                    Data_max = max(max(Data));
                    Data_min = min(min(Data));
                    Data = (Data-Data_min)./(Data_max-Data_min);
                           
                            
                  end
                  if project_to_default_anat
                      
                      
                      % Calculate interpolation matrix
                      % Interpolation is calculated here
                      srcSurfFile = surface_file;
                      isInteractive = 0;
                      isStopWarped = [];
                      try
                      [Wmat, sSrcSubj, sDestSubj, srcSurfMat, destSurfMat, isStopWarped] = ...
                            tess_interp_tess2tess( srcSurfFile, destSurfFile, isInteractive, isStopWarped );
                      catch
                          error('BRAINSTORM SHOULD BE RUNNING IN THE BACKGROUND, DO NOT USE CLEAR ALL')
                      end
                      % Apply interpolation to the Data/ImageGridAmp matrix
                      % A local copy of the interpolation function exists in the
                      % code folder
                      disp('DEFAULT PROJECTION')
                      Data = muliplyInterp_local...
                        (Wmat, double(Data), Kernel.nComponents);
                   end
               catch
                   continue
             end 
            end

    %% Atlas based signal extraction
          atlas_name = 'Mindboggle';
          if ~strcmpi(atlas_name,atlas_name_check)
              error('YOU HAVE NOT CHECKED THE ATLAS CHECK ROIs LINE 40 above too')
          end

          if project_to_default_anat
               surface_file = destSurfFile;
               surface_file_data = load([database_anat 'anat' filesep surface_file]);
          end
          brk = 0;  
          for atname = 1:1:length(surface_file_data.Atlas)
              if strcmpi(surface_file_data.Atlas(atname).Name,atlas_name)
                     brk = 1;
                     Scouts = surface_file_data.Atlas(atname).Scouts;
                     PCA_scout = [];
                     for sc = 1:1:length(ROIs)
                        
                        V = Scouts(ROIs(sc)).Vertices;
                        dt = Data(V,:);
                        pca_opts.stack_data='yes';
                        pca_opts.storage='full';
                        pca_opts.precision='double';
                        pca_opts.eig_solver='selective';
                        
                        try
                            [V, Lambda] = icatb_calculate_pca(dt, min(25,size(dt,1)-1),...
                                'type', 'standard', 'whiten', 0, 'verbose', 1, ...
                                'preproc_type', 'remove mean per timepoint', ...
                            'pca_options', pca_opts, 'varToLoad',...
                            [], 'dims', size(dt));
                        catch
                            try 
                            [V, Lambda] = icatb_calculate_pca(dt, 20,...
                                'type', 'standard', 'whiten', 0, 'verbose', 1, ...
                                'preproc_type', 'remove mean per timepoint', ...
                            'pca_options', pca_opts, 'varToLoad',...
                            [], 'dims', size(dt));
                            catch
                                try
                                    [V, Lambda] = icatb_calculate_pca(dt, 4,...
                                        'type', 'standard', 'whiten', 0, 'verbose', 1, ...
                                        'preproc_type', 'remove mean per timepoint', ...
                                    'pca_options', pca_opts, 'varToLoad',...
                                    [], 'dims', size(dt));
                                catch
                                    disp('UNSTABLE SOURCE RECONSTRUCTION..maybe due to low number of trials')
                                end
                            end


                        end

                        V = V';
                        V_last = V(end,:); 
                        PCA_scout = [PCA_scout;V_last];

                     end
              end
             if brk==1
                break
             end
          end
          
            if rank(PCA_scout) == 42
            
                icasig = icatb_remove_mean(PCA_scout)';

                mean_PCA_scout = mean(PCA_scout,2);
                variance_PCA_scout = std(PCA_scout,[],2);

                mean_PCA_scout_1 = repmat(mean_PCA_scout,[1,size(PCA_scout,2)]);
                variance_PCA_scout_1 = repmat(variance_PCA_scout,[1,size(PCA_scout,2)]);

                PCA_scout_mean_sub = PCA_scout - mean_PCA_scout_1;
                PCA_scout_var_norm = PCA_scout_mean_sub ./ variance_PCA_scout_1;


                [L]= symmetric_orthogonalise_custom(PCA_scout_var_norm);

                L = L - repmat(mean(L,2),[1,size(L,2)]);
                L = L./repmat(std(L,[],2),[1,size(L,2)]);


                MEG_signals_per_trial = L';      
                data_subject{1,iSubject}{1,good_trials_count} = MEG_signals_per_trial;
                T{1,iSubject}(1,good_trials_count) = size(MEG_signals_per_trial,1);
                BEHAVIOR{iSubject,2}(1,good_trials_count) = Subj_behavior;
                BEHAVIOR{iSubject,1} = current_subject;
                combined_across_trials = [combined_across_trials;MEG_signals_per_trial];

            else
                bad_trials_count = bad_trials_count + 1;
                good_trials_count = good_trials_count - 1;
                msgbox('SHORT TRIAL DETECTED, was trashed as a bad trial',['Short trial' current_subject])
            end




      
        end    
    end
    
    
    subject_total_trials(iSubject) = total_trials_count;
    subject_bad_trials(iSubject) = bad_trials_count;
    
end


%% Save data
clck = char(datetime);
clck = strrep(clck,':','_');
clck = strrep(clck,' ','_');
clck = strrep(clck,'-','_');


save(['Time_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'T','-v7.3')
save(['Behavior_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'BEHAVIOR','-v7.3')
save(['data_subject_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'data_subject','-v7.3')
save(['Total_trials_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'subject_total_trials','-v7.3')
save(['Bad_trials_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'subject_bad_trials','-v7.3')
save(['sampling_frequency_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'downsample_rate','-v7.3')


%% RANDOM SHUFFLE DATASET

for shuffle_sub = 1:size(data_subject,2)
    
    % Shuffle response1 and response2 trials within a subject so that the HMM is not
    % stuck in some local minima for a specific response type
    s = RandStream('dsfmt19937');
    trial_shuffle = randperm(s,size(BEHAVIOR{shuffle_sub,2},2));
    shuffled_behavior = BEHAVIOR{shuffle_sub,2};
    shuffled_behavior = shuffled_behavior(trial_shuffle);
    BEHAVIOR{shuffle_sub,2} = shuffled_behavior;
    shuffled_data = data_subject{1,shuffle_sub};
    shuffled_data = shuffled_data(trial_shuffle);
    data_subject{1,shuffle_sub} = shuffled_data';
    shuffled_T =  T{1,shuffle_sub};
    shuffled_T = shuffled_T(trial_shuffle);
    T{1,shuffle_sub} = shuffled_T;

    DATA{shuffle_sub,1} = cell2mat(data_subject{1,shuffle_sub});
    
end

save(['shuffled_DATA_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'DATA','-v7.3')
save(['shuffled_Time_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'T','-v7.3')
save(['shuffled_Behavior_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'BEHAVIOR','-v7.3')
save(['shuffled_data_subject_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'data_subject','-v7.3')
save(['shuffled_Total_trials_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'subject_total_trials','-v7.3')
save(['shuffled_Bad_trials_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'subject_bad_trials','-v7.3')
save(['shuffled_sampling_frequency_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'downsample_rate','-v7.3')
