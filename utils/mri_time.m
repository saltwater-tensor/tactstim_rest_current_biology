% Read mri timing

mri_path = 'Q:\tactstim_resting\MRI\';
subjs = dir (mri_path);
id_count = 1;
dates = [];
subj_id = [];
for s = 7:38

    if subjs(s).isdir

        slices = fullfile(subjs(s).folder, subjs(s).name);
        slices = dir([slices '\slices']);
        if ~isempty(slices)
            slice_name = slices(6).name;
            slice_name = strsplit(slice_name,'.');
            mri_t = slice_name{end};
            mri_t = mri_t(1:8);
            measurement_date = datestr(datetime(mri_t,"Format","uuuuMMdd"));
            dates{id_count,1} = measurement_date;
            subj_id{id_count,1} = ['TR_' subjs(s).name];
            id_count = id_count + 1;
        end
    end
end
t = table(subj_id,dates, 'VariableNames',{'Name', 'date'});
writetable(t,[mri_path 'final_list.csv'])

