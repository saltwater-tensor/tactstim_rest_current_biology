function p_perm_vec = temporal_test(lftimes_mat_resp1, lftimes_mat_resp2,p_value_to_test)

    for st = 1:6
        if iscell(lftimes_mat_resp1)
            vec1 = cell2mat(lftimes_mat_resp1(1:end,st));
            vec2 = cell2mat(lftimes_mat_resp2(1:end,st));
        else
            vec1 = (lftimes_mat_resp1(1:end,st));
            vec2 = (lftimes_mat_resp2(1:end,st));
        end
        [p_perm, observeddifference, effectsize] = permutationTest(vec1, vec2, 1000);
        % Sided test
        [p_perm_larger, observeddifference, effectsize] = permutationTest(vec1, vec2, 1000,'sidedness','larger');
        [p_perm_smaller, observeddifference, effectsize] = permutationTest(vec1, vec2, 1000,'sidedness','smaller');
        p_tbl_permtest_larger(st) = p_perm_larger;
        p_tbl_permtest_smaller(st) = p_perm_smaller;
        p_perm_vec(st) = p_perm;
    end
    p_perm_vec = mafdr(p_perm_vec,'BHFDR',true);
    p_perm_vec_disp = p_perm_vec;
    p_perm_vec_disp(p_perm_vec_disp >0.05) = NaN;

    % Correcting for multiple comparisons one sided tests 
    p_permtest_mafdr_larger = mafdr(p_tbl_permtest_larger(:),'BHFDR',true);
        
    p_permtest_mafdr_smaller = mafdr(p_tbl_permtest_smaller(:),'BHFDR',true);

    % Sidedness
    p_disp_permtest_mafdr_larger = p_permtest_mafdr_larger;
    p_disp_permtest_mafdr_larger(p_disp_permtest_mafdr_larger > p_value_to_test) = NaN;
    p_disp_permtest_mafdr_smaller = p_permtest_mafdr_smaller;
    p_disp_permtest_mafdr_smaller(p_disp_permtest_mafdr_smaller > p_value_to_test) = NaN;

    
    figure('Name','Permutation test two sided','NumberTitle','off') 
    heatmap(p_perm_vec_disp,'Colormap',jet)

    figure('Name','Permutation test smaller','NumberTitle','off') 
    heatmap(p_disp_permtest_mafdr_smaller,'Colormap',jet)

    figure('Name','Permutation test larger','NumberTitle','off') 
    heatmap(p_disp_permtest_mafdr_larger,'Colormap',jet)
    


end