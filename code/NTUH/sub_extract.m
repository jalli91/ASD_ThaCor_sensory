function beta_corrected_p_matrix = sub_extract(ROI_index1, ROI_index2, ROI_name1, ROI_name2, test_matrix)

    hm_matrix_mean = [];
    hm_matrix_int = [];
    hm_matrix_p = [];

    for i = 1:length(ROI_index1)
        for j = 1:length(ROI_index2)
            FC_test_in = test_matrix{ROI_index1(i), ROI_index2(j)};
            beta_p_matrix{i, j} = FC_test_in;
            if ~isnan(FC_test_in)
                hm_matrix_mean(i, j) = FC_test_in(1);
            end    
        end
    end 

    p_values = []; % Store all p-values for FDR correction
    p_indices = []; % Store indices of ROI pairs for FDR correction
    % get adjusted p value matrix to control false discovery rate, q = 0.05
    % Step 1: collect p-values from h_matrix_ex
    for i = 1:size(beta_p_matrix, 1)
        for j = 1:(size(beta_p_matrix, 2)) 
            if ~isnan(beta_p_matrix{i, j})
                p_values = [p_values; beta_p_matrix{i, j}(3)]; % Collect p-value for FDR
                p_indices = [p_indices; i, j]; % Track ROI indices
            end
        end
    end
    % Step 2: Apply FDR correction using Benjamini-Hochberg procedure
    q = 0.05; % Set your FDR level
    [adj_p_values, ~, ~, ~] = fdr_BH(p_values, q); % Apply FDR
    
    % Step 3: Update significance matrix based on FDR-adjusted p-values
    for k = 1:length(p_values)
        i = p_indices(k, 1);
        j = p_indices(k, 2);
        adj_p = adj_p_values(k);
        hm_matrix_p(i, j) = adj_p; % Update with FDR-corrected p-value
        beta_p_matrix{i, j}(4) = adj_p;
    end

    for i = 1:size(beta_p_matrix, 1)
        beta_corrected_p_matrix{i+1, 1} = ROI_name1(i);
        for j = 1:(size(beta_p_matrix, 2)) 
            beta_corrected_p_matrix{1, 4*j-2} = ROI_name2(j)+'_b';
            beta_corrected_p_matrix{1, 4*j-1} = ROI_name2(j)+'_int';
            beta_corrected_p_matrix{1, 4*j} = ROI_name2(j)+'_p';
            beta_corrected_p_matrix{1, 4*j+1} = ROI_name2(j)+'_corrected_p';
            beta_corrected_p_matrix{i+1, 4*j-2} = beta_p_matrix{i, j}(1);
            beta_corrected_p_matrix{i+1, 4*j-1} = beta_p_matrix{i, j}(2);
            beta_corrected_p_matrix{i+1, 4*j} = beta_p_matrix{i, j}(3);
            beta_corrected_p_matrix{i+1, 4*j+1} = beta_p_matrix{i, j}(4);
        end
    end 
    
    % Plot heatmap
    % Create a custom colormap: negative values are dark blue, positive values are dark red
    % Near-zero values will be lighter in color.
    nColors = 100;  % Number of colors in the colormap
    blueGradient = [linspace(0.5, 1, nColors/2)', linspace(0.5, 1, nColors/2)', linspace(1, 1, nColors/2)'];  % Light to dark blue
    redGradient  = [linspace(1, 1, nColors/2)', linspace(0.5, 1, nColors/2)', linspace(0.5, 1, nColors/2)'];  % Light to dark red
    customColormap = [blueGradient; flipud(redGradient)];  % Combine blue and red gradients

    hm = figure();
    imagesc(hm_matrix_mean);
    colorbar;
    colormap(customColormap);
    caxis([-0.5, 0.5]);  % Set color axis to ensure 0 is white, -1.5 is blue, and 1.5 is red
    axis equal tight;  % Make grid cells square
    xticks(1:length(ROI_index2));
    xticklabels(ROI_name2);
    yticks(1:length(ROI_index1));
    yticklabels(ROI_name1);
    add_significance_markers(hm_matrix_p);

end

function ROI_sig_index = add_significance_markers(p_value_matrix)
    % Helper function to add significance markers based on p-values
    ROI_sig_index = [];
    for row = 1:size(p_value_matrix, 1)
        for col = 1:size(p_value_matrix, 2)
            if any(p_value_matrix(row, col))
                p_value = p_value_matrix(row, col);
            % Add significance markers based on p-value thresholds
                if p_value <= 0.05
                    ROI_sig_index = [ROI_sig_index; [row col]];
                    markerText = get_significance_marker(p_value);
                    text(col, row, markerText, 'Color', 'black', 'FontSize', 25, 'HorizontalAlignment', 'center');
                end
            end    
        end
    end
end

function markerText = get_significance_marker(p_value)
    % Returns appropriate marker text based on p-value thresholds
    if p_value > 0.01
        markerText = '*';
    elseif p_value > 0.001
        markerText = '**';
    elseif p_value > 0.0001
        markerText = '***';
    else
        markerText = '****';
    end
end