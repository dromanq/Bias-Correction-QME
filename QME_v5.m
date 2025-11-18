%This file applies bias-correction using the Quantile Method for Extremes
%(QME). V5 applies the QME to 3-hourly data. The data is preprocessed on a
%monthly basis. Bias-correction is applied on a monthly basis as well.

function Model_COSMO_complete = QME_v5(Observed, Model_COSMO, Model_COSMO_complete, Date, Date_comp, stations)
%% SECTION 1: PREPROCESSING RAINFALL DATA
disp("Pre-processing data")
tic

% Adjust small observed values
Observed(Observed == 0.1) = 0.2;

% Define seasons and years
seasons = { [9,10,11], [12,1,2], [3,4,5], [6,7,8] }; % SON, DJF, MAM, JJA
years = 2006:2023;
numStations = size(Observed, 2);
numSeasons = numel(seasons);
minH_range = 0:0.000001:0.2; % Possible minH values

% Initialize storage matrices
optMinH = zeros(12, numStations); 
min_errorRD_month = inf(12, numStations); 

% Precompute time-based masks
monthIndices = month(Date);
yearIndices = year(Date);

% --- Compute Rainy Days Before Correction ---
numRD_Osea1 = arrayfun(@(s) sum(ismember(monthIndices, seasons{s}) & Observed > 0), 1:numSeasons, 'UniformOutput', false);
numRD_Msea1 = arrayfun(@(s) sum(ismember(monthIndices, seasons{s}) & Model_COSMO > 0), 1:numSeasons, 'UniformOutput', false);
errorRD_seasonal1 = (cell2mat(numRD_Osea1) - cell2mat(numRD_Msea1)) ./ cell2mat(numRD_Osea1);

numRD_Oyea1 = arrayfun(@(y) sum(yearIndices == years(y) & Observed > 0), 1:numel(years), 'UniformOutput', false);
numRD_Myea1 = arrayfun(@(y) sum(yearIndices == years(y) & Model_COSMO > 0), 1:numel(years), 'UniformOutput', false);
errorRD_yea1 = cell2mat(numRD_Oyea1) - cell2mat(numRD_Myea1);

numRD_Omo1 = arrayfun(@(m) sum(monthIndices == m & Observed > 0), 1:12, 'UniformOutput', false);
numRD_Mmo1 = arrayfun(@(m) sum(monthIndices == m & Model_COSMO > 0), 1:12, 'UniformOutput', false);
errorRD_mo1 = cell2mat(numRD_Omo1) - cell2mat(numRD_Mmo1);

% --- Find Optimal minH Per Month ---

calculated = true; % Define manually if you want to compute the optimalMinH

if ~exist('optMinH.mat') || ~calculated
    optMinH = computeOptimalMinH(Observed, Model_COSMO, monthIndices, minH_range);
else
    load('optMinH.mat')
end
toc

% --- Apply Optimal minH ---
for i = 1:numStations
    for m = 1:12
        minH_m = optMinH(m, i);
        
        % Ensure monthMask is based on Date_comp
        monthMask = (month(Date_comp) == m); 
        
        % Ensure Model_COSMO_complete has the correct size
        Model_COSMO_complete = Model_COSMO_complete(1:length(Date_comp), :);

        % Apply threshold correction safely
        idx_Mcomp = find(monthMask & Model_COSMO_complete(:, i) < minH_m);
        Model_COSMO_complete(idx_Mcomp, i) = 0;
    end
end

Model_COSMO = Model_COSMO_complete(1:length(Model_COSMO),:);

% --- Compute Rainy Days After Correction ---
numRD_Osea2 = arrayfun(@(s) sum(ismember(monthIndices, seasons{s}) & Observed > 0), 1:numSeasons, 'UniformOutput', false);
numRD_Msea2 = arrayfun(@(s) sum(ismember(monthIndices, seasons{s}) & Model_COSMO > 0), 1:numSeasons, 'UniformOutput', false);
errorRD_seasonal2 = (cell2mat(numRD_Osea2) - cell2mat(numRD_Msea2));

numRD_Oyea2 = arrayfun(@(y) sum(yearIndices == years(y) & Observed > 0), 1:numel(years), 'UniformOutput', false);
numRD_Myea2 = arrayfun(@(y) sum(yearIndices == years(y) & Model_COSMO > 0), 1:numel(years), 'UniformOutput', false);
errorRD_yea2 = cell2mat(numRD_Oyea2) - cell2mat(numRD_Myea2);

numRD_Omo2 = arrayfun(@(m) sum(monthIndices == m & Observed > 0), 1:12, 'UniformOutput', false);
numRD_Mmo2 = arrayfun(@(m) sum(monthIndices == m & Model_COSMO > 0), 1:12, 'UniformOutput', false);
errorRD_mo2 = cell2mat(numRD_Omo2) - cell2mat(numRD_Mmo2);

%% SECTION 2: CALCULATE DATA DISTRIBUTIONS TO SCALE DATA AND APPLY LIMITS
disp("Calculating data distributions")
tic

% Parameters
min_hRF = 0;
bin_W = 0.1;
xtr = 5;
cal_smth = 0.001; % Smooth factor applied after calibration factors

var_name = 'pr'; 
cal_trh = 'mnth'; % If calibration is monthly, use 'mnth'. If calibration is seasonal, use 'sea'.

tmp = strcmp(cal_trh, 'sea') * 4 + strcmp(cal_trh, 'mnth') * 12;

% Initialize storage structures
dist_obs_ssn = struct();
dist_model_ssn = struct();

numStations = length(stations);
numSeasons = tmp;

MaxAbs = zeros(1, numStations);
reso = zeros(1, numStations);
sf = zeros(1, numStations);

% Precompute time-based masks
monthIndices = month(Date);

% Loop through each station
for i = 1:numStations
    
    % Extract observed and modeled data
    OBS = Observed(:, i);
    MOD = Model_COSMO(:, i);

    % Define upper limit
    MaxAbs(i) = max([OBS; MOD]) + 1;
    
    % Define resolution for histogram bins
    reso(i) = round(MaxAbs(i) / bin_W, 0);
    edges = (0:reso(i) + 1)';

    % Initialize histograms
    dist_obs_ssn_i = zeros(reso(i) + 1, numSeasons);
    dist_model_ssn_i = zeros(reso(i) + 1, numSeasons);

    % Limit data to ensure it's within valid range
    OBS = max(OBS, 0);
    OBS = min(OBS, MaxAbs(i));
    MOD = max(MOD, 0);
    MOD = min(MOD, MaxAbs(i));

    % Compute scale factor
    sf(i) = s_f(MaxAbs(i), reso(i), var_name);
    
    % Scale and round data
    OBS = round(scale_data(OBS, var_name, sf(i)));
    MOD = round(scale_data(MOD, var_name, sf(i)));

    % Compute histograms by season
    for j = 1:numSeasons
        
        if strcmp(cal_trh, 'mnth')
            switch j
                case 1, seasonMask = (monthIndices == 1);
                case 2, seasonMask = (monthIndices == 2);
                case 3, seasonMask = (monthIndices == 3);
                case 4, seasonMask = (monthIndices == 4);
                case 5, seasonMask = (monthIndices == 5);
                case 6, seasonMask = (monthIndices == 6);
                case 7, seasonMask = (monthIndices == 7);
                case 8, seasonMask = (monthIndices == 8);
                case 9, seasonMask = (monthIndices == 9);
                case 10, seasonMask = (monthIndices == 10);
                case 11, seasonMask = (monthIndices == 11);
                case 12, seasonMask = (monthIndices == 12);
            end
        else
            switch j
                case 1, seasonMask = (monthIndices == 9 | monthIndices == 10 | monthIndices == 11); sea = 'Aut';
                case 2, seasonMask = (monthIndices == 12 | monthIndices == 1 | monthIndices == 2); sea = 'Win';
                case 3, seasonMask = (monthIndices == 3 | monthIndices == 4 | monthIndices == 5); sea = 'Spr';
                case 4, seasonMask = (monthIndices == 6 | monthIndices == 7 | monthIndices == 8); sea = 'Sum';
            end
        end       
        
        OBS_sea = OBS(seasonMask);
        MOD_sea = MOD(seasonMask);

        dist_obs_ssn_i(:, j) = histcounts(OBS_sea(OBS_sea >= min_hRF), edges)';
        dist_model_ssn_i(:, j) = histcounts(MOD_sea(MOD_sea >= min_hRF), edges)';
    end

    % Store histograms in struct
    dist_obs_ssn.(stations{i}) = dist_obs_ssn_i;
    dist_model_ssn.(stations{i}) = dist_model_ssn_i;
end

toc
%% SECTION 3: CALCULATE BIAS-CORRECTION FACTORS
disp('Calculating calibration factors');
tic

% Unscaled resolution array for each station
unscaled_reso = struct();
for i = 1:numStations
    unscaled_reso.(stations{i}) = unscale_data((0:reso(i))', var_name, sf(i));
end

xtr = double(xtr); % Ensure it's a double type
bias_arr_ssn = struct();

% Loop through each station
for i = 1:numStations
    
    % Retrieve station-specific histograms
    dist_obs_ssn_i = dist_obs_ssn.(stations{i});
    dist_model_loc_i = dist_model_ssn.(stations{i});
    unscaled_reso_i = unscaled_reso.(stations{i});
    
    % Preallocate bias array
    bias_arr_ssn_i = zeros(reso(i) + 1, numSeasons);
    
    % Loop through each season
    for ssn = 1:numSeasons
        dist_obs_loc = dist_obs_ssn_i(:, ssn);
        dist_model_loc = dist_model_loc_i(:, ssn);

        % Normalize histogram counts
        tots_obs = sum(dist_obs_loc);
        tots_model = sum(dist_model_loc);
        
        if tots_obs ~= tots_model
            scale_factor = tots_model / tots_obs;
            dist_obs_loc = dist_obs_loc * scale_factor;
            if round(sum(dist_obs_loc)) ~= round(sum(dist_model_loc))
                error('Sample size mismatch after scaling');
            end
        end

        % Identify nonzero quantile positions
        vals_obs = find(dist_obs_loc > 0);
        vals_model = find(dist_model_loc > 0);
        if isempty(vals_obs) || isempty(vals_model)
            error('Empty histogram detected');
        end

        % Initialize bias correction array
        bias_arr_loc = zeros(reso(i) + 1, 1);

        % Quantile matching (increasing)
        obs_pos = vals_obs(1);
        cum_obs = cumsum(dist_obs_loc);
        cum_model = cumsum(dist_model_loc);
        for mdl_pos = vals_model(1):vals_model(end)
            while cum_obs(obs_pos) < cum_model(mdl_pos)
                obs_pos = obs_pos + 1;
            end
            bias_arr_loc(mdl_pos) = obs_pos;
        end

        % Quantile matching (decreasing)
        obs_pos2 = vals_obs(end);
        for mdl_pos = vals_model(end):-1:vals_model(1)
            while sum(dist_obs_loc(obs_pos2:end)) < sum(dist_model_loc(mdl_pos:end))
                obs_pos2 = obs_pos2 - 1;
            end
            bias_arr_loc(mdl_pos) = (bias_arr_loc(mdl_pos) + obs_pos2) / 2; % Average both matches
        end

        % Adjust tails (if extreme values need correction)
        if xtr > 1
            % Lower tail
            av_obs_low = compute_tail_avg(dist_obs_loc, unscaled_reso_i, xtr, 'lower');
            [av_model_low, i_p_low] = compute_tail_avg(dist_model_loc, unscaled_reso_i, xtr, 'lower'); 
            i_p_low = i_p_low - 1;
            bias_arr_loc(1:i_p_low) = min(unscaled_reso_i(1:i_p_low) + av_obs_low - av_model_low, bias_arr_loc(2:i_p_low + 1));

            % Upper tail
            av_obs_high = compute_tail_avg(dist_obs_loc, unscaled_reso_i, xtr, 'upper');
            [av_model_high, i_p_high] = compute_tail_avg(dist_model_loc, unscaled_reso_i, xtr, 'upper'); 
            i_p_high = i_p_high + 1;

            cal_vals = strcmp(var_name, 'pr') * (unscaled_reso_i * (av_obs_high / av_model_high)) + ...
                       (~strcmp(var_name, 'pr')) * (unscaled_reso_i + av_obs_high - av_model_high);

            bias_arr_loc(i_p_high:end) = max(cal_vals(i_p_high:end), bias_arr_loc(i_p_high - 1:end - 1));
        end

        % Linear interpolation over adjacent values that are the same
        for w = 2:reso(i)
            if bias_arr_loc(w) == bias_arr_loc(w - 1)
                ww = w;
                while ww < reso(i) && bias_arr_loc(ww) == bias_arr_loc(ww - 1)
                    ww = ww + 1;
                end
                if ww < reso(i)
                    increment = (bias_arr_loc(ww) - bias_arr_loc(w - 1)) / (ww - (w - 1));
                    for www = w:(ww - 1)
                        bias_arr_loc(www) = bias_arr_loc(w - 1) + (www - (w - 1)) * increment;
                    end
                end
            end
        end

        % Fill data gap with data interpolation
        j = 2;
        ri_max = -45;
        ri = (bias_arr_loc(j) - bias_arr_loc(j - 1))/(1);
        while ri >= ri_max && j < length(bias_arr_loc)
            j = j + 1;
            ri = (bias_arr_loc(j) - bias_arr_loc(j - 1))/(1);
        end

        xi = j-1:length(bias_arr_loc);
        mi = (bias_arr_loc(end) - bias_arr_loc(j - 1)) / (xi(end) - xi(1));
        bi = bias_arr_loc(j - 1) - mi * xi(1);
        k = 1;
        for n = xi
            bias_arr_loc(n) = mi * xi(k) + bi;
            k = k + 1;
        end

        % Unscale, smooth, and rescale bias correction factors
        ef = 10;
        bias_arr_loc = smooth(unscale_data(bias_arr_loc, var_name, sf(i)), cal_smth, 'rloess');
%         bias_arr_loc = max(0, min(reso(i), bias_arr_loc.*ef)); % Ensure within valid range
        bias_arr_loc = scale_data(bias_arr_loc, var_name, sf(i));

        % Store results
        bias_arr_ssn_i(:, ssn) = bias_arr_loc;
    end

    bias_arr_ssn.(stations{i}) = bias_arr_ssn_i;
end

toc
%% SECTION 4: APPLYING BIAS-CORRECTION FACTORS TO THE COMPLETE MODEL DATASET

disp('Applying calibration factors');
tic

% --- Define season indices ---
m_com = month(Date_comp);
if any(isnan(m_com)) % Check integrity
    error('NaN values detected in Date_comp, check data integrity!');
end

if strcmp(cal_trh, 'mnth')
    sea_c = m_com;
else
    sea_c = zeros(length(m_com), 1);
    sea_c(ismember(m_com, [9, 10, 11])) = 1; % Autumn (SON)
    sea_c(ismember(m_com, [12, 1, 2])) = 2;  % Winter (DJF)
    sea_c(ismember(m_com, [3, 4, 5])) = 3;   % Spring (MAM)
    sea_c(ismember(m_com, [6, 7, 8])) = 4;   % Summer (JJA)
end

% --- Create output directory ---
mkdir("Bias_corrections");

% --- Convert bias correction factors to anomalies ---
for i = 1:numStations
    bias_arr_ssn_i = bias_arr_ssn.(stations{i}) - (0:reso(i))';

    % Unscale bias correction
    u_bias_i = unscale_data(bias_arr_ssn_i, var_name, sf(i));

    % Generate plot
    figure("Visible", "off");
    plot(unscale_data(1:length(bias_arr_ssn_i), var_name, sf(i)), u_bias_i);
    xlabel("Rainfall (mm(3h)^{-1})");
    ylabel("Bias corr. (mm)");
    grid on;
    title(stations{i});

    if strcmp(cal_trh, 'mnth')
        legend({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}, ...
            'Location','best','NumColumns',2)
    else
        legend('SON','DJF','MAM','JJA','Location','best')
    end
    
    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 16, 'FontName', 'Times');
    
    % Save figure
    exportgraphics(gcf, fullfile('Bias_corrections', strcat("Bias-Corr_", stations{i}, ".png")), 'Resolution', 300);

    % Update bias array
    bias_arr_ssn.(stations{i}) = bias_arr_ssn_i;
end

% --- Prepare model data for correction ---
for i = 1:numStations
    % Limit data range
    Model_COSMO_complete(:, i) = min(max(Model_COSMO_complete(:, i), 0), MaxAbs(i));

    % Scale data
    Model_COSMO_complete(:, i) = scale_data(Model_COSMO_complete(:, i), var_name, sf(i));
end

% --- Apply bias correction ---
for station = 1:length(stations)
    
    % Select the bias array for this particular station
    bias_arr_ssn_i = bias_arr_ssn.(stations{station});
    
    for i = 1:length(Model_COSMO_complete)        
        val = round(Model_COSMO_complete(i, station));
        Model_COSMO_complete(i, station) = Model_COSMO_complete(i, station) + bias_arr_ssn_i(val + 1, sea_c(i));
    end
end

% --- Unscale corrected data & store results ---
BC_Model_complete = struct();
for i = 1:numStations
    Model_COSMO_complete(:, i) = unscale_data(Model_COSMO_complete(:, i), var_name, sf(i));
    BC_Model_complete.(stations{i}).data = Model_COSMO_complete(:, i);
end

toc

%% Internal functions

function cmb_hist = comb_hist(obs,mod_unc,mod_corr,min_hRF,statID,sea,bin_W)
    figure('WindowState', 'maximized') 
    
    tiledlayout(1,2)
    
    nexttile
    
    histogram(obs(obs>=min_hRF),"Normalization","probability","BinWidth",bin_W)        
    
    hold on
    
    histogram(mod_unc(mod_unc>=min_hRF),"Normalization","probability","BinWidth",bin_W)        
    
    title(strcat(statID,'-',sea))
    legend("obs","mod-unc")
    xlabel("Rainfall (mm)")
    ylabel("Frequency (-)")
    set(gca,'FontSize',24,'FontName', 'Times')
    set(gcf,'color','w');
    grid on
    xlim([10 100])
    ylim([0 0.005])
    
    %Plot histograms from corrected data
    
    nexttile
    
    histogram(obs(obs>=min_hRF),"Normalization","probability","BinWidth",bin_W)        
    
    hold on
    
    histogram(mod_corr(mod_corr>=min_hRF),"Normalization","probability","BinWidth",bin_W)        
    
    title(strcat(statID,'-',sea))
    legend("obs","mod-corr")
    xlabel("Rainfall (mm)")
    ylabel("Frequency (-)")
    set(gca,'FontSize',24,'FontName', 'Times')
    set(gcf,'color','w');
    grid on
    xlim([10 100])
    ylim([0 0.005])
end

% Helper function for tail averaging
function [avg, i] = compute_tail_avg(dist_loc, unscaled_reso, xtr, tail)
    if strcmp(tail, 'lower')
        idx = 1:length(dist_loc);
    elseif strcmp(tail, 'upper')
        idx = length(dist_loc):-1:1;
    else
        error('Invalid tail type');
    end
    
    count = 0;
    sum_vals = 0;
    i = 1;
    while count < xtr && i <= length(idx)
        sum_vals = sum_vals + dist_loc(idx(i)) * unscaled_reso(idx(i));
        count = count + dist_loc(idx(i));
        i = i + 1;
    end
    
    i = idx(i);
    avg = sum_vals / max(1, count);
end

function sf_opt = s_f(MaxAbs,reso,var_name)
    if strcmp(var_name, 'pr')
        % Define the objective function to minimize the absolute error
        objFun = @(sf) abs(reso - (log(MaxAbs * sf + 1) * sf));

        % Set initial guess and bounds for sf
        sf0 = 1; % Initial guess
        lb = 0; % Lower bound (scale factor should be non-negative)
        ub = 100; % Upper bound (adjust as needed)
        
        % Optimize using fmincon
        options = optimoptions('fmincon', 'Display', 'off'); % Suppress output
        sf_opt = fmincon(objFun, sf0, [], [], [], [], lb, ub, [], options);
    else
        error('Variable name not recognized or scaling not implemented for this variable.');
    end
end

function scaled_data = scale_data(data, var_name, sf)
    % Scale the data array 'data' for better representation of extremes in histograms.
    % Inputs:
    %   - data: Input data to scale
    %   - var_name: Variable type ('pr' in this case)
    % Output:
    %   - scaled_data: Scaled data array

%     sf = 14.54279172;

    if strcmp(var_name, 'pr')
        scaled_data = log(data * sf + 1) * sf;
    else
        error('Variable name not recognized or scaling not implemented for this variable.');
    end
end

function unscaled_data = unscale_data(data, var_name, sf)
    % Revert scaled data back to its original values.
    % Inputs:
    %   - data: Scaled data to unscale
    %   - var_name: Variable type ('pr' in this case)
    % Output:
    %   - unscaled_data: Original unscaled data array

%     sf = 14.54279172;

    if strcmp(var_name, 'pr')
        unscaled_data = (exp(data / sf) - 1) / sf;
    else
        error('Variable name not recognized or unscaling not implemented for this variable.');
    end
end

function [optMinH, errRD] = computeOptimalMinH(Observed, Model_COSMO, monthIndices, minH_range)

% Initialize output arrays
errRD = zeros(12, size(Observed,2));
optMinH = zeros(12, size(Observed,2));

% Loop over stations
for i = 1:size(Observed,2)
    % Loop over months
    for m = 1:12
        % Create mask for the current month
        monthMask = (monthIndices == m);
        
        % Count observed rainy days
        numRD_Omo = sum(Observed(monthMask, i) > 0);
        
        % Extract model data for the current month and station
        modelData = Model_COSMO(monthMask, i);
        
        % Compare each model value to all minH thresholds (N×H matrix)
        exceedMatrix = modelData > minH_range;  % Implicit broadcasting (N×H)
        
        % Count exceedances for each minH (1×H vector)
        numRD_Mmo = sum(exceedMatrix, 1);
        
        % Find the minH that minimizes the absolute error
        [errRD(m, i), bestIdx] = min(abs(numRD_Omo - numRD_Mmo));
        optMinH(m, i) = minH_range(bestIdx);
    end
end

save('optMinH.mat','optMinH')

end

end