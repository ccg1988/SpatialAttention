% by CCG @ 2025-12-23

clear; clc; close all
load('SRF_shift_test_analysis.mat')
AvP_units_id = SRF_shift_test_analysis.unitsHvP;
N_units = size(AvP_units_id, 1); %>>>>>>>>>>>>>>208
rP_trial_mean = nan(N_units, 1);
rA_trial_mean = nan(N_units, 1);
rP_trial_sig = [] ;
rP_trial_nonsig = [] ;
rA_trial_sig = [] ;
rA_trial_nonsig = [] ;
AvP_units_raw = SRF_shift_test_analysis.hits_vs_passive;

%% 1. Prepare and Reshape the Data
% We need to convert the data from "Wide" (Col 1 vs Col 2) to "Long"
% (One single column for firing rate, with a new column defining "Condition")

% A. Extract the columns
passive_data = AvP_units_raw(:, 1);  % Passive Firing Rates
active_data  = AvP_units_raw(:, 2);  % Active Firing Rates
location     = AvP_units_raw(:, 28); % Location
monkey_id    = AvP_units_raw(:, 30); % Monkey ID
unit_id_raw  = AvP_units_raw(:, 31); % Unit ID

% B. Create a Unique Unit Identifier
% Important: If Monkey 2 has "Unit 1" and Monkey 3 has "Unit 1", MATLAB
% will think they are the same cell. We combine MonkeyID and UnitID to fix this.
unique_unit_id = (monkey_id * 1000) + unit_id_raw; 

% C. Stack the data into "Long" format
% We stack Passive on top of Active
y_firing_rate = [passive_data; active_data];

% Create the Condition label (0 = Passive, 1 = Active)
% The first half is Passive (0), the second half is Active (1)
n_obs = length(passive_data);
condition = [zeros(n_obs, 1); ones(n_obs, 1)]; 

% Duplicate the grouping variables to match the new length
loc_long  = [location; location];
unit_long = [unique_unit_id; unique_unit_id];

%% 2. Create the Table
% Linear Mixed Models in MATLAB require a 'table' data structure.
% We cast Condition and Location as categorical variables.
tbl = table(y_firing_rate, condition, loc_long, unit_long, ...
    'VariableNames', {'FiringRate', 'Condition', 'Location', 'UnitID'});

% Convert to categorical factors for the stats model
tbl.Condition = categorical(tbl.Condition, [0 1], {'Passive', 'Active'});
tbl.Location  = categorical(tbl.Location); 
tbl.UnitID    = categorical(tbl.UnitID);

%% 3. Run the Linear Mixed-Effects Model (LMM)
% Formula: FiringRate ~ Condition + Location + (1 | UnitID)
% Translation: "Model FiringRate based on Condition and Location (Fixed Effects),
% while acknowledging that data points with the same UnitID have a shared
% random baseline (Random Intercept)."

disp('Running Linear Mixed-Effects Model...');
lme = fitlme(tbl, 'FiringRate ~ Condition + Location + (1|UnitID)');

%% 4. Display Results
disp('--------------------------------------------------');
disp('Model Coefficients (Estimate of the effect size):');
disp(lme.Coefficients);

disp('--------------------------------------------------');
disp('ANOVA Table (Use these p-values for the paper):');
% This tests if "Condition" (Active vs Passive) is statistically significant
% accounting for the repeated measures.
stats = anova(lme); 
disp(stats);

% Optional: Check for normality of residuals (Reviewer #2 "Analytically Rigorous")
figure;plotResiduals(lme, 'histogram');
title('Residuals of the LMM (Should look roughly Bell-Shaped)');