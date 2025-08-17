%% run_CINDERellA.m
% Example script for running CINDERellA Bayesian Network Learning
% 
% This script demonstrates how to:
% 1. Set up the CINDERellA environment
% 2. Load gene expression data
% 3. Run network inference
% 4. Evaluate results against known networks
%
% Author: Shinya Tasaki, Ph.D.
% Compatible with MATLAB up to R2016a

%% Clear workspace and setup
clear all; close all; clc;

%% Step 1: Setup CINDERellA environment
fprintf('Setting up CINDERellA environment...\n');
CINDERellA_PATH = './functions';
addpath(CINDERellA_PATH);

%% Step 2: Load expression data
fprintf('Loading expression data...\n');
% Load your gene expression data (genes as rows, samples as columns)
expdata = read_exp('test_data/exp.txt');

% Display data information
fprintf('Data loaded successfully!\n');
fprintf('Number of genes: %d\n', size(expdata.data, 1));
fprintf('Number of samples: %d\n', size(expdata.data, 2));

%% Step 3: Run CINDERellA with basic parameters
fprintf('\nRunning CINDERellA network inference...\n');

% Basic run with default parameters
CINDERellA(expdata.data);

% Alternative: Run with custom parameters
% CINDERellA(expdata.data, ...
%     'output_dir', 'my_results', ...
%     'max_parents', 3, ...
%     'runtime_minutes', 30, ...
%     'num_samples', 1000, ...
%     'edge_threshold', 0.33, ...
%     'layout', 'force');

%% Step 4: Load results and evaluate (if true network is available)
fprintf('\nLoading and evaluating results...\n');

try
    % Load true network for comparison (if available)
    network = read_network('test_data/network.txt', size(expdata.data, 1));
    
    % Load learned edge frequencies
    edgefrq_data = dlmread('./CINDERellA_results/edgefrq.txt');
    
    % Convert to sparse matrix format
    nGenes = size(expdata.data, 1);
    edgefrq = sparse(edgefrq_data(:,1), edgefrq_data(:,2), edgefrq_data(:,3), nGenes, nGenes);
    
    % Perform evaluation
    fprintf('Evaluating network performance...\n');
    [AUCPR, AUCROC] = evaluation(edgefrq, network.data, 'plot', 1);
    
    % Display results
    fprintf('\n=== Evaluation Results ===\n');
    fprintf('Area Under Precision-Recall Curve (AUCPR): %.4f\n', AUCPR);
    fprintf('Area Under ROC Curve (AUCROC): %.4f\n', AUCROC);
    fprintf('==========================\n');
    
catch ME
    fprintf('Evaluation skipped: %s\n', ME.message);
    fprintf('This is normal if no true network file is available.\n');
end

%% Step 5: Display output files
fprintf('\n=== Output Files Generated ===\n');
output_dir = './CINDERellA_results/';
output_files = {
    'edgefrq.txt', 'Edge frequencies (main result)';
    'Mcmc.mat', 'Sampled networks';
    'Param.mat', 'Analysis parameters';
    'LS.mat', 'Local scores';
    'mcmc_diagnostics.png', 'MCMC convergence plot';
    'network_visualization.png', 'Network visualization'
};

for i = 1:size(output_files, 1)
    filename = [output_dir, output_files{i, 1}];
    if exist(filename, 'file')
        fprintf('✓ %s - %s\n', output_files{i, 1}, output_files{i, 2});
    else
        fprintf('✗ %s - Not found\n', output_files{i, 1});
    end
end

fprintf('\nResults saved in: %s\n', output_dir);
fprintf('\nAnalysis complete! Check the output directory for results.\n');

%% Optional: Advanced usage examples (commented out)
%{
%% Example 1: Using prior knowledge to constrain edges
fprintf('\n=== Example: Using Prior Knowledge ===\n');
nGenes = size(expdata.data, 1);
prior = ones(nGenes, nGenes);  % Start with all edges allowed

% Example: Disallow specific edges based on biological knowledge
% prior(1, 2) = 0;  % Disallow edge from gene 1 to gene 2
% prior(3, 4) = 0;  % Disallow edge from gene 3 to gene 4

% Run with prior constraints
CINDERellA(expdata.data, ...
    'prior_matrix', prior, ...
    'output_dir', 'results_with_prior', ...
    'runtime_minutes', 1);

%% Example 2: Parameter exploration for different settings
fprintf('\n=== Example: Parameter Exploration ===\n');

% Test different numbers of maximum parents
max_parents_list = [2, 3, 4];
runtime_minutes = 1;

for i = 1:length(max_parents_list)
    mp = max_parents_list(i);
    output_dir = sprintf('results_mp%d', mp);
    
    fprintf('Running with max_parents = %d...\n', mp);
    CINDERellA(expdata.data, ...
        'max_parents', mp, ...
        'runtime_minutes', runtime_minutes, ...
        'output_dir', output_dir);
end

%% Example 3: Different MCMC samplers comparison
fprintf('\n=== Example: MCMC Sampler Comparison ===\n');

samplers = {'M.c2PB', 'M.REV50', 'M.1PB'};
for i = 1:length(samplers)
    sampler = samplers{i};
    output_dir = sprintf('results_%s', strrep(sampler, '.', '_'));
    
    fprintf('Running with sampler = %s...\n', sampler);
    CINDERellA(expdata.data, ...
        'sampler', sampler, ...
        'runtime_minutes', 1, ...
        'output_dir', output_dir);
end
%}