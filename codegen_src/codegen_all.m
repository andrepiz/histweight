close all
clear
clc

%% SCRIPT NAME
% codegen_all
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Main codegen script to start code generation of ABRAM toolbox core modules
% -------------------------------------------------------------------------------------------------------------
%% NEEDED FROM BASE WORKSPACE
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUT TO BASE WORKSPACE
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-11-2024          Pietro Califano        Script initialized, added placeholders
% 14-11-2024          Pietro Califano        Added codegen setup for histweight (test passed)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------



%% GENERAL CODER CONFIG SPECIFICATION
coder_config = coder.config('mex');
coder_config.TargetLang = 'C++';
coder_config.GenerateReport = true;
coder_config.LaunchReport = true;

if strcmpi(coder_config.Name, 'MexCodeConfig')
    coder_config.EnableJIT = false;
elseif strcmpi(coder_config.Name, 'Embedded')
end

coder_config.MATLABSourceComments = true;

% DEVNOTE: setup
% ----------------------------- Optimizations ---------------------------
% 
%                  CacheDynamicArrayDataPointer: true
%                     EnableAutoParallelization: false
%                                  EnableMemcpy: true
%                                  EnableOpenMP: true
%                               MemcpyThreshold: 64
%                            NumberOfCpuThreads: 0
%                            OptimizeReductions: false
%                              SIMDAcceleration: 'Portable'
% 
% -----------------------------------------------------------------------
coder_config.EnableAutoParallelization = true;
coder_config.EnableAutoParallelizationReporting = true;
coder_config.EnableOpenMP = true;
coder_config.OptimizeReductions = true;
coder_config.NumberOfCpuThreads = 6;


%% Codegen execution: histweight
histweight_codegen_setup;

targetFcnName = 'histweight';
makeCodegen(targetFcnName, args_cell, coder_config)


clearvars -except coder_config
return
%% Codegen execution: parhistweight
% targetFcnName = 'parhistweight';
% makeCodegen(targetFcnName, args_cell, coder_config)

% clearvars -except coder_config







