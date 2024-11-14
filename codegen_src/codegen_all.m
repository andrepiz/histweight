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







