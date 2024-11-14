function [] = makeCodegen(targetFcnName, args_cell, coder_config)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Automatic code generation makers
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-04-2022        Pietro Califano         First version. Very basic codegen call.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% assert(strcmpi(coder_config.Name, 'MexCodeConfig'), 'Only MEX config is currently supported!')

%% Coder settings
if nargin < 3
    fprintf("No coder configuration object specified. Using default configuration...\n")
    coder_config = coder.config('mex');
    coder_config.TargetLang = 'C++';
    coder_config.GenerateReport = true;
    coder_config.LaunchReport = true;
    coder_config.EnableJIT = false;
    coder_config.MATLABSourceComments = true;
end

% IF PARALLEL REQUIRED
% coder_config.EnableAutoParallelization = true;
% coder_config.EnableOpenMP = true;

%% Target function details
% Get number of outputs
numOutputs = nargout(targetFcnName);
% Extract filename and add MEX indication
[~, targetFcnName, ~] = fileparts(fullfile(targetFcnName));
outputFcnName = strcat(targetFcnName, '_MEX');

% numOfInputs; % ADD ASSERT to size of args_cell from specification functions

% ENTRY-POINT FUNCTION ARGUMENTS DEFINITION
% NOTE: third option argument variable_dimensions is an array of bools, one
% for each dimension of the array

% coder.typeof(example_value, size_vector, variable_dims);

% coder.getArgTypes % Function call to automatically specify input
% arguments. Note that this takes the specific sizes used in the function
% call


% Defining structures 
% struct1.fieldname1 = coder.typeof(0,[3 5],1);
% struct1.fieldname2 = magic(3);
% coder.typeof(struct1); % Defines

% Defining nested structures
% S = struct('fieldname1', double(0),'fieldname2',single(0)); % Inner structure
% SuperS.x = coder.typeof(S);
% SuperS.y = single(0); 
% coder.typeof(SuperS) % Outer structure 

% IMPORTANT: in recent MATLAB versions (>2021), arguments can be speficied
% directly in function --> no need to define them outside in -args
% arguments
%     u (1,4) double
%     v (1,1) double

%% CODEGEN CALL
fprintf("\nCODE GENERATION EXECUTION: STARTED\n")
% Execute code generation
codegenCommands = {strcat(targetFcnName,'.m'), "-config", coder_config,...
    "-args", args_cell, "-nargout", numOutputs, "-o", outputFcnName};
codegen(codegenCommands{:});
fprintf("\nCODE GENERATION EXECUTION: COMPLETED\n")
end
