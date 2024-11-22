%% SCRIPT NAME
% codegen_all
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Codegen setup script for histweight functional module
% -------------------------------------------------------------------------------------------------------------
%% NEEDED FROM BASE WORKSPACE
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUT TO BASE WORKSPACE
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-11-2024       Pietro Califano     Script intialized
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------


%% Function input specification
% Include and Source paths

% dCoords        (2, :) double {ismatrix}
% dValues        (1, :) double {isvector}
% dLimits        (:, 2) double {isvector} = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
% dGranularity   (1, 1) double {isscalar} = 1 % uint32?
% charMethod     (1, :) string            = 'area'
% bFlagProgress  (1, 1) logical           = false;
% bVECTORIZED    (1, 1) logical           = false;
% bDEBUG_MODE    (1, 1) logical           = false;

MAX_SIZE = 1e6;

dCoords        = coder.typeof(0,        [2, MAX_SIZE], [false, true]);
dValues        = coder.typeof(0,        [1, MAX_SIZE], [false, true]);
dLimits        = coder.typeof(0,        [2, 2], [false, false]);
dGranularity   = coder.typeof(0,        [1, 1], [false, false]);
i8MethodId     = coder.typeof(int8(1),  [1, 1], [false, false]);
bFlagProgress  = coder.typeof(false,    [1, 1], [false, false]);
bVECTORIZED    = coder.typeof(false,    [1, 1], [false, false]);
bDEBUG_MODE    = coder.typeof(false,    [1, 1], [false, false]);
dGaussianSigma = coder.typeof(0,        [1, 1], [false, false]);
dWindowSize    = coder.typeof(0,        [1, 1], [false, false]);


args_cell{1} = dCoords      ;
args_cell{2} = dValues      ;
args_cell{3} = dLimits      ;
args_cell{4} = dGranularity ;
args_cell{5} = i8MethodId   ;
args_cell{6} = bFlagProgress;
args_cell{7} = bVECTORIZED  ;
args_cell{8} = bDEBUG_MODE  ;
args_cell{9} = dGaussianSigma  ;
args_cell{10} = dWindowSize  ;





