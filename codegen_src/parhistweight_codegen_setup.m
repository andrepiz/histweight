% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Codegen setup script for parhistweight functional module
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



SIZE_SOLVED_FOR    = 15;
CHBV_MAXDEG        = 20;
MAX_CHBV_SWITCHES  = 5;
ATM_MODEL_LUT_SIZE = 28;

i_dxState       = coder.typeof(0, [2*SIZE_SOLVED_FOR, 1], [false, false]); 
i_dPxState      = coder.typeof(0, [SIZE_SOLVED_FOR, SIZE_SOLVED_FOR], [false, false]);

i_dDeltaTime    = coder.typeof(0, [1, 1], [false, false]);
i_dStateTimetag = coder.typeof(0, [1, 1], [false, false]);

% Dynamics parameter struct
strMoonEPHdata.dChbvPolycoeffs = coder.typeof(0,        [3*CHBV_MAXDEG, 1], [true, false]);
strMoonEPHdata.ui8PolyDeg      = coder.typeof(uint8(0), [1,1],            [false,false]);
strMoonEPHdata.dTimeLowBound   = coder.typeof(0,        [1,1],            [false,false]); 
strMoonEPHdata.dTimeUpBound    = coder.typeof(0,        [1,1],            [false,false]);

i_strDynParams.strMoonEPHdata = coder.typeof(orderfields(strMoonEPHdata));

% DEVNOTE: remember to orderfields before codegen

% DEFINITION OF ARGUMENTS CELLS 
% DEVNOTE: Achtung: the struct inputs MUST be defined in the same way. Otherwise the C struct has to be
% defined twice, and this is not possible with the same name, i.e. duplicated memory storage needed. 
% The functions won't accept the same struct as input because of missing or additional fields.


args_cell{1} 









