%----------------------------------------------------------------------------------------------------------
%-- Path_Defs.m: Path Definitions for various common and protocol-specific analysis routines.  You can
%--	simply modify BASE_PATH to change the base directory for all routines, e.g., when moving code
%--	to a different machine.  GCD, 1/3/2000
%----------------------------------------------------------------------------------------------------------

%Base path specification for protocol-specific analysis routines
BASE_PATH = 'Z:\LabTools\Matlab\TEMPO_Analysis\';
global g_path_defs     %Jing 10/15/2012

%add path for common tools
if isempty(g_path_defs) || g_path_defs ~= 7       %Jing 10/15/2012
    
junk_str = [BASE_PATH 'CommonTools\'];
addpath(junk_str);


g_path_defs = 7;      %Jing 10/15/2012
end