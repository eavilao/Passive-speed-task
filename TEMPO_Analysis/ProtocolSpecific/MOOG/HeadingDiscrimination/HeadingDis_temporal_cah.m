%-----------------------------------------------------------------------------------------------------------------------
%heading discrimination task: computing CP over the 1s window around the 'peak' of each cell's PSTH.
%--	02/22/08 AC
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_temporal_cah(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

OutputPath=['Z:\Users\Aihua\z_tempOutputs\'];
% OutputPath=['C:\Aihua\z_TempOutputs\'];
% OutputPath=['C:\output\Aihua\'];

[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Z:\Users\Aihua\HeadingDiscriAihua');
if Protocol == 100 | Protocol == 107  | Protocol == 104  
    Tuning_Mid1s(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);         
else
    % CP_PeakTime(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
    % CP_MaxVar(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
    % CP_MaxRate(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
%     CP_Mid1s(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
%    CP_Mid1s_fewReps(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
   CP_Forward_Choice(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analayse the eye data
% CP_EyeMovement(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analayse the Reaction time data
% CP_EyeMovement(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
% Psychometric_RT_cah(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compute the time course of choice probability
% addpath('Z:\Users\Aihua\HeadingDiscriSheng');
% if Protocol == 100 | Protocol == 107  | Protocol == 104  
%     [StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);    
%     WindowInterval=400;
%     f3dsimple_temporal(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath, WindowInterval);         
% else
%     HeadingDis_cum_origin_temporal(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);
% end


return;


