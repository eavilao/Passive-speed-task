function MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

switch(Analysis{1})
    case  'turn off'     
        system('shutdown -s') ;
    case 'Plot Tuning Surface'
        DirectionTuningPlot_3D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface (Aihua)'
            DirectionTuningPlot_3D_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Plot Tuning Surface (Aihua)1'
%         DirectionTuningPlot_3D_Aihua1(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Mean Variance Firing Output (Tunde)',
        Mean_Var_Firing_Output(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface (Tanya)',
        DirectionTuningPlot_3D_Tanya(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_yong'
        DirectionTuningPlot_3D_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning SurfaceLFP_yong'
        DirectionTuningPlot_3DLFP_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface (Sheng)'
        DirectionTuningPlot_3D_sheng(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_LFP_xiaodong'
        DirectionTuningPlot_3D_LFP_xiaodong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Plot Lambert Tuning Surface'
        Lambert_DirectionTuningPlot_3D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'TuningPlot'
        TuningPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin,PATH, FILE);
    case 'Plot Tuning Azimuth'
        DirectionTuningPlot_1D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Azimuth_Adhira'
        DirectionTuningPlot_1D_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);   
    case 'Direction/Disparity Tuning (Aihua)'
        Direction_Disparity_Tuning_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);         
    case 'Plot Tuning Azimuth PSTH',
        Azimuth_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'DirectionTuningPlot_3D_pairwiseunits_yong'
        DirectionTuningPlot_3D_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_3D_pairwiseunits_sam'
        DirectionTuningPlot_3D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_eyetrace'
        DirectionTuningPlot_1D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_eyetrace_noSpikes'
        DirectionTuningPlot_1D_eyetrace_noSpikes(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'DirectionTuningPlot_1D_firingrate'
        DirectionTuningPlot_1Dfiringrate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_firingrate (Tunde)'
        DirectionTuningPlot_1Dfiringrate_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'DirectionTuningPlot_1D_pairwiseunits_yong'
        DirectionTuningPlot_1D_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_pairwiseunits_sam'
        DirectionTuningPlot_1D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_ves'
        Directiontuningplot_ves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Event Times'
        EventTimes(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot PSTH'
        MOOG_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH (Tanya)'
        MOOG_PSTH_Tanya(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_yong'
        PSTH_Translation_yong(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'MOOG_PSTH_xiongjie'
        MOOG_PSTH_xiongjie(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_Anuk'
        PSTH_Anuk(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Corr_Analysis'
        PSTH_Corr_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Corr_Analysis_fix'
        PSTH_Corr_Analysis_fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_Anuk_fix'
        PSTH_Anuk_fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run Frequency Analysis'
        Frequency_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run Frequency Analysis alternate'
        Frequency_Analysis_alternate(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_Anuk_HTI_fix'
        PSTH_Anuk_HTI_fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Spike Rasters'
        PlotRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, Protocol);
    case 'Plot Psycho_neuro_cum'
        HeadingDis_cum(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cum_Adhira'
        HeadingDis_cum_Adhira(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cum (Sheng)'
        HeadingDis_cum_Sheng(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Heading_tuning'
        HeadingDis_fixation(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric'
        Psychometric(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Conflict)'
        Psychometric_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Psychometric (Reaction Time Task)'
        Psychometric_RT(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Reaction Time Analysis (LIP_Adhira)'
        Reaction_Time_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Psychometric Conflict'
        Psychometric_1I_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Adaptation)'
        Psychometric_Adaptation(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);                
    case 'Plot Psychometric (Adam)'
        Psychometric_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);  
    case 'Plot Psychometric (Kalpana)'
        Psychometric_Kalpana(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE); 
    case 'Plot Psychometric (Adaptation, Adam)'
        Psychometric_Adaptation_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
     case 'Plot Psychometric (Ari)'
        Psychometric_Ari(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);          
    case 'Plot Microstim'
        Heading_microstim(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Psychometric_pursuit'
        Psychometric_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Psychometric_eye_pursuit'
        addpath('Z:\Users\sheng\training')
     
        Psychometric_eye_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
   
    case 'Plot eye_pursuit_HD_eyemovement'
        addpath('Z:\Users\sheng\training')
        try,
         Psychometric_eye_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
        catch,
        end
    case 'Plot PURSUIT_HEADING_eyemovement'
        PURSUIT_HEADING_eyemovement(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Accelerometer_cum'
        Accelerometer_cum(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Accelerometer_cum (Tunde)'
        Accelerometer_cum(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot CP Distribution'
        Heading_CP_Distrib(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot HeadingDiscimination_PSTH'
        HeadingDis_cum_PSTH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Calculate HeadingDiscimination_DFT'
        HeadingDis_cumDFT(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Calculate HeadingDiscimination_Corr'
        HeadingDis_cumCorr(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'HeadingDis_cum_eyetrace'
        HeadingDis_cum_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'HeadingDis_cum_eyetrace_staircase'
        HeadingDis_cum_eyetrace_staircase(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'shiftwindow'
        HeadingDis_cum_shiftwindow(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cumVarience'
        HeadingDis_cumVarience(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cumTimecourse'
        HeadingDis_cumTimecourse(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot session_yong'
        HeadingDis_cum_session(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'HeadingDis_cum_pairwiseunits_yong'
        HeadingDis_cum_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'HeadingDis_cum_pairwiseunits_sam'
        HeadingDis_cum_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_Fix_sam'
        DirectionTuningPlot_Fix_SUSU(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_Az_VF_sam'
        DirectionTuningPlot_1D_Az_VF_SUSU(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'LessSampleFit'
        LessSampleFit(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Eye trace'
        DirectionTuning3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Eye trace (tunde)'
        DirectionTuning3D_eyetrace_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);            
    case 'Plot Fixation Tuning Surface'
        DirectionTuningPlot_Fix(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Fixation Tuning Surface (Tunde)'
        DirectionTuningPlot_Fix_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Fixation Tuning Surface (Aihua)'
        DirectionTuningPlot_Fix_aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);     
    case 'Wrapped Gaussian Fit (Tunde)'
        WrappedGaussianFit_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Lambert Fixation Tuning Surface',
        Lambert_DirectionTuningPlot_3D_fix(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Fixation PSTH'
        MOOG_PSTH_Fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Output Data for Curve Fitting'
        DirectionTuningPlot_Curvefit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_conflict'
        Direction2d_cue_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_conflict (Aihua)'
        Direction2d_cue_conflict_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Fit_cue_conflict_2D'
        Direction2d_cue_conflict_fit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Info_cue_conflict_2D'
        Direction2d_cue_conflict_info(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Basic_cue_conflict_2D'
        Direction2d_cue_conflict_basic(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning_yong'
        Direction2d_cue_conflict_basic_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Fit Optic Flow Tuning (Zack)'
        Heading_CurveFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Fit Optic Flow Tuning_Rot (Zack)' % for rotation
        Heading_CurveFit_Rot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation Tuning Curves_Azimuth'
        DirectionTuningPlot_1D_Az_VF(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation Tuning Curves_Azimuth (Tunde)'
        DirectionTuningPlot_1D_Az_VF_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation Tuning Curves_Elevation'
        DirectionTuningPlot_1D_El_VF(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation PSTH'
        MOOG_PSTH_Fix_1D(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot 1D Fixation Tuning Curves_Azimuth Head-Eye'
        DirectionTuningPlot_1D_Az_VF_headeye(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Rotation Tuning 3D'
        Rotation3Dtuning(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Plot Rotation Tuning 3D (Aihua)'
        Rotation3Dtuning_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Plot Rotation Tuning 3D (Adhira)'
        Rotation3Dtuning_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);  
    case 'Plot Rotation Tuning 3D (Sheng)'
        Rotation3Dtuning_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);                   
    case 'Plot Lambert Rotation Tuning 3D'
        Lambert_Rotation3Dtuning(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case '2D Frequency_Analysis'
        Frequency_Analysis_2d(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Rotation Frequency Analysis'
        Rotation_Frequency_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation PSTH'
        Rotation_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation PSTH (Tanya)'
        Rotation_PSTH_Tanya(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);    
     case 'Plot Rotation PSTH_xiongjie'
        Rotation_PSTH_xiongjie(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation PSTH_yong'
        PSTH_yong(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation Eye Trace (Katsu)'
        Rotation3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Rotation Eye Trace (Tunde)'
        Rotation3D_eyetrace_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Output Firing rate'
        FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Output Firing rate (Tunde)'
        FiringRate_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Output ROT_Firing rate'
        Rotation_FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot  TILT_TRANSLATION PSTH'
        TiltTrans_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);  %added by AHC
    case 'Plot Tilt/Trans Direction tuning'
        DirectionTuning2D_TiltTrans(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Pursuit Direction Tuning'
        DirectionTuning2D_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Pursuit PSTH'
        DirectionTuning2D_pursuit_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);%added by KT
    case 'MU activity'
        %DirectionTuningPlot_3D_munit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        MunitAnalysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'MU activity (Xiaodong)'
        %DirectionTuningPlot_3D_munit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        MunitAnalysis_Xiaodong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Plot PSTH_Anuk_rotation'
        PSTH_Anuk_rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Corr_Analysis_rotation'
        PSTH_Corr_Analysis_rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Gaussfit'
        %PSTH_Anuk_rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
        PSTH_Gaussfit(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Psychometric 2I'
        psychometric_2I(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric 2I Conflict'
        Psychometric_2I_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Sinusoid Analysis'
        Sinusoid_Analysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Accelerlometer study'
        Accelerometer_study(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tilt/Trans Eye Trace'
        TiltTrans_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Horizontal only (Katsu)'
        Horizontal_only(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Pursuit Eye Trace'
        Pursuit_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Rotation Horizontal only (Katsu)'
        Rotation_Horizontal_only(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Rotation_Frontal (Katsu)'
        Rotation_Frontal(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Save 4direct Pursuit Eye Trace'
        Pursuit_eachCell_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot  First1sec_PSTH'
        First1sec_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);  %added by Katsu
    case 'FanoFactor (Katsu)'
        FanoFactor(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Rotation3D_DDIperm (Katsu)'
        Rotation3D_DDIperm(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Direction3D_DDIperm (Katsu)'
        Direction3D_DDIperm(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Temporal_Analysis_Translation (Katsu)'
        Temporal_Analysis_Translation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Temporal_Analysis_Rotation (Katsu)'
        Temporal_Analysis_Rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Compute Vertical Slice'
        VerticalSlice(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Lambert Rotation Tuning 3D (Katsu)'
        Katsu_Lambert_plotter(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Output Rotation PSTH (Katsu)'
        Rotation_PSTH_output(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Analysis);
    case 'Output Translation PSTH (Katsu)'
        Translation_PSTH_output(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Analysis);
    case 'Frontal parralel plane'
        Frontal_parralel_plane(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Temporal Analysis (cah)'
        HeadingDis_temporal_cah(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot_PSTH_Michael'
        MOOG_PSTH_Michael(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'DirectionTuningPlot_1D (Tunde)'
        DirectionTuningPlot_1D_Az_VF_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Plot 1Dtuning_trapezoid (Yong)'
        AZIMUTH_TUNING_1D_TRAP_Tuning_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Plot 1Dtuningcoherence_trapezoid (Yong)'
        AZIMUTH_TUNING_1D_TRAP_Coherence_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Plot 1Dtuningtimecourse_trapezoid (Yong)'
        AZIMUTH_TUNING_1D_TRAP_timecourse_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Accelerometer_cum_1D (Jing)'
        Accelerometer_cum_1D_Jing(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Delayed Saccade Analysis (Adhira)'
        Delayed_Saccade_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);               
    case 'Delayed Saccade Analysis (Yong)'
        Delayed_Saccade_Analysis_Yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Delayed Saccade Analysis (Adam)'
        Delayed_Saccade_Analysis_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Memory Saccade Analysis (Adhira)'
        Memory_Saccade_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Galvo Pursuit'
        Delayed_Saccade_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Pursuit Visual (Adhira)'
        Pursuit_Visual_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);   
    case 'Pursuit-Combine (Adhira)'
        Pursuit_Combine_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Pursuit -Visual+Combine (Adhira)'
        Pursuit_VisualCombine_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);     
    case 'Plot Psychometric (Rot)'  %Jing 10/18/2012
        Psychometric_Rot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Rot, Two-Target)'  %Jing 12/17/2012
        Psychometric_Rot_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot YAW_ROT_FIXATION'
        YAW_ROT_FIXATION(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'FineRotDiscrm (Psycho)'
         Psychometric_FineRotDiscrm_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'FineRotDiscrm (Neuron)'    
             cindex = strfind(upper(FILE),  'C');
                    rindex = strfind(upper(FILE),  'R');
        if  str2num(FILE(cindex+1:rindex-1)) <= 19  
                        
                         HeadingRaw_cum_finehd_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                         Heading_PSTH_fine_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  10)
         else
                         HeadingRaw_cum_finehd(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        
                         Heading_PSTH_fine(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  10)
         end
    case 'FineRotDiscrm (Two-target)'
%         if sum(sum(data.spike_data(1,:,:))) < 100
         if isempty(data.spike_data)
             if SpikeChan > 1
                if SpikeChan == 1 
                   SpikeChan = 1;
                   Psychometric_FineRotDiscrm_2targ_CI(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                else
                   SpikeChan = 1;
                   seperatefile(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%                    Psychometric_FineRotDiscrm_2targ_deg(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                end 
             else
                Psychometric_FineRotDiscrm_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
             end
         else
            if  SpikeChan == 2
                Psychometric_FineRotDiscrm_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
            else
                if  SpikeChan == 3
                seperatefile(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%                    Psychometric_FineRotDiscrm_2targ_deg(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

                else
                    cindex = strfind(upper(FILE),  'C');
                    rindex = strfind(upper(FILE),  'R');
                    if  str2num(FILE(cindex+1:rindex-1)) <= 19  
                        
                         HeadingRaw_cum_finehd_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                         Heading_PSTH_fine_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  10)

                    else
                         HeadingRaw_cum_finehd(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        
                         Heading_PSTH_fine(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  10)
                    end
%                   Heading_PSTH2(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5)
                end
            end
         end
     case 'FineRotDiscrm  afferent (Two-target)'
%         if sum(sum(data.spike_data(1,:,:))) < 100
         if isempty(data.spike_data)
            Psychometric_FineRotDiscrm_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
         else
             % now afferent use middle 400,  vn use middle 500 
                 fix = choose(FILE);
                 if fix == 1
%                     HeadingRaw_cum_finehd_aff_opp(data, Protocol,  Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                    HeadingRaw_cum_finehd_aff_sv_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);      
%                     Heading_PSTH2_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5)
                 else
%                      HeadingRaw_cum_finehd_aff(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                     HeadingRaw_cum_finehd_aff_sv(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                    
%                      Heading_PSTH2(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5);
                 end
%                  HeadingRaw_cum_finehd_aff_ISI_2(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%                  HeadingRaw_cum_finehd_aff_ISI(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%                   HeadingRaw_cum_finehd_aff_ISI_c(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%                   HeadingRaw_cum_finehd_aff_ISI_c2(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%             HeadingRaw_cum_finehd_aff_opp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
         end
         
    case 'Plot Tuning Surface variance'
        addpath('Z:\Users\sheng\program\3dsimple')
         m_pval = 0.01;
         Outputpath = 'D:\rotation\';
         WindowInterval = 400;
%         if Protocol ~= 101
%              m_time = 400;
             % get spike count information, for spikes stastics
%              f3dsimple_spikes(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Outputpath, m_time, m_pval);
%          % for 1d   
%            f3dsimple_temporal(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Outputpath,  WindowInterval);
           %  for 3d
           try
                     if sum(sum(data.spike_data(SpikeChan,:,:))) > 100
                             f3dsimple_temporal_3d(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Outputpath,  WindowInterval);
                     end
                         
           catch
               
               if sum(sum(data.spike_data(1,:,:))) > 100
                  f3dsimple_temporal_3d(data, 1, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Outputpath,  WindowInterval);
               end
               
           end
           
          
         
          
             
    case  'Plot PURSUIT_TRAP'
             Yaw_pursuit(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
             Heading_PSTH_pursuit(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5);
      case  'Plot PURSUIT Eye Track'
             Yaw_pursuit_eyetrack(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
             Heading_PSTH_pursuit(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5);
             
    case 'Plot Psychometric (Varying FP)'  %Jing 11/05/2012
        Psychometric_VaryingFP(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Varying FP, Two-Target)'  %Jing 11/05/2012
        Psychometric_VaryingFP_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Two-target)'
                Psychometric_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE); 
    case 'Plot Psychometric (Rot, Courtney)'
        if isempty(data.spike_data)
             if SpikeChan == 3
                 seperatefile2(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%                  Psychometric_courtney_3(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
             else
                  Psychometric_courtney(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
             end
        else
              if SpikeChan == 2  
                   Psychometric_courtney(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
              else
                  if  SpikeChan  == 4
                        HeadingRaw_cum_old(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                        Heading_PSTH_yaw_old(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5)
                  else
                      if SpikeChan == 3
                        seperatefile2(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);  
%                         Psychometric_courtney_3(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                      else
                        HeadingRaw_cum(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                        Heading_PSTH_yaw(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,  5);
                      end
                  end
              end
        end

    case 'Plot Psychometric (Gabor)'  %Jing 10/19/2016
        Psychometric_Gabor(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Gabor,Two-target)'  %Jing 10/19/2016
        Psychometric_Gabor_2targ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
end
    
return;