function prs = default_prs(exp_name,monk_id)

switch exp_name
    case 'linearspeed'
        prs.InfoFile = 'monkeyInfoFile_linearSpeed';
        prs.tempo(1).label = 'OUTCOME';
        prs.tempo(1).col = 9;
        prs.tempo(1).val = '0'; % 0 for successful trial
        prs.tempo(2).label = 'AMPLITU';
        prs.tempo(2).col = 17:21;
        prs.tempo(3).label = 'STIM_TYPE';
        prs.tempo(3).col = 11:15;
        switch monk_id
            case 41
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.5 -0.5];
            case 44
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.5 -0.5];
            case 45
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.5 -0.5];
        end
        prs.filepath = ['Z:\Users\Kaushik\Plexon files\Speed experiment\linearspeed\'];
        prs.tstim_on = 0;
        prs.tbeg_acc = 0.4;
        prs.tbeg_dec = 1;
        prs.tstim_off = 1.5;
    case 'angularspeed'
        prs.InfoFile = 'monkeyInfoFile_angularSpeed';
        prs.tempo(1).label = 'OUTCOME';
        prs.tempo(1).col = 9;
        prs.tempo(1).val = '0'; % 0 for successful trial
        prs.tempo(2).label = 'ROT_AMPLI';
        prs.tempo(2).col = 21:25;
        prs.tempo(3).label = 'STIM_TYPE';
        prs.tempo(3).col = 11:15;
        switch monk_id
            case 41
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.5 -0.5];
            case 44
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.5 -0.8];
            case 45
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.5 -0.8];
        end
        prs.filepath = ['Z:\Users\Kaushik\Plexon files\Speed experiment\angularspeed\'];
        prs.tstim_on = 0;
        prs.tbeg_acc = 0.4;
        prs.tbeg_dec = 1;
        prs.tstim_off = 1.5;
    case '1DAzi'
        prs.InfoFile = 'monkeyInfoFile_1DAzi';
        prs.tempo(1).label = 'OUTCOME';
        prs.tempo(1).col = 9;
        prs.tempo(1).val = '0'; % 0 for successful trial
        prs.tempo(2).label = 'AZIMUTH';
        prs.tempo(2).col = 9:13;
        prs.tempo(3).label = 'STIM_TYPE';
        prs.tempo(3).col = 11:15;
        switch monk_id
            case 41
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.9 -0.9];
            case 44
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.9 -0.9];
            case 45
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.9 -0.9];
        end
        prs.filepath = ['/Users/eavilao/Documents/Temp_data/1DAziHD/'];
        prs.tstim_on = 0;
        prs.tbeg_acc = 0.4;
        prs.tbeg_dec = 2.4;
        prs.tstim_off = 2.8;
    case 'HD'
        prs.InfoFile = 'monkeyInfoFile_HD';
        prs.tempo(1).label = 'OUTCOME';
        prs.tempo(1).col = 9;
        prs.tempo(1).val_corr = '0'; % 0 for successful trial
        prs.tempo(1).val_err = '5'; % 5 for incorrect trial
        prs.tempo(2).label = 'STIM_TYPE';
        prs.tempo(2).col = 13:15;
        prs.tempo(3).label = 'HEADING';
        prs.tempo(3).col = 9:13;
        switch monk_id
            case 44
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.25 -0.25]; % 0.25 after stim on / 0.25 before stim ends
            case 45
                prs.tspk = [-0.5 +0.5];
                prs.nspk = [0.25 -0.25]; % 0.25 after stim on / 0.25 before stim ends
        end
        prs.filepath = ['/Users/eavilao/Documents/Temp_data/1DAziHD/'];
        prs.tstim_on = 0;
        prs.tbeg_acc = 0.1;
        prs.tbeg_dec = 1;
        prs.tstim_off = 1;
end
prs.dt = 0.01; % temporal resolution (seconds)
prs.tsmooth = 0.1; %0.05 speed paper; %0.025; % width of gaussian (seconds)for speed protocols %0.1 for 1DAzi;  
prs.tsmooth_pop = 0.1; %0.04 speed paper; % width of gaussian (seconds) for population averaged psth %0.1 for 1DAzi;
prs.alpha = 0.05; % significance level
prs.Nd = 100; % LFP downsample factor
prs.spkwf_fs = 40000; % sampling rate of spike waveform
 % time-windows for transient detection
prs.twin = 0.1;  
prs.twin_current = prs.twin/2;
prs.twin_previous = prs.tbeg_acc - prs.twin/2;
prs.win_after_motion = 0.25;
prs.doBootstrap = 1; % bootstrapping to calc std dev in visual responses
prs.bootnum = 100; % bootstrap iterations