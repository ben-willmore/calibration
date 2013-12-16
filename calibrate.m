%% new calibration code BW Nov 2013
tdt50k = 48828.125;

%% parameters that can be changed
sampleRates = [tdt50k tdt50k*2];
cutoffs{1} = [500 tdt50k/2];
cutoffs{2} = [1000 32000];

zBusNum = 1;
deviceName = 'RX6';
channelNames = {'Left-hand earphone (BLUE)'; 'Right-hand earphone (RED)'};
golay_rms = 0.05*10.^(-12/20);
recording_highpass_f = 150; % Hz; filter out frequencies below this to deal with low-frequency noise

%%
addpath('./general');
addpath('./jan');

%%
fprintf('\n\n\n');
fprintf('======================================\n');
subdir = input('Subdirectory to save in: \n   ','s');
dirname = fixpath(['./' subdir]);
fprintf('======================================\n\n');
mkdir_nowarning(dirname);

%% don't change below this
n_sampleRates = length(sampleRates);
n_channels = length(channelNames);

if ~exist('calibs', 'var')
    calib = struct;
    for samplerate_idx = 1:n_sampleRates
        calib.sampleRate = sampleRates(samplerate_idx);
        calibs{samplerate_idx} = calib;
    end
end

%% 1. get the level of the reference tone
fprintf('== Reference tone: measuring mic response: \n');
for samplerate_idx = 1:n_sampleRates
    calib = calibs{samplerate_idx};
    fprintf('= %0.0fHz: ', calib.sampleRate);
    if isfield(calib, 'recorded_rms_volts_per_pascal')
        fprintf('already done, skipping\n');
    else
        happy = false;
        while ~happy
            if exist('fakeData', 'var') && fakeData
              fprintf('Using fake data\n');
              f = load('fakeCalib.mat');
              calib.refTone = f.fakeCalib.refTone;
            else
              fprintf('Plug in reference sound and press a key...\n');
              pause;
              % record tone and calculate RMS
              % (this plays silence)
              calib.refTone.irf = play_and_analyse_golay(calib.sampleRate, zBusNum, deviceName, 11, 0, 5, 0, [], recording_highpass_f);
            end

            rec = calib.refTone.irf.input_buffer.chan1 - mean(calib.refTone.irf.input_buffer.chan1);
            rec = rec(round(25/1000*calib.sampleRate):end-round(25/1000*calib.sampleRate)); % chop off ends which get messed up by filtering
            t = ((1:length(rec))-1)/calib.sampleRate;
            figure(1);
            subplot(2,2,1);
            plot(t,rec);
            subplot(2,2,2);
            n_plot = round(calib.sampleRate/100);
            plot(t(1:300),rec(1:300));

            subplot(2,2,[3 4]);
            [Pxx, F] = pwelch(rec, [], [], [], calib.sampleRate);
            diff = abs(F-1000);
            f = find(diff<100);
            dF = F(2)-F(1);
            rms_volts = sqrt(trapz(Pxx(f))*dF);
            %plot(F,Pxx)

            loglog(F, sqrt(Pxx));
            hold all
            loglog(F(f), sqrt(Pxx(f)), 'ro-');
            hold off;
        
            s = '';
            while ~strcmpi(s, 'y') && ~strcmpi(s, 'n')
              s = input('  - happy [y/n]', 's');
              if strcmpi(s, 'y')
                  happy = true;
              end
            end
        end
        calib.recorded_rms_volts_per_pascal = rms_volts;

        calibs{samplerate_idx} = calib;
        fprintf('-> Done\n');
    end
end

%% 2. relative calibration for each channel in turn
fprintf('\n== Relative calibration: measuring frequency response:\n');
for channel = 1:n_channels
    channelName = channelNames{channel};
    fprintf('%s:\n', channelName);
    plugged_in = false;

    for samplerate_idx = 1:n_sampleRates
        calib = calibs{samplerate_idx};
                    
        fprintf('= %0.0fHz: ', calib.sampleRate);

        if ~isfield(calib, 'relCalibs')
            calib.relCalibs = {};
        end
        
        if length(calib.relCalibs)>=channel
            fprintf('already done, skipping\n');
            continue;
        end
        if ~plugged_in
            fprintf('Plug in %s and press a key...\n', channelName);
            pause;
            plugged_in = true;
        end
        calib.relCalibs{channel} = relative_calib(calib.sampleRate, zBusNum, deviceName, golay_rms, channelName(1), ...
            cutoffs{samplerate_idx}, calibs{samplerate_idx}.recorded_rms_volts_per_pascal, recording_highpass_f, dirname);
        calib.relCalibs{channel}.channelIdx = channel;
        calib.relCalibs{channel}.channelName = channelName;
        calibs{samplerate_idx} = calib;
        fprintf('-> Done');
    end
end

%% 3. save
save(sprintf('%s/calibration.mat', dirname), 'calibs');

full_calibs = calibs;
calibs = {};
for samplerate_idx = 1:length(full_calibs)
    full_calib = full_calibs{samplerate_idx};
    calib = struct;
    calib.sampleRate = full_calib.sampleRate;
    calib.recorded_rms_volts_per_pascal = full_calib.recorded_rms_volts_per_pascal;
    
    relCalibs = {};
    for channel_idx = 1:length(full_calib.relCalibs)
        full_relcalib = full_calib.relCalibs{channel_idx};
        relCalib = struct;
        relCalib.channelIdx = full_relcalib.channelIdx;
        relCalib.channelName = full_relcalib.channelName;
        relCalib.ADrate = full_relcalib.ADrate;
        relCalib.filter = full_relcalib.filter;
        relCalibs{channel_idx} = relCalib;
    end
    calib.relCalibs = [relCalibs{:}];
    calibs{samplerate_idx} = calib;
end
calibs = [calibs{:}];

save(sprintf('%s/compensation_filters.mat', dirname), 'calibs');
