%% new calibration code BW Nov 2013
tdt50k = 48828.125;

%% parameters that can be changed
sampleRates = [tdt50k tdt50k*2];
zBusNum = 1;
deviceName = 'RX6';
channelNames = {'Left-hand earphone (BLUE)'; 'Right-hand earphone (RED)'};
golay_rms = 0.05;

%%
addpath('./general');

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

%% 1. get the rms recorded voltage corresponding to 1Pa using the reference tone
fprintf('== Reference tone: measuring mic response: \n');
for samplerate_idx = 1:n_sampleRates
    calib = calibs{samplerate_idx};
    fprintf('= %0.0fHz: ', calib.sampleRate);
    if isfield(calib, 'recorded_rms_volts_per_pascal')
        fprintf('already done, skipping\n');
    else
        if exist('fakeData', 'var') && fakeData
          fprintf('Using fake data\n');
          f = load('fakeCalib.mat');
          calib.refTone = f.fakeCalib.refTone;
        else
          fprintf('Plug in reference sound and press a key\n');
          pause;
          % record tone
          calib.refTone.irf = play_and_analyse_golay(calib.sampleRate, zBusNum, deviceName, 11, 0, 5, 0, []);
        end
        
        % subtract mean
        rec = calib.refTone.irf.input_buffer.chan1 - mean(calib.refTone.irf.input_buffer.chan1);

        % tests
        if false
          % to test code, try a pure tone
          t = (1:length(rec))*1/calib.sampleRate;
          rec = sin(2*pi*1000*t);
        elseif false
          % try a filtered version of rec
          [B,A] = butter(3,[900 1100]/(calib.sampleRate/2));
          rec = filtfilt(B,A,rec);
        end

        refTone.total_variance = var(rec);

        % plot
        t = ((1:length(rec))-1)/calib.sampleRate;
        figure(1);
        subplot(2,2,1);
        plot(t,rec);
        
        % plot zoomed in
        subplot(2,2,2);
        n_plot = round(calib.sampleRate/100);
        plot(t(1:300),rec(1:300));

        % plot PSD
        subplot(2,2,[3 4]);
        [refTone.powerspec, refTone.freqs] = pwelch(rec, [], [], [], calib.sampleRate);
  
        % find region of power spectrum from 900:1100Hz to allow (amply) for discreteness of FFT
        d = abs(refTone.freqs-1000);
        f = find(d<100);
        d_freq = refTone.freqs(2)-refTone.freqs(1);

        % integrate over spectrum, multiplied by delta_frequency
        refTone.total_power = trapz(refTone.powerspec)*d_freq;
        % total_power should be equal to var(rec), but it isn't quite (2% difference for current data)
        %assert(abs(total_variance-var(rec))<eps(total_variance)); % doesn't quite work
        refTone.tone_power = trapz(refTone.powerspec(min(f):max(f)))*d_freq; % variance 
        refTone.tone_percent_power = refTone.tone_power/refTone.total_power*100;
        refTone.tone_rms = sqrt(refTone.tone_power);
        loglog(refTone.freqs, sqrt(refTone.powerspec))
        hold all
        loglog(refTone.freqs(f), sqrt(refTone.powerspec(f)), 'ro-');
        hold off;
        title(sprintf('Overall variance = %0.2e; total power = %0.2e; tone power = %0.2e (%0.1f%%); tone RMS = %0.2e', ...
            refTone.total_variance, refTone.total_power, refTone.tone_power, refTone.tone_percent_power, refTone.tone_rms))
  
        refTone.level = 94; % assuming reference tone is 94dB
        refTone.pressure_pascals = 10.^((refTone.level-94)/20); % 94dB is 1 Pascal
        refTone.rms_volts_recorded_per_pascal = refTone.tone_rms/refTone.pressure_pascals;
        calib.refTone = refTone;
        calibs{samplerate_idx} = calib;
        fprintf('-> Done\n');
    end
end

%% 2. relative calibration for each channel in turn
fprintf('\n== Relative calibration: measuring frequency response:\n');
for channel = 1:n_channels
    channelName = channelNames{channel};
    fprintf('Plug in %s...\n', channelName);
    pause;
    
    for samplerate_idx = 1:n_sampleRates
        calib = calibs{samplerate_idx};
        calib.relCalib = {};
        calib.relCalib{channel} = struct;
        fprintf('= %0.0fHz: ', calib.sampleRate);
        if isfield(calib, 'compensationFilters')
            fprintf('already done, skipping\n');
        else
            relCalibs = {};
            calib = relative_calib(calib.sampleRate, zBusNum, deviceName, golay_rms, channelName(1), cutoffs, dirname);
            fprintf('-> Done');
        end
    end
end



