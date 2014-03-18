%% new calibration code BW Nov 2013
tdt50k = 48828.125;

%% parameters that can be changed
sampleRates = [tdt50k tdt50k*2]; % sample rates at which to calibrate
cutoffs{1} = [200 tdt50k/2]; % low and high freq cutoffs for sampleRates(1)
cutoffs{2} = [200 32000]; % low and high freq cutoffs for sampleRates(2)

zBusNum = 1; % zBus number (usually 1)
deviceName = 'RX6'; % device to use (there must be a corresponding .rcx file)
channelNames = {'Left-hand earphone (BLUE)'; 'Right-hand earphone (RED)'}; % you must have one name per channel
golay_rms = 0.05*10.^(-12/20)*10; % RMS voltage to play golay codes at -- adjust to a reasonable level
recording_highpass_f = 150; % Hz; filter out frequencies below this to deal with low-frequency noise

%% don't change below this
addpath('./general');
addpath('./jan');

%% initial setup
fprintf('\n\n\n');
fprintf('======================================\n');
subdir = input('Subdirectory to save in: \n   ','s');
dirname = fixpath(['./' subdir]);
fprintf('======================================\n\n');
mkdir_nowarning(dirname);


n_sampleRates = length(sampleRates);
n_channels = length(channelNames);

if ~exist('calibs', 'var')
    calib = struct;
    for samplerate_idx = 1:n_sampleRates
        calib.sampleRate = sampleRates(samplerate_idx);
        calib.cutoffs = cutoffs{samplerate_idx};
        calibs{samplerate_idx} = calib;
    end
end

%% 1. get the level of the reference tone
fprintf('== Reference tone: measuring mic response: \n');
for samplerate_idx = 1:n_sampleRates
    calib = calibs{samplerate_idx};
    fprintf('= %0.0fHz: ', calib.sampleRate);
    if isfield(calib, 'refTone')
        fprintf('already done, skipping\n');
    else
        happy = false;
        while ~happy
            if exist('fakeData', 'var') && fakeData
              fprintf('Using fake data\n');
              f = load('fakeCalib.mat');
              calib.refTone = f.fakeCalib.refTone;
            else
              fprintf('Plug in reference sound set to 94dB and press a key...\n');
              pause;
              % record tone and calculate RMS
              % (this plays silence)
              calib.refTone.irf = play_and_analyse_golay(calib.sampleRate, zBusNum, deviceName, 11, 0, 5, 0, [], recording_highpass_f);
            end

            % get RMS of signal at 1000Hz
            % assuming ref tone is at 94dB (RMS=1Pa), then this gives us
            % the volts/pascal of the mic
            recorded_rms_volts_per_pascal = getToneRMS(calib.sampleRate, calib.refTone.irf.input_buffer.chan1, 1000);

            s = '';
            while ~strcmpi(s, 'y') && ~strcmpi(s, 'n')
              s = input('  - happy [y/n]', 's');
              if strcmpi(s, 'y')
                  happy = true;
              end
            end
        end
        calib.refTone.recorded_rms_volts_per_pascal = recorded_rms_volts_per_pascal;
        calibs{samplerate_idx} = calib;
        fprintf('-> Done\n');
    end
end

%% 2. calibration for each channel in turn
fprintf('\n== Calibration:\n');
for channel = 1:n_channels
    channelName = channelNames{channel};
    fprintf('%s:\n', channelName);
    plugged_in = false;

    for samplerate_idx = 1:n_sampleRates
        calib = calibs{samplerate_idx};
                   
        if ~isfield(calib, 'channel')
            calib.channel = {};
        end
        
        if length(calib.channel)>=channel
            fprintf('already done, skipping\n');
            continue;
        end
        if ~plugged_in
            fprintf('Plug in %s and press a key...\n', channelName);
            pause;
            plugged_in = true;
        end
        fprintf('= %0.0fHz: ', calib.sampleRate);
        calib.channel{channel} = relative_calib(calib.sampleRate, zBusNum, deviceName, golay_rms, channelName(1), ...
            cutoffs{samplerate_idx}, calibs{samplerate_idx}.refTone.recorded_rms_volts_per_pascal, recording_highpass_f, dirname);
        calib.channel{channel}.channelIdx = channel;
        calib.channel{channel}.channelName = channelName;
        % now play 1000Hz compensated, record RMS etc.
        calib.channel{channel}.absCalib = struct;
        t = ((1:calib.sampleRate)-1)/calib.sampleRate;
        snd = sin(2*pi*1000*t);
        snd = snd/rms(snd); % 1Pa RMS = 94dB
        
        fudgeGain = 30;
        calib.channel{channel}.absCalib.irf = play_and_analyse_sound(calib.sampleRate, zBusNum, deviceName, snd, calib.channel{channel}.filter/fudgeGain, recording_highpass_f, calib.refTone.recorded_rms_volts_per_pascal);
        recorded_rms_volts_per_pascal = getToneRMS(calib.sampleRate, calib.channel{channel}.absCalib.irf.input_buffer.chan1, 1000);
        filter_correction_factor = calib.refTone.recorded_rms_volts_per_pascal/fudgeGain/recorded_rms_volts_per_pascal;
        calib.channel{channel}.filter = calib.channel{channel}.filter * filter_correction_factor;
        
%         for freq = [500 1000 5000 10000]
%             snd = sin(2*pi*freq*t);
%             snd = snd/rms(snd); % 1Pa RMS = 94dB
%             calib.channel{channel}.absCalib.irf2 = play_and_analyse_sound(calib.sampleRate, zBusNum, deviceName, snd, calib.channel{channel}.filter, recording_highpass_f);
%             recorded_rms_volts_per_pascal = getToneRMS(calib.sampleRate, calib.channel{channel}.absCalib.irf2.input_buffer.chan1, freq);
%             recorded_pa = recorded_rms_volts_per_pascal/calib.refTone.recorded_rms_volts_per_pascal;
%         end
        %recorded_rms_volts_per_pascal = getToneRMS(calib.channel{channel}.sampleRate, calib.channel{channel}.refTone.irf.input_buffer.chan1, 1000);
        
        calibs{samplerate_idx} = calib;
        fprintf('-> Done\n');
    end
end

%% 3. save
save(sprintf('%s/calibration.mat', dirname), 'calibs');


%% 4. save minimal information (inc filters) in a second file
full_calibs = calibs;
calibs = {};
for samplerate_idx = 1:length(full_calibs)
    full_calib = full_calibs{samplerate_idx};
    calib = struct;
    calib.sampleRate = full_calib.sampleRate;
    calib.reftone_rms_volts_per_pascal = full_calib.refTone.recorded_rms_volts_per_pascal;
    calib.cutoffs = full_calib.cutoffs;
    channels = {};
    for channel_idx = 1:length(full_calib.channel)
        full_channel = full_calib.channel{channel_idx};
        channel = struct;
        channel.channelIdx = full_channel.channelIdx;
        channel.channelName = full_channel.channelName;
        channel.ADrate = full_channel.ADrate;
        channel.filter = full_channel.filter;
        channels{channel_idx} = channel;
    end
    calib.channels = [channels{:}];
    calibs{samplerate_idx} = calib;
end
calibs = [calibs{:}];

save(sprintf('%s/compensation_filters.mat', dirname), 'calibs');
