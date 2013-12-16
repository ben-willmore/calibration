%% test calibration

%calibdir = input('Which calibration directory? ', 's');
calibdir = 'calib.2012.12.13';
load(sprintf('%s/compensation_filters.mat', calibdir));

%% needs to loop through sample rates and channels, but for now...
samplerate_idx = 1;
channel_idx = 1;

calib = calibs(samplerate_idx);
relCalib = calib.relCalibs(channel_idx);

fs = calib.sampleRate;

%% fm sweep
f_min = 200;
f_max = 20000;
level = 80; % desired level in dB

f_c = f_max/2;
fmsweep = fmmod(logspace(log10(f_min), log10(f_max), fs)/f_c-1, f_c, fs, f_c); % 
fmsweep = fmsweep/rms(fmsweep); % now 94dB (RMS=1)
fmsweep = fmsweep * 10.^((level-94)/20); % now level dB

play_and_analyse_sound(fs, 1, 'RX6', fmsweep, calib.recorded_rms_volts_per_pascal, relCalib.filter, 150);

%%
freqs = logspace(log10(500),log10(32000), 14);
freqs = freqs(1:end-1);
level = 90; % dB

tones = zeros(1, 2*ceil(fs));
gap = floor(length(tones)/length(freqs));
tonedur = floor(0.5*gap);
ramplen = floor(5/1000*fs);
ramp = (1:ramplen)/ramplen;
env = [ramp ones(1, tonedur-2*ramplen) fliplr(ramp)];
t_tone = ((1:tonedur)-1)/fs;
    
for freq_idx = 1:length(freqs)
    freq = freqs(freq_idx);
    offset = (freq_idx-1)*gap;
    carrier = sin(2*pi*freq*t_tone);
    tones(offset+1:offset+tonedur) = env.*carrier;
end

% set RMS=1; 94dB
tones = tones*sqrt(2);    

% set desired level
tones = tones * 10^((level-94)/20); % now level dB

play_and_analyse_sound(fs, 1, 'RX6', tones, calib.recorded_rms_volts_per_pascal, relCalib.filter, 150);
subplot(2,2,2);
ylim([30 100]);

%% band-limited noise
freqs = logspace(log10(500),log10(16000), 14);
freqs = freqs(1:end-1);
level = 90; % dB

blnoise = zeros(1, 2*ceil(fs));
gap = floor(length(blnoise)/length(freqs));
tonedur = floor(0.5*gap);
ramplen = floor(5/1000*fs);
ramp = (1:ramplen)/ramplen;
env = [ramp ones(1, tonedur-2*ramplen) fliplr(ramp)];
t_burst = ((1:tonedur)-1)/fs;
    
for freq_idx = 1:length(freqs)
    freq = freqs(freq_idx);
    offset = (freq_idx-1)*gap;
    carrier = randn(1, tonedur);
    [N, Wn] = buttord([freq freq*1.5]/(fs/2), [freq*.8 freq*1.5/.8]/(fs/2), 2, 20);
    [B, A] = butter(N, Wn);
    carrier = filtfilt(B, A, carrier);
    carrier = carrier/rms(carrier);
    blnoise(offset+1:offset+tonedur) = env.*carrier;
end

% RMS = 1 already ; 94 dB 
% set desired level
blnoise = blnoise * 10.^((level-94)/20); % now level dB

play_and_analyse_sound(fs, 1, 'RX6', blnoise, calib.recorded_rms_volts_per_pascal, relCalib.filter, 150);
subplot(2,2,2);
ylim([30 100]);

%% wide(r) band noise
f_min = 500;
f_max = 600;
level = 90; % dB
freqs = f_min;

blnoise = zeros(1, 2*ceil(fs));
gap = floor(length(blnoise)/length(freqs));
tonedur = floor(0.5*gap);
ramplen = floor(5/1000*fs);
ramp = (1:ramplen)/ramplen;
env = [ramp ones(1, tonedur-2*ramplen) fliplr(ramp)];
t_burst = ((1:tonedur)-1)/fs;
    
for freq_idx = 1:length(freqs)
    freq = freqs(freq_idx);
    offset = (freq_idx-1)*gap;
    carrier = randn(1, tonedur);
    [N, Wn] = buttord([f_min f_max]/(fs/2), [f_min*.8 f_max/.8]/(fs/2), 2, 20);
    [B, A] = butter(N, Wn);
    carrier = filtfilt(B, A, carrier);
    carrier = carrier/rms(carrier);
    blnoise(offset+1:offset+tonedur) = env.*carrier;
end

% RMS = 1 already ; 94 dB 
% set desired level
blnoise = blnoise * 10.^((level-94)/20); % now level dB

play_and_analyse_sound(fs, 1, 'RX6', blnoise, calib.recorded_rms_volts_per_pascal, relCalib.filter, 150);
subplot(2,2,2);
ylim([30 100]);


%% click
click = zeros(1, round(0.1*fs));
click(100) = 10;
play_and_analyse_sound(fs, 1, 'RX6', click, calib.recorded_rms_volts_per_pascal, [], 150);
