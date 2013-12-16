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
f_min = 200;
f_max = 20000;
level = 80; % desired level in dB

f_c = f_max/2;
fmsweep = fmmod(logspace(log10(f_min), log10(f_max), fs)/f_c-1, f_c, fs, f_c); % 
fmsweep = fmsweep/rms(fmsweep); % now 94dB (RMS=1)


fmsweep = fmsweep * 10.^((level-94)/20);
play_and_analyse_sound(fs, 1, 'RX6', fmsweep, calib.recorded_rms_volts_per_pascal, relCalib.filter, 150);
