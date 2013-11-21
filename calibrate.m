%% new calibration code BW Nov 2013
load fakeCalib;

tdt50k = 48828.125;
zBusNum = 1;
deviceName = 'RX6';

calib = struct;
calib.sampleRate = tdt50k;

%% 1. get the level of the reference tone
if ~isfield(calib, 'reference_rms')
  fprintf('Plug in reference sound and press a key\n');
  pause;
  % record tone and calculate RMS
  %calib.refTone.irf = play_and_analyse_golay(calib.sampleRate, zBusNum, deviceName, 11, 0, 5, 0, []);
  calib.refTone = fakeCalib.refTone;
  rec = calib.refTone.irf.input_buffer.chan1 - mean(calib.refTone.irf.input_buffer.chan1);
  figure(1);
  subplot(2,1,1);
  plot(rec);
  hold all;
  [B,A] = butter(3,[900 1100]/(calib.sampleRate/2));
  flt = filtfilt(B,A,rec);
  plot(flt);
  hold off;
  subplot(2,1,2);
  [Pxx, F] = pwelch(rec, [], [], [], calib.sampleRate);
  diff = abs(F-1000);
  f = find(diff<100);
  val = sum(sqrt(Pxx(f)))
  %plot(F,Pxx)
  loglog(F, sqrt(Pxx))
  hold all
  loglog(F(f), sqrt(Pxx(f)), 'ro-');
  hold off;
  
  
  ft = fft(calib.refTone.irf.input_buffer.chan1); % take first half of FFT only
  amp_spec = ft(1:(length(ft+1)/2));
  freqs = linspace(0,lengthcalib.sampleRate/2
  pwelch(calib.refTone.irf.input_buffer.chan1,5000);
  
end