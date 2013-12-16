function irf=play_and_analyse_sound(samplerate, zbusnum, devicename, snd, rms_voltage_per_pascal, compensationFilter, highpass)
%   *** adapted from o4golayrec ***
%
% irf = o4golayrec(golaylen, trim, gap);
%   records fourier spectra of impulse response functions to "irf" from 
%   golay codes of length golaylen presented once through tdt system 3.
%     - trim (in ms) specifies how much of the recorded signal delay to trim off
%     - gap (in ms) specifies how much of a time delay to leave between presentations 
%         of each code in the pair

%% prepare
% =========

%% parse sample rate
tdt50k = 48828.125;

% available sample rates
sampleRates =   [0.125 0.25 0.5 1 2 4 8]*tdt50k;
sampleRateIDs = [    0    1   2 3 4 5 6];

diff = abs(sampleRates-samplerate);
f = find(diff==min(diff), 1);

if length(f)==1
  samplerate = sampleRates(f);
  sampleRateID = sampleRateIDs(f);
  %fprintf('Sample rate is %0.0f Hz\n',samplerate); 
else
  error('Unknown sample rate');
end

%% connnect RP
global RP;
if ~isempty(RP)
  irf.ADrate=double(invoke(RP,'GetSFreq'));
  if round(irf.ADrate)~=round(samplerate)
    RP = [];
  end
end

if isempty(RP)
  RP=actxcontrol('RPco.x',[5 5 26 26]);
  if invoke(RP,['Connect' devicename],'GB',zbusnum) == 0
    error(['Cannot connect to ' devicename ' on GB 1 or USB 1']);
  end;
  if invoke(RP,'LoadCOFsf',['playSig_' devicename '.rcx'], sampleRateID) == 0
    error(['Cannot load playSig_' devicename '.rcx']);
  end;
  if invoke(RP,'Run') == 0
    error('RCOx Circuit failed to run.');
  end;
end;

% check sample rate
irf.ADrate=double(invoke(RP,'GetSFreq'));
if irf.ADrate==0 
    error('GetSFreq failed!');
end;

if round(irf.ADrate)~=round(samplerate)
  error(sprintf('Wrong sample rate on device (%0.0f vs %0.0f)', irf.ADrate, samplerate));
end


%% build signal to play
% ===================================================

padding = round(25/1000*samplerate);
outbuf=[zeros(1, padding) snd zeros(1, padding)];
outbuf = outbuf * rms_voltage_per_pascal;

if exist('compensationFilter') && ~isempty(compensationFilter)
  outbuf=conv(compensationFilter, outbuf);
end;

if max(abs(outbuf))>10
    error('too loud!!');
end

timestep=1/irf.ADrate;
taxis=timestep:timestep:timestep*length(outbuf);
t_last=taxis(length(outbuf));

%% record
% ==========

pause(0.1);
%inbuf=RX6_play(outbuf*.1)';

% find out sample conversion rate

% work out how long sample will be in milliseconds
sigDur=(length(outbuf)/irf.ADrate)*1000;
% pump signal to play into buffer
if invoke(RP,'WriteTagVex','SigData',0,'F32',outbuf) == 0,
  error('WriteTagVEX SigData failed');
end;
% reset record / playback buffer indeces
invoke(RP,'SoftTrg',2);

if invoke(RP,'SetTagVal','SigLen',sigDur+5) == 0,
  error('SetTagVal SigLen failed');
end;

% trigger DA conversion
invoke(RP,'SoftTrg',1);
pause(sigDur/1000);
inbuf.chan1_unfiltered=invoke(RP,'ReadTagVex','RecOut',0,length(outbuf),'F32','F64',1);
%inbuf.chan2=invoke(RP,'ReadTagVex','RecOut2',0,length(outbuf),'F32','F64',1);

t = ((1:length(snd))-1)/samplerate;

if exist('highpass', 'var') && highpass~=0 && isfinite(highpass)
   [N, Wn] = buttord(highpass/(samplerate/2), 0.6*highpass/(samplerate/2), 1, 20);
   [B, A] = butter(N, Wn, 'high');
   inbuf.chan1 = filtfilt(B, A, inbuf.chan1_unfiltered);
else
    inbuf.chan1 = inbuf.chan1_unfiltered;
end

if ~isempty(compensationFilter)
  inbuf.chan1 = inbuf.chan1(padding+length(compensationFilter)/2:end-padding-length(compensationFilter)/2);
else
  inbuf.chan1 = inbuf.chan1(padding+1:end-padding);
end

figure(9);
subplot(2,2,1);
plot(t, snd/rms_voltage_per_pascal);
subplot(2,2,3);
plot(t, inbuf.chan1);
subplot(2,2,2);
plot(t, 94+20*log10(movingstd(snd, round(10/1000*samplerate))));
subplot(2,2,2);
hold all;
plot(t, 94+20*log10(movingstd(inbuf.chan1, round(10/1000*samplerate))/rms_voltage_per_pascal));
hold off;



%% analyse channel 1
% ====================
%keyboard
%keyboard
