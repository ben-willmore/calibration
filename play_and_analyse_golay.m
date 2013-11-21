function irf=play_and_analyse_golay(samplerate, zbusnum, devicename, golaylen, trim, gap, rmsvalue, firFlt)
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

% prepare golay
[ga,gb]=golay(golaylen);
ga=[ga zeros(1,3*length(ga))];
gb=[gb zeros(1,3*length(gb))];
initrms = rms([ga gb]);
ga = ga/initrms*rmsvalue;
gb = gb/initrms*rmsvalue;
trimtaps=0;


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
  fprintf('Sample rate is %0.0f Hz\n',samplerate); 
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


% parse input
if exist('trim')
  trimtaps=ceil((trim/1000)*irf.ADrate);
end;
delaytaps=0;
if exist('gap')
  delaytaps=ceil((gap/1000)*irf.ADrate);
end;


%% build signal to play from golay pair and delays
% ===================================================

outbuf=[ga zeros(1,trimtaps) zeros(1,delaytaps) gb zeros(1,trimtaps)];
if exist('firFlt') & ~isempty(firFlt)
  outbuf=conv(firFlt,outbuf);
end;

timestep=1/irf.ADrate;
taxis=timestep:timestep:timestep*length(outbuf);
t_last=taxis(length(outbuf));

%% record
% ==========

pause(0.5);
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
inbuf.chan1=invoke(RP,'ReadTagVex','RecOut',0,length(outbuf),'F32','F64',1);
%inbuf.chan2=invoke(RP,'ReadTagVex','RecOut2',0,length(outbuf),'F32','F64',1);

%% analyse channel 1
% ====================

suma=inbuf.chan1(trimtaps+1:trimtaps+length(ga));
bstart=length(ga)+2*trimtaps+delaytaps;
sumb=inbuf.chan1(bstart+1:bstart+length(gb));

irf.chan1=golanal(suma(:)',sumb(:)',ga(:)',gb(:)');
irf.input_buffer = inbuf;
irf.input_buffer.t_axis = taxis;

irf.spectrum.f_axis = 0:(irf.ADrate/2)/(length(irf.chan1)/2-1):irf.ADrate/2;
irf.spectrum.power  = 20*log10(abs(irf.chan1(1:length(irf.chan1)/2)));

