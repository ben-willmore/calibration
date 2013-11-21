function irf=play_and_analyse_golay(samplerate, zbusnum, devicename, golaylen, trim, gap, firFlt, silent)
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
trimtaps=0;

% connnect RP
global RP;
if isempty(RP),
  RP=actxcontrol('RPco.x',[5 5 26 26]);
  if invoke(RP,['Connect' devicename],'GB',1) == 0
    if invoke(RP,['Connect' devicename],'USB',1) == 0
      error(['Cannot connect to ' devicename ' on GB 1 or USB 1');
    end;
  end;
  if invoke(RP,'LoadCOF',['playSig' devicename '.rcx']) == 0
    error(['Cannot load playSig' devicename '.rcx']);
  end;
  if invoke(RP,'Run') == 0
    error('RCOx Circuit failed to run.');
  end;
end;
irf.ADrate=double(invoke(RP,'GetSFreq'));

keyboard


% parse input
if exist('trim')
  trimtaps=ceil((trim/1000)*irf.ADrate);
end;
delaytaps=0;
if exist('gap')
  delaytaps=ceil((gap/1000)*irf.ADrate);
end;


if ~(irf.ADrate==48828.125)
  error('TDT:error',['wrong sample rate (' n2s(irf.ADrate) ' instead of 48828.125)']);
end

%% build signal to play from golay pair and delays
% ===================================================

outbuf=[ga zeros(1,trimtaps) zeros(1,delaytaps) gb zeros(1,trimtaps)];
if exist('firFlt')
  outbuf=conv(firFlt,outbuf);
end;

timestep=1/irf.ADrate;
taxis=timestep:timestep:timestep*length(outbuf);
t_last=taxis(length(outbuf));

%% record
% ==========

pause(0.5);
inbuf=RX6_play(outbuf*.1)';

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

