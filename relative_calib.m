function calib = relative_calib(sampleRate, zBusNum, deviceName, golay_rms, channelID, cutoffs, dirname)
%% record calibration curve
% ==================================================

cont=1;
while cont==1
  
  % record silence
  fprintf('  - recording silence...');
  calib.silence = cell(1,10);
  for ii=1:10
    calib.silence{ii} = play_and_analyse_golay(sampleRate, zBusNum, deviceName, 11, 0, 5, 0, []);
  end
  calib.silence = [calib.silence{:}];
  fprintf(' [done]\n');
  
  % record impulse
  fprintf('  - recording impulse response...');
  calib.irf = cell(1,10);
  for ii=1:10
    calib.irf{ii} = play_and_analyse_golay(sampleRate, zBusNum, deviceName, 11, 0, 5, golay_rms, []);
  end
  calib.irf = [calib.irf{:}];
  fprintf(' [done]\n');
  
  % plot
  figure;
  clf;
  sp1 = subplot(2,2,1);
  hold all;
  p1 = plot_input_buffer(calib.silence);
  set(p1,'color',[1 0 0]);
  title(sprintf('%s: silence','fontsize', channelID),14,'fontweight','bold');
  sp2 = subplot(2,2,2); hold all;
  p2 = plot_input_buffer(calib.irf);
  set(p2,'color',[0 0 1]);
  title(sprintf('%s: impulse response', channelID), 'fontsize',14,'fontweight','bold');
  yl = max(abs([get(sp1,'ylim') get(sp2,'ylim')]))*[-1 1];
  set([sp1 sp2], 'ylim', yl);
  sp3 = subplot(2,1,2); hold on;
  p3a = plot_spectrum(calib.silence);
  set(p3a(1:(end-2)),'color',[1 0 0]);
  p3b = plot_spectrum(calib.irf);
  set(p3b(1:(end-2)),'color',[0 0 1]);
  title(sprintf('%s: fft', channelID), 'fontsize',14,'fontweight','bold');
  set(gca,'xtick',[1e2 1e3 1e4 4e4],'xticklabel',{'100','1k','10k','40k'});
  
  % are we happy
  s=input('  - happy? [y/n]  ','s');
  if ~isempty(s) && (s=='y' || s=='Y')
    cont=0;
  end
end


%% save figure
% ===============

% save figure
fprintf('  - saving figure...');
savefigure(dirname,sprintf('calibration.fig.01 - %s at 50kHz', channelID));
fprintf(' [done]\n');

  
%% make compensation filters
% ============================

fprintf('  - constructing compensation filter...\n');

% put results into calib
calib.ADrate = calib.irf(1).ADrate;
calib.mean_power = mean(20*log10(abs(reach(calib.irf,'chan1''')')));    
calib.mean_intensity = 10 .^ (calib.mean_power/20);
  
% inverse filter
len.L = L(calib.mean_intensity);
norm1 = ones(1,len.L) ./ calib.mean_intensity;

% flatten above cutoffs(2)
fromBin = round(len.L*cutoffs(2)/calib.ADrate);
segLen = len.L/2 - fromBin;
norm1(fromBin+1:fromBin+2*segLen) = ones(1,2*segLen);

% flatten below cutoffs(1)
fromBinLowFreqs=round(len.L*cutoffs(1)/calib.ADrate);
norm1(1:fromBinLowFreqs) = ones(1,fromBinLowFreqs) *10^(0/20);
norm1((end-fromBinLowFreqs+2):end) = ones(1,fromBinLowFreqs-1) *10^(0/20);

% zero DC term
norm1(1)=0;

% construct filter
n=abs(jdecimate(jdecimate(norm1')))';
s=length(n)/2; % a slope s on the phases centers the fir filter
calib.filter = real(ifft( n .* exp(j*[-s*pi:pi:(s-1)*pi])));   %#ok<*IJCL>
calib.filter = calib.filter/norm(calib.filter); % length 1 -- no effect on amplitude
    
   
%% test calibration curve
% ==========================

cont=1;
while cont==1
  
  % play calibrated sound
  fprintf('  - recording compensated impulse response...');
  calib.irfc = cell(1,10);
  for ii=1:10
    calib.irfc{ii} = play_and_analyse_golay(sampleRate, zBusNum, deviceName, 11, 0, 5, golay_rms, []);
  end
  calib.irfc = [calib.irfc{:}];
  fprintf(' [done]\n');
  
  % plot
  figure(2); clf;
  sp1 = subplot(2,1,1); hold all;
  p1 = plot_input_buffer(calib.irfc);
  set(p1,'color',[0 0 1]);
  title(sprintf('%s: compensated impulse response', channelID),'fontsize',14,'fontweight','bold');
  yl = max(get(sp1,'ylim'))*[-1 1];
  %set(sp1, 'ylim', yl);
  sp2 = subplot(2,1,2); hold on;
  p2 = plot_spectrum(calib.irfc);
  set(p2(1:(end-2)),'color',[0 0 1]);
  plot([cutoffs(1) cutoffs(2)],[0 0],'r--','linewidth',2);
  title(sprintf('%s: fft', channelID), 'fontsize',14,'fontweight','bold');
  set(gca,'xtick',[1e2 1e3 1e4 4e4],'xticklabel',{'100','1k','10k','40k'});
  
  % happy?
  s=input('  - happy? [y/n]  ','s');
  if ~isempty(s) && (s=='y' || s=='Y')
    cont=0;
  end
end

% save figure
fprintf('  - saving figure...');
savefigure(dirname,sprintf('calibration.fig.02 - %s at 50kHz compensated', channelID));
fprintf(' [done]\n');
