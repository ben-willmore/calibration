function p = plot_spectrum(calib, smoothed, rms_volts_per_pascal)

%keyboard

  meanpower = 0*calib(1).spectrum.power;

  p = [];
  for ii=1:length(calib)
    p(ii) = plot(calib(ii).spectrum.f_axis, calib(ii).spectrum.power);
    meanpower = meanpower + calib(ii).spectrum.power/length(calib);
  end
    
  p(ii+1) = plot( calib(1).spectrum.f_axis, meanpower, 'y', 'linewidth', 4);
  p(ii+2) = plot( calib(1).spectrum.f_axis, meanpower, 'k', 'linewidth', 3);
  
% aesthetics
  set(gca,'fontsize',10,'fontweight','bold');
  set(gca,'xlim',[0,calib(1).ADrate/2]);
  %ylim([-80 20]);
  set(gca,'xscale','log');
  ylabel('dB','fontsize',12,'fontweight','bold');
  xlabel('Hz','fontsize',12,'fontweight','bold');
  grid on;