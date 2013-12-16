function p = plot_input_buffer(calib)
  
  p = [];
  for ii=1:length(calib)
    p(ii) = plot(calib(ii).input_buffer.t_axis*1000, calib(ii).input_buffer.chan1);
  end
  
  xlim([0 max(calib(1).input_buffer.t_axis)]*1000);
  xlabel('time (ms)','fontsize',12,'fontweight','bold');
  ylabel('Volt','fontsize',12,'fontweight','bold');

end