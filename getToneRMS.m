function rms_volts = getToneRMS(fs, rawRecording, toneFreq)

rec = rawRecording - mean(rawRecording);
rec = rec(round(25/1000*fs):end-round(25/1000*fs)); % chop off ends which get messed up by filtering
t = ((1:length(rec))-1)/fs;
figure(1);
subplot(2,2,1);
plot(t,rec);
subplot(2,2,2);
n_plot = round(fs/100);
plot(t(1:300),rec(1:300));

subplot(2,2,[3 4]);
[Pxx, F] = pwelch(rec, [], [], [], fs);
diff = abs(F-toneFreq);
f = find(diff<(toneFreq/10));
dF = F(2)-F(1);
rms_volts = sqrt(trapz(Pxx(f))*dF);

loglog(F, sqrt(Pxx));
hold all
loglog(F(f), sqrt(Pxx(f)), 'ro-');
hold off;
