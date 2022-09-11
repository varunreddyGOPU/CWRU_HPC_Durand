function fx = butt_band(x,FLOW,FHIGH, FS)

FCUT_HIGH = FHIGH/(FS/2);
FCUT_LOW = FLOW/(FS/2);
[bhigh,ahigh] = butter(7,FCUT_LOW, 'high');
[blow,alow] = butter(7,FCUT_HIGH, 'low');

fx = filtfilt(blow, alow, filtfilt(bhigh,ahigh, x'));
fx = fx';
end