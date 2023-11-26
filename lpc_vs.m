clc
clearvars
close all
%%
[x, Fs] = audioread("aaa.wav");
w = 4096;
h = 2048;
win = hanning(w);
win_s = hanning(w);
pad_length = mod(length(x), w);
x = x(:,1);
a_lpc = zeros(100, 1);
x = vertcat(x, zeros(pad_length,1));
out = zeros(size(x));
n_frames = length(x)/h - 1;
%% go over the signal in frames
for i=1:n_frames
    frame = x((i-1)*h+1:(i-1)*h+w);
    %% 2) Perform ptch detection using auto-correlation method.
    %  Consider only frequencies between 60 Hz and 500 Hz
    
    
    % compute autocorrelation (hint: use two output arguments of xcorr() and
    % pay attention to MATLAB normalization)
    [r,rlags] = xcorr(frame, frame, 'coeff');


    % consider correlation only for positive lags, including 0
    rpos = r(rlags >=0);

    % find maximum peak within the accepted range
    Fmin = 60;
    Fmax = 500;
    lagmin = floor(Fs/Fmax);
    lagmax = min(ceil(Fs/Fmin), w-1);
    
    rc = rpos(lagmin+1:lagmax+1);
    %+1 needed to index coherently with matlab
    
    [rc_maxv,rc_maxi] = max(rc);
    %rc_maxi is indexed in matlab way (starting from 1)
    
    r_maxi = lagmin+rc_maxi-1;
    %rc_maxi-1 compensates for Matlab indexing, so that rc_maxi is 0 based indexed
    
    r_maxv = rc_maxv;
    
    ptch_lag = r_maxi;
    %ptch_lag = pitch_period;
    ptch = Fs/ptch_lag;
  
    p = 100;
    R = toeplitz(rpos(1:p)); %coefficients from 0 to p-1
    a = R\rpos(2:p+1); %rpos(2:p+1) : coefficients from 1 to p
    a_lpc = [a_lpc, a];
    % Alternatively, you can use the lpc() function, be aware of what is returned
    %a = lpc(frame,p);
    %a = -a(2:end);
    %a=a';
    fprintf('Frame %i  -  Pitch lag: %d samples - freq: %.2f Hz\n',i, ptch_lag,ptch);
    
    %% 4) Plot the prediction error and its magnitude spectrum
    e = filter([1; -a],1,frame);

    u = zeros(length(frame),1);
    u(1:ptch_lag:end) = 1;
    
    % normalize the energy: force the energy of the
    % impulse train to be equal to that of the residual
    u = u / std(u) * std(e); % energy normalization

    frame_synth = filter(1,[1; -a],u);
%     frame_synth = frame_synth/max(abs(frame_synth));
    frame_synth = frame_synth.*win_s;
    out((i-1)*h + 1:(i-1)*h + w, :) = out((i-1)*h + 1:(i-1)*h + w, :) + frame_synth;
end

%%
out = out./(max(abs(out)));
% sound(s,Fs);
% pause(1);
player = audioplayer(out,Fs);
play(player);
%% aaa
%now open aaa.wav 
[s, Fs] = audioread("aaa.wav");
s = s(:,1);
pad_length = mod(length(s), w);
s = vertcat(s, zeros(pad_length,1));
out = zeros(size(s));

%extract a frame from the middle
frame = s(floor(length(s)/2)+1:floor(length(s)/2)+w);
%compute autocorrelation
[r_full,rlags] = xcorr(frame, frame, 'coeff');
p = 24;
%consider correlation only for positive lags, including 0
r = r_full(rlags >=0);
R = toeplitz(r(1:p));
a = R\r(2:p+1);

u = zeros(length(frame),1);
u(1:ptch_lag:end) = 1;

out = zeros(30*w,1);
for i = 1:30
    frame_synth = filter(1,[1; -a],u);
    frame_synth = frame_synth.*win_s;
    out((i-1)*h + 1:(i-1)*h + w, :) = out((i-1)*h + 1:(i-1)*h + w, :) + frame_synth;
end


pl = audioplayer(out,Fs);
play(pl);


%%


% find maximum peak within the accepted range
Fmin = 60;
Fmax = 500;
lagmin = floor(Fs/Fmax);
lagmax = min(ceil(Fs/Fmin), w-1);
    
rc = rpos(lagmin+1:lagmax+1);
%+1 needed to index coherently with matlab
    
[rc_maxv,rc_maxi] = max(rc);
%rc_maxi is indexed in matlab way (starting from 1)
    
r_maxi = lagmin+rc_maxi-1;
%rc_maxi-1 compensates for Matlab indexing, so that rc_maxi is 0 based indexed
    
r_maxv = rc_maxv;
    
ptch_lag = r_maxi;
%ptch_lag = pitch_period;
ptch = Fs/ptch_lag;
ptch_lag = Fs/98;
p = 10;
R = toeplitz(rpos(1:p)); %coefficients from 0 to p-1
a = R\rpos(2:p+1); %rpos(2:p+1) : coefficients from 1 to p
a = a_lpc(:,40)
e = filter([1; -a],1,frame);
%%
u = zeros(length(frame),1);
u(1:0.5*ptch_lag:end) = 1;
    
    % normalize the energy: force the energy of the
    % impulse train to be equal to that of the residual
u = u / std(u) * std(e); % energy normalization
%%
out = zeros(30*w,1);
for i = 1:30
    frame_synth = filter(1,[1; -a],u);
    frame_synth = frame_synth.*win_s;
    out((i-1)*h + 1:(i-1)*h + w, :) = out((i-1)*h + 1:(i-1)*h + w, :) + frame_synth;
end

%%
pl = audioplayer(out,Fs) ;
play(pl);