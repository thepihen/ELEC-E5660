% This file contains the (failed) attempt to implement Serra's sines +
% noise model I talk about in the end of the report. I am including this
% because a non-negligible number of hours went into this

clc
clearvars
close all
%%
[x,Fs] = audioread("<YOUR_AUDIO_FILE_HERE>");%i was using copyrighted vocals
%that i can't upload...
%% extract frame from audio signal
x1 = x(15*Fs:30*Fs,1); %consider a mono signal for simplicity
L = length(x1);
h = 1024;
w = 2048;
win = hann(w,'periodic');
%% extract frames
n_frames =  ceil((length(x1) - w)/h) + 1;
frames = zeros(w, n_frames);

for i=1:n_frames
    if i==n_frames
        %zero pad
        frame = x1((i-1)*h+1:end);
        frame = [frame; zeros(w - length(frame),1)];
        frames(:,i) = frame.*win;
    else
        frames(:,i) = x1((i-1)*h+1:(i-1)*h+w).*win;
    end
end

%% take their FFT
frames_uncut_f = zeros(size(frames));
for i=1:n_frames
    frames_uncut_f(:,i) = (fft(frames(:,i)));
end
frames_f = abs(frames_uncut_f(1:w/2, :));
%%

% Compute the spectrogram
[S, F, T] = spectrogram(x1, hann(w), h, w, Fs);

% S: The spectrogram matrix
% F: The frequency values
% T: The time instants

% Plot the spectrogram
figure;
imagesc(T, F, 10 * log10(abs(S))); % Use log-scale for better visualization
axis xy; % Set the y-axis to be frequency in Hz
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram');
colorbar; % Add a color bar to show intensity

%% plot it
figure();
t = 0:h/Fs:(n_frames*h - h)/Fs;
f = 1:2048;
ff = 1:1024;
f_rescale_factor = Fs/(w);
t_rescale_factor = h/Fs;
f_plot = ff*f_rescale_factor;
axis xy; % Set the y-axis to be frequency in Hz
xlabel('Time (s)');
ylabel('Frequency (Hz)');
imagesc(t,f_plot,10 * log10(frames_f));
colorbar;
set(gca,'YDir','normal')

%% peak detection
peaks = zeros(size(frames_f));
%this will basically be a binary mask
%%
for i=1:n_frames 
    peaks(:,i) = peak_detect(frames_f(:,i), 10, 10, 0.7, 10);
end
%%
figure();
plot(f_plot,frames_f(:,50));
hold on;
plot(f_plot,frames_f(:,50).*peaks(:,50), "*", "Color","red");
%% remove some outliers and non-useful peaks
%deltas = frames_uncut_f.*peaks;
deltas = frames_f.*peaks;
for i = 1:n_frames
    deltas(size(deltas,1)/2:end,i) = 0;
end
%TODO: add more outlier removal techniques
%% ifft
frames_rec = zeros(size(frames));
for i=1:n_frames
    fvec = deltas(:,i);%frames_uncut_f(:,i);
    fvec = vertcat(fvec, 0, flipud(fvec(2:end)));
    frames_rec(:,i) = ifft(fvec);
end
%% overlap and add
x_out = zeros((h+2)*n_frames,1);
for i=1:n_frames
    x_out((h*(i - 1))+1:(h*(i - 1))+w) = x_out((h*(i - 1))+1:(h*(i - 1))+w)+ frames_rec(:,i).*win;
end
%% redivide
frames_resyn = zeros(w, n_frames);
for i=1:n_frames
    if i==n_frames && length(x_out((i-1)*h+1:end))<w
        %zero pad
        frame = x_out((i-1)*h+1:end);
        frame = [frame; zeros(w - length(frame),1)];
        frames_resyn(:,i) = frame.*win;
    else
        frames_resyn(:,i) = x_out((i-1)*h+1:(i-1)*h+w).*win;
    end
end

%%
frames_resyn_f = zeros(size(frames_resyn));
for i=1:n_frames
    frames_resyn_f(:,i) = abs(fft(frames_resyn(:,i)));
end
%%
err_final = abs(frames_uncut_f) - frames_resyn_f;
%% overlap and add
err_out = zeros((h+2)*n_frames,1);
for i=1:n_frames
    err_out((h*(i - 1))+1:(h*(i - 1))+w) = err_out((h*(i - 1))+1:(h*(i - 1))+w)+ ifft(err_final(:,i)).*win;
end
%% listen to the error
pl_err=audioplayer(err_out, Fs);
play(pl_err)

%% take their FFT
frames_uncut_f = zeros(size(frames));
for i=1:n_frames
    frames_uncut_f(:,i) = (fft(frames(:,i)));
end
frames_f = abs(frames_uncut_f(1:w/2, :));

%% ifft
frames_rec_f = zeros(size(frames));
for i=1:n_frames
    fvec = deltas(:,i);%frames_uncut_f(:,i);
    fvec = vertcat(fvec, 0, flipud(fvec(2:end)));
    frames_rec_f(:,i) = fvec;
end
frames_err = (abs(frames_uncut_f) - frames_rec_f);
%% overlap and add
err_out = zeros((h+2)*n_frames,1);
for i=1:n_frames
    err_out((h*(i - 1))+1:(h*(i - 1))+w) = err_out((h*(i - 1))+1:(h*(i - 1))+w)+ ifft(frames_err(:,i)).*win;
end
%% listen to the error
pl_err=audioplayer(err_out, Fs);
play(pl_err)
%% overlap and add
x_out = zeros((h+2)*n_frames,1);
for i=1:n_frames
    x_out((h*(i - 1))+1:(h*(i - 1))+w) = x_out((h*(i - 1))+1:(h*(i - 1))+w)+ frames_rec(:,i).*win;
end
%%

% Compute the spectrogram
[S2, F2, T2] = spectrogram(x_out, hann(w), h, w, Fs);

% S: The spectrogram matrix
% F: The frequency values
% T: The time instants

% Plot the spectrogram
figure;
imagesc(T2, F2, 10 * log10(abs(S2))); % Use log-scale for better visualization
axis xy; % Set the y-axis to be frequency in Hz
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram');
colorbar; % Add a color bar to show intensity

%%
pl=audioplayer(x_out, Fs);
play(pl)
%% compute error
err = x1 - x_out(1:length(x1));
%% listen to the error
pl_err=audioplayer(err, Fs);
play(pl_err)

%%
function pp = peak_detect(x, pre, post, delta, wait)
%x: input signal
%pre: number of samples to take before the current one
%post: number of samples to take after the current one
%delta: threshold
%wait: after it finds a peak, how many samples should the algorithm wait
%before looking for another one

    %if x is not mono make it
    x = x(:,1);
    n = length(x);
    pp = zeros(n,1);
    waiting = 0;
    
    for i=1:n
        vec = zeros(pre+post+1,1);
        if(waiting>0)
            waiting = waiting - 1;
        else
            if(i+post>n)
                diff = i+post-n;
                vec(1:length(x)-(i-pre)+1) = x(i-pre:end); 
            elseif(i-pre<1)
                diff = 1 - (i-pre);
                vec(diff+1:end) = x(1:i+post);
            else
                vec = x(i-pre:i+post);
            end
            avg = mean(vec);
            [maxv, maxi] = max(vec);
            sampl = vec(1+pre);
            if(maxi==1+pre)
                if(sampl>(1+delta)*avg)
                    pp(i) = 1;
                    waiting = wait;
                end
            end
        end
    end
end