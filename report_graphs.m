% this file was used to create the graphs included in the report

clearvars
close all
clc
%% load data and setup
source = "aaa.wav";
load('aaa_lpc_coeffs.mat', 'a_lpc');
[x, Fs] = audioread(source);
lpc_coeff = a_lpc(:,45)';

n = 0:2048-1;
f0 = 139;
pitch = Fs/f0;
P = pitch;
M = 2*floor(P/2) + 1;
rat = Fs/f0;
%% build inputs

input_BLIT = (M/P)*sin(pi* M*n/P)./(M*sin(pi*n/P));
input_BLIT(mod(n,P)==0) = 1; 
input_BLIT(1) = 1;
input_BLIT(isnan(input_BLIT)) = 1;
input_BLIT = input_BLIT';
input_IMPT = zeros(size(input_BLIT));
input_IMPT(1:P:end) = 1;
%% plot inputs

t = (0:length(input_IMPT)-1)/Fs;
figure;
subplot(211)
plot(t,input_BLIT)
grid on
xlabel("Time [s]")
title("BLIT")

subplot(212)
plot(t, input_IMPT)
grid on
xlabel("Time [s]")
title("Impulse Train")

sgtitle(['BLIT and Impulse Train for f= ',num2str(f0),' Hz, F_s = 44100 Hz']) 
%% compute the output of the LPC shaping filter 
block_BLIT = filter(1,[1; -lpc_coeff'], input_BLIT);
block_IMPT = filter(1,[1; -lpc_coeff'], input_IMPT);
%% plot results
figure;
subplot(311)
plot(t,block_BLIT)
grid on
xlabel("Time [s]")
title("BLIT")

subplot(312)
plot(t, block_IMPT)
grid on
xlabel("Time [s]")
title("Impulse Train")

subplot(313)
plot(t, block_IMPT - block_BLIT)
grid on
xlabel("Time [s]")
title("block_IMPT - block_BLIT")


sgtitle(['LPC output for BLIT and Impulse Train inputs for f= ',num2str(f0),' Hz, F_s = 44100 Hz']) 

%% plot results in frequency
block_BLIT_padded = [block_BLIT; zeros(44100-2048, 1)];
block_IMPT_padded = [block_IMPT; zeros(44100-2048, 1)];
fft_BLIT  = abs((fft(block_BLIT_padded)));
fft_IMPT = abs((fft(block_IMPT_padded)));

figure;
subplot(311)
plot(fft_BLIT(1:floor(length(fft_BLIT)/2)))
grid on
xlabel("Time [s]")
title("BLIT")

subplot(312)
plot(fft_IMPT(1:floor(length(fft_BLIT)/2)))
grid on
xlabel("Time [s]")
title("Impulse Train")

subplot(313)
plot(fft_BLIT(1:floor(length(fft_BLIT)/2)) - fft_IMPT(1:floor(length(fft_BLIT)/2)))
grid on
xlabel("Time [s]")
title("Difference")

sgtitle(['FFT of LPC output for BLIT and Impulse Train inputs for f= ',num2str(f0),' Hz, F_s = 44100 Hz']) 

