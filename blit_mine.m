clearvars
close all
clc

source = "aaa.wav";
load('aaa_lpc_coeffs.mat', 'a_lpc');
[x, Fs] = audioread(source);
lpc_coeff = a_lpc(1,:)';
input = zeros(length(x),1);
% n = 0:2048-1;
n = 0:2048-1;
f0 = 110;
Fs = 44100;
pitch = Fs/f0;
P = pitch;
M = 2*floor(P/2) + 1;
w0T = 2*pi*f0/Fs;
input(1:pitch:end) = 1;


output = 0;
for i=1:5 
    if(mod(i,2)~=0 && i == 1)
        f0 = f0 * sqrt(2);
    elseif (mod(i,2)~=0 && i>1)
        f0 = f0 * 2^(5/12);
    else
        f0 = f0 * 2^(7/12);
    end
    pitch = Fs/f0;
    P = pitch;
    M = 2*floor(P/2) + 1;
    y = (M/P)*sincM((M/P)*n,M);
    y = y';
    %This is not enough yet. If you feed this as input it will not work as
    %there are still many peaks due to zeros in the denominator.
    %We will set all those values to one
    y(mod(n,P)==0) = 1; %this is a certified cantinaro moment though
    y(1) = 1;
    y(isnan(y)) = 1;
    amp=0.5;
    block = filter(1,[1; -lpc_coeff], y);
    block = block/max(abs(block)) * amp;
    block = block.*hamming(length(block));
    block_final = zeros(length(block)*21,1);
    w = length(block);
    h = w/2;
    
    for i=1:40
        block_final(1+((i-1)*h):w+((i-1)*h)) = block_final(1+((i-1)*h):w+((i-1)*h))+block;
    end
    % block = repmat(block, 30,1);
    output = [output;block_final];
end
output = [output, output];
audiowrite("AAAs-demo.wav",output,Fs);



%%
a = sin((M/P)*n);
b = (M/P)*sin(w0T*pi*n)./(M*sin(w0T*pi*n/M));
b = (M/P)*sincM(w0T*n, M);
b(isnan(b)) = 1;
figure;
plot(b)
title("b")
%% build a train of pulses
%https://www.music.mcgill.ca/~gary/307/week5/node14.html
%https://ccrma.stanford.edu/~jos/SpecEnv/LPC_Envelope_Example_Speech.html
%^^these are the inspirations
% y = (M/P)*sin(pi* M*n/P)./(M*sin(pi*n/P));
y = (M/P)*sincM((M/P)*n,M);
y = y';
%This is not enough yet. If you feed this as input it will not work as
%there are still many peaks due to zeros in the denominator.
%We will set all those values to one
y(mod(n,P)==0) = 1; %this is a certified cantinaro moment though
y(1) = 1;
y(isnan(y)) = 1;
plot(y)
title("y")
%%
input = (M/P)*sincM((M/P)*n, M);
input(isnan(input))=1;
IN = fftshift(fft(input));
IN = abs(IN);
figure;
plot(input);
title("Input")
figure;
plot(abs(fft(input)))
title("FFT input")
%%
amp=0.5;
%%

block_in = filter(1,[1; -lpc_coeff], input);
block_in = block_in/max(abs(block_in)) * amp;
%%
block = filter(1,[1; -lpc_coeff], y);
block = block/max(abs(block)) * amp;
%%
figure;
plot(abs(fft(block_in)));
hold on
plot(abs(fft(block)))
figure
plot(abs(fft(block_in)) - abs(fft(block)))

%%
block = repmat(block, 50,1);
soundsc(block,Fs)
% input(1:pitch:end) = 1;
% block = filter(1,[1; -a], input);
% block = block/max(abs(block)) * amp;


function y = sincM(x, M)
    y = sin(pi*x)./(M*sin(pi*x/M));
end