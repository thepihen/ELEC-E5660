clearvars
close all
clc

%% load audio
[x, Fs] = audioread("COME NON MAI (Italy Demo).mp3");
%% setup device writer
deviceWriter = audioDeviceWriter('SampleRate',Fs, ...
    "Driver", "ASIO", "Device", "FL Studio ASIO","SupportVariableSizeInput", true,...
    "BitDepth", "24-bit integer");
% asiosettings(deviceWriter.Device);
deviceWriter.BufferSize = 512;

%%
N = 4096;
n_frames = ceil(length(x)/N);
for i = 1:n_frames
    frame = x((i-1) * N + 1:i*N, :);
    for j = 1:N/512
        sframe = frame(1:512, :);
        frame = circshift(frame, -length(sframe));
        deviceWriter(sframe);
    end
end