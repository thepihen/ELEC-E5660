%% main, testing version
% main version made for testing without having to re-select the audio file
% and the MIDI device each time

%% startup
clearvars
close all
clc

%% important variables
Fs = 44100;

%% load previously computed LPC coefficients
load('coeffs/aa_lpc_coeffs.mat');
load('coeffs/oo_lpc_coeffs.mat');
load('coeffs/uu_lpc_coeffs.mat');
% a small note on these LPC coefficients. They are order 100 LPC coefficients
% extracted from 7 recordings per vowel. The recordings cover C,E,G of 2
% octaves, plus the C from a third one
% this is, unfortunately, the vocal range I, a completely untrained singer,
% was able to (more or less) consistently hit notes on
%% play file
% player = audioplayer(x, Fs);
% play(player);

%% choose MIDI device (possibly make an intelligent choice)
userIn = "d";
while(userIn ~= "y" && userIn ~= "n")
    userIn = input("Do you wish to play from your external MIDI keyboard? (y/n)","s");
end
playFromMIDI = (userIn == "y");
mididevinfo;
selectedDevice = "KOMPLETE KONTROL M32 MIDI";

%% audioDeviceWriter object to play audio to sound card
%needed for lower latency for real-time audio
%buffer of 192 samples should guarantee a pretty low latency
% deviceWriter = audioDeviceWriter('SampleRate',44100, ...
%     "Driver", "ASIO", "Device", "Focusrite USB ASIO","SupportVariableSizeInput", true,...
%     "BitDepth", "24-bit integer");
%deviceWriter.BufferSize = 192;
deviceWriter = audioDeviceWriter('SampleRate',44100, ...
    "Driver", "ASIO", "Device", "FL Studio ASIO","SupportVariableSizeInput", true,...
    "BitDepth", "24-bit integer");
% asiosettings(deviceWriter.Device);
deviceWriter.BufferSize = 512;
if playFromMIDI
    simplesynth(selectedDevice, deviceWriter, aa_lpc_coeffs, Fs);
else
    simplesynthNoMIDI(deviceWriter, oo_lpc_coeffs, Fs);
end

function simplesynthNoMIDI(deviceWriter,a_lpc, Fs)
deviceWriter.SupportVariableSizeInput = true;
% a = a_lpc(:,40);
a = a_lpc(:,1);
amp = 0.6;
freq = 100;
prevFreq = 0;
remainder = zeros(deviceWriter.BufferSize,1);
N = deviceWriter.BufferSize;
N = 512;
h = deviceWriter.BufferSize/2;
notes = [59,61,64,70,76,80,85];
currNote = 1;
delayBetweenNotes = 1;
noteTic = tic;
hasFinishedProcessingPreviousBlock = true;
currBlock = 1;
while true

    freq = note2freq(notes(currNote));
    
    if(toc(noteTic)>delayBetweenNotes)
        currNote = currNote + 1;
        if currNote > 7
            currNote = 1;
        end
        noteTic = tic;
    end
    
    if(freq~=prevFreq)
        ll = 8192;
        input = zeros(ll,1);
        n = 0:ll-1;
        n = n';
        w0T = 2*pi*freq/Fs;
        pitch = (Fs/freq);
        P = pitch;
        M = 2*floor(P/2) + 1;
        input = (M/P)*sin(pi* M*n/P)./(M*sin(pi*n/P));
        %This is not enough yet. If you feed this as input it will not work as
        %there are still many peaks due to zeros in the denominator.
        %We will set all those values to one
        input(mod(n,P)==0) = 1; %this is a certified cantinaro moment though
        input(1) = 1;
        input(isnan(input)) = 1;

        block = filter(1,[1; -a], input);
        block = block/max(abs(block)) * amp;
        
        prevFreq = freq;
        currBlock = 1;
    end
    sblock = block(1:N,:);
    if(currBlock == 1)
        sblock2 = sblock.*hamming(length(sblock), 'periodic');
        sblock(1:length(sblock)/2, :) = sblock2(1:length(sblock2)/2, :);
    end
    if(mod(currBlock, ll/N)==0)
        block = circshift(block, -length(sblock));
        nextBlock = block(1:N,:);
        fakeBlock = [nextBlock; sblock].*hamming(size(sblock,1)*2);
        sblock = fakeBlock(length(fakeBlock)/2 + 1:end) + fakeBlock(1:length(fakeBlock)/2);
        currBlock = currBlock + 1; %you're skipping one
    end
    % if(currBlock = ll/N)
    %     sblock2 = sblock.*hamming(length(sblock), 'periodic');
    %     sblock(1:length(sblock)/2, :) = sblock2(1:length(sblock2)/2, :);
    % end
    block = circshift(block, -length(sblock));
    deviceWriter([sblock, sblock]);
    currBlock = currBlock + 1;
    % block = circshift(block, -N);
end
end


function yes = isNoteOn(msg)
yes = msg.Type == midimsgtype.NoteOn ...
    && msg.Velocity > 0;
end

function yes = isNoteOff(msg)
yes = msg.Type == midimsgtype.NoteOff ...
    || (msg.Type == midimsgtype.NoteOn && msg.Velocity == 0);
end

function freq = note2freq(note)
freqA = 440;
noteA = 69;
freq = freqA * 2.^((note-noteA)/12);
end

function y = sincM(x, M)
y = sin(pi*x)./(M*sin(pi*x/M));
end

function y = sinc(x)
y = sin(pi * x) ./ (pi * x);
y(x == 0) = 1;
end



%https://www.music.mcgill.ca/~gary/307/week5/node14.html
%https://ccrma.stanford.edu/~jos/SpecEnv/LPC_Envelope_Example_Speech.html
function simplesynth(midiDeviceName, deviceWriter,a_lpc, Fs)
midiInput = mididevice(midiDeviceName);
deviceWriter.SupportVariableSizeInput = true;
a = a_lpc(:,40);
amp = 0.6;
freq = 100;
prevFreq = 0;
remainder = zeros(deviceWriter.BufferSize,1);
N = deviceWriter.BufferSize;
h = deviceWriter.BufferSize/2;
while true
    msgs = midireceive(midiInput);
    for i = 1:numel(msgs)
        msg = msgs(i);
        if isNoteOn(msg)
            freq = note2freq(msg.Note);
            amp = msg.Velocity/127;
        elseif isNoteOff(msg)
            if msg.Note == msg.Note
                amp = 0;
            end
        end
    end
    if(freq~=prevFreq)
        ll = 4096;
        input = zeros(ll,1);
        n = 0:ll-1;
        n = n';
        w0T = 2*pi*freq/Fs;
        pitch = (Fs/freq);
        disp(pitch);
        P = pitch;
        M = 2*floor(P/2) + 1;
        % input = (M/P)*sincM((M/P * w0T)*n, M);


        input = (M/P)*sin(pi* M*n/P)./(M*sin(pi*n/P));
        %This is not enough yet. If you feed this as input it will not work as
        %there are still many peaks due to zeros in the denominator.
        %We will set all those values to one
        input(mod(n,P)==0) = 1; %this is a certified cantinaro moment though
        input(1) = 1;
        input(isnan(input)) = 1;

        % t = 0:1/Fs:(512-1)/Fs;
        % T = Fs/freq;
        % % Compute the BLIT signal
        % N = length(length(a));
        % blit = sinc(t / T) .* sum(sinc((-N:N) / T), 2);

        % Normalize the signal
        % blit = blit / max(abs(blit));
        % input = blit';


        % input = (M/P)*sincM((M/P)*n, M);
        % w0T = 2*pi*freq/Fs;
        % nharm = floor((Fs/2)/freq); % number of harmonics
        % sig = zeros(1,512);
        % % Synthesize bandlimited impulse train
        % for i=1:nharm
        %     sig = sig + cos(i*w0T*n);
        % end
        % input = sig;
        % input = input.*hamming(size(input,1));
        %input(1:pitch:end) = 1;
        block = filter(1,[1; -a], input);
        block = block/max(abs(block)) * amp;
        %we must window to reduce unwanted artifacts (due to a implicit
        %boxcar windowing)
        % block = block .* hann(length(block));
        % blockNew = circshift(block, -length(block)/2);
        % block = block+blockNew;

        %another idea here would be to take into account how long the block
        %is, how long the buffer is and take the difference as hop

        % block = block(1:512).*hamming(512);
    end
    %some OLA magic is needed here now
    % block2 = block(1:N);
    % block2 = block2 .* hann(N);
    % block2 = block2+remainder;
    % remainder= block(1+h:N+h).*hann(N);
    % block2(1+h:N) = block2(1+h:N)+remainder(1:N-h);
    % remainder(1:h) = remainder(N-h+1:end);
    % remainder(h+1:end) = 0;

    % block = block.*hamming(size(block,1));
    deviceWriter([block(1:N), block(1:N)]);
    block = circshift(block, -N);
end
end
