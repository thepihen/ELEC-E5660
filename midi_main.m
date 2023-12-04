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
playFromMIDI = (userIn == "y" || userIn == "Y");
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
deviceWriter.BufferSize = 256;
if playFromMIDI
    simplesynthOptionalMIDI(deviceWriter, aa_lpc_coeffs, Fs, true, selectedDevice);
else
    simplesynthOptionalMIDI(deviceWriter, aa_lpc_coeffs, Fs, false, "none");
end

function simplesynthOptionalMIDI(deviceWriter,a_lpc, Fs, playingFromMidiFlag, selectedDevice)
if(selectedDevice~="none")
    midiInput = mididevice(selectedDevice);
end
deviceWriter.SupportVariableSizeInput = true;
% a = a_lpc(:,40);
% a = a_lpc(:,1);
amp = 0.6;
freq = 0;
prevFreq = -1;
% remainder = zeros(deviceWriter.BufferSize,1);
N = deviceWriter.BufferSize;
% N = 512;
h = deviceWriter.BufferSize/2;
notes = [59,61,63,70,76,80,85];
currNote = 1; %for playing with a MIDI keyboard
arrCurrNote = 1; %for playing without a MIDI keyboard
delayBetweenNotes = 1;
noteTic = tic;
% hasFinishedProcessingPreviousBlock = true;
currBlock = 1;
%buildEnvelope(Fs, levels, times)
%levels = [aLev, dLev, sLev]
%times = [a-d, d-s, s, r]
%the sustain one is not really needed, is for when you're playing without a
%keyboard
env_aLev = 0;
env_dLev = 1;
env_sLev = 0.2;
env_ad = 1;
env_ds = 0.5;
env_s = 3;
env_r = 0.2;

[envelope,envTimes] = buildEnvelope(Fs, [env_aLev, env_dLev, env_sLev], [env_ad, env_ds, env_s, env_r]);
envelope = envelope';
% disp(size(envelope));
ll = 4096;%block length
env_AD = envelope(1:envTimes(1)+envTimes(2)); %this has to be applied automatically
env_R = envelope(end-envTimes(3)+1:end);
remainder_AD = mod(length(env_AD), ll);
%envelope = [envelope, zeros(1, remainder)]';
times = (length(env_AD)/ll) *(ll/N);
envelope = env_AD;
currTime = 1;

isPlayingNote = false;
while true

    if (playingFromMidiFlag)
        msgs = midireceive(midiInput);

        for i = 1:numel(msgs)
            msg = msgs(i);
            switch (msg.Type)
                case midimsgtype.NoteOn
                    freq = note2freq(msg.Note);
                    amp = msg.Velocity/127;
                    envelope = env_AD;
                    times = (length(env_AD)/ll) *(ll/N);
                    currNote = msg.Note;
                    isPlayingNote = true;

                    %buildEnvelope(Fs, levels, times)
                    %levels = [aLev, dLev, sLev]
                    %times = [a-d, d-s, s, r]
                    %the sustain one is not really needed, is for when you're playing without a
                    %keyboard
                    [envelope,envTimes] = buildEnvelope(Fs, [env_aLev, ...
                        env_dLev, env_sLev], [env_ad, env_ds, env_s, env_r]);
                    envelope = envelope';
                    % disp(size(envelope));
                    ll = 4096;%block length
                    env_AD = envelope(1:envTimes(1)+envTimes(2)); %this has to be applied automatically
                    env_R = envelope(end-envTimes(3)+1:end);
                    remainder_AD = mod(length(env_AD), ll);
                    %envelope = [envelope, zeros(1, remainder)]';
                    times = (length(env_AD)/ll) *(ll/N);
                    envelope = env_AD;

                    currTime = 1;

                case midimsgtype.NoteOff
                    % if msg.Note == msg.Note
                    if(currNote == msg.Note)
                        envelope = env_R;
                        times = (length(envelope)/ll) *(ll/N);
                        isPlayingNote = false;
                        % end
                        currTime = 1;
                    end
                case midimsgtype.ControlChange
                    switch msg.CCNumber
                        case 18
                            %let's set a maximum of 5 seconds, and a
                            %minimum of 10 milliseconds
                            env_ad = (msg.CCValue/127 * 0.49) +0.01;
                            disp("attack changed")
                        case 19
                            disp("decay chaged")
                            %same here
                            env_ds = (msg.CCValue/127 * 4.99) +0.01;
                        case 20
                            %this only changes the sustain level
                            env_sLev = (msg.CCValue/127);%from 0 to 1
                            disp("sustain changed")
                        case 21
                            disp("release changed")
                            env_r= (msg.CCValue/127 * 3) ;
                    end
            end
        end

    else
        currNote = notes(arrCurrNote);
        freq = note2freq(notes(arrCurrNote));
        if(toc(noteTic)>delayBetweenNotes)
            arrCurrNote = arrCurrNote + 1;
            if arrCurrNote > 7
                arrCurrNote = 1;
            end
            noteTic = tic;
        end
    end



    if(freq~=prevFreq)
        currTime = 1;
        lpcInd = matchClosestLPCCoeff(currNote);
        a = a_lpc(:,lpcInd);
        % input = zeros(ll,1);
        n = 0:ll-1;
        n = n';
        % w0T = 2*pi*freq/Fs;
        blockFreq = freq + randn(1);
        pitch = (Fs/blockFreq);
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
    %there is a discontinuity between the beginning and the end (actually
    %between the end and the beinning). A better algorithm would loop an
    %area in the middle, but right now I don't have the time to follow
    %that approach, hence this tries to address said discontinuity
    %with some overlap-and-add between the last frame and the first one.
    %This means one frame is "skipped"
    if(mod(currBlock, ll/N)==0)
        block = circshift(block, -length(sblock));
        nextBlock = block(1:N,:);
        fakeBlock = [nextBlock; sblock].*hamming(size(sblock,1)*2);
        sblock = fakeBlock(length(fakeBlock)/2 + 1:end) + fakeBlock(1:length(fakeBlock)/2);
        currBlock = currBlock + 1; %we're "skipping" one, we need to account for that
        %or the algorithm will only work on the first round
    end
    % if(currBlock = ll/N)
    %     sblock2 = sblock.*hamming(length(sblock), 'periodic');
    %     sblock(1:length(sblock)/2, :) = sblock2(1:length(sblock2)/2, :);
    % end
    disp(currTime+ " / " + times);
    if(currTime <= times )
        if(currTime == times )
            %remainder
            coveredEnv = [envelope(end-remainder+1:end); ones(length(envelope)-remainder, 1)*env_sLev*isPlayingNote];
        else
            coveredEnv = envelope((currTime-1)*N +1:currTime*N);
        end
    else
        %the envelope will be incomplete
        coveredEnv = ones(size(coveredEnv))*env_sLev*isPlayingNote;
    end
    disp(coveredEnv(10));
    disp(env_s)
    sblock = sblock.*coveredEnv*amp;
    currTime = currTime +1;


    block = circshift(block, -length(sblock));
    deviceWriter([sblock, sblock]);
    currBlock = currBlock + 1;
    % block = circshift(block, -N);
    % disp(currTime + " / "+times);


end
end


function block = blockTransition(prevBlock, currBlock)
windowed = currBlock .* hanning(length(currBlock));
windowed_prev = prevBlock .* hanning(length(prevBlock));
block = [windowed(1:length(currBlock)/2)+windowed_prev(length(windowed_prev)/2 + 1:end); currBlock(length(currBlock)/2+1:end)];
end

function c = matchClosestLPCCoeff(midiNoteNumber)
arr = [60,64,67,72,76,79,84];
[~,c] = min(abs(arr-midiNoteNumber));
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

function [env,lengths] = buildEnvelope(Fs, levels, times)
%levels = [aLev, dLev, sLev] %the sustain level is implied to be the
%same as the release level
%times = [a-d, d-s, s, r]
%a,d,s,r = attack, decay , sustain, release in seconds
aEnv = linspace(levels(1),levels(2), round(Fs*times(1)));
dEnv = linspace(levels(2), levels(3), round(Fs*times(2)));
sEnv = linspace(levels(3), levels(3), round(Fs*times(3)));
rEnv = linspace(levels(3), 0, round(Fs*times(4)));
env = [aEnv, dEnv, sEnv, rEnv];
lengths = [length(aEnv), length(dEnv), length(rEnv)];
end



%https://www.music.mcgill.ca/~gary/307/week5/node14.html
%https://ccrma.stanford.edu/~jos/SpecEnv/LPC_Envelope_Example_Speech.html
function simplesynth(midiDeviceName, deviceWriter,a_lpc, Fs)
midiInput = mididevice(midiDeviceName);
deviceWriter.SupportVariableSizeInput = true;
a = a_lpc(:,40);
amp = 0.6;
freq = 0;
prevFreq = -1;
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


%% TODOS
%-if you are moving from one note to the next (noteOff after a noteOn) then
%the pitch should "slide" to the new one (using the blockTransition
%function)
%
%-make the script work for any buffer size, always keeping 256-samples
%sub-buffers you will use to build larger ones
%
%-move to midi
%
%-midi noteOn, noteOff detection, and maybe CC controls
