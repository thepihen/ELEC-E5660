% ELEC-E5660 - LPC-based vocal synthesizer
% main file
% student: Francesco Colotti (id 101613215)

%% startup
clearvars
close all
clc

%% SETUP
%1) Setup midi device for play
mididevinfo; %use mididevinfo to retrieve list of available midi devices
selectedDevice = "KOMPLETE KONTROL M32 MIDI"; %example, used in development

%2) Set sampling frequency
Fs = 44100;

%3) Setup audioDeviceWriter object to write buffers
%A buffer size of 192 samples should guarantee a pretty low latency
%but may not be supported by all devices
%A device with an ASIO driver is preferable

%deviceWriter = audioDeviceWriter('SampleRate',44100, ...
% "Driver", "ASIO", "Device", "Focusrite USB ASIO","SupportVariableSizeInput", true,...
% "BitDepth", "24-bit integer"); % example of input
%deviceWriter.BufferSize = 192; %set manually just in case 

%deviceWriter = audioDeviceWriter('SampleRate',Fs, ...
%    "Driver", "ASIO", "Device", "FL Studio ASIO","SupportVariableSizeInput", true,...
%    "BitDepth", "24-bit integer"); %used in development
% asiosettings(deviceWriter.Device);
%deviceWriter.BufferSize = 256; % check this is supported by the actual
%device

%if all this information is not available one can simply call
deviceWriter = audioDeviceWriter('SampleRate', Fs);


%% load previously computed LPC coefficients
load('coeffs/aa_lpc_coeffs.mat');
load('coeffs/oo_lpc_coeffs.mat');
load('coeffs/uu_lpc_coeffs.mat');
% a small note on these LPC coefficients. They are order 10 to 100 LPC coefficients
% extracted from 7 recordings per vowel. The recordings cover C,E,G of 2
% octaves, plus the C from a third one
% this is, unfortunately, the vocal range I, a completely untrained singer,
% was able to (more or less) consistently hit notes on

%% Choose whether to play with MIDI device
userIn = "d";
while(userIn ~= "y" && userIn ~= "n")
    userIn = input("Do you wish to play from your external MIDI keyboard? (y/n)","s");
end
playFromMIDI = (userIn == "y" || userIn == "Y");

if playFromMIDI
    simplesynthOptionalMIDI(deviceWriter, oo_lpc_coeffs, Fs, true, selectedDevice);
else
    simplesynthOptionalMIDI(deviceWriter, aa_lpc_coeffs, Fs, false, "none");
end

%--------------------------------------------------------------------------



function simplesynthOptionalMIDI(deviceWriter,aa_lpc, Fs, playingFromMidiFlag, selectedDevice)
if(selectedDevice~="none")
    midiInput = mididevice(selectedDevice);
end
deviceWriter.SupportVariableSizeInput = true;

amp = 0.6;
freq = 0;
prevFreq = -1;
N = deviceWriter.BufferSize;
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
reverb = reverberator;
lpcOrder = 100;
a_lpc = aa_lpc(1:lpcOrder, :, lpcOrder-9); %-9 because of course indexing starts at 1
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
                        case 14
                            lpcOrder = round((msg.CCValue/127) *(size(aa_lpc,3) - 1))+10;
                            disp(lpcOrder)
                        case 15
                            reverb.WetDryMix = msg.CCValue/127;
                        case 18
                            %let's set a maximum of 5 seconds, and a
                            %minimum of 10 milliseconds
                            env_ad = (msg.CCValue/127 * 0.49) +0.01;
                            disp("attack changed")
                        case 19
                            disp("decay chaged")
                            %same here
                            env_ds = (msg.CCValue/127 * 0.49) +0.01;
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
        a = aa_lpc(1:lpcOrder, lpcInd, lpcOrder-9);
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
    % disp(currTime+ " / " + times); %DEBUG
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
    sblock = sblock.*coveredEnv*amp;
    currTime = currTime +1;


    block = circshift(block, -length(sblock));
    deviceWriter(reverb([sblock, sblock]));
    currBlock = currBlock + 1;

end
end


function c = matchClosestLPCCoeff(midiNoteNumber)
arr = [60,64,67,72,76,79,84];
[~,c] = min(abs(arr-midiNoteNumber));
end

function freq = note2freq(note)
freqA = 440;
noteA = 69;
freq = freqA * 2.^((note-noteA)/12);
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
% uncomment if you want to plot the envelope
% figure;
% t = (0:length(env)-1)/Fs;
% plot(t, env);
% grid on
% xlabel("Time [s]");
% ylabel("Env");
% title("Envelope");
lengths = [length(aEnv), length(dEnv), length(rEnv)];
end