clc
clearvars
close all
%%
%go in directory "samples". There are 3 subdirectories "aa", "oo", "uu"
%each of them contains samples with a name in the following format:
%<vowel>_##.wav , where <vowel> is aa, oo or uu, and ## corresponds to the midi note

Fs = 44100; %just in case
p = 100;
w = 4096;
h = 2048;
win = hann(w);

aa_lpc_coeffs = zeros(p, 7);
oo_lpc_coeffs = zeros(p, 7);
uu_lpc_coeffs = zeros(p, 7);

%open the directory
dirinfo = dir('samples');
%cycle through the subdirectories
for i = 3:length(dirinfo) %the first two directories are "." and ".."
    %open the current directory
    subdirinfo = dir(['samples/' dirinfo(i).name]);
    for j = 3:length(subdirinfo)
        %get each file
        [frame, Fs] = audioread(['samples/' dirinfo(i).name '/' subdirinfo(j).name]);
        %extract a frame from the middle of the file
        frame = frame(floor(length(frame)/2)-floor(w/2):floor(length(frame)/2)+floor(w/2)-1, 1);
        %apply the window
        % frame = frame.*win;
        [r,rlags] = xcorr(frame, frame, 'coeff');
        % consider correlation only for positive lags, including 0
        rpos = r(rlags >=0);
        R = toeplitz(rpos(1:p));
        a = R\rpos(2:p+1);
        switch i
            %there are only 3 vowels...
            case 3
                aa_lpc_coeffs(:, j-2) = a;
            case 4
                oo_lpc_coeffs(:, j-2) = a;
            case 5
                uu_lpc_coeffs(:, j-2) = a;
        end
    end
end
%%
save("aa_lpc_coeffs", 'aa_lpc_coeffs');
save("oo_lpc_coeffs", 'oo_lpc_coeffs');
save("uu_lpc_coeffs", 'uu_lpc_coeffs');
%%
%if you want to test
% freq = 220;
% a = uu_lpc_coeffs(:,3);
% imp = zeros(200000, 1);
% pitch = Fs/freq;
% imp(1:pitch:end) = 1;
% out = filter(1, [1; -a ], imp);
% out = out./max(out);
% pl = audioplayer(out, Fs);
% pl.play()
