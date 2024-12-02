function [ simin,nbsecs,fs ] = initparams_stereo( toplay, fs )
% Creates simin from toplay that will be supplied to recplay.mdl
%
% INPUT:
% toplay    T1X2    Signal supplied as a 2-column vector of T1 samples.
% fs        1X1     Sampling frequency [Hz].
%
% OUTPUT:
% simin     T2X2   Signal to be supplied to recplay.mdl consisting of a
%                  signal to be played in column 1 and another signal to be played in column 2.
% nbsecs    1X1    Length of simin [s].
% fs        1X1    Sampling frequency [Hz].

%% Add a synchronisation pulse
f_sync = 5000;
t = 0:1/fs:1;
sync_pulse = sin(2*pi*f_sync*t)'; 

% Append the synchronisation pulse and zeros at least equal to the length
% of the IR before toplay and Insert 2s of silence at the beginning of toplay and 1s of silcence at the end of toplay
toplay = [zeros(2*fs, 2); sync_pulse, sync_pulse; zeros(fs/10, 2); toplay; zeros(fs, 2)];

% Scale toplay to be between [-1,1] 
toplay = rescale(toplay, -1, 1);

%% Create simin
simin = toplay;

%% Calculate number of seconds simin is
nbsecs = size(simin, 1)/fs;

end

