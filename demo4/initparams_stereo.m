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
sync_pulse = ; 

% Append the synchronisation pulse and zeros at least equal to the length
% of the IR before toplay and Insert 2s of silence at the beginning of toplay and 1s of silcence at the end of toplay
toplay = ;

% Scale toplay to be between [-1,1] 
toplay = ;

%% Create simin
simin = ;

%% Calculate number of seconds simin is
nbsecs = ;

end

