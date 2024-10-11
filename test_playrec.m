% Tests the functionality of initparams() and recplay.mdl

%% Generate a sinewave
fs = ; % Sampling rate [Hz]
f0 = ; % Frequency of sinewave [Hz]
t_max = ; % Length of the signal [s]

% Generate sinewave
sinewave = ;

% Call to initparams()
[] = initparams();

%% Play and record signal
% Call to recplay.mdl to play simin and record simout
sim('recplay'); 

% Retrieve recorded output
out=simout.signals.values; 

%% Playback
% Play the sound using sound() or soundsc()