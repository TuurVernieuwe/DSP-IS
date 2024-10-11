function [ simin,nbsecs,fs ] = initparams( toplay, fs, varargin )
% Creates simin from toplay that will be supplied to recplay.mdl
%
% INPUT:
% toplay    T1X1    Signal supplied as a column vector of T1 samples.
% fs        1X1     Sampling frequency [Hz].
% varargin:
% sync_puls T2X1    Synchronization pulse (you can ignore this until
%                   exercise 5.2.1)
%
% OUTPUT:
% simin     T2X2   Signal to be supplied to recplay.mdl consisting of a
%                  signal to be played in column 1 and zeros in column 2.
% nbsecs    1X1    Length of simin [s].
% fs        1X1    Sampling frequency [Hz].

%% Extract inputs
if nargin == 3
    sync_pulse = varargin{1};
end

%% Asserting assumption of toplay being a column vector
[r,c] = size(toplay);
assert(c==1,'toplay is not a column vector');

%% Add a synchronisation pulse
% you can ignore this until exercise 5.2.1
if nargin == 3
    % Append the synchronisation pulse and zeros at least equal to the length
    % of the IR before toplay
    toplay;
end

%% Append silence to toplay
% Insert 2s of silence at the beginning of toplay and 1s of silcence at the end of toplay
toplay = [];

% Scale toplay to be between [-1,1] -
% -> You can leave this commented until exercise 1.2.7
% toplay = ;

%% Create simin
% Two column matrix: one column is toplay, one column is silence
simin = [];

%% Calculate number of seconds in simin
nbsecs = ;

%% Asserting correctness of the dimensions of simin
% nbsecs should be a scalar
assert(max(size(nbsecs))==1,'nbsecs is the wrong dimension');
% simin should be at least the length of toplay and 3s of silence 
assert(length(simin(1,:))==2 || length(simin) >= (r+3*fs),'simin has the wrong dimensions');
end