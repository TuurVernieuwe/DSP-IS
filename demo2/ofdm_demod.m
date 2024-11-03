function [ data_seq, CHANNELS ] = ofdm_demod(OFDM_seq, N, Lcp, varargin )
% OFDM demodulation
%
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T1 samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% varargin  
% For Session 3 
% empty
%
% For Session 4
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% channel       T2X1            Impulse response of channel.
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
%                               (you can ignore this until exercise 4.3)
% equalization  1X1             If 1 channel equalization is performed, if 0 no
%                               channel equalization is performed.
%
% For session 5
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
%                               (you can ignore this until exercise 4.3)
% trainblock    T3X1            Training block of T3 QAM symbols (you can
%                               ignore this until exercise 6.1.1)
%
% For Session 6
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
%                               (you can ignore this until exercise 4.3)
% trainblock    T3X1            Training block of T3 QAM symbols (you can
%                               ignore this until exercise 6.1.1)
% Lt            1X1             Number of training frames (you can ignore this
%                               until exercise 6.1.1)
% Ld            1X1             Number of data frames (you can ignore this
%                               until exercise 6.1.1)
% nbPackets     1X1             Number of packets, where one packet consist of
%                               a training and data frame. (you can ignore 
%                               this until exercise 6.1.1)
% Session 7
% Lt            1X1             Number of training frames
% M             1X1             QAM-ary constellation size
% trainblock    T3X1            Training block of T3 QAM symbols
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
% mu            1X1             Adaptive filter stepsize
% alpha         1X1             Adaptive filter regularisation factor
% type          1X1             Adaptive filter type (supported: 'NLMS')
%
% OUTPUT:
% data_seq      T4X1            QAM sequence of T4 symbols.
% CHANNELS      N/2-1XP         Frequency domain estiamted channel with P either nbPackets or 1. (
%                               you can ignore this until exercise 5.1.4)    

%% Extract input arguments
if nargin == 4 % Session 3
    streamLength = varargin{1};
elseif nargin == 5 % Session 4
    streamLength = varargin{1};
    channel = varargin{2};
    % ON_OFF_mask = varargin{3};
    % equalization = varargin{4};
elseif nargin == 6 % Session 5
    streamLength = varargin{1};
    ON_OFF_mask = varargin{2};
    trainblock = varargin{3};
elseif nargin == 9 % Session 6
    streamLength = varargin{1};
    ON_OFF_mask = varargin{2};
    trainblock = varargin{3};
    Lt = varargin{4};
    Ld = varargin{5};
    nbPackets = varargin{6};
elseif nargin == 10 % Session 7
    Lt = varargin{1};
    M = varargin{2};
    trainblock = varargin{3};
    ON_OFF_mask = varargin{4};
    mu = varargin{5};
    alpha = varargin{6};
    type = varargin{7};
end

%% Perform OFDM demodulation
% Reshape the received OFDM sequence (serial to parallel conversion)
OFDM_matrix = reshape(OFDM_seq, N+Lcp, []);

% Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
OFDM_matrix = OFDM_matrix(Lcp+1:end, :);

% Apply fft operation
QAM_matrix = fft(OFDM_matrix);

% Remove the redundant parts of QAM_matrix
QAM_matrix = QAM_matrix(2:size(QAM_matrix, 1)/2, :);

% Apply channel equalisation (you can ignore this until exercise 4.2.3)
if exist('channel', 'var')
    CHANNELS = fft(channel, N/2-1); % Compute CFR with N points

    % Equalize: scale each subcarrier by the inverse of the channel response
    QAM_matrix = QAM_matrix ./ CHANNELS; % Equalization
else
    CHANNELS = []; % No equalization applied
end
% Apply on-off mask (you can ignore this until exercise 4.3)
% QAM_matrix = ;

% Supply streamLength number of symbols (you can ignore this until exercise 4.2)
QAM_seq = reshape(QAM_matrix, [], 1);
data_seq = QAM_seq(1:streamLength);

end