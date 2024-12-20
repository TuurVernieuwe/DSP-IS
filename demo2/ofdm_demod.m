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
% streamLength  1X1             Lengtha of QAM sequence after call to qammod [samples].
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
% CHANNELS      N/2-1XP         Frequency domain estimated channel with P either nbPackets or 1. (
%                               you can ignore this until exercise 5.1.4)    

%% Extract input arguments
if nargin == 4 % Session 3
    streamLength = varargin{1};
elseif nargin == 7 % Session 4
    streamLength = varargin{1};
    channel = varargin{2};
    ON_OFF_mask = varargin{3};
    equalization = varargin{4};
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
if ~(nargin == 9 || nargin == 10)
    % Reshape the received OFDM sequence (serial to parallel conversion)
    OFDM_matrix = reshape(OFDM_seq, N+Lcp, []);

    % Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
    OFDM_matrix = OFDM_matrix(Lcp+1:end, :);

    % Apply fft operation
    QAM_matrix = fft(OFDM_matrix);

    % Remove the redundant parts of QAM_matrix
    QAM_matrix = QAM_matrix(2:size(QAM_matrix, 1)/2, :);

    % Apply channel equalisation (you can ignore this until exercise 4.2.3)
    if exist("equalization", "var")
        if equalization
            CHANNELS = fft(channel, N);
            QAM_matrix = QAM_matrix ./ CHANNELS(2:N/2);
        end
    elseif exist("trainblock", "var")
        CHANNELS = zeros(N/2-1, 1);
        width = size(QAM_matrix, 2);
        index = 1;
        for i=1:N/2-1
            if ON_OFF_mask(i)
                CHANNELS(i) = (trainblock(index)*ones(width, 1)) \ QAM_matrix(i,:).';
                index = index + 1;
            end
        end
        QAM_matrix = QAM_matrix ./ CHANNELS;
    end

    % Apply on-off mask (you can ignore this until exercise 4.3)
    QAM_matrix = QAM_matrix(logical(ON_OFF_mask),:);

    % Supply streamLength number of symbols (you can ignore this until exercise 4.2)
    QAM_seq = reshape(QAM_matrix, [], 1);
    data_seq = QAM_seq(1:streamLength);
elseif (nargin == 9)
    bins = sum(ON_OFF_mask);
    data_seq = zeros(bins*Ld*nbPackets,1);
    CHANNELS = zeros(N/2-1, nbPackets);
    for i = 1:nbPackets
        % Reshape the received OFDM sequence (serial to parallel conversion)
        OFDM_matrix = reshape(OFDM_seq(1+(i-1)*(Lt+Ld)*(N+Lcp):i*(Lt+Ld)*(N+Lcp)), N+Lcp, []);
        
        % Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
        OFDM_matrix = OFDM_matrix(Lcp+1:end, :);
        
        % Apply fft operation
        QAM_matrix = fft(OFDM_matrix);
        
        % Remove the redundant parts of QAM_matrix
        QAM_matrix = QAM_matrix(2:size(QAM_matrix, 1)/2, :);
        
        % Split matrix
        trainpacket = QAM_matrix(:, 1:Lt);
        QAM_matrix = QAM_matrix(:, Lt+1:end);
        
        % Calculate channel frequency response
        CHANNEL = zeros(N/2-1, 1);
        index = 1;
        for j=1:N/2-1
            if ON_OFF_mask(j)
                CHANNEL(j) = (trainblock(index)*ones(length(trainpacket(j,:).'), 1)) \ trainpacket(j,:).';
                index = index + 1;
            end
        end
        QAM_matrix = QAM_matrix ./ CHANNEL;
        CHANNELS(:, i) = CHANNEL;
        
        % Apply on-off mask (you can ignore this until exercise 4.3)
        QAM_matrix = QAM_matrix(logical(ON_OFF_mask),:);
        
        % Supply streamLength number of symbols (you can ignore this until exercise 4.2)
        QAM_seq = reshape(QAM_matrix, [], 1);
        data_seq(1+(i-1)*bins*Ld:i*bins*Ld) = QAM_seq;
    end
    data_seq = data_seq(1:streamLength);
else % nargin = 10
    % Reshape the received OFDM sequence (serial to parallel conversion)
    OFDM_matrix = reshape(OFDM_seq, N+Lcp, []);

    % Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
    OFDM_matrix = OFDM_matrix(Lcp+1:end, :);

    % Apply fft operation
    QAM_matrix = fft(OFDM_matrix);

    % Remove the redundant parts of QAM_matrix
    QAM_matrix = QAM_matrix(2:size(QAM_matrix, 1)/2, :);

    % Split matrix
    trainpacket = QAM_matrix(:, 1:Lt);
    QAM_matrix = QAM_matrix(:, Lt+1:end);
    
    % Estimate channel frequency response
    % Make initial channel estimation based on training frames
    CHANNEL = zeros(N/2-1, 1);
    index = 1;
    for j=1:N/2-1
        if ON_OFF_mask(j)
            CHANNEL(j) = (trainblock(index)*ones(length(trainpacket(j,:).'), 1)) \ trainpacket(j,:).';
            index = index + 1;
        end
    end
    W = 1./conj(CHANNEL);
    QAM_matrix(:,1) = conj(W).*QAM_matrix(:,1);
    
    % DD equalization
    Ld = size(QAM_matrix, 2);
    CHANNELS = zeros(N/2-1, Ld);
    CHANNELS(:,1) = 1./conj(W);
    for i = 1:Ld-1 % Columns
        % NMLS filter implementation.
        % Initialise filters, reconstructed transmitted signal and error
        rec_Xk = zeros(N/2-1,1); 
        Ek = zeros(N/2-1,1);
        for n = 1:N/2-1 % Rows
            if ON_OFF_mask(n)
                % Apply filter.
                estXk = W(n)' * QAM_matrix(n, i+1);
                % Reconstruct transmitted signal.
                rec_bits = qam_demod(estXk, M, log2(M));
                rec_Xk(n) = qam_mod(rec_bits, M);
                % Calculate error signal.
                Ek(n) = rec_Xk(n) - estXk;
                % Update filter.
                W(n) = W(n) + mu/(alpha + QAM_matrix(n,i+1)'.*QAM_matrix(n,i+1)) .* QAM_matrix(n,i+1) .* conj(Ek(n));
            end
            CHANNELS(:, n+1) = 1./conj(W);
        end
        % Apply equalisation with the updated values
        QAM_matrix(:,i+1) = conj(W).*QAM_matrix(:,i+1);
    end

    % Apply on-off mask (you can ignore this until exercise 4.3)
    QAM_matrix = QAM_matrix(logical(ON_OFF_mask),:);

    % Supply streamLength number of symbols (you can ignore this until exercise 4.2)
    data_seq = reshape(QAM_matrix, [], 1);
end