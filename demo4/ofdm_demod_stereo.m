function [ data_seq, CHANNELS] = ofdm_demod_stereo( OFDM_seq, N, Lcp, trainblock, Lt, Ld, M, nbPackets, Equalization, mu, alpha )
% Performs stereo OFDM demodulation
% 
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples].
% trainblock    T2X1            Training block of T2 QAM symbols 
% M             1X1             QAM-ary constellation size
% Lt            1X1             Number of training frames
% Ld            1X1             Number of data frames
% nbPackets     1X1             Number of packets, where one packet consist of
%                               a training and data frame.
% Equalization  String          Equalization mode: If "fixed" (for transmit_pic_stereo_a), training frames are
%                               transmitted first to estimate the channel and
%                               no additional updating of the channel is
%                               performed thereafter as only data frames
%                               follow. If "packet" packet-based equalization should be
%                               used. (for transmit_pic_stereo_b)
% mu            1X1             NLMS stepsize
% alpha         1X1             NLMS regularization factor
% 
% OUTPUT:
% data_seq      T3X1            QAM sequence of T3 symbols.
% CHANNELS      N/2-1XP         Frequency domain estimated combined channel for each 1..P.

%% Perform OFDM demodulation
if Equalization == "fixed"
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
    
    % Calculate channel frequency response
    CHANNELS = zeros(N/2-1, 1);
    for j=1:N/2-1
        CHANNELS(j) = (trainblock(j)*ones(length(trainpacket(j,:).'), 1)) \ trainpacket(j,:).';
    end
    figure
    plot(abs(CHANNELS));
    QAM_matrix = QAM_matrix ./ CHANNELS;
    
    % Supply streamLength number of symbols (you can ignore this until exercise 4.2)
    data_seq = reshape(QAM_matrix, [], 1);
elseif Equalization == "packet"
    bins = sum(ON_OFF_mask);
    data_seq = zeros(bins*Ld*nbPackets,1);
    CHANNELS = zeros(N/2-1, nbPackets);
    for i = 1:nbPackets
        % Reshape the received OFDM sequence (serial to parallel conversion)
        OFDM_matrix = reshape(OFDM_seq((i-1)*(Lt+Ld)*(N+Lcp)+1:i*(Lt+Ld)*(N+Lcp)), N+Lcp, []);
        
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
end

