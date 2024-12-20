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
    QAM_matrix = QAM_matrix ./ CHANNELS;
    
    % Supply streamLength number of symbols (you can ignore this until exercise 4.2)
    data_seq = reshape(QAM_matrix, [], 1);
elseif Equalization == "packet"
    bins = N/2-1;
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
        
        % Estimate channel frequency responses
        % Make initial channel estimation based on training frames
        CHANNEL = zeros(N/2-1, 1);
        for j=1:N/2-1
            CHANNEL(j) = (trainblock(j)*ones(length(trainpacket(j,:).'), 1)) \ trainpacket(j,:).';
        end
        W = 1./conj(CHANNEL);
        QAM_matrix(:,1) = conj(W).*QAM_matrix(:,1);
        
        % DD equalization
        Ld = size(QAM_matrix, 2);
        for j = 1:Ld-1 % Columns
            % NMLS filter implementation.
            % Initialise filters, reconstructed transmitted signal and error
            rec_Xk = zeros(N/2-1,1); 
            Ek = zeros(N/2-1,1);
            for n = 1:N/2-1 % Rows
                % Apply filter.
                estXk = W(n)' * QAM_matrix(n, j+1);
                % Reconstruct transmitted signal.
                rec_bits = qam_demod(estXk, M, log2(M));
                rec_Xk(n) = qam_mod(rec_bits, M);
                % Calculate error signal.
                Ek(n) = rec_Xk(n) - estXk;
                % Update filter.
                W(n) = W(n) + mu/(alpha + QAM_matrix(n,j+1)'.*QAM_matrix(n,j+1)) .* QAM_matrix(n,j+1) .* conj(Ek(n));
            end
            % Apply equalisation with the updated values
            QAM_matrix(:,j+1) = conj(W).*QAM_matrix(:,j+1);
        end
        CHANNELS(:,i) = CHANNEL;%1./conj(W);
        
        % Supply streamLength number of symbols (you can ignore this until exercise 4.2)
        QAM_seq = reshape(QAM_matrix, [], 1);
        data_seq(1+(i-1)*bins*Ld:i*bins*Ld) = QAM_seq;
    end
end
end