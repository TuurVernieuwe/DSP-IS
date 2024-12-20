function [OFDM_seq, a, b, nbPackets] = ofdm_mod_stereo( QAM_seq, N, Lcp, H, Lt, Ld, trainblock, Equalization)
% Stereo OFDM modulation
% 
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples].
% H             N/2-1X2     Two frequency domain channels for a DFT size of N.
% Lt            1X1         Number of training frames.
% Ld            1X1         Number of data frames.
% trainblock    T2X1        Training block of T2 QAM symbols
% Equalization  String      Equalization mode: If "fixed" (for transmit_pic_stereo_a), training frames are
%                           transmitted first to estimate the channel and
%                           no additional updating of the channel is
%                           performed thereafter as only data frames
%                           follow. If "packet" packet-based equalization should be
%                           used. (for transmit_pic_stereo_b)
% 
% OUTPUT:
% OFDM_seq      T3X2        Two channel time domain OFDM sequence of length T3 samples.
% a             N/2-1X1     Per-bin factor associated with the first stream.
% b             N/2-1X1     Per-bin factor associated with the second stream.
% nbPackets     1X1         Number of packets, where one packet consist of
%                           a training and data frame.

%% Construct the OFDM sequence
% Padding to make the QAM_seq a multiple of used bins
bins = N/2-1;
padlength = bins - mod(length(QAM_seq), bins);
if padlength ~= bins
    QAM_seq = [QAM_seq; zeros(padlength, 1)];
end

% Construct QAM-matrix
QAM_matrix_raw = reshape(QAM_seq, N/2-1, []);
trainblock_raw = trainblock;

% Calculate per-bin factors
[a, b] = fixed_transmitter_side_beamformer(H(:,1), H(:,2));

if Equalization == "fixed"
    nbPackets = 1;
    for i = 1:2
        % Apply factors
        if i == 1
            QAM_matrix = a.*QAM_matrix_raw;
            trainblock = a.*trainblock_raw;
        else
            QAM_matrix = b.*QAM_matrix_raw;
            trainblock = b.*trainblock_raw;
        end
        
        % Construct the OFDM frames according to Figure 2 in session 3
        fOFDM_frame = [zeros(1, size(QAM_matrix, 2)); QAM_matrix; zeros(1, size(QAM_matrix, 2)); conj(flip(QAM_matrix))];
        fOFDM_frame_trainblock = [0; trainblock; 0; conj(flip(trainblock))];
        
        % Apply the inverse Fourier transform (IFFT)
        OFDM_frame = ifft(fOFDM_frame);
        OFDM_frame_trainblock = ifft(fOFDM_frame_trainblock);
        
        % Add in the cyclic prefix
        OFDM_frame = [OFDM_frame(end-Lcp+1:end, :); OFDM_frame];
        OFDM_frame_trainblock = [OFDM_frame_trainblock(end-Lcp+1:end, :); OFDM_frame_trainblock];
        
        % Construct OFDM sequence
        OFDM_seq(:,i) = [repmat(OFDM_frame_trainblock, Lt, 1); reshape(OFDM_frame, [],1);];
    end
elseif Equalization == "packet"
    % Padding to make the amount of columns a multiple of Ld
    fOFDM_pad = Ld - mod(size(QAM_matrix_raw, 2), Ld);
    if ~(Ld == fOFDM_pad)
        QAM_matrix_raw = [QAM_matrix_raw, zeros(size(QAM_matrix_raw, 1), fOFDM_pad)];
    end
    nbPackets = size(QAM_matrix_raw, 2)/Ld;
    len = (N+Lcp)*(Lt+Ld); % Number of elements in one packet
    OFDM_seq = zeros(len*nbPackets, 2);
    for i = 1:2
        % Apply factors
        if i == 1
            QAM_matrix = a.*QAM_matrix_raw;
            trainblock = a.*trainblock_raw;
        else
            QAM_matrix = b.*QAM_matrix_raw;
            trainblock = b.*trainblock_raw;
        end

        % Construct the OFDM frames according to Figure 2 in session 3
        fOFDM_frame = [zeros(1, size(QAM_matrix, 2)); QAM_matrix; zeros(1, size(QAM_matrix, 2)); conj(flip(QAM_matrix))];
        fOFDM_frame_trainblock = [0; trainblock; 0; conj(flip(trainblock))];
        
        % Apply the inverse Fourier transform (IFFT)
        OFDM_frame = ifft(fOFDM_frame);
        OFDM_frame_trainblock = ifft(fOFDM_frame_trainblock);
        
        % Add in the cyclic prefix
        OFDM_frame = [OFDM_frame(end-Lcp+1:end, :); OFDM_frame];
        OFDM_frame_trainblock = [OFDM_frame_trainblock(end-Lcp+1:end, :); OFDM_frame_trainblock];
        
        % Construct serialized training packet
        trainpacket = repmat(OFDM_frame_trainblock, Lt, 1);
        
         % Construct serialized packets
        for j = 1:nbPackets
            OFDM_seq(1+len*(j-1):len*j, i) = [trainpacket; reshape(OFDM_frame(:, 1+Ld*(j-1):Ld*j), [], 1)];
        end
    end
end
end