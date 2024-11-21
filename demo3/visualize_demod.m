t_packet = 0.2; % Duration of one packet [s]
idx = 1;
while 1
    try
        figure(1);
        subplot(2,2,1);
        CHANNELS = CHANNELS.*ON_OFF_mask;
        est_h = ifft([0; CHANNELS; 0; flip(conj(CHANNELS))]);
        plot(real(est_h()));
        xlim([0 200])
        ylim([-0.1 0.1]); % plot the channel impulse response
        title('Channel in time domain')
        
        subplot(2,2,3);
        T0 = fs/N;
        plot(T0:T0:T0*(N/2-1), pow2db(abs(CHANNELS).^2).'); % plot the channel frequency response
        ylim([-50 0])
        xlabel('Frequency [Hz]')
        ylabel('dB')
        title('Channel in frequency domain (no DC)')

        subplot(2,2,2);
        colormap(colorMap); image(imageData); axis image; title('Original image'); % the original image
        
        subplot(2,2,4);
        bins = sum(ON_OFF_mask);
        demoded_bitstream = rx_bits(1:idx*bins*Ld); % the bit stream that is received by the receiver so far
        imageRx = bitstreamtoimage(demoded_bitstream, imageSize, bitsPerPixel);
        colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
        drawnow;
    catch e
        fprintf(1, "%s\n\n", e.message)
        continue
    end
    
    drawnow;
    pause(t_packet);
    idx = idx + 1;
end
