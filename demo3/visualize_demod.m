
t_packet = 0.2; % Duration of one packet [s]

idx = 1;
while 1
    try
        figure(1);
        subplot(2,2,1);
        est_h = ifft([0; CHANNELS; 0; flip(conj(CHANNELS))], N);
        plot(0:t_packet*1000, est_h(1:t_packet*1000+1).');
        ylim([-0.1 0.1]); % plot the channel impulse response
        xlabel('Time [ms]')
        ylabel('Impulse response')
        
        subplot(2,2,3);
        plot(0:N/2-2, pow2db(abs(CHANNELS).^2).'); % plot the channel frequency response
        ylim([-50 0])
        xlabel('dB')
        ylabel('Frequency [Hz]')

        subplot(2,2,2);
        colormap(colorMap); image(imageData); axis image; title('Original image'); % the original image
        
        subplot(2,2,4);
        demoded_bitstream = ; % the bit stream that is received by the receiver so far
        imageRx = bitstreamtoimage(demoded_bitstream, imageSize, bitsPerPixel);
        colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
        drawnow;
    catch e
        fprintf(1, "%s\n\n", e.message)
        continue
    end
    
    drawnow;
    pause(t_packet);
    idx = idx + 1;
end
