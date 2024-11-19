global CHANNELS data_seq imageData imageSize bitsPerPixel;
t_packet = 0.2; % Duration of one packet [s]

idx = 1;
while 1
    try
        figure(1);
        subplot(2,2,1);
        est_h = ifft([0; CHANNELS; 0; flip(conj(CHANNELS))], N);
        plot(0:t_packet*10, est_h(1:length(h)));
        ylim([-1 1]); % plot the channel impulse response
        xlabel('Time [ms]')
        ylabel('Impulse response')
        
        subplot(2,2,3);
        plot(0:12000, pow2db(abs(CHANNELS).^2)); % plot the channel frequency response
        xlabel('dB')
        ylabel('Frequency [Hz]')

        subplot(2,2,2);
        colormap(colorMap); image(imageData); axis image; title('Original image'); % the original image
        
        subplot(2,2,4);
        demoded_bitstream = data_seq; % the bit stream that is received by the receiver so far
        imageRx = bitstreamtoimage(demoded_bitstream, imageSize, bitsPerPixel);
        colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
        drawnow;
    catch
        break;
    end
    
    pause(t_packet);
    idx = idx + 1;
end