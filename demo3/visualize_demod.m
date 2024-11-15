t_packet = ; % Duration of one packet [s]

idx = 1;
while 1
    try
        figure(1);
        subplot(2,2,1);
        plot();
        ylim([-1 1]); % plot the channel impulse response
        xlabel('')
        ylabel('')
        
        subplot(2,2,3);
        plot(); % plot the channel frequency response
        xlabel('')
        ylabel('')

        subplot(2,2,2);
        colormap(colorMap); image(imageData); axis image; title('Original image'); % the original image
        
        subplot(2,2,4);
        demoded_bitstream = ; % the bit stream that is received by the receiver so far
        imageRx = bitstreamtoimage(demoded_bitstream, imageSize, bitsPerPixel);
        colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
        drawnow;
    catch
        break;
    end
    
    pause();
    idx = idx + 1;
end