function Rx = simulate_channel(Tx, Nswitch,filepath,smoothing_factor)
% Simulates a slowly changing channel, with a channel change every Nswitch
% number of samples.
% 
% INPUT:
% Tx                TX1     Stream to transmit of length T samples.
% Nswitch           1X1     The channel changes every Nswitch number of samples.
% filepath          String  Filepath to channel to use.
% smoothing_factor  1X1     Smoothing factor to model a slowly changing
%                           channel. If 0 a violently changing channel is modelled, if 1 the channel does not change. 
% 
% OUTPUT:
% Rx        TX1     Received stream of length T samples.

%% Pad the Transmitted stream with zeros for the length to be a multiple of Nswitch
Tx_length = length(Tx);
P = ceil(length(Tx)/Nswitch); % Number of channel changes
padLength = Nswitch*P - length(Tx);
Tx = [Tx; zeros(padLength,1)]; 
Rx = zeros(size(Tx));

%% Apply channel
load(filepath); % Initial channel
h_init = h;
for p=1:P
    out = fftfilt(h, Tx); % Apply channel
    Rx( (p-1)*Nswitch+1 : p*Nswitch ) = out( (p-1)*Nswitch+1 : p*Nswitch );
    h = smoothing_factor * h + (1-smoothing_factor) * (max(abs(h_init))*randn(size(h,1),1)); % channel is changing slowly 
end

%% Return received stream of same length as transmitted stream
Rx = Rx(1:Tx_length);

end
