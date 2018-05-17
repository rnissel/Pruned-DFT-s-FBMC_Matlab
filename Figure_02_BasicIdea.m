% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================       
% This script illustrates the basic idea of pruned DFT spread FBMC by 
% showing the underlying basis pulses. 

% Reproduces Figure 2 of "Pruned DFT Spread FBMC: Low PAPR, Low Latency, 
% High Spectral Efficiency", R. Nissel and M. Rupp, IEEE Transactions on 
% Communications

clear; close all;

%% Parameters
Nfft = 512;          % FFT size
L    = 16;           % Number of subcarrier


%% DFT matrices and RRC window
% Full DFT matrices
W_Nfft = fft(eye(Nfft))/sqrt(Nfft);
W_L    = fft(eye(L))/sqrt(L);

% Pruned DFT matrix
W_LxLo2 = W_L(:,L/4+(1:L/2));

% Root Raised Cosine window
RRCwindow = sqrt( ( cos( 2*pi*( (0:Nfft-1)-Nfft/2 ) .' / Nfft ) + 1 ) ) / sqrt(2);


%% Transmit Matrices
% OFDM
G_OFDM = W_Nfft'*[ eye(L)*sqrt(Nfft/L); zeros(Nfft-L,L) ];

% SC-FDMA, i.e., DFT spread OFDM
G_SCFDMA = G_OFDM*W_L;

% Pruned DFT spread OFDM
G_pDFTsOFDM = G_OFDM*W_LxLo2;

% Pruned DFT spread FBMC (unscaled)
G_pDFTsFBMCunscaled = repmat(RRCwindow,1,L/2).*G_pDFTsOFDM;

% ################## Pruned DFT spread FBMC #########################
G_pDFTsFBMC = G_pDFTsFBMCunscaled*diag(1./sqrt(diag(G_pDFTsFBMCunscaled'*G_pDFTsFBMCunscaled)))*sqrt(Nfft/L);
% ###################################################################


% Add zeros in time for illustration purpose
G_OFDM              = [ zeros(Nfft/4,L);   G_OFDM;              zeros(Nfft/4,L) ];
G_SCFDMA            = [ zeros(Nfft/4,L);   G_SCFDMA;            zeros(Nfft/4,L) ];
G_pDFTsOFDM         = [ zeros(Nfft/4,L/2); G_pDFTsOFDM;         zeros(Nfft/4,L/2) ];
G_pDFTsFBMCunscaled = [ zeros(Nfft/4,L/2); G_pDFTsFBMCunscaled; zeros(Nfft/4,L/2) ];
G_pDFTsFBMC         = [ zeros(Nfft/4,L/2); G_pDFTsFBMC;         zeros(Nfft/4,L/2) ];


%% Plot Results
t = ( -Nfft/4+1 : Nfft+Nfft/4 ) / Nfft;

figure(2);
% OFDM
subplot(5,1,1);
plot(t, abs(G_OFDM(:,1)).^2 ); hold on;   
plot(t, sum(abs(G_OFDM).^2,2) , 'black' ); 
ylabel('Power');
ylim([0 1.1]); xlim([-0.2 1.2]);

% SC-FDMA
subplot(5,1,2);
plot(t, abs(G_SCFDMA).^2 );
ylabel('Power');
ylim([0 1.1]); xlim([-0.2 1.2]);

% Pruned DFT spread OFDM
subplot(5,1,3);
plot(t, abs(G_pDFTsOFDM).^2 );
ylabel('Power');
ylim([0 1.1]); xlim([-0.2 1.2]);

% Pruned DFT spread FBMC (unscaled)
subplot(5,1,4);
plot(t, abs(G_pDFTsFBMCunscaled).^2 ); hold on;
plot(t,abs( [ zeros(Nfft/4,1); RRCwindow ; zeros(Nfft/4,1)] ).^2,'black');
ylabel('Power');
ylim([0 1.1]); xlim([-0.2 1.2]);

% ################## Pruned DFT spread FBMC #########################
subplot(5,1,5);
plot(t, abs(G_pDFTsFBMC).^2 );
ylabel('Power');
xlabel('Normalized Time, tF');
ylim([0 1.1]); xlim([-0.2 1.2]);
% ###################################################################


%% Calculate the SIR for Pruned DFT spread FBMC 
SignalPower       = sum(diag(G_pDFTsFBMC'*G_pDFTsFBMC).^2);
InterferencePower = sum(sum(abs(G_pDFTsFBMC'*G_pDFTsFBMC).^2)) - SignalPower;

SIR_dB = 10*log10(SignalPower./InterferencePower);
disp(['Signal-to-Interference Ratio: ' int2str(SIR_dB) 'dB']);




