% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================   
% This script simulates the Bit Error Ratio (BER) of pruned DFT spread
% FBMC, SC-FDMA, OFDM and FBMC. 
% The Bit Error Probability (BEP) is also calculated if "CalculateTheory"
% is set to true.

% Allows to reproduce Figure 10 11 12 of "Pruned DFT Spread FBMC: Low PAPR,
% Low Latency, High Spectral Efficiency", R. Nissel and M. Rupp, IEEE 
% Transactions on Communications

clear; close all;
addpath('./Theory');

%% Parameters
% Simulation
M_SNR_dB            = [0:5:30];                     % Signal-to-Noise Ratio in dB
NrRepetitions       = 30;                           % Number of Monte Carlo repetitions
CalculateTheory     = false;                        % If set to true, calculate the BEP and the SINR. Set to "false" because the theoretical calculations require large matrix multiplications, increasing the simulation/calculation time significantly.

% FBMC and OFDM Parameters
NrSubcarriers       = 256;                          % Number of subcarriers
QAM_ModulationOrder = 16;                           % Modulation order, 4, 16, 64,...
SubcarrierSpacing   = 15e3;                         % Subcarrier spacing (15kHz, same as LTE)
CarrierFrequency    = 2.5e9;                        % Carrier Frequency
K_FBMC              = 30;                           % Number of FBMC symbols in time
K_OFDMnoCP          = 15;                           % Number of OFDM symbols in time (no CP)
K_OFDM              = 14;                           % Number of OFDM symbols in time (same as in LTE)
CP_Length           = 1/SubcarrierSpacing/14;       % LTE CP Length in seconds
CP_Length_FBMC_DFT  = 0;                            % CP in the frequency domain for the DFT spreading aproach. Multiple of two: 0, 2, 4... Can usually be set to zero

SamplingRate        = 15e3*14*12*2;                 % Sampling rate, should approximatly match the power-delay profile of the channel. "*14" due to the CP

% Channel
PowerDelayProfile   = 'TDL-A_300ns';                % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
Velocity_kmh        = 200;                          % Velocity in km/h


% #########################################################################
% % In the paper:
% M_SNR_dB        = [0:1.5:30];
% SamplingRate    = 15e3*14*12*8; 
% NrRepetitions   = 1000;
% CalculateTheory = true;
% #########################################################################


% ######################### For Figure 10 #################################
% CalculateTheory = true; 
% #########################################################################


% ######################### For Figure 11 #################################
% Velocity_kmh = 0; 
% #########################################################################


%% FBMC Object
FBMC = Modulation.FBMC(...
    NrSubcarriers,...                               % Number of subcarriers
    K_FBMC,...                                      % Number of FBMC symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    'Hermite-OQAM',...                              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    4, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                          % Initial phase shift
    true ...                                        % Polyphase implementation
    );
% The only difference between DFT_FBMC and FBMC is the prototype filter, which is slightly reduced in DFT_FBMC (improves the SIR a litte bit and reduces the complexity)
FBMC_DFT = Modulation.FBMC(...          
    NrSubcarriers,...                               % Number of subcarriers
    K_FBMC,...                                      % Number of FBMC symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    'HermiteCut-OQAM',...                           % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman
    4, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                          % Initial phase shift
    true ...                                        % Polyphase implementation
    );


%% OFDM Object
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+0*SamplingRate)*K_OFDMnoCP)/2)/SamplingRate;
OFDMnoCP = Modulation.OFDM(...
    NrSubcarriers,...                               % Number of  active subcarriers
    K_OFDMnoCP,...                                  % Number of OFDM Symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    0, ...                                          % Cyclic prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                         % Zero guard length (s)
    );
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round(SamplingRate/SubcarrierSpacing)+CP_Length*SamplingRate)*K_OFDM)/2)/SamplingRate;
OFDM = Modulation.OFDM(...
    NrSubcarriers,...                               % Number of  active subcarriers
    K_OFDM,...                                      % Number of OFDM Symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    CP_Length, ...                                  % Cyclic prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                         % Zero guard length (s)
    );

%% Check Number of Samples
if  OFDM.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal || OFDMnoCP.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = OFDM.Nr.SamplesTotal;

%% Modulation Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');


%% Channel Object
ChannelModel = Channel.FastFading(...
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    PowerDelayProfile,...                                               % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    N,...                                                               % Number of total samples
    Velocity_kmh/3.6*CarrierFrequency/2.998e8,...                       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                                                         % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    1,...                                                               % Number of transmit antennas
    1,...                                                               % Number of receive antennas
    true ...                                                            % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

%% DFT Matrix
DFTMatrix = fft(eye(NrSubcarriers))/sqrt(NrSubcarriers);

%% Generate coding matrix for pruned DFT spread FBMC
TrueNrMCSymbols = FBMC_DFT.Nr.MCSymbols;
FBMC_DFT.SetNrMCSymbols(1);
D_temp = FBMC_DFT.GetFBMCMatrix;
FBMC_DFT.SetNrMCSymbols(TrueNrMCSymbols);

% Note that, if CP_Length==0, then T_CP and R_CP are identity matrices
T_CP                                                        = zeros(NrSubcarriers,NrSubcarriers-CP_Length_FBMC_DFT);
T_CP(1:CP_Length_FBMC_DFT/2,end-CP_Length_FBMC_DFT/2+1:end) = eye(CP_Length_FBMC_DFT/2);
T_CP(CP_Length_FBMC_DFT/2+1:end-CP_Length_FBMC_DFT/2,:)     = eye(NrSubcarriers-CP_Length_FBMC_DFT);
T_CP(end-CP_Length_FBMC_DFT/2+1:end,1:CP_Length_FBMC_DFT/2) = eye(CP_Length_FBMC_DFT/2);
R_CP                                                        = zeros(NrSubcarriers,NrSubcarriers-CP_Length_FBMC_DFT);
R_CP(CP_Length_FBMC_DFT/2+1:end-CP_Length_FBMC_DFT/2,:)     = eye(NrSubcarriers-CP_Length_FBMC_DFT);

% DFT matrix for the coding process
W = fft( eye(NrSubcarriers-CP_Length_FBMC_DFT) ) / sqrt( NrSubcarriers-CP_Length_FBMC_DFT );

% Diagonal elements of the FBMC transmission matrix after DFT spreading despreading
a = abs(diag(W'*R_CP'*D_temp*T_CP*W));
a = a+randn(size(a))*10^-12; %  randn so that sorting is unique

% Sort a
a_Tilde = sort(a,'descend');

% Get index representing the largest values of a
alpha       = a_Tilde((NrSubcarriers-CP_Length_FBMC_DFT)/2);
Index_Tilde = (a>=alpha);

% Pruned DFT matrix
W_Tilde = W(:,Index_Tilde) ;

% One-tap scaling of the data symbols
b_Tilde = sqrt(1./(a(Index_Tilde)));

% Final coding matrix for one FBMC symbol
Cf_DFTspread_TX = T_CP*W_Tilde*diag(b_Tilde);
Cf_DFTspread_RX = R_CP*W_Tilde*diag(b_Tilde);

C_DFTspread_TX = kron(sparse(eye(K_FBMC)),Cf_DFTspread_TX);
C_DFTspread_RX = kron(sparse(eye(K_FBMC)),Cf_DFTspread_RX);



%% Get OFDM and FBMC Transmit and Receive Matrices
GTX_OFDM = sparse(OFDM.GetTXMatrix);
GRX_OFDM = sparse(OFDM.GetRXMatrix');

GTX_OFDMnoCP = sparse(OFDMnoCP.GetTXMatrix);
GRX_OFDMnoCP = sparse(OFDMnoCP.GetRXMatrix');

G_FBMC = sparse(FBMC.GetTXMatrix);

G_FBMC_DFT = sparse(FBMC_DFT.GetTXMatrix);


%% Normalize OFDM and FBMC (the default matrices are normalized to have unit transmit power for unit power data symbols)
NormalizationOFDM     = sqrt((GRX_OFDM(:,1)'*GRX_OFDM(:,1)));
NormalizationOFDMnoCP = sqrt((GRX_OFDMnoCP(:,1)'*GRX_OFDMnoCP(:,1)));
NormalizationFBMC     = 1/sqrt((G_FBMC(:,1)'*G_FBMC(:,1)));
NormalizationFBMC_DFT = 1/sqrt((G_FBMC_DFT(:,1)'*G_FBMC_DFT(:,1)));

GTX_OFDM     = GTX_OFDM*NormalizationOFDM;
GRX_OFDM     = GRX_OFDM/NormalizationOFDM;

GTX_OFDMnoCP = GTX_OFDMnoCP*NormalizationOFDMnoCP;
GRX_OFDMnoCP = GRX_OFDMnoCP/NormalizationOFDMnoCP;

G_FBMC       = G_FBMC*NormalizationFBMC;

G_FBMC_DFT   = G_FBMC_DFT*NormalizationFBMC;


%% Preallocate for parfor
BER_OFDM         = nan(length(M_SNR_dB),NrRepetitions);
BER_OFDMnoCP     = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC         = nan(length(M_SNR_dB),NrRepetitions);

BER_DFT_OFDM     = nan(length(M_SNR_dB),NrRepetitions);
BER_DFT_OFDMnoCP = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_DFT     = nan(length(M_SNR_dB),NrRepetitions);

SINR_DFT_OFDM_dB     = nan(NrSubcarriers,K_OFDM,length(M_SNR_dB),NrRepetitions);
SINR_DFT_OFDMnoCP_dB = nan(NrSubcarriers,K_OFDMnoCP,length(M_SNR_dB),NrRepetitions);
SINR_FBMC_DFT_dB     = nan((NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC,length(M_SNR_dB),NrRepetitions);

SINR_Approx_DFT_OFDM_dB     = nan(NrSubcarriers,K_OFDM,length(M_SNR_dB),NrRepetitions);
SINR_Approx_DFT_OFDMnoCP_dB = nan(NrSubcarriers,K_OFDMnoCP,length(M_SNR_dB),NrRepetitions);
SINR_Approx_FBMC_DFT_dB     = nan((NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC,length(M_SNR_dB),NrRepetitions);


%% Start Simulation (Calculation)
tic;
for i_rep = 1:NrRepetitions   
    % Generate Bit Stream
    BinaryDataStream          = randi([0 1],K_OFDMnoCP*NrSubcarriers*log2(QAM.ModulationOrder),1);
    BinaryDataStream_OFDM     = randi([0 1],K_OFDM*NrSubcarriers*log2(QAM.ModulationOrder),1);          % Reduced bit rate due to the CP
    BinaryDataStream_FBMC_DFT = randi([0 1],K_OFDMnoCP*(NrSubcarriers-CP_Length_FBMC_DFT)*log2(QAM.ModulationOrder),1);  % Maybe a reduced bit rate due to the frequency CP (but not necesarry most of the time!)

    % Map Bit Stream to Symbols
    x_OFDM     = reshape(QAM.Bit2Symbol(BinaryDataStream_OFDM),NrSubcarriers,K_OFDM);
    x_OFDMnoCP = reshape(QAM.Bit2Symbol(BinaryDataStream),NrSubcarriers,K_OFDMnoCP);
    x_FBMC_DFT = reshape(QAM.Bit2Symbol(BinaryDataStream_FBMC_DFT),(NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC);
    x_FBMC     = reshape(PAM.Bit2Symbol(BinaryDataStream),NrSubcarriers,K_FBMC)/sqrt(2);  % 1/sqrt(2) => same TX power for OFDM and FBMC

    % Generate Transmit Signal in the Time Domain
    s_OFDM         = OFDM.Modulation(x_OFDM)*NormalizationOFDM;
    s_OFDMnoCP     = OFDMnoCP.Modulation(x_OFDMnoCP)*NormalizationOFDMnoCP;
    s_FBMC         = FBMC.Modulation(x_FBMC)*NormalizationFBMC;

    s_DFT_OFDM     = OFDM.Modulation(DFTMatrix*x_OFDM)*NormalizationOFDM;
    s_DFT_OFDMnoCP = OFDMnoCP.Modulation(DFTMatrix*x_OFDMnoCP)*NormalizationOFDMnoCP;
    s_FBMC_DFT     = FBMC_DFT.Modulation(Cf_DFTspread_TX*x_FBMC_DFT)*NormalizationFBMC_DFT;

    % Channel
    ChannelModel.NewRealization;
    H = ChannelModel.GetConvolutionMatrix{1};
    
    r_OFDM_noNoise         = H * s_OFDM;
    r_OFDMnoCP_noNoise     = H * s_OFDMnoCP;
    r_FBMC_noNoise         = H * s_FBMC;

    r_DFT_OFDM_noNoise     = H * s_DFT_OFDM;
    r_DFT_OFDMnoCP_noNoise = H * s_DFT_OFDMnoCP;
    r_FBMC_DFT_noNoise     = H * s_FBMC_DFT;
    
    noise_unitPower = sqrt(1/2)*( randn(N,1) + 1j * randn(N,1) );

    % Calculate One Tap Channels
    % Note that diag(G'*H*G)==sum((G'*H).*G.',2)
    h_OFDM      = full( reshape( sum((GRX_OFDM'*H).*GTX_OFDM.',2), NrSubcarriers, [] ) ); 
    h_OFDMnoCP  = full( reshape( sum((GRX_OFDMnoCP'*H).*GTX_OFDMnoCP.',2), NrSubcarriers, [] ) ); 
    h_FBMC      = full( reshape( sum((G_FBMC'*H).*G_FBMC.',2), NrSubcarriers, [] ) ); 
    h_FBMC_DFT  = full( reshape( sum((G_FBMC_DFT'*H).*G_FBMC_DFT.',2), NrSubcarriers, [] ) ); 

    % Precalculate Stuff to Improve the Computation Time
    if CalculateTheory
        Precalc_OFDM     = GRX_OFDM' * H *GTX_OFDM * kron(sparse(eye(K_OFDM)),DFTMatrix);
        Precalc_OFDMnoCP = GRX_OFDMnoCP' * H *GTX_OFDMnoCP * kron(sparse(eye(K_OFDMnoCP)),DFTMatrix);
        Precalc_FBMC     = G_FBMC_DFT' * H * G_FBMC_DFT * C_DFTspread_TX;
    end

    for i_SNR = 1:length(M_SNR_dB)
        SNR_dB  = M_SNR_dB(i_SNR); 
        Pn      = 10^(-SNR_dB/10); 

        %% Simulation
        % Add Noise
        noise = sqrt(Pn) * noise_unitPower;

        r_OFDM         = r_OFDM_noNoise     + noise;
        r_OFDMnoCP     = r_OFDMnoCP_noNoise + noise;
        r_FBMC         = r_FBMC_noNoise     + noise;

        r_DFT_OFDM     = r_DFT_OFDM_noNoise     + noise;
        r_DFT_OFDMnoCP = r_DFT_OFDMnoCP_noNoise + noise;
        r_FBMC_DFT     = r_FBMC_DFT_noNoise     + noise;

        % Received Symbols (Demodulation)
        y_OFDM         = OFDM.Demodulation(r_OFDM) / NormalizationOFDM;
        y_OFDMnoCP     = OFDMnoCP.Demodulation(r_OFDMnoCP) / NormalizationOFDMnoCP;
        y_FBMC         = FBMC.Demodulation(r_FBMC) / NormalizationFBMC;

        y_DFT_OFDM     = OFDM.Demodulation(r_DFT_OFDM) / NormalizationOFDM;
        y_DFT_OFDMnoCP = OFDMnoCP.Demodulation(r_DFT_OFDMnoCP) / NormalizationOFDMnoCP;
        y_FBMC_DFT     = FBMC_DFT.Demodulation(r_FBMC_DFT) / NormalizationFBMC_DFT;

        % ZF Equalizer for OFDM and FBMC
        x_est_OFDM      = y_OFDM ./ h_OFDM;
        x_est_OFDMnoCP  = y_OFDMnoCP ./ h_OFDMnoCP;
        x_est_FBMC      = real( y_FBMC ./ h_FBMC ) * sqrt(2);

        % One-tap (scaled) MMSE followed by despreading
        Scaling_DFT_OFDM     =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_OFDM ).^2 ),1)), NrSubcarriers,1);
        Scaling_DFT_OFDMnoCP =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_OFDMnoCP ).^2 ),1)), NrSubcarriers,1);
        Scaling_FBMC_DFT     =  repmat(1./(mean(1 ./( 1 + Pn./abs( h_FBMC(CP_Length_FBMC_DFT/2+1:end-CP_Length_FBMC_DFT/2,:) ).^2 ),1)), NrSubcarriers,1);
        
        e_DFT_OFDM     = Scaling_DFT_OFDM     .*conj(h_OFDM)    ./( abs(h_OFDM).^2     + Pn );
        e_DFT_OFDMnoCP = Scaling_DFT_OFDMnoCP .*conj(h_OFDMnoCP)./( abs(h_OFDMnoCP).^2 + Pn );
        e_FBMC_DFT     = Scaling_FBMC_DFT     .*conj(h_FBMC)    ./( abs(h_FBMC).^2     + Pn );
        
        x_est_DFT_OFDM     = DFTMatrix'       * (y_DFT_OFDM     .* e_DFT_OFDM);
        x_est_DFT_OFDMnoCP = DFTMatrix'       * (y_DFT_OFDMnoCP .* e_DFT_OFDMnoCP);
        x_est_FBMC_DFT     = Cf_DFTspread_RX' * (y_FBMC_DFT     .* e_FBMC_DFT);
        
        % Symbol to Bit
        DetectedBitStream_OFDM         = QAM.Symbol2Bit(x_est_OFDM);
        DetectedBitStream_OFDMnoCP     = QAM.Symbol2Bit(x_est_OFDMnoCP);
        DetectedBitStream_FBMC         = PAM.Symbol2Bit(x_est_FBMC);

        DetectedBitStream_DFT_OFDM     = QAM.Symbol2Bit(x_est_DFT_OFDM);
        DetectedBitStream_DFT_OFDMnoCP = QAM.Symbol2Bit(x_est_DFT_OFDMnoCP);
        DetectedBitStream_FBMC_DFT     = QAM.Symbol2Bit(x_est_FBMC_DFT);

        % Bit Error Ratio
        BER_OFDM(i_SNR,i_rep)         = mean( BinaryDataStream_OFDM ~= DetectedBitStream_OFDM );
        BER_OFDMnoCP(i_SNR,i_rep)     = mean( BinaryDataStream      ~= DetectedBitStream_OFDMnoCP );
        BER_FBMC(i_SNR,i_rep)         = mean( BinaryDataStream      ~= DetectedBitStream_FBMC );

        BER_DFT_OFDM(i_SNR,i_rep)     = mean( BinaryDataStream_OFDM     ~= DetectedBitStream_DFT_OFDM );
        BER_DFT_OFDMnoCP(i_SNR,i_rep) = mean( BinaryDataStream          ~= DetectedBitStream_DFT_OFDMnoCP );
        BER_FBMC_DFT(i_SNR,i_rep)     = mean( BinaryDataStream_FBMC_DFT ~= DetectedBitStream_FBMC_DFT );
        
        %% Theoretical calculations for pruned DFT spread FBMC and SC-FDMA
        if CalculateTheory
            % Equalizer
            E_DFT_OFDM      =  sparse( diag( e_DFT_OFDM(:)     ) );
            E_DFT_OFDMnoCP  =  sparse( diag( e_DFT_OFDMnoCP(:) ) );
            E_FBMC_DFT      =  sparse( diag( e_FBMC_DFT(:)     ) );

            % Full transmission matrix
            D_DFT_OFDM     = kron(sparse(eye(K_OFDM)),DFTMatrix)' * E_DFT_OFDM * Precalc_OFDM;
            Gamma_DFT_OFDM = kron(sparse(eye(K_OFDM)),DFTMatrix)' * E_DFT_OFDM * GRX_OFDM';

            D_DFT_OFDMnoCP     = kron(sparse(eye(K_OFDMnoCP)),DFTMatrix)' * E_DFT_OFDMnoCP * Precalc_OFDMnoCP;
            Gamma_DFT_OFDMnoCP = kron(sparse(eye(K_OFDMnoCP)),DFTMatrix)' * E_DFT_OFDMnoCP * GRX_OFDMnoCP';    

            D_FBMC_DFT     = C_DFTspread_RX' * E_FBMC_DFT * Precalc_FBMC;
            Gamma_FBMC_DFT = C_DFTspread_RX' * E_FBMC_DFT * G_FBMC_DFT';

            % Calculate the approximated SINR
            SINR_Approx_DFT_OFDM_dB(:,:,i_SNR,i_rep)     = 10*log10(repmat(1./(1./mean(1./(1+Pn./abs(h_OFDM).^2)) - 1),NrSubcarriers,1));
            SINR_Approx_DFT_OFDMnoCP_dB(:,:,i_SNR,i_rep) = 10*log10(repmat(1./(1./mean(1./(1+Pn./abs(h_OFDMnoCP).^2)) - 1),NrSubcarriers,1));
            SINR_Approx_FBMC_DFT_dB(:,:,i_SNR,i_rep)     = 10*log10(repmat(1./(1./mean(1./(1+Pn./abs(h_FBMC_DFT(CP_Length_FBMC_DFT/2+1:end-CP_Length_FBMC_DFT/2,:)).^2)) - 1),(NrSubcarriers-CP_Length_FBMC_DFT)/2 , 1 ));

            % Calculate the true SINR
            SINR_DFT_OFDM_dB(:,:,i_SNR,i_rep)     = 10*log10(1./(reshape( sum(abs(D_DFT_OFDM-eye(size(D_DFT_OFDM))).^2,2) + sum(abs(Gamma_DFT_OFDM).^2,2)*Pn , size(x_OFDM))));
            SINR_DFT_OFDMnoCP_dB(:,:,i_SNR,i_rep) = 10*log10(1./(reshape( sum(abs(D_DFT_OFDMnoCP-eye(size(D_DFT_OFDMnoCP))).^2,2) + sum(abs(Gamma_DFT_OFDMnoCP).^2,2)*Pn , size(x_OFDMnoCP))));
            SINR_FBMC_DFT_dB(:,:,i_SNR,i_rep)     = 10*log10(1./(reshape( sum(abs(D_FBMC_DFT-eye(size(D_FBMC_DFT))).^2,2) + sum(abs(Gamma_FBMC_DFT).^2,2)*Pn , size(x_FBMC_DFT))));
        end
    end
    TimePassed = toc;
    if mod(i_rep,1)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes = ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60/60) 'hours']);
    end
end


if CalculateTheory
    %% Calculate the Bit Error Probability (BEP)
    M_SNR_dB_morePoints = -30:0.1:50;
    ReferenceBitErrorProability = BitErrorProbabilityAWGN( M_SNR_dB_morePoints, QAM.SymbolMapping , QAM.BitMapping );

    % Use interpolation because it is faster than directly using the function "BitErrorProbabilityAWGN(.)"
    [a,b,c,d] = size(SINR_DFT_OFDM_dB);
    BEP_DFT_OFDM_all = reshape( interp1( M_SNR_dB_morePoints , ReferenceBitErrorProability , SINR_DFT_OFDM_dB(:) ,'spline'), a,b,c,d);
    BEP_DFT_OFDM = squeeze(mean(mean(mean(BEP_DFT_OFDM_all,1),2),4));

    [a,b,c,d] = size(SINR_DFT_OFDMnoCP_dB);
    BEP_DFT_OFDMnoCP_all = reshape( interp1( M_SNR_dB_morePoints , ReferenceBitErrorProability , SINR_DFT_OFDMnoCP_dB(:) ,'spline'), a,b,c,d);
    BEP_DFT_OFDMnoCP = squeeze(mean(mean(mean(BEP_DFT_OFDMnoCP_all,1),2),4));

    [a,b,c,d] = size(SINR_FBMC_DFT_dB);
    BEP_FBMC_DFT_all = reshape( interp1( M_SNR_dB_morePoints , ReferenceBitErrorProability , SINR_FBMC_DFT_dB(:) ,'spline'), a,b,c,d);
    BEP_FBMC_DFT = squeeze(mean(mean(mean(BEP_FBMC_DFT_all,1),2),4));

    
    %% Error between the true SINR and the approximated SINR
    Error_SINR_DFT_OFDM_dB     = 10*log10(1./(squeeze(mean(mean(mean(10.^(-SINR_DFT_OFDM_dB./10),1),2),4)))) - 10*log10(1./(squeeze(mean(mean(mean(10.^(-SINR_Approx_DFT_OFDM_dB./10),1),2),4))));
    Error_SINR_DFT_OFDMnoCP_dB = 10*log10(1./(squeeze(mean(mean(mean(10.^(-SINR_DFT_OFDMnoCP_dB./10),1),2),4)))) - 10*log10(1./(squeeze(mean(mean(mean(10.^(-SINR_Approx_DFT_OFDMnoCP_dB./10),1),2),4))));
    Error_SINR_FBMC_DFT_dB     = 10*log10(1./(squeeze(mean(mean(mean(10.^(-SINR_FBMC_DFT_dB./10),1),2),4)))) - 10*log10(1./(squeeze(mean(mean(mean(10.^(-SINR_Approx_FBMC_DFT_dB./10),1),2),4))));
end


if CalculateTheory
    %% Plot the Approximation Error
    figure(10);
    plot(M_SNR_dB,-Error_SINR_DFT_OFDM_dB,'black');
    hold on;
    plot(M_SNR_dB,-Error_SINR_DFT_OFDMnoCP_dB,'red');
    plot(M_SNR_dB,-Error_SINR_FBMC_DFT_dB,'blue');
    ylim([-2 8]);
    xlabel('Signal-to-Noise Ratio');
    ylabel('SINRapprox - SINR [dB]');
    legend({'SC-FDMA','SC-FDMD (noCP)','Pruned DFT-s FBMC'})   
    
    %% Plot the BEP and BER
    figure(12);
    markersize =4;      
    semilogy(M_SNR_dB,BEP_DFT_OFDM,'- black','Markersize',markersize);
    hold on;
    semilogy(M_SNR_dB,BEP_DFT_OFDMnoCP,'- red','Markersize',markersize);
    semilogy(M_SNR_dB,BEP_FBMC_DFT,'- blue','Markersize',markersize);
    semilogy(M_SNR_dB,mean(BER_DFT_OFDM,2),'x black','Markersize',markersize);
    semilogy(M_SNR_dB,mean(BER_DFT_OFDMnoCP,2),'s red','Markersize',markersize);
    semilogy(M_SNR_dB,mean(BER_FBMC_DFT,2),'o blue','Markersize',markersize);
    ylim([10^-4 1]);
    xlabel('Signal-to-Noise Ratio');
    ylabel('Bit Error Probability, Bit Error Ratio');    
    legend({'SC-FDMA','SC-FDMA (noCP)','Pruned DFT-s FBMC'})       
else
    figure(12);
    markersize =4;      
    semilogy(M_SNR_dB,mean(BER_DFT_OFDM,2),'-x black','Markersize',markersize);
    hold on;
    semilogy(M_SNR_dB,mean(BER_DFT_OFDMnoCP,2),'-s red','Markersize',markersize);
    semilogy(M_SNR_dB,mean(BER_FBMC_DFT,2),'-o blue','Markersize',markersize);
    ylim([10^-4 1]);
    xlabel('Signal-to-Noise Ratio');
    ylabel('Bit Error Ratio');    
    legend({'SC-FDMA','SC-FDMA (noCP)','Pruned DFT-s FBMC'})            
end


%% Plot BER, including multicarrier modulation
figure(100)
semilogy(M_SNR_dB,mean(BER_DFT_OFDM,2),'black');
hold on;
semilogy(M_SNR_dB,mean(BER_DFT_OFDMnoCP,2),'red');
semilogy(M_SNR_dB,mean(BER_FBMC_DFT,2),'blue');
semilogy(M_SNR_dB,mean(BER_OFDM,2),'-- black');
semilogy(M_SNR_dB,mean(BER_OFDMnoCP,2),'-- red');
semilogy(M_SNR_dB,mean(BER_FBMC,2),'-- blue');
xlabel('Signal-to-Noise Ratio');
ylabel('Bit Error Ratio');   
legend({'SC-FDMA','SC-FDMD (noCP)','Pruned DFT-s FBMC','OFDM','OFDM (noCP)','FBMC'})      

figure(12);

%% Show the Bitrates
BitRate_OFDM     = length(BinaryDataStream_OFDM)/(OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols);
BitRate_OFDMnoCP = length(BinaryDataStream)/(OFDMnoCP.PHY.TimeSpacing*OFDMnoCP.Nr.MCSymbols);
BitRate_FBMC     = length(BinaryDataStream)/(FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols);
BitRate_FBMC_DFT = length(BinaryDataStream_FBMC_DFT)/(FBMC_DFT.PHY.TimeSpacing*FBMC_DFT.Nr.MCSymbols);
fprintf('=================================================================\n');
fprintf('                  |  CP-OFDM | OFDM(noCP) |    FBMC  | pDFTsFBMC |\n');
fprintf('Bit Rate [Bits/s] |%10.0f|%12.0f|%10.0f|%11.0f|\n',BitRate_OFDM,BitRate_OFDMnoCP,BitRate_FBMC,BitRate_FBMC_DFT);
fprintf('=================================================================\n');


%% Save Results
SaveStuff = false;
if SaveStuff 
    Name = ['.\Results\SISO_BER_' PowerDelayProfile '_v' int2str(Velocity_kmh) '_' int2str(QAM_ModulationOrder) '_' int2str(NrSubcarriers) '.mat'];
    meanBER_DFT_OFDM     = mean(BER_DFT_OFDM,2);
    meanBER_DFT_OFDMnoCP = mean(BER_DFT_OFDMnoCP,2);
    meanBER_FBMC_DFT     = mean(BER_FBMC_DFT,2);
    
    meanBEP_DFT_OFDM     = mean(BEP_DFT_OFDM,2);
    meanBEP_DFT_OFDMnoCP = mean(BEP_DFT_OFDMnoCP,2);
    meanBEP_FBMC_DFT     = mean(BEP_FBMC_DFT,2); 
    
    meanBER_OFDM     = mean(BER_OFDM,2);
    meanBER_OFDMnoCP = mean(BER_OFDMnoCP,2);
    meanBER_FBMC     = mean(BER_FBMC,2);    
 
    save(Name,...
        'meanBER_DFT_OFDM',...
        'meanBER_DFT_OFDMnoCP',...
        'meanBER_FBMC_DFT',...        
        'meanBER_OFDM',...
        'meanBER_OFDMnoCP',...
        'meanBER_FBMC',...        
        'M_SNR_dB',...
        'NrRepetitions',...
        'meanBEP_DFT_OFDM',...
        'meanBEP_DFT_OFDMnoCP',...
        'meanBEP_FBMC_DFT',...
        'Error_SINR_DFT_OFDM_dB',...
        'Error_SINR_DFT_OFDMnoCP_dB',...
        'Error_SINR_FBMC_DFT_dB' ...        
        );
end


