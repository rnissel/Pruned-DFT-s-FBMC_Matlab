% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================      
% This script shows the simulated peak-to-average-power ratio (PAPR) and
% additionally the simulated power in time and in frequency. 

% Reproduces Figure 9 of "Pruned DFT Spread FBMC: Low PAPR, Low Latency, 
% High Spectral Efficiency", R. Nissel and M. Rupp, IEEE Transactions on 
% Communications 

clear; close all;

%% Parameters
% Simulation
NrRepetitions           = 2000;                         % Number of Monte Carlo repetitions.

% FBMC and OFDM Parameters
NrSubcarriers           = 256;                          % Number of subcarriers
QAM_ModulationOrder     = 4;                            % Modulation order, 4, 16, 64,...
SubcarrierSpacing       = 15e3;                         % Subcarrier spacing (15kHz, same as LTE)
CarrierFrequency        = 2.5e9;                        % Carrier Frequency
K_FBMC                  = 30;                           % Number of FBMC symbols in time
K_OFDMnoCP              = 15;                           % Number of OFDM symbols in time (no CP)
K_OFDM                  = 14;                           % Number of OFDM symbols in time (same as in LTE)
CP_Length               = 1/SubcarrierSpacing/14;       % LTE CP Length in seconds

CP_Length_FBMC_DFT      = 0;                            % Frequency domain CP for pruned DFT spread FBMC. Multiple of two: 0, 2, 4... Can usually be set to zero
CP_Length_FBMC_DFT_LCP  = 26;                           % Same as "CP_Length_FBMC_DFT" but different value to check the effect of the frequency CP on the PAPR
CP_Length_FFT_FBMC      = 0;                            % "Time domain" CP for the FFT-FBMC scheme. 

SamplingRate            = SubcarrierSpacing*14*12*2;    % Sampling rate, should approximatly match the power-delay profile of the channel. "14", so that the CP fits the sampling rate

PseudoOverlappingFactor = 4;                            % Pseudo overlapping factor to keep the implementation simple. The true overlapping factor is 1.56. 


% #########################################################################
% % In the paper:
% SamplingRate  = 15e3*14*12*8*2; 
% NrRepetitions = 200000;
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
    PseudoOverlappingFactor, ...                    % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                          % Initial phase shift
    true ...                                        % Polyphase implementation
    );
FBMC_PhaseCondition = Modulation.FBMCphaseCondition(...
    NrSubcarriers,...                               % Number of subcarriers
    K_FBMC,...                                      % Number of FBMC symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    'Hermite-OQAM',...                              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    PseudoOverlappingFactor, ...                    % Overlapping factor (also determines oversampling in the frequency domain)
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
    PseudoOverlappingFactor, ...                    % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                          % Initial phase shift
    true ...                                        % Polyphase implementation
    );
FBMC_TimeSpread = Modulation.FBMC(...          
    NrSubcarriers,...                               % Number of subcarriers
    K_FBMC,...                                      % Number of FBMC symbols
    SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                                % Sampling rate (Samples/s)
    0,...                                           % Intermediate frequency first subcarrier (Hz)
    false,...                                       % Transmit real valued signal
    'PHYDYAS-OQAM',...                              % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman
    PseudoOverlappingFactor, ...                    % Overlapping factor (also determines oversampling in the frequency domain)
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


%% DFT Matrix
DFTMatrix               = fft(eye(NrSubcarriers))/sqrt(NrSubcarriers);

%% Generate coding matrix for pruned DFT spread FBMC
TrueNrMCSymbols     = FBMC_DFT.Nr.MCSymbols;
FBMC_DFT.SetNrMCSymbols(1);
D_temp              = FBMC_DFT.GetFBMCMatrix;
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

% Reduced DFT matrix
W_Tilde = W(:,Index_Tilde) ;

% One-tap scaling of the data symbols
b_Tilde = sqrt(2./(a(Index_Tilde)));

% Final coding matrix for one FBMC symbol
C_DFTspread_TX = T_CP*W_Tilde*diag(b_Tilde);
C_DFTspread_RX = R_CP*W_Tilde*diag(b_Tilde);


%% Larger Frequency CP for the PAPR Plot (essentially the same as above)
T_CP2                                                                = zeros(NrSubcarriers,NrSubcarriers-CP_Length_FBMC_DFT_LCP);
T_CP2(1:CP_Length_FBMC_DFT_LCP/2,end-CP_Length_FBMC_DFT_LCP/2+1:end) = eye(CP_Length_FBMC_DFT_LCP/2);
T_CP2(CP_Length_FBMC_DFT_LCP/2+1:end-CP_Length_FBMC_DFT_LCP/2,:)     = eye(NrSubcarriers-CP_Length_FBMC_DFT_LCP);
T_CP2(end-CP_Length_FBMC_DFT_LCP/2+1:end,1:CP_Length_FBMC_DFT_LCP/2) = eye(CP_Length_FBMC_DFT_LCP/2);
R_CP2                                                                            = zeros(NrSubcarriers,NrSubcarriers-CP_Length_FBMC_DFT_LCP);
R_CP2(CP_Length_FBMC_DFT_LCP/2+1:end-CP_Length_FBMC_DFT_LCP/2,:)     = eye(NrSubcarriers-CP_Length_FBMC_DFT_LCP);

W2              = fft( eye(NrSubcarriers-CP_Length_FBMC_DFT_LCP) ) / sqrt( NrSubcarriers-CP_Length_FBMC_DFT_LCP );
a2              = abs(diag(W2'*R_CP2'*D_temp*T_CP2*W2));
a2              = a2+randn(size(a2))*10^-12; %  randn so that sorting is unique
a_Tilde2        = sort(a2,'descend');
alpha2          = a_Tilde2((NrSubcarriers-CP_Length_FBMC_DFT_LCP)/2);
Index_Tilde2    = (a2>=alpha2);
W_Tilde2        = W2(:,Index_Tilde2) ;
b_Tilde2        = sqrt(2./(a2(Index_Tilde2)));
C_DFTspread_TX2 = T_CP2*W_Tilde2*diag(b_Tilde2);
C_DFTspread_RX2 = R_CP2*W_Tilde2*diag(b_Tilde2);


%% Spreading in Time, FFT-FBMC
TrueNrSubcarriers   = FBMC_TimeSpread.Nr.Subcarriers;
FBMC_TimeSpread.SetNrSubcarriers(1);
D_temp              = FBMC_TimeSpread.GetFBMCMatrix;
FBMC_TimeSpread.SetNrSubcarriers(TrueNrSubcarriers);

T_CP3                                                        = zeros(K_FBMC,K_FBMC-CP_Length_FFT_FBMC);
T_CP3(1:CP_Length_FFT_FBMC/2,end-CP_Length_FFT_FBMC/2+1:end) = eye(CP_Length_FFT_FBMC/2);
T_CP3(CP_Length_FFT_FBMC/2+1:end-CP_Length_FFT_FBMC/2,:)     = eye(K_FBMC-CP_Length_FFT_FBMC);
T_CP3(end-CP_Length_FFT_FBMC/2+1:end,1:CP_Length_FFT_FBMC/2) = eye(CP_Length_FFT_FBMC/2);
R_CP3                                                        = zeros(K_FBMC,K_FBMC-CP_Length_FFT_FBMC);
R_CP3(CP_Length_FFT_FBMC/2+1:end-CP_Length_FFT_FBMC/2,:)     = eye(K_FBMC-CP_Length_FFT_FBMC);

W3              = fft( eye(K_FBMC-CP_Length_FFT_FBMC) ) / sqrt( K_FBMC-CP_Length_FFT_FBMC );
a3              = abs(diag(W3'*R_CP3'*D_temp*T_CP3*W3));
a3              = a3+randn(size(a3))*10^-12; %  randn so that sorting is unique
a_Tilde3        = sort(a3,'descend');
alpha3          = a_Tilde3((K_FBMC-CP_Length_FFT_FBMC)/2);
Index_Tilde3    = (a3>=alpha3);
W_Tilde3_1      = W3(:,Index_Tilde3);
W_Tilde3_2      = W3(:,not(Index_Tilde3));

% W_Tilde3_1      = W_Tilde3_1*diag(sqrt(2./(a3(Index_Tilde3))));
% W_Tilde3_2      = W_Tilde3_2*diag(sqrt(2./(a3(not(Index_Tilde3)))));

C_FFT_FBMC_TX = kron(kron(T_CP3*W_Tilde3_1,sparse(eye(FBMC_TimeSpread.Nr.Subcarriers/2))),[1,0;0,0]) + kron(kron(T_CP3*W_Tilde3_2,sparse(eye(FBMC_TimeSpread.Nr.Subcarriers/2))),[0,0;0,1]);
C_FFT_FBMC_RX = kron(kron(R_CP3*W_Tilde3_1,sparse(eye(FBMC_TimeSpread.Nr.Subcarriers/2))),[1,0;0,0]) + kron(kron(R_CP3*W_Tilde3_2,sparse(eye(FBMC_TimeSpread.Nr.Subcarriers/2))),[0,0;0,1]);

G_temp = FBMC_TimeSpread.GetTXMatrix;
NormalizationFactor = sum(sum(abs(G_temp).^2))/sum(sum(abs(G_temp*C_FFT_FBMC_TX).^2));

C_FFT_FBMC_TX = C_FFT_FBMC_TX*sqrt(NormalizationFactor);
C_FFT_FBMC_RX = C_FFT_FBMC_RX/sqrt(NormalizationFactor);


%% Preallocate
Simulated_TransmitPowerOverTime_OFDM         = zeros(N,1);
Simulated_TransmitPowerOverTime_OFDMnoCP     = zeros(N,1);
Simulated_TransmitPowerOverTime_FBMC         = zeros(N,1);
Simulated_TransmitPowerOverTime_DFT_OFDM     = zeros(N,1);
Simulated_TransmitPowerOverTime_DFT_OFDMnoCP = zeros(N,1);
Simulated_TransmitPowerOverTime_FBMC_DFT     = zeros(N,1);
dt                                           = 1/SamplingRate; % time-spacing between samples

Simulated_PowerSpectralDensity_OFDM         = zeros(N,1);
Simulated_PowerSpectralDensity_OFDMnoCP     = zeros(N,1);
Simulated_PowerSpectralDensity_FBMC         = zeros(N,1);
Simulated_PowerSpectralDensity_DFT_OFDM     = zeros(N,1);
Simulated_PowerSpectralDensity_DFT_OFDMnoCP = zeros(N,1);
Simulated_PowerSpectralDensity_FBMC_DFT     = zeros(N,1);
df                                          = 1./(N*dt); % Frequency-spacing between samples


%% Start Simulation
tic;
for i_rep = 1:NrRepetitions   
    %% Generate Bit Stream
    BinaryDataStream            = randi( [0 1] , K_OFDMnoCP*NrSubcarriers*log2(QAM.ModulationOrder) ,1);
    BinaryDataStream_OFDM       = randi( [0 1] , K_OFDM*NrSubcarriers*log2(QAM.ModulationOrder) ,1);                            % Reduced bit rate due to the CP
    BinaryDataStream_FBMC_DFT   = randi( [0 1] , K_OFDMnoCP*(NrSubcarriers-CP_Length_FBMC_DFT)*log2(QAM.ModulationOrder) ,1);   % Maybe a reduced bit rate due to the CP which, however, is not necesarry most of the time => same bit rate

    BinaryDataStream_FBMC_DFT_LCP = randi( [0 1] , K_OFDMnoCP*(NrSubcarriers-CP_Length_FBMC_DFT_LCP)*log2(QAM.ModulationOrder) ,1);   % Maybe a reduced bit rate due to the CP which, however, is not necesarry most of the time => same bit rate
    BinaryDataStream_FFT_FBMC     = randi( [0 1] , K_OFDMnoCP*(NrSubcarriers-CP_Length_FFT_FBMC)*log2(QAM.ModulationOrder) ,1);   % Maybe a reduced bit rate due to the CP which, however, is not necesarry most of the time => same bit rate


    %% Map Bit Stream to Symbols
    x_OFDM      = reshape( QAM.Bit2Symbol(BinaryDataStream_OFDM) , NrSubcarriers, K_OFDM );
    x_OFDMnoCP  = reshape( QAM.Bit2Symbol(BinaryDataStream) , NrSubcarriers , K_OFDMnoCP );
    x_FBMC_DFT  = reshape( QAM.Bit2Symbol(BinaryDataStream_FBMC_DFT), (NrSubcarriers-CP_Length_FBMC_DFT)/2 , K_FBMC );
    x_FBMC      = reshape( PAM.Bit2Symbol(BinaryDataStream), NrSubcarriers , K_FBMC );

    x_FBMC_DFT_LCP  = reshape( QAM.Bit2Symbol(BinaryDataStream_FBMC_DFT_LCP), (NrSubcarriers - CP_Length_FBMC_DFT_LCP)/2 , K_FBMC );
    x_FBMC_FFT_FBMC = QAM.Bit2Symbol(BinaryDataStream_FFT_FBMC);


    %% Generate Transmit Signal in the Time Domain
    s_OFDM              = OFDM.Modulation( x_OFDM );
    s_OFDMnoCP          = OFDMnoCP.Modulation( x_OFDMnoCP );
    s_FBMC              = FBMC.Modulation( x_FBMC );

    s_DFT_OFDM          = OFDM.Modulation( DFTMatrix*x_OFDM );
    s_DFT_OFDMnoCP      = OFDMnoCP.Modulation( DFTMatrix*x_OFDMnoCP );
    s_FBMC_DFT          = FBMC_DFT.Modulation( C_DFTspread_TX*x_FBMC_DFT );
    s_FBMC_DFT_LCP      = FBMC_DFT.Modulation( C_DFTspread_TX2*x_FBMC_DFT_LCP ); % Large frequency CP

    s_FBMC_FFT_FBMC     = FBMC_TimeSpread.Modulation( reshape( C_FFT_FBMC_TX*x_FBMC_FFT_FBMC(:) , NrSubcarriers, K_FBMC )  ); 

    % Simple DFT spreading which keeps the OQAM structure
    x_FBMC_complex              = x_FBMC(:,1:2:end) + 1j*x_FBMC(:,2:2:end); % OQAM destaggering
    x_FBMC_complex_DFT          = DFTMatrix*x_FBMC_complex;              
    x_FBMC_simpleton            = nan(size(x_FBMC));
    x_FBMC_simpleton(:,1:2:end) = real(x_FBMC_complex_DFT);                 % OQAM staggering
    x_FBMC_simpleton(:,2:2:end) = imag(x_FBMC_complex_DFT);                 % OQAM staggering
    s_FBMC_SimpletonDFT         = FBMC.Modulation( x_FBMC_simpleton );

    s_FBMC_SimpletonDFT_PhaseCon = FBMC_PhaseCondition.Modulation( x_FBMC_simpleton );


    %% Calculate Transmit Power
    Simulated_TransmitPowerOverTime_OFDM         = abs(s_OFDM).^2         + Simulated_TransmitPowerOverTime_OFDM;
    Simulated_TransmitPowerOverTime_OFDMnoCP     = abs(s_OFDMnoCP).^2     + Simulated_TransmitPowerOverTime_OFDMnoCP;
    Simulated_TransmitPowerOverTime_FBMC         = abs(s_FBMC).^2         + Simulated_TransmitPowerOverTime_FBMC;
    Simulated_TransmitPowerOverTime_DFT_OFDM     = abs(s_DFT_OFDM).^2     + Simulated_TransmitPowerOverTime_DFT_OFDM;
    Simulated_TransmitPowerOverTime_DFT_OFDMnoCP = abs(s_DFT_OFDMnoCP).^2 + Simulated_TransmitPowerOverTime_DFT_OFDMnoCP;
    Simulated_TransmitPowerOverTime_FBMC_DFT     = abs(s_FBMC_DFT).^2     + Simulated_TransmitPowerOverTime_FBMC_DFT;

    
    %% Calculate Power Spectral Density
    Simulated_PowerSpectralDensity_OFDM         = abs(fft(s_OFDM)).^2         + Simulated_PowerSpectralDensity_OFDM;
    Simulated_PowerSpectralDensity_OFDMnoCP     = abs(fft(s_OFDMnoCP)).^2     + Simulated_PowerSpectralDensity_OFDMnoCP;
    Simulated_PowerSpectralDensity_FBMC         = abs(fft(s_FBMC)).^2         + Simulated_PowerSpectralDensity_FBMC;
    Simulated_PowerSpectralDensity_DFT_OFDM     = abs(fft(s_DFT_OFDM)).^2     + Simulated_PowerSpectralDensity_DFT_OFDM;
    Simulated_PowerSpectralDensity_DFT_OFDMnoCP = abs(fft(s_DFT_OFDMnoCP)).^2 + Simulated_PowerSpectralDensity_DFT_OFDMnoCP;
    Simulated_PowerSpectralDensity_FBMC_DFT     = abs(fft(s_FBMC_DFT)).^2     + Simulated_PowerSpectralDensity_FBMC_DFT;

    
    %% For the PAPR calculation we chop the signal in time-periods of length 1/F
    s_OFDM_reshaped                 = reshape( s_OFDM((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)) , OFDM.Implementation.FFTSize+OFDM.Implementation.CyclicPrefix , OFDM.Nr.MCSymbols );
    s_OFDM_reshaped                 = s_OFDM_reshaped(OFDM.Implementation.CyclicPrefix+1:end ,:); % Not necessary

    s_OFDMnoCP_reshaped             = reshape( s_OFDMnoCP((OFDMnoCP.Implementation.ZeroGuardSamples+1):(end-OFDMnoCP.Implementation.ZeroGuardSamples)) , OFDMnoCP.Implementation.FFTSize+OFDMnoCP.Implementation.CyclicPrefix , OFDMnoCP.Nr.MCSymbols );
    s_OFDMnoCP_reshaped             = s_OFDMnoCP_reshaped(OFDMnoCP.Implementation.CyclicPrefix+1:end,:); % Not necessary

    s_FBMC_reshaped                 = reshape( s_FBMC((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)) , OFDM.Implementation.FFTSize , FBMC.Nr.MCSymbols/2 );
    s_FBMC_reshaped(:,[1 end])      = [];

    % DFT
    s_DFT_OFDM_reshaped             = reshape( s_DFT_OFDM((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)) , OFDM.Implementation.FFTSize+OFDM.Implementation.CyclicPrefix , OFDM.Nr.MCSymbols );
    s_DFT_OFDM_reshaped             = s_DFT_OFDM_reshaped(OFDM.Implementation.CyclicPrefix+1:end,:); % Not necessary

    s_DFT_OFDMnoCP_reshaped         = reshape( s_DFT_OFDMnoCP((OFDMnoCP.Implementation.ZeroGuardSamples+1):(end-OFDMnoCP.Implementation.ZeroGuardSamples)) , OFDMnoCP.Implementation.FFTSize+OFDMnoCP.Implementation.CyclicPrefix , OFDMnoCP.Nr.MCSymbols);
    s_DFT_OFDMnoCP_reshaped         = s_DFT_OFDMnoCP_reshaped(OFDMnoCP.Implementation.CyclicPrefix+1:end,:); % Not necessary

    s_FBMC_DFT_reshaped             = reshape(s_FBMC_DFT((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)),OFDM.Implementation.FFTSize,FBMC.Nr.MCSymbols/2);
    s_FBMC_DFT_reshaped(:,[1 end])  = [];

    s_FBMC_SimpletonDFT_reshaped            = reshape(s_FBMC_SimpletonDFT((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)),OFDM.Implementation.FFTSize,FBMC.Nr.MCSymbols/2);
    s_FBMC_SimpletonDFT_reshaped(:,[1 end]) = [];

    s_FBMC_SimpletonDFT_PhaseCon_reshaped            = reshape(s_FBMC_SimpletonDFT_PhaseCon((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)),OFDM.Implementation.FFTSize,FBMC.Nr.MCSymbols/2);
    s_FBMC_SimpletonDFT_PhaseCon_reshaped(:,[1 end]) = [];

    s_FBMC_DFT_LCP_reshaped                 = reshape(s_FBMC_DFT_LCP((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)),OFDM.Implementation.FFTSize,FBMC.Nr.MCSymbols/2);
    s_FBMC_DFT_LCP_reshaped(:,[1 end])      = [];

    s_FBMC_FFT_FBMC_reshaped                 = reshape(s_FBMC_FFT_FBMC((OFDM.Implementation.ZeroGuardSamples+1):(end-OFDM.Implementation.ZeroGuardSamples)),OFDM.Implementation.FFTSize,FBMC.Nr.MCSymbols/2);
    s_FBMC_FFT_FBMC_reshaped(:,[1 end])      = [];

    
    %% Signal Power, Should be 1
    AveragePower_OFDM(i_rep)         = mean(mean( abs(s_OFDM_reshaped).^2 ));
    AveragePower_OFDMnoCP(i_rep)     = mean(mean( abs(s_OFDMnoCP_reshaped).^2 ));
    AveragePower_FBMC(i_rep)         = mean(mean( abs(s_FBMC_reshaped).^2 ));

    AveragePower_DFT_OFDM(i_rep)     = mean(mean( abs(s_DFT_OFDM_reshaped).^2 ));
    AveragePower_DFT_OFDMnoCP(i_rep) = mean(mean( abs(s_DFT_OFDMnoCP_reshaped).^2 ));
    AveragePower_FBMC_DFT(i_rep)     = mean(mean( abs(s_FBMC_DFT_reshaped).^2 ));

    AveragePower_FBMC_SimpletonDFT(i_rep)           = mean(mean( abs(s_FBMC_SimpletonDFT_reshaped).^2 ));
    AveragePower_FBMC_SimpletonDFT_PhaseCon(i_rep)  = mean(mean( abs(s_FBMC_SimpletonDFT_PhaseCon_reshaped).^2 ));
    AveragePower_FBMC_DFT_LCP(i_rep)                = mean(mean( abs(s_FBMC_DFT_LCP_reshaped).^2 ));

    AveragePower_FBMC_FFT_FBMC(i_rep) = mean(mean( abs(s_FBMC_FFT_FBMC_reshaped).^2 ));

    
    %% Peak power values for each time intervall 1/F (NFFT Samples)
    PeakPower_OFDM(i_rep,:)          = max( abs(s_OFDM_reshaped).^2 ,[],1);
    PeakPower_OFDMnoCP(i_rep,:)      = max( abs(s_OFDMnoCP_reshaped).^2 ,[],1);
    PeakPower_FBMC(i_rep,:)          = max( abs(s_FBMC_reshaped).^2 ,[],1);

    PeakPower_DFT_OFDM(i_rep,:)      = max( abs(s_DFT_OFDM_reshaped).^2 ,[],1);
    PeakPower_DFT_OFDMnoCP(i_rep,:)  = max( abs(s_DFT_OFDMnoCP_reshaped).^2 ,[],1);
    PeakPower_FBMC_DFT(i_rep,:)      = max( abs(s_FBMC_DFT_reshaped).^2 ,[],1);

    PeakPower_FBMC_SimpletonDFT(i_rep,:)          = max( abs(s_FBMC_SimpletonDFT_reshaped).^2 ,[],1);
    PeakPower_FBMC_SimpletonDFT_PhaseCon(i_rep,:) = max( abs(s_FBMC_SimpletonDFT_PhaseCon_reshaped).^2 ,[],1);
    PeakPower_FBMC_DFT_LCP(i_rep,:)               = max( abs(s_FBMC_DFT_LCP_reshaped).^2 ,[],1);

    PeakPower_FBMC_FFT_FBMC(i_rep,:)      = max( abs(s_FBMC_FFT_FBMC_reshaped).^2 ,[],1);

    TimePassed = toc;
    if mod(i_rep,10)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
    end
end


%% Plot the Peak-to-average power ratio
PAPR_OFDM         = 10*log10(PeakPower_OFDM(:));
PAPR_OFDMnoCP     = 10*log10(PeakPower_OFDMnoCP(:));
PAPR_FBMC         = 10*log10(PeakPower_FBMC(:));

PAPR_DFT_OFDM     = 10*log10(PeakPower_DFT_OFDM(:));
PAPR_DFT_OFDMnoCP = 10*log10(PeakPower_DFT_OFDMnoCP(:));
PAPR_FBMC_DFT     = 10*log10(PeakPower_FBMC_DFT(:));

PAPR_FBMC_SimpletonDFT = 10*log10(PeakPower_FBMC_SimpletonDFT(:));
PAPR_FBMC_SimpletonDFT_PhaseCon = 10*log10(PeakPower_FBMC_SimpletonDFT_PhaseCon(:));
PAPR_FBMC_DFT_LCP      = 10*log10(PeakPower_FBMC_DFT_LCP(:));
PAPR_FBMC_FFT_FBMC     = 10*log10(PeakPower_FBMC_FFT_FBMC(:));

[CCDF_OFDM , CCDF_OFDM_xAxis]                 = ecdf(PAPR_OFDM); CCDF_OFDM=1-CCDF_OFDM;
[CCDF_OFDMnoCP , CCDF_OFDMnoCP_xAxis]         = ecdf(PAPR_OFDMnoCP); CCDF_OFDMnoCP=1-CCDF_OFDMnoCP;
[CCDF_FBMC , CCDF_FBMC_xAxis]                 = ecdf(PAPR_FBMC); CCDF_FBMC=1-CCDF_FBMC;

[CCDF_DFT_OFDM , CCDF_DFT_OFDM_xAxis]         = ecdf(PAPR_DFT_OFDM); CCDF_DFT_OFDM=1-CCDF_DFT_OFDM;
[CCDF_DFT_OFDMnoCP , CCDF_DFT_OFDMnoCP_xAxis] = ecdf(PAPR_DFT_OFDMnoCP); CCDF_DFT_OFDMnoCP=1-CCDF_DFT_OFDMnoCP;
[CCDF_FBMC_DFT , CCDF_FBMC_DFT_xAxis]         = ecdf(PAPR_FBMC_DFT); CCDF_FBMC_DFT=1-CCDF_FBMC_DFT;

[CCDF_FBMC_SimpletonDFT , CCDF_FBMC_SimpletonDFT_xAxis] = ecdf(PAPR_FBMC_SimpletonDFT); CCDF_FBMC_SimpletonDFT=1-CCDF_FBMC_SimpletonDFT;
[CCDF_FBMC_SimpletonDFT_PhaseCon , CCDF_FBMC_SimpletonDFT_PhaseCon_xAxis] = ecdf(PAPR_FBMC_SimpletonDFT_PhaseCon); CCDF_FBMC_SimpletonDFT_PhaseCon=1-CCDF_FBMC_SimpletonDFT_PhaseCon;
[CCDF_FBMC_DFT_LCP , CCDF_FBMC_DFT_LCP_xAxis]         = ecdf(PAPR_FBMC_DFT_LCP); CCDF_FBMC_DFT_LCP=1-CCDF_FBMC_DFT_LCP;
[CCDF_FBMC_FFT_FBMC , CCDF_FBMC_FFT_FBMC_xAxis]         = ecdf(PAPR_FBMC_FFT_FBMC); CCDF_FBMC_FFT_FBMC=1-CCDF_FBMC_FFT_FBMC;

MaxMinAll_noDFT = max([min(CCDF_OFDM_xAxis),min(CCDF_OFDMnoCP_xAxis),min(CCDF_FBMC_xAxis),min(CCDF_FBMC_FFT_FBMC_xAxis)]);
MinMaxAll_noDFT = min([max(CCDF_OFDM_xAxis),max(CCDF_OFDMnoCP_xAxis),max(CCDF_FBMC_xAxis),max(CCDF_FBMC_FFT_FBMC_xAxis)]);

CCDF_OFDM_xAxis_LessPoints = linspace(MaxMinAll_noDFT,MinMaxAll_noDFT,100);
[x, index] = unique(CCDF_OFDM_xAxis);
CCDF_OFDM_LessPoints = 10.^interp1(x,log10(CCDF_OFDM(index)),CCDF_OFDM_xAxis_LessPoints);

CCDF_OFDMnoCP_xAxis_LessPoints = linspace(MaxMinAll_noDFT,MinMaxAll_noDFT,100);
[x, index] = unique(CCDF_OFDMnoCP_xAxis);
CCDF_OFDMnoCP_LessPoints = 10.^interp1(x,log10(CCDF_OFDMnoCP(index)),CCDF_OFDMnoCP_xAxis_LessPoints);

CCDF_FBMC_xAxis_LessPoints = linspace(MaxMinAll_noDFT,MinMaxAll_noDFT,100);
[x, index] = unique(CCDF_FBMC_xAxis);
CCDF_FBMC_LessPoints = 10.^interp1(x,log10(CCDF_FBMC(index)),CCDF_FBMC_xAxis_LessPoints);

CCDF_FBMC_FFT_FBMC_xAxis_LessPoints = linspace(MaxMinAll_noDFT,MinMaxAll_noDFT,100);
[x, index] = unique(CCDF_FBMC_FFT_FBMC_xAxis);
CCDF_FBMC_FFT_FBMC_LessPoints = 10.^interp1(x,log10(CCDF_FBMC_FFT_FBMC(index)),CCDF_FBMC_FFT_FBMC_xAxis_LessPoints);       


MaxMinAll_DFT = max([min(CCDF_DFT_OFDM_xAxis),min(CCDF_DFT_OFDMnoCP_xAxis),min(CCDF_FBMC_DFT_xAxis)]);
MinMaxAll_DFT = min([max(CCDF_DFT_OFDM_xAxis),max(CCDF_DFT_OFDMnoCP_xAxis),max(CCDF_FBMC_DFT_xAxis)]);

CCDF_DFT_OFDM_xAxis_LessPoints = linspace(MaxMinAll_DFT,MinMaxAll_DFT,100);
[x, index] = unique(CCDF_DFT_OFDM_xAxis);
CCDF_DFT_OFDM_LessPoints = 10.^interp1(x,log10(CCDF_DFT_OFDM(index)),CCDF_DFT_OFDM_xAxis_LessPoints);

CCDF_DFT_OFDMnoCP_xAxis_LessPoints = linspace(MaxMinAll_DFT,MinMaxAll_DFT,100);
[x, index] = unique(CCDF_DFT_OFDMnoCP_xAxis);
CCDF_DFT_OFDMnoCP_LessPoints = 10.^interp1(x,log10(CCDF_DFT_OFDMnoCP(index)),CCDF_DFT_OFDMnoCP_xAxis_LessPoints);

CCDF_FBMC_DFT_xAxis_LessPoints = linspace(MaxMinAll_DFT,MinMaxAll_DFT,100);
[x, index] = unique(CCDF_FBMC_DFT_xAxis);
CCDF_FBMC_DFT_LessPoints = 10.^interp1(x,log10(CCDF_FBMC_DFT(index)),CCDF_FBMC_DFT_xAxis_LessPoints);

MaxMinAll_SimpletonDFT = max([min(CCDF_FBMC_SimpletonDFT_xAxis)]);
MinMaxAll_SimpletonDFT = min([max(CCDF_FBMC_SimpletonDFT_xAxis)]);

CCDF_FBMC_SimpletonDFT_xAxis_LessPoints = linspace(MaxMinAll_SimpletonDFT,MinMaxAll_SimpletonDFT,100);
[x, index] = unique(CCDF_FBMC_SimpletonDFT_xAxis);
CCDF_FBMC_SimpletonDFT_LessPoints = 10.^interp1(x,log10(CCDF_FBMC_SimpletonDFT(index)),CCDF_FBMC_SimpletonDFT_xAxis_LessPoints);


MaxMinAll_SimpletonDFT_PhaseCon = max([min(CCDF_FBMC_SimpletonDFT_PhaseCon_xAxis)]);
MinMaxAll_SimpletonDFT_PhaseCon = min([max(CCDF_FBMC_SimpletonDFT_PhaseCon_xAxis)]);

CCDF_FBMC_SimpletonDFT_xAxis_LessPoints_PhaseCon = linspace(MaxMinAll_SimpletonDFT_PhaseCon,MinMaxAll_SimpletonDFT_PhaseCon,100);
[x, index] = unique(CCDF_FBMC_SimpletonDFT_PhaseCon_xAxis);
CCDF_FBMC_SimpletonDFT_LessPoints_PhaseCon = 10.^interp1(x,log10(CCDF_FBMC_SimpletonDFT_PhaseCon(index)),CCDF_FBMC_SimpletonDFT_xAxis_LessPoints_PhaseCon);    

MaxMinAll_DFT_LCP = max(min(CCDF_FBMC_DFT_LCP_xAxis));
MinMaxAll_DFT_LCP = min(max(CCDF_FBMC_DFT_LCP_xAxis));

CCDF_FBMC_DFT_LCP_xAxis_LessPoints = linspace(MaxMinAll_DFT_LCP,MinMaxAll_DFT_LCP,100);
[x, index] = unique(CCDF_FBMC_DFT_LCP_xAxis);
CCDF_FBMC_DFT_LCP_LessPoints = 10.^interp1(x,log10(CCDF_FBMC_DFT_LCP(index)),CCDF_FBMC_DFT_LCP_xAxis_LessPoints);   


figure(9);
markersize = 3;
s1 = semilogy(CCDF_OFDM_xAxis_LessPoints,CCDF_OFDM_LessPoints,'Color',[127 127 0]/255); hold on;
s2 = semilogy(CCDF_FBMC_xAxis_LessPoints,CCDF_FBMC_LessPoints,'magenta'); hold on;
semilogy(CCDF_OFDM_xAxis_LessPoints(30:10:100),CCDF_OFDM_LessPoints(30:10:100),'s','Markersize',markersize,'Color',[127 127 0]/255); hold on;    
semilogy(CCDF_FBMC_xAxis_LessPoints(35:10:100),CCDF_FBMC_LessPoints(35:10:100),'magenta  o','Markersize',markersize); hold on;

s3 = semilogy(CCDF_FBMC_SimpletonDFT_xAxis_LessPoints,CCDF_FBMC_SimpletonDFT_LessPoints,'Color',[127 127 127]/255); hold on;
semilogy(CCDF_FBMC_SimpletonDFT_xAxis_LessPoints(35:10:100),CCDF_FBMC_SimpletonDFT_LessPoints(35:10:100),'o','Markersize',markersize,'Color',[127 127 127]/255); hold on;

s4 = semilogy(CCDF_FBMC_SimpletonDFT_xAxis_LessPoints_PhaseCon,CCDF_FBMC_SimpletonDFT_LessPoints_PhaseCon,'Color',[80 80 80]/255); hold on;
semilogy(CCDF_FBMC_SimpletonDFT_xAxis_LessPoints_PhaseCon(35:10:100),CCDF_FBMC_SimpletonDFT_LessPoints_PhaseCon(35:10:100),'o','Markersize',markersize,'Color',[80 80 80]/255); hold on;


s5 = semilogy(CCDF_DFT_OFDM_xAxis_LessPoints,CCDF_DFT_OFDM_LessPoints,'black'); hold on;
s6 = semilogy(CCDF_FBMC_DFT_xAxis_LessPoints,CCDF_FBMC_DFT_LessPoints,'blue'); hold on;
semilogy(CCDF_DFT_OFDM_xAxis_LessPoints(30:10:100),CCDF_DFT_OFDM_LessPoints(30:10:100),'black s','Markersize',markersize); hold on;
semilogy(CCDF_FBMC_DFT_xAxis_LessPoints(35:10:100),CCDF_FBMC_DFT_LessPoints(35:10:100),'blue o','Markersize',markersize); hold on;

s7 = semilogy(CCDF_FBMC_DFT_LCP_xAxis_LessPoints,CCDF_FBMC_DFT_LCP_LessPoints,'Color',[0 0 127/255]); hold on;
semilogy(CCDF_FBMC_DFT_LCP_xAxis_LessPoints(35:10:100),CCDF_FBMC_DFT_LCP_LessPoints(35:10:100),'o','Markersize',markersize,'Color',[0 0 127/255]); hold on;

s8 = semilogy(CCDF_FBMC_FFT_FBMC_xAxis_LessPoints,CCDF_FBMC_FFT_FBMC_LessPoints,'Color',[0 0 180/255]); hold on;
semilogy(CCDF_FBMC_FFT_FBMC_xAxis_LessPoints(33:10:100),CCDF_FBMC_FFT_FBMC_LessPoints(33:10:100),'o','Markersize',markersize,'Color',[0 0 180/255]); hold on;

ylim([10^-3 1]);
xlim([4 11]);

ylabel('CCDF');
xlabel('Peak-to-Average Power Ratio [dB]');
legend([s1 s2 s3 s4 s5 s6 s7 s8],{'OFDM','FBMC','Simple DFT-s FBMC','Simple DFT-s FBMC (phase)','SC-FDMA','Pruned DFT-s FBMC','Pruned DFT-s FBMC (L_CP = 26)','FFT-FBMC'});


%% Plot the Simulated Transmit Power
Simulated_TransmitPowerOverTime_OFDM         = Simulated_TransmitPowerOverTime_OFDM/i_rep;
Simulated_TransmitPowerOverTime_OFDMnoCP     = Simulated_TransmitPowerOverTime_OFDMnoCP/i_rep;
Simulated_TransmitPowerOverTime_FBMC         = Simulated_TransmitPowerOverTime_FBMC/i_rep;
Simulated_TransmitPowerOverTime_DFT_OFDM     = Simulated_TransmitPowerOverTime_DFT_OFDM/i_rep;
Simulated_TransmitPowerOverTime_DFT_OFDMnoCP = Simulated_TransmitPowerOverTime_DFT_OFDMnoCP/i_rep;
Simulated_TransmitPowerOverTime_FBMC_DFT     = Simulated_TransmitPowerOverTime_FBMC_DFT/i_rep;

figure(100);
t = (0:N-1)*dt-OFDM.PHY.ZeroGuardTimeLength;
plot(t/1e-3,Simulated_TransmitPowerOverTime_OFDM,': black'); hold on;
plot(t/1e-3,Simulated_TransmitPowerOverTime_FBMC,': blue'); hold on;
plot(t/1e-3,Simulated_TransmitPowerOverTime_DFT_OFDM,'black'); hold on;
plot(t/1e-3,Simulated_TransmitPowerOverTime_FBMC_DFT,'blue'); hold on;
legend({'CP-OFDM','FBMC','DFT-OFDM','DFT-FBMC'})
ylabel('Simulated Transmit Power');
xlabel('Time [ms]');


%% Plot the Simulated Power Spectral Density
figure(101);
f = (0:N-1)*SamplingRate/N;
IndexInBand = f<(FBMC.Nr.Subcarriers-1)*FBMC.PHY.SubcarrierSpacing; % In band minus 1
NormalizeFactorPSD = mean([Simulated_PowerSpectralDensity_OFDM(IndexInBand);Simulated_PowerSpectralDensity_OFDMnoCP(IndexInBand);Simulated_PowerSpectralDensity_FBMC(IndexInBand);Simulated_PowerSpectralDensity_DFT_OFDM(IndexInBand);Simulated_PowerSpectralDensity_DFT_OFDMnoCP(IndexInBand);Simulated_PowerSpectralDensity_FBMC_DFT(IndexInBand)]);

f = f+0.5*FBMC.PHY.SubcarrierSpacing; % -0.5 because the first subcarrier starts at -0.5*FBMC.PHY.SubcarrierSpacing => Shift by 0.5 so that the bandwidth starts at exactly 0
plot(f/1e6,10*log10(Simulated_PowerSpectralDensity_OFDM/NormalizeFactorPSD),': black'); hold on;
plot(f/1e6,10*log10(Simulated_PowerSpectralDensity_FBMC/NormalizeFactorPSD),': blue'); hold on;
plot(f/1e6,10*log10(Simulated_PowerSpectralDensity_DFT_OFDM/NormalizeFactorPSD),'black'); hold on;
plot(f/1e6,10*log10(Simulated_PowerSpectralDensity_FBMC_DFT/NormalizeFactorPSD),'blue'); hold on;
legend({'CP-OFDM','FBMC','DFT-OFDM','DFT-FBMC'})

UpperTransmitFrequency  = (FBMC.Nr.Subcarriers)*FBMC.PHY.SubcarrierSpacing;
LowerPlotFrequency      = (FBMC.Nr.Subcarriers-5)*FBMC.PHY.SubcarrierSpacing;
UpperPlotFrequency      = (FBMC.Nr.Subcarriers+20)*FBMC.PHY.SubcarrierSpacing;

plot([UpperTransmitFrequency/1e6,UpperTransmitFrequency/1e6],[-100,20],'-','Color',[1 1 1]*0.5);
xlim([LowerPlotFrequency/1e6 UpperPlotFrequency/1e6]);
ylim([-60 2]);

xlabel('Frequency [MHz]');
ylabel('Simulated Power Spectral Density');


%% Show the Bitrates
BitRate_OFDM     = length(BinaryDataStream_OFDM)/(OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols);
BitRate_OFDMnoCP = length(BinaryDataStream)/(OFDMnoCP.PHY.TimeSpacing*OFDMnoCP.Nr.MCSymbols);
BitRate_FBMC     = length(BinaryDataStream)/(FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols);
BitRate_FBMC_DFT = length(BinaryDataStream_FBMC_DFT)/(FBMC_DFT.PHY.TimeSpacing*FBMC_DFT.Nr.MCSymbols);
fprintf('=================================================================\n');
fprintf('                  |  CP-OFDM | OFDM(noCP) |    FBMC  | pDFTsFBMC |\n');
fprintf('Bit Rate [Bits/s] |%10.0f|%12.0f|%10.0f|%11.0f|\n',BitRate_OFDM,BitRate_OFDMnoCP,BitRate_FBMC,BitRate_FBMC_DFT);
fprintf('=================================================================\n');

