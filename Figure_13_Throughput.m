% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================   
% This script simulates the throughput of pruned DFT spread FBMC, SC-FDMA,
% OFDM and FBMC. Note that the achievable rate can also be calculated, but
% the parameter "CalculateTheory" must be set to true. 

% Allows to reproduce Figure 13 of "Pruned DFT Spread FBMC: Low PAPR, Low 
% Latency, High Spectral Efficiency", R. Nissel and M. Rupp, IEEE
% Transactions on Communications
% #### We recommend using parfor, that is, commenting out line 264-266 ####

clear; close all;
addpath('./Theory');

%% Parameters
% Simulation
M_SNR_dB            = [-10:4:30];                   % Signal-to-Noise Ratio in dB
NrRepetitions       = 4;                            % Number of Monte Carlo repetitions
CalculateTheory     = false;                        % If set to true, calculate the achievable rate, an upper bound of the throughput. To keep the simulation time short we set it to false. 

% FBMC and OFDM Parameters
NrSubcarriers       = 256;                          % Number of subcarriers
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
% M_SNR_dB        = [-10:1.25:30];
% SamplingRate    = 15e3*14*12*8; 
% NrRepetitions   = 1000;
% CalculateTheory = true;
% #########################################################################


%% Adaptive Modulation and Coding (CQI Table)
% The first column represents the modulation order: 4, 16, 64, 256, 1024...
% The second column represents the code rate (must be between zero and one)
% Currently, values are chosen according to the (old) LTE standard:
M_CQI = [4  ,  78/1024;...
         4  , 120/1024;...
         4  , 193/1024;...
         4  , 308/1024;...
         4  , 449/1024;...
         4  , 602/1024;...
         16 , 378/1024;...
         16 , 490/1024;...
         16 , 616/1024;...
         64 , 466/1024;...
         64 , 567/1024;...
         64 , 666/1024;...
         64 , 772/1024;...
         64 , 873/1024;...
         64 , 948/1024]; % page 48 of http://www.etsi.org/deliver/etsi_ts/136200_136299/136213/08.08.00_60/ts_136213v080800p.pdf 

if not(strcmp(mexext,'mexw64'))  
    % We use a win64 mexfile for code rates smaller than 1/3 => only works 
    % in 64-bit Windows
    IndexCodeRateSmallerOneThird =  find(M_CQI(:,2)<1/3);
    if  numel(IndexCodeRateSmallerOneThird)>0
        M_CQI(IndexCodeRateSmallerOneThird,:) = [];
        warning('A code rate smaller than 1/3 is only supported in Windows 64-bit => CQI values which contain a code rate smaller than 1/3 are discarded!');
    end    
end
    
%% FBMC Objects
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
% We use the short notation "FBMC_DFT" to indicate pruned DFT spread FBMC.
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


%% OFDM Objects
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


%% Pre-initialize CQI: Turbo Coder and QAM
for i_cqi = 1:size(M_CQI,1)
    QAMModulationOrder  = M_CQI(i_cqi,1);
    PAMModulationOrder  = sqrt(QAMModulationOrder);
    CodeRate            = M_CQI(i_cqi,2);
    
    QAM{i_cqi}          = Modulation.SignalConstellation(QAMModulationOrder,'QAM');
    PAM{i_cqi}          = Modulation.SignalConstellation(PAMModulationOrder,'PAM');

    NrTransmittedBits_OFDM     = NrSubcarriers*K_OFDM*log2(QAMModulationOrder);
    NrTransmittedBits_OFDMnoCP = NrSubcarriers*K_OFDMnoCP*log2(QAMModulationOrder);
    NrTransmittedBits_FBMC     = NrSubcarriers*K_FBMC*log2(PAMModulationOrder);
    NrTransmittedBits_FBMC_DFT = (NrSubcarriers-CP_Length_FBMC_DFT)/2*K_FBMC*log2(QAMModulationOrder);
  
    TurboCoding_OFDM{i_cqi}     = Coding.TurboCoding( NrTransmittedBits_OFDM     , round(CodeRate*NrTransmittedBits_OFDM));
    TurboCoding_OFDMnoCP{i_cqi} = Coding.TurboCoding( NrTransmittedBits_OFDMnoCP , round(CodeRate*NrTransmittedBits_OFDMnoCP));
    TurboCoding_FBMC{i_cqi}     = Coding.TurboCoding( NrTransmittedBits_FBMC     , round(CodeRate*NrTransmittedBits_FBMC));
    TurboCoding_FBMC_DFT{i_cqi} = Coding.TurboCoding( NrTransmittedBits_FBMC_DFT , round(CodeRate*NrTransmittedBits_FBMC_DFT));
end


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

% Diagonal elements of the FBMC transmission matrix after DFT spreading and despreading
a = abs(diag(W'*R_CP'*D_temp*T_CP*W));
a = a+randn(size(a))*10^-12; %  randn so that sorting is unique

% Sort a
a_Tilde = sort(a,'descend');

% Get index which represents the largest values of a
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
SINR_FBMC_dB     = nan(NrSubcarriers,K_FBMC,length(M_SNR_dB),NrRepetitions);
SINR_FBMC_DFT_dB = nan((NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC,length(M_SNR_dB),NrRepetitions);

% Preallocate simulation results (needed for parfor)
M_Througput_OFDM            = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_OFDMnoCP        = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_FBMC            = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );

M_Througput_DFT_OFDM        = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_DFT_OFDMnoCP    = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );
M_Througput_FBMC_DFT        = nan( length(M_SNR_dB) , NrRepetitions , size(M_CQI,1) );

disp('The simulation may take a while ... ');
%% Start Simulation (Calculation)
NrWorkers = 1;                                       % Conventional FOR loop                   
for i_Rep = 1:NrRepetitions                          % Conventional FOR loop  
% cluster     = parcluster('local');                 % PARFOR
% NrWorkers   = cluster.NumWorkers;                  % PARFOR 
% parfor i_Rep = 1:NrRepetitions                     % PARFOR 
    tic;
    % Channel 
    ChannelModel.NewRealization;
    H = ChannelModel.GetConvolutionMatrix{1};
    noise_unitPower = sqrt(1/2)*( randn(N,1) + 1j * randn(N,1) );
      
    % Calculate One Tap Channels
    % Note that diag(G'*H*G)==sum((G'*H).*G.',2)
    h_OFDM      = full( reshape( sum((GRX_OFDM'*H).*GTX_OFDM.',2), NrSubcarriers, [] ) ); 
    h_OFDMnoCP  = full( reshape( sum((GRX_OFDMnoCP'*H).*GTX_OFDMnoCP.',2), NrSubcarriers, [] ) ); 
    h_FBMC      = full( reshape( sum((G_FBMC'*H).*G_FBMC.',2), NrSubcarriers, [] ) ); 
    h_FBMC_DFT  = full( reshape( sum((G_FBMC_DFT'*H).*G_FBMC_DFT.',2), NrSubcarriers, [] ) ); 
 
    % Preallocate simulation result for one realization (needed for parfor)
    M_Througput_OFDM_OneRealization       = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_OFDMnoCP_OneRealization   = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_FBMC_OneRealization       = nan( length(M_SNR_dB) , size(M_CQI,1) );
    
    M_Througput_DFT_OFDM_OneRealization     = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_DFT_OFDMnoCP_OneRealization = nan( length(M_SNR_dB) , size(M_CQI,1) );
    M_Througput_FBMC_DFT_OneRealization     = nan( length(M_SNR_dB) , size(M_CQI,1) );    
   
    % Simulate over different modulation orders and code rates
    for i_cqi = 1:size(M_CQI,1)
        % Generate Data Bit Stream
        BinaryDataStream_OFDM     = randi( [0 1] , TurboCoding_OFDM{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_OFDMnoCP = randi( [0 1] , TurboCoding_OFDMnoCP{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_FBMC     = randi( [0 1] , TurboCoding_FBMC{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_FBMC_DFT = randi( [0 1] , TurboCoding_FBMC_DFT{i_cqi}.NrDataBits , 1 );
        
        % Update Interleaving of the Turbo Coder
        TurboCoding_OFDM{i_cqi}.UpdateInterleaving;
        TurboCoding_OFDMnoCP{i_cqi}.UpdateInterleaving;        
        TurboCoding_FBMC{i_cqi}.UpdateInterleaving;        
        TurboCoding_FBMC_DFT{i_cqi}.UpdateInterleaving;        
       
        % Turbo Coding of the Data Bits
        CodedBits_OFDM      = TurboCoding_OFDM{i_cqi}.TurboEncoder( BinaryDataStream_OFDM );
        CodedBits_OFDMnoCP  = TurboCoding_OFDMnoCP{i_cqi}.TurboEncoder( BinaryDataStream_OFDMnoCP );
        CodedBits_FBMC      = TurboCoding_FBMC{i_cqi}.TurboEncoder( BinaryDataStream_FBMC );
        CodedBits_FBMC_DFT  = TurboCoding_FBMC_DFT{i_cqi}.TurboEncoder( BinaryDataStream_FBMC_DFT );
               
        % Bit Interleaving
        BitInterleaving_OFDM     = randperm( TurboCoding_OFDM{i_cqi}.NrCodedBits );
        BitInterleaving_OFDMnoCP = randperm( TurboCoding_OFDMnoCP{i_cqi}.NrCodedBits );
        BitInterleaving_FBMC     = randperm( TurboCoding_FBMC{i_cqi}.NrCodedBits );
        BitInterleaving_FBMC_DFT = randperm( TurboCoding_FBMC_DFT{i_cqi}.NrCodedBits );
        
        CodedBits_OFDM     = CodedBits_OFDM(     BitInterleaving_OFDM );
        CodedBits_OFDMnoCP = CodedBits_OFDMnoCP( BitInterleaving_OFDMnoCP );
        CodedBits_FBMC     = CodedBits_FBMC(     BitInterleaving_FBMC );
        CodedBits_FBMC_DFT = CodedBits_FBMC_DFT( BitInterleaving_FBMC_DFT );
             
        % Map Bit Stream to Symbols
        x_OFDM     = reshape(QAM{i_cqi}.Bit2Symbol(CodedBits_OFDM)    , NrSubcarriers,K_OFDM);
        x_OFDMnoCP = reshape(QAM{i_cqi}.Bit2Symbol(CodedBits_OFDMnoCP), NrSubcarriers,K_OFDMnoCP);
        x_FBMC     = reshape(PAM{i_cqi}.Bit2Symbol(CodedBits_FBMC)    , NrSubcarriers,K_FBMC)/sqrt(2);  % 1/sqrt(2) => same TX power for OFDM and FBMC
       
        x_FBMC_DFT = reshape(QAM{i_cqi}.Bit2Symbol(CodedBits_FBMC_DFT),(NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC);

        % Generate Transmit Signal in the Time Domain
        s_OFDM         = OFDM.Modulation(x_OFDM)*NormalizationOFDM;
        s_OFDMnoCP     = OFDMnoCP.Modulation(x_OFDMnoCP)*NormalizationOFDMnoCP;
        s_FBMC         = FBMC.Modulation(x_FBMC)*NormalizationFBMC;

        s_DFT_OFDM     = OFDM.Modulation(DFTMatrix*x_OFDM)*NormalizationOFDM;
        s_DFT_OFDMnoCP = OFDMnoCP.Modulation(DFTMatrix*x_OFDMnoCP)*NormalizationOFDMnoCP;
        s_FBMC_DFT     = FBMC_DFT.Modulation(Cf_DFTspread_TX*x_FBMC_DFT)*NormalizationFBMC_DFT;

        % Channel    
        r_OFDM_noNoise         = H * s_OFDM;
        r_OFDMnoCP_noNoise     = H * s_OFDMnoCP;
        r_FBMC_noNoise         = H * s_FBMC;

        r_DFT_OFDM_noNoise     = H * s_DFT_OFDM;
        r_DFT_OFDMnoCP_noNoise = H * s_DFT_OFDMnoCP;
        r_FBMC_DFT_noNoise     = H * s_FBMC_DFT;

        % Simulate over different SNR values
        for i_SNR = 1:length(M_SNR_dB)
            SNR_dB  = M_SNR_dB(i_SNR); 
            Pn      = 10^(-SNR_dB/10); 

            % Add Noise
            noise = sqrt(Pn) * noise_unitPower;

            r_OFDM         = r_OFDM_noNoise     + noise;
            r_OFDMnoCP     = r_OFDMnoCP_noNoise + noise;
            r_FBMC         = r_FBMC_noNoise     + noise;

            r_DFT_OFDM     = r_DFT_OFDM_noNoise     + noise;
            r_DFT_OFDMnoCP = r_DFT_OFDMnoCP_noNoise + noise;
            r_FBMC_DFT     = r_FBMC_DFT_noNoise     + noise;

            % Received Symbols (Demodulation)
            y_OFDM         = OFDM.Demodulation(r_OFDM)/NormalizationOFDM;
            y_OFDMnoCP     = OFDMnoCP.Demodulation(r_OFDMnoCP)/NormalizationOFDMnoCP;
            y_FBMC         = FBMC.Demodulation(r_FBMC)/NormalizationFBMC;

            y_DFT_OFDM     = OFDM.Demodulation(r_DFT_OFDM)/NormalizationOFDM;
            y_DFT_OFDMnoCP = OFDMnoCP.Demodulation(r_DFT_OFDMnoCP)/NormalizationOFDMnoCP;
            y_FBMC_DFT     = FBMC_DFT.Demodulation(r_FBMC_DFT)/NormalizationFBMC_DFT;

            % ZF Equalizer for OFDM and FBMC
            x_est_OFDM      = y_OFDM./h_OFDM;
            x_est_OFDMnoCP  = y_OFDMnoCP./h_OFDMnoCP;
            x_est_FBMC      = real(y_FBMC./h_FBMC)*sqrt(2);

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

            % Calculate LLR Values
            LLR_OFDM     = QAM{i_cqi}.LLR_AWGN( x_est_OFDM(:)     , Pn .* 1./abs(h_OFDM(:)).^2);
            LLR_OFDMnoCP = QAM{i_cqi}.LLR_AWGN( x_est_OFDMnoCP(:) , Pn .* 1./abs(h_OFDMnoCP(:)).^2);
            LLR_FBMC     = PAM{i_cqi}.LLR_AWGN( x_est_FBMC(:)     , 2*Pn .* 1./abs(h_FBMC(:)).^2);
        
            AWGNequivalentNoise_DFT_OFDM     = repmat(1./mean(1./(1+Pn./abs(h_OFDM).^2),1)-1 ,NrSubcarriers,1);
            AWGNequivalentNoise_DFT_OFDMnoCP = repmat(1./mean(1./(1+Pn./abs(h_OFDMnoCP).^2),1)-1,NrSubcarriers,1);
            AWGNequivalentNoise_FBMC_DFT     = repmat(1./mean(1./(1+Pn./abs(h_FBMC_DFT(CP_Length_FBMC_DFT/2+1:end-CP_Length_FBMC_DFT/2,:)).^2),1)-1, (NrSubcarriers-CP_Length_FBMC_DFT)/2 , 1 );
                     
            LLR_DFT_OFDM     = QAM{i_cqi}.LLR_AWGN( x_est_DFT_OFDM(:)     , AWGNequivalentNoise_DFT_OFDM(:));
            LLR_DFT_OFDMnoCP = QAM{i_cqi}.LLR_AWGN( x_est_DFT_OFDMnoCP(:) , AWGNequivalentNoise_DFT_OFDMnoCP(:));
            LLR_FBMC_DFT     = QAM{i_cqi}.LLR_AWGN( x_est_FBMC_DFT(:)     , AWGNequivalentNoise_FBMC_DFT(:));
                       
            % Bitdeinterleaving
            LLR_OFDM(BitInterleaving_OFDM)         = LLR_OFDM;
            LLR_OFDMnoCP(BitInterleaving_OFDMnoCP) = LLR_OFDMnoCP;
            LLR_FBMC(BitInterleaving_FBMC)         = LLR_FBMC;
            
            LLR_DFT_OFDM(BitInterleaving_OFDM)         = LLR_DFT_OFDM;
            LLR_DFT_OFDMnoCP(BitInterleaving_OFDMnoCP) = LLR_DFT_OFDMnoCP;
            LLR_FBMC_DFT(BitInterleaving_FBMC_DFT)     = LLR_FBMC_DFT;           

            % Decode Bits
            DecodedBits_OFDM     = TurboCoding_OFDM{i_cqi}.TurboDecoder( LLR_OFDM );
            DecodedBits_OFDMnoCP = TurboCoding_OFDMnoCP{i_cqi}.TurboDecoder( LLR_OFDMnoCP );            
            DecodedBits_FBMC     = TurboCoding_FBMC{i_cqi}.TurboDecoder( LLR_FBMC );
            
            DecodedBits_DFT_OFDM     = TurboCoding_OFDM{i_cqi}.TurboDecoder( LLR_DFT_OFDM );
            DecodedBits_DFT_OFDMnoCP = TurboCoding_OFDMnoCP{i_cqi}.TurboDecoder( LLR_DFT_OFDMnoCP );            
            DecodedBits_FBMC_DFT     = TurboCoding_FBMC_DFT{i_cqi}.TurboDecoder( LLR_FBMC_DFT );            
              
            % Simulated throughput after decoding (all bits must be correctly detected. If one bit is wrong, the throughput is zero)
            M_Througput_OFDM_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_OFDM  == BinaryDataStream_OFDM ) * length(BinaryDataStream_OFDM)/(OFDM.PHY.TimeSpacing*(OFDM.Nr.MCSymbols));
            M_Througput_OFDMnoCP_OneRealization(i_SNR,i_cqi) = all( DecodedBits_OFDMnoCP  == BinaryDataStream_OFDMnoCP ) * length(BinaryDataStream_OFDMnoCP)/(OFDMnoCP.PHY.TimeSpacing*(OFDMnoCP.Nr.MCSymbols));
            M_Througput_FBMC_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_FBMC  == BinaryDataStream_FBMC ) * length(BinaryDataStream_FBMC)/(FBMC.PHY.TimeSpacing*(FBMC.Nr.MCSymbols));
 
            M_Througput_DFT_OFDM_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_DFT_OFDM  == BinaryDataStream_OFDM ) * length(BinaryDataStream_OFDM)/(OFDM.PHY.TimeSpacing*(OFDM.Nr.MCSymbols));
            M_Througput_DFT_OFDMnoCP_OneRealization(i_SNR,i_cqi) = all( DecodedBits_DFT_OFDMnoCP  == BinaryDataStream_OFDMnoCP ) * length(BinaryDataStream_OFDMnoCP)/(OFDMnoCP.PHY.TimeSpacing*(OFDMnoCP.Nr.MCSymbols));
            M_Througput_FBMC_DFT_OneRealization(i_SNR,i_cqi)     = all( DecodedBits_FBMC_DFT  == BinaryDataStream_FBMC_DFT ) * length(BinaryDataStream_FBMC_DFT)/(FBMC_DFT.PHY.TimeSpacing*(FBMC_DFT.Nr.MCSymbols));

        end
    end

    M_Througput_OFDM(:,i_Rep,:)      = M_Througput_OFDM_OneRealization;
    M_Througput_OFDMnoCP(:,i_Rep,:)  = M_Througput_OFDMnoCP_OneRealization;
    M_Througput_FBMC(:,i_Rep,:)      = M_Througput_FBMC_OneRealization;
    
    M_Througput_DFT_OFDM(:,i_Rep,:)      = M_Througput_DFT_OFDM_OneRealization;
    M_Througput_DFT_OFDMnoCP(:,i_Rep,:)  = M_Througput_DFT_OFDMnoCP_OneRealization;
    M_Througput_FBMC_DFT(:,i_Rep,:)      = M_Througput_FBMC_DFT_OneRealization;   
    
    
    if CalculateTheory
        %% Calculate Theoretical SINR
        % Precalculate Stuff to Improve Computation time
        Precalc_FBMC     = G_FBMC' * H * G_FBMC;      
        Precalc_FBMC_DFT = G_FBMC_DFT' * H * G_FBMC_DFT * C_DFTspread_TX;  

        SINR_FBMC_dB_OneRealization = nan(NrSubcarriers, K_FBMC, length(M_SNR_dB) );    
        SINR_FBMC_DFT_dB_OneRealization = nan((NrSubcarriers-CP_Length_FBMC_DFT)/2, K_FBMC, length(M_SNR_dB) );

        for i_SNR = 1:length(M_SNR_dB)
            SNR_dB  = M_SNR_dB(i_SNR); 
            Pn      = 10^(-SNR_dB/10); 

            % Equalizer
            E_FBMC = sparse(diag(1./h_FBMC(:)));

            Scaling_FBMC_DFT = repmat(1./(mean(1 ./ ( 1 + Pn./abs( h_FBMC(CP_Length_FBMC_DFT/2+1:end-CP_Length_FBMC_DFT/2,:) ).^2 ),1)), NrSubcarriers,1);
            e_FBMC_DFT       = Scaling_FBMC_DFT .* conj(h_FBMC) ./ ( abs(h_FBMC).^2 + Pn );
            E_FBMC_DFT       = sparse( diag( e_FBMC_DFT(:) ) );

            % Transmission Matrix
            D_FMBC = real(E_FBMC*Precalc_FBMC);
            Gamma_FBMC = 1./abs(h_FBMC(:)).^2;
            SINR_FBMC_dB_OneRealization(:,:,i_SNR) = 10*log10(1./(reshape( sum(abs(D_FMBC-eye(size(D_FMBC))).^2,2) + Gamma_FBMC * Pn, [NrSubcarriers, K_FBMC])));

            D_FBMC_DFT     = C_DFTspread_RX' * E_FBMC_DFT * Precalc_FBMC_DFT;
            Gamma_FBMC_DFT = C_DFTspread_RX' * E_FBMC_DFT * G_FBMC_DFT'; 
            SINR_FBMC_DFT_dB_OneRealization(:,:,i_SNR) = 10*log10(1./(reshape( sum(abs(D_FBMC_DFT-eye(size(D_FBMC_DFT))).^2,2) + sum(abs(Gamma_FBMC_DFT).^2,2)*Pn , [(NrSubcarriers-CP_Length_FBMC_DFT)/2, K_FBMC])));
        end

        SINR_FBMC_dB(:,:,:,i_Rep)     = SINR_FBMC_dB_OneRealization;
        SINR_FBMC_DFT_dB(:,:,:,i_Rep) = SINR_FBMC_DFT_dB_OneRealization;

    end
    
    TimePassed = toc;
    disp(['Realization ' int2str(i_Rep) ' of ' int2str(NrRepetitions) ' needed ' int2str(TimePassed) 's. Total simulation time:' int2str(TimePassed*NrRepetitions/NrWorkers/60) 'minutes']);

end


if CalculateTheory
    %% Calculate Achievable Rate
    Load_BICM_QAM = load('Theory\BICM_Capacity_AWGN_4_16_64_QAM');
    Load_BICM_PAM = load('Theory\BICM_Capacity_AWGN_2_4_8_PAM');
    
    BICM_QAM = Load_BICM_QAM.C_max;
    BICM_PAM = Load_BICM_PAM.C_max;
    
    SNR_Ref_QAM = Load_BICM_QAM.SNR_dB;
    SNR_Ref_PAM = Load_BICM_PAM.SNR_dB;

    [a,b,c,d] = size(SINR_FBMC_dB);
    OneTapRate_FBMC = reshape( interp1( SNR_Ref_PAM , BICM_PAM , SINR_FBMC_dB(:) ,'spline'), a,b,c,d);
    R_FBMC = squeeze(mean(sum(sum(OneTapRate_FBMC,1),2),4))/(FBMC.PHY.TimeSpacing*(FBMC.Nr.MCSymbols));

    [a,b,c,d] = size(SINR_FBMC_DFT_dB);
    OneTapRate_FBMC_DFT = reshape( interp1( SNR_Ref_QAM , BICM_QAM , SINR_FBMC_DFT_dB(:) ,'spline'), a,b,c,d);
    R_FBMC_DFT = squeeze(mean(sum(sum(OneTapRate_FBMC_DFT,1),2),4))/(FBMC_DFT.PHY.TimeSpacing*(FBMC_DFT.Nr.MCSymbols));
end


%% Maximize over CQI => perfect feedback
Througput_OFDM      = max(M_Througput_OFDM,[],3);
Througput_OFDMnoCP  = max(M_Througput_OFDMnoCP,[],3);
Througput_FBMC      = max(M_Througput_FBMC,[],3);

Througput_DFT_OFDM     = max(M_Througput_DFT_OFDM,[],3);
Througput_DFT_OFDMnoCP = max(M_Througput_DFT_OFDMnoCP,[],3);
Througput_FBMC_DFT     = max(M_Througput_FBMC_DFT,[],3);


%% Plot Throughput
if CalculateTheory
    figure(13);
    plot(M_SNR_dB,mean(Througput_DFT_OFDM,2)/1e6,'- black'); hold on;
    plot(M_SNR_dB,mean(Througput_DFT_OFDMnoCP,2)/1e6,'- red'); hold on;
    plot(M_SNR_dB,mean(Througput_FBMC_DFT,2)/1e6,'- blue'); hold on;
    plot(M_SNR_dB,mean(Througput_FBMC,2)/1e6,'- magenta'); hold on;
    plot(M_SNR_dB,R_FBMC_DFT/1e6,': blue');
    plot(M_SNR_dB,R_FBMC/1e6,': magenta');  
    
    xlabel('Signal-to-Noise Ratio [dB]');
    ylabel('Achievable Rate, Throughput [Mbit/s]');
    legend({'SC-FDMA (with CP)','SC-FDMA (no CP)','Pruned DFT-s FBMC','FBMC-OQAM','Rate p-DFT-s FBMC','Rate FBMC-OQAM'},'Location','NorthWest');
else
    figure(13);
    plot(M_SNR_dB,mean(Througput_DFT_OFDM,2)/1e6,'- black'); hold on;
    plot(M_SNR_dB,mean(Througput_DFT_OFDMnoCP,2)/1e6,'- red'); hold on;
    plot(M_SNR_dB,mean(Througput_FBMC_DFT,2)/1e6,'- blue'); hold on;
    plot(M_SNR_dB,mean(Througput_FBMC,2)/1e6,'- magenta'); hold on;
    
    xlabel('Signal-to-Noise Ratio [dB]');
    ylabel('Throughput [Mbit/s]');
    legend({'SC-FDMA (with CP)','SC-FDMA (no CP)','Pruned DFT-s FBMC','FBMC-OQAM'},'Location','NorthWest');    
end

figure(100);
plot(M_SNR_dB,mean(Througput_OFDM,2)/1e6,'- black'); hold on;
plot(M_SNR_dB,mean(Througput_OFDMnoCP,2)/1e6,'- red'); hold on;
plot(M_SNR_dB,mean(Througput_FBMC,2)/1e6,'- magenta'); hold on;
xlabel('Signal-to-Noise Ratio [dB]');
ylabel('Throughput [Mbit/s]');
legend({'OFDM (with CP)','OFDM (no CP)','FBMC-OQAM'},'Location','NorthWest');   
title('Comparision of Multicarrier Schemes (OFDM vs FBMC-OQAM)');

figure(13);

%% Save Results
SaveStuff = false;
if SaveStuff 
    Name = ['.\Results\SISO_Throughput_' PowerDelayProfile '_v' int2str(Velocity_kmh) '_' int2str(NrSubcarriers) '.mat'];
    
    meanThrougput_OFDM     = mean(Througput_OFDM,2);
    meanThrougput_OFDMnoCP = mean(Througput_OFDMnoCP,2);
    meanThrougput_FBMC     = mean(Througput_FBMC,2);
    
    meanThrougput_DFT_OFDM     = mean(Througput_DFT_OFDM,2);
    meanThrougput_DFT_OFDMnoCP = mean(Througput_DFT_OFDMnoCP,2);
    meanThrougput_FBMC_DFT     = mean(Througput_FBMC_DFT,2);   
  
    meanThrougput_DFT_OFDM = mean(Througput_DFT_OFDM,2);

    save(Name,...
        'R_FBMC', ...
        'R_FBMC_DFT', ...
        'C_FBMC', ...
        'C_FBMC_DFT', ...        
        'meanThrougput_OFDM', ...
        'meanThrougput_OFDMnoCP', ...
        'meanThrougput_FBMC', ...    
        'meanThrougput_DFT_OFDM', ...
        'meanThrougput_DFT_OFDMnoCP', ...
        'meanThrougput_FBMC_DFT', ... 
        'M_SNR_dB',...
        'NrRepetitions');

end

