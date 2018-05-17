% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================   
% This script simulates the MIMO throughput of pruned DFT spread FBMC, 
% OFDM and FBMC-OQAM. We consider the special case of a low delay spread in
% combination with a low bandwidth, allowing to despread before 
% equalization in pruned DFT spread FBMC. This allows to straightforwardly
% use all known MIMO OFDM methods in FBMC, including ML symbol detection. 

% Allows to reproduce Figure 14 of "Pruned DFT Spread FBMC: Low PAPR, Low 
% Latency, High Spectral Efficiency", R. Nissel and M. Rupp, IEEE
% Transactions on Communications

clear; close all;

%% Parameters
% Simulation
M_SNR_dB            = [0:4:20];                     % Signal-to-Noise Ratio in dB
NrRepetitions       = 4;                            % Number of Monte Carlo repetitions, in the paper: 1000

% FBMC and OFDM Parameters
NrSubcarriers       = 64;                           % Number of subcarriers
SubcarrierSpacing   = 15e3;                         % Subcarrier spacing (15kHz, same as LTE)
CarrierFrequency    = 2.5e9;                        % Carrier Frequency
K_FBMC              = 30;                           % Number of FBMC symbols in time
K_OFDMnoCP          = 15;                           % Number of OFDM symbols in time (no CP)
K_OFDM              = 14;                           % Number of OFDM symbols in time (same as in LTE)
CP_Length           = 1/SubcarrierSpacing/14;       % LTE CP Length in seconds
CP_Length_FBMC_DFT  = 2;                            % CP in the frequency domain for the DFT spreading aproach. Multiple of two: 0, 2, 4... Can usually be set to zero
ICIpowerAware       = false;                        % MEMORY problems => false. If set to true, calculate the channel induced interference power which is later used as addional "noise" power for the LLR calculation => Interference aware receiver.

UseSphereDecoder    = true;                         % Sphere decoder is faster than ML detection                     

SamplingRate        = 15e3*14*96*4;                 % Sampling rate, should approximatly match the power-delay profile of the channel. "*14" due to the CP

% Channel
PowerDelayProfile   = 'TDL-A_9ns'; % 9ns => 10ns after sampling. Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
Velocity_kmh        = 3;                            % Velocity in km/h


% #########################################################################
% % In the paper:
% M_SNR_dB         = [0:1.25:20];
% NrRepetitions    = 1000;
% UseSphereDecoder = false;  
% #########################################################################


%% LTE CQI
% CQI Table: The first column represents the modulation order: 4, 16, 64, 256. The second column represents the code rate (must be between zero and one). Currently, values are chosen according to the (old) LTE standard:
M_CQI = [4 ,  78/1024;...
         4 , 120/1024;...
         4 , 193/1024;...
         4 , 308/1024;...
         4 , 449/1024;...
         4 , 602/1024;...
         16, 378/1024;...
         16, 490/1024;...
         16, 616/1024;...
         64, 466/1024;...
         64, 567/1024;...
         64, 666/1024;...
         64, 772/1024;...
         64, 873/1024;...
         64, 948/1024]; % page 48 of http://www.etsi.org/deliver/etsi_ts/136200_136299/136213/08.08.00_60/ts_136213v080800p.pdf 
% Considered CQI combinations for two layers. We assume that both layers use the same CQI to keep the evaluation time low

if not(strcmp(mexext,'mexw64'))  
    % We use a win64 mexfile for code rates smaller than 1/3 => only works 
    % in 64-bit Windows
    IndexCodeRateSmallerOneThird =  find(M_CQI(:,2)<1/3);
    if  numel(IndexCodeRateSmallerOneThird)>0
        M_CQI(IndexCodeRateSmallerOneThird,:) = [];
        warning('A code rate smaller than 1/3 is only supported in Windows 64-bit => CQI values which contain a code rate smaller than 1/3 are discarded!');
    end    
end
M_Index_CQI_MIMO = [(1:size(M_CQI,1))' (1:size(M_CQI,1))']; 


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


%% Check Number of Samples
if  OFDMnoCP.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = OFDMnoCP.Nr.SamplesTotal;


%% Channel Object
ChannelModel = Channel.FastFading(...
    SamplingRate,...                                % Sampling rate (Samples/s)
    PowerDelayProfile,...                           % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    N,...                                           % Number of total samples
    Velocity_kmh/3.6*CarrierFrequency/2.998e8,...   % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                                     % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    50, ...                                        	% Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    2,...                                          	% Number of transmit antennas
    2,...                                          	% Number of receive antennas
    true ...                                       	% Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );


%% Coding and QAM Objects
for i_cqi = 1:size(M_Index_CQI_MIMO,1)
    QAMModulationOrder_Antenna1 = M_CQI(M_Index_CQI_MIMO(i_cqi,1),1); % Stream 1
    QAMModulationOrder_Antenna2 = M_CQI(M_Index_CQI_MIMO(i_cqi,2),1); % Stream 2
    PAMModulationOrder_Antenna1 = sqrt(QAMModulationOrder_Antenna1);  % Stream 1
    PAMModulationOrder_Antenna2 = sqrt(QAMModulationOrder_Antenna2);  % Stream 2
    
    CodeRate_Antenna1           = M_CQI(M_Index_CQI_MIMO(i_cqi,1),2);
    CodeRate_Anteanna2          = M_CQI(M_Index_CQI_MIMO(i_cqi,2),2);
        
    QAM_Antenna1{i_cqi}         = Modulation.SignalConstellation(QAMModulationOrder_Antenna1,'QAM');
    QAM_Antenna2{i_cqi}         = Modulation.SignalConstellation(QAMModulationOrder_Antenna2,'QAM');    
    PAM_Antenna1{i_cqi}         = Modulation.SignalConstellation(PAMModulationOrder_Antenna1,'PAM');
    PAM_Antenna2{i_cqi}         = Modulation.SignalConstellation(PAMModulationOrder_Antenna2,'PAM');
    
    TurboCoding_Antenna1{i_cqi}  = Coding.TurboCoding(...
        log2(QAMModulationOrder_Antenna1)*NrSubcarriers*K_OFDMnoCP,...                              % Number transmitted bits
        round(CodeRate_Antenna1*log2(QAMModulationOrder_Antenna1)*NrSubcarriers*K_OFDMnoCP)...      % Number data bits
        );  
    TurboCoding_Antenna2{i_cqi}  = Coding.TurboCoding(...
        log2(QAMModulationOrder_Antenna2)*NrSubcarriers*K_OFDMnoCP,...                              % Number transmitted bits
        round(CodeRate_Anteanna2*log2(QAMModulationOrder_Antenna2)*NrSubcarriers*K_OFDMnoCP)...     % Number data bits
        );    
    FBMC_DFT_TurboCoding_Antenna1{i_cqi}  = Coding.TurboCoding(...
        log2(QAMModulationOrder_Antenna1)*(NrSubcarriers-CP_Length_FBMC_DFT)/2*K_FBMC,...                           % Number transmitted bits
        round(CodeRate_Antenna1*log2(QAMModulationOrder_Antenna1)*(NrSubcarriers-CP_Length_FBMC_DFT)/2*K_FBMC)...   % Number data bits
        );  
    FBMC_DFT_TurboCoding_Antenna2{i_cqi}  = Coding.TurboCoding(...
        log2(QAMModulationOrder_Antenna2)*(NrSubcarriers-CP_Length_FBMC_DFT)/2*K_FBMC,...                           % Number transmitted bits
        round(CodeRate_Anteanna2*log2(QAMModulationOrder_Antenna2)*(NrSubcarriers-CP_Length_FBMC_DFT)/2*K_FBMC)...  % Number data bits
        );        
end


%% DFT Matrix
DFTMatrix               = fft(eye(NrSubcarriers))/sqrt(NrSubcarriers);

%% Generate coding matrix for the novel DFT spreading concept
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


%% DFT Spread no CP
% DFT matrix for the coding process
W = fft( eye(NrSubcarriers) ) / sqrt( NrSubcarriers );

% Diagonal elements of the FBMC transmission matrix after DFT spreading despreading
a = abs(diag(W'*D_temp*W));
a = a+randn(size(a))*10^-12; %  randn so that sorting is unique

% Sort a
a_Tilde = sort(a,'descend');

% Get index representing the largest values of a
alpha       = a_Tilde((NrSubcarriers)/2);
Index_Tilde = (a>=alpha);

% Reduced DFT matrix
W_Tilde = W(:,Index_Tilde) ;

% One-tap scaling of the data symbols
b_Tilde = sqrt(2./(a(Index_Tilde)));

% Final coding matrix for one FBMC symbol
C_DFTspread_TX_noCP = W_Tilde*diag(b_Tilde);
C_DFTspread_RX_noCP = W_Tilde*diag(b_Tilde);



%% ICI Power Aware 
if not(ICIpowerAware)
    PI_OFDMnoCP = 0;
    PI_OFDM     = 0;
    PI_FBMC     = 0;
else
    % we reduce the sampling rate and number of subcarriers so that the vectorized channel correlation matrix can be calculated
    SamplingRate_Smaller    = SamplingRate/4;
    NrSubcarriersTemp       = floor(SamplingRate_Smaller/SubcarrierSpacing);
    OFDMnoCP_Temp = Modulation.OFDM(...
        NrSubcarriersTemp,...                           % Number of  active subcarriers
        3,...                                           % Number of OFDM Symbols
        SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
        SamplingRate_Smaller,...                        % Sampling rate (Samples/s)
        0,...                                           % Intermediate frequency first subcarrier (Hz)
        false,...                                       % Transmit real valued signal
        0, ...                                          % Cyclic prefix length (s) 1/SubcarrierSpacing/(K/2-1)
        0 ...                                           % Zero guard length (s)
        );
    ChannelModel_Temp = Channel.FastFading(...
        SamplingRate_Smaller,...                        % Sampling rate (Samples/s)
        PowerDelayProfile,...                           % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
        OFDMnoCP_Temp.Nr.SamplesTotal,...               % Number of total samples
        Velocity_kmh/3.6*CarrierFrequency/2.998e8,...   % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
        'Jakes',...                                     % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
        50, ...                                         % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
        1,...                                           % Number of transmit antennas
        1,...                                           % Number of receive antennas
        false ...                                       % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
        );
    [PS_OFDMnoCP,PI_OFDMnoCP] = OFDMnoCP_Temp.GetSignalAndInterferencePowerQAM(...
        ChannelModel_Temp.GetCorrelationMatrix,eye(NrSubcarriersTemp*3),0, round(NrSubcarriersTemp/2),2);
    disp(['OFDM(noCP) SIR: ' int2str(10*log10(PS_OFDMnoCP/PI_OFDMnoCP)) 'dB']);

    FBMC_Temp = Modulation.FBMC(...
        7,...                                           % Number of subcarriers
        7,...                                           % Number of FBMC symbols
        SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
        SamplingRate_Smaller,...                        % Sampling rate (Samples/s)
        0,...                                           % Intermediate frequency first subcarrier (Hz)
        false,...                                       % Transmit real valued signal
        'Hermite-OQAM',...                              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
        4, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
        0, ...                                          % Initial phase shift
        true ...                                        % Polyphase implementation
        );
    ChannelModel_Temp = Channel.FastFading(...
        SamplingRate_Smaller,...                                            % Sampling rate (Samples/s)
        PowerDelayProfile,...                                               % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
        FBMC_Temp.Nr.SamplesTotal,...                                       % Number of total samples
        Velocity_kmh/3.6*CarrierFrequency/2.998e8,...                       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
        'Jakes',...                                                         % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
        50, ...                                                             % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
        1,...                                                               % Number of transmit antennas
        1,...                                                               % Number of receive antennas
        false ...                                                           % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
        );
    [PS_FBMC,PI_FBMC] = FBMC_Temp.GetSignalAndInterferencePowerOQAM(...
        ChannelModel_Temp.GetCorrelationMatrix,eye(7*7),0, 4,4);    
    disp(['FBMC-OQAM SIR:  ' int2str(10*log10(PS_FBMC/PI_FBMC)) 'dB']);

end

InterferenceMatrix = kron(eye(FBMC.Nr.MCSymbols),C_DFTspread_RX)'*FBMC_DFT.GetFBMCMatrix*kron(eye(FBMC.Nr.MCSymbols),C_DFTspread_TX)/2;
PI_DFT_Intrinsic = mean(sum(abs(InterferenceMatrix-eye(size(InterferenceMatrix,1))).^2,2));

disp(['Intrinsic interference FBMC DFT Spread: ' int2str(10*log10(1/PI_DFT_Intrinsic)) 'dB']);

    
%% Preallocate
Throughput_OFDMnoCP     = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));
Throughput_OFDMnoCP_DFT = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));
Throughput_FBMC         = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));
Throughput_FBMC_DFT     = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));

Throughput_OFDMnoCP_ML      = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));
Throughput_FBMC_DFT_ML      = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));
Throughput_FBMC_DFTnoCP_ML  = nan(length(M_SNR_dB),NrRepetitions,size(M_Index_CQI_MIMO,1));

disp('The simulation may take a while ... ');
% parfor i_rep = 1:NrRepetitions                        % PARFOR
for i_rep = 1:NrRepetitions                             % Conventional FOR loop
    tic;                                            
    %% New Channel Realization
    ChannelModel.NewRealization;

    %% Perfect Channel Knowledge
    h_OFDMnoCP = 1/sqrt(2) * ChannelModel.GetTransferFunction( OFDMnoCP.GetTimeIndexMidPos , OFDMnoCP.Implementation.FFTSize , OFDMnoCP.Implementation.IntermediateFrequency+(1:OFDMnoCP.Nr.Subcarriers) );
    H_OFDMnoCP = permute(h_OFDMnoCP,[3 4 1 2]);

    h_FBMC     = 1/sqrt(2) * ChannelModel.GetTransferFunction( FBMC.GetTimeIndexMidPos , FBMC.Implementation.FFTSize , FBMC.Implementation.IntermediateFrequency+(1:FBMC.Nr.Subcarriers) );
    H_FBMC     = permute(h_FBMC,[3 4 1 2]);

    H_FBMC_FrequencyMean     = repmat(mean(H_FBMC,3),[1 1 (NrSubcarriers-CP_Length_FBMC_DFT)/2 1]);
    H_FBMCnoCP_FrequencyMean = repmat(mean(H_FBMC,3),[1 1 NrSubcarriers/2 1]);

    %% Preallocate2
    Throughput_OFDMnoCP_OneRealization     = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));
    Throughput_OFDMnoCP_DFT_OneRealization = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));
    Throughput_FBMC_OneRealization         = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));
    Throughput_FBMC_DFT_OneRealization     = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));

    Throughput_OFDMnoCP_ML_OneRealization      = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));
    Throughput_FBMC_DFT_ML_OneRealization      = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));
    Throughput_FBMC_DFTnoCP_ML_OneRealization  = nan(length(M_SNR_dB),size(M_Index_CQI_MIMO,1));

    % Simulate over different modulation orders and code rates
    for i_cqi = 1:size(M_Index_CQI_MIMO,1);

        % Generate Bit Stream
        BinaryDataStream_Antenna1          = randi([0 1],TurboCoding_Antenna1{i_cqi}.NrDataBits,1);
        BinaryDataStream_Antenna2          = randi([0 1],TurboCoding_Antenna2{i_cqi}.NrDataBits,1);

        BinaryDataStream_FBMC_DFT_Antenna1 = randi([0 1],FBMC_DFT_TurboCoding_Antenna1{i_cqi}.NrDataBits,1);  % Maybe a reduced bit rate due to the CP (but not necesarry most of the time!)
        BinaryDataStream_FBMC_DFT_Antenna2 = randi([0 1],FBMC_DFT_TurboCoding_Antenna2{i_cqi}.NrDataBits,1);  % Maybe a reduced bit rate due to the CP (but not necesarry most of the time!)

        % Channel Coding
        TurboCoding_Antenna1{i_cqi}.UpdateInterleaving;
        TurboCoding_Antenna2{i_cqi}.UpdateInterleaving;
        FBMC_DFT_TurboCoding_Antenna1{i_cqi}.UpdateInterleaving;
        FBMC_DFT_TurboCoding_Antenna2{i_cqi}.UpdateInterleaving;

        CodedBits_Antenna1          = TurboCoding_Antenna1{i_cqi}.TurboEncoder(BinaryDataStream_Antenna1);
        CodedBits_Antenna2          = TurboCoding_Antenna2{i_cqi}.TurboEncoder(BinaryDataStream_Antenna2);
        CodedBits_FBMC_DFT_Antenna1 = FBMC_DFT_TurboCoding_Antenna1{i_cqi}.TurboEncoder(BinaryDataStream_FBMC_DFT_Antenna1);
        CodedBits_FBMC_DFT_Antenna2 = FBMC_DFT_TurboCoding_Antenna2{i_cqi}.TurboEncoder(BinaryDataStream_FBMC_DFT_Antenna2);


        % Bit Interleaving
        BitInterleaving_Antenna1          = randperm(TurboCoding_Antenna1{i_cqi}.NrCodedBits);
        BitInterleaving_Antenna2          = randperm(TurboCoding_Antenna2{i_cqi}.NrCodedBits);
        BitInterleaving_FBMC_DFT_Antenna1 = randperm(FBMC_DFT_TurboCoding_Antenna1{i_cqi}.NrCodedBits);
        BitInterleaving_FBMC_DFT_Antenna2 = randperm(FBMC_DFT_TurboCoding_Antenna2{i_cqi}.NrCodedBits);

        CodedBits_Antenna1          = CodedBits_Antenna1(BitInterleaving_Antenna1);
        CodedBits_Antenna2          = CodedBits_Antenna2(BitInterleaving_Antenna2);
        CodedBits_FBMC_DFT_Antenna1 = CodedBits_FBMC_DFT_Antenna1(BitInterleaving_FBMC_DFT_Antenna1);
        CodedBits_FBMC_DFT_Antenna2 = CodedBits_FBMC_DFT_Antenna2(BitInterleaving_FBMC_DFT_Antenna2);


        % Map Bit Stream to Symbols
        x_OFDMnoCP_Antenna1 = reshape(QAM_Antenna1{i_cqi}.Bit2Symbol(CodedBits_Antenna1),NrSubcarriers,K_OFDMnoCP);
        x_OFDMnoCP_Antenna2 = reshape(QAM_Antenna2{i_cqi}.Bit2Symbol(CodedBits_Antenna2),NrSubcarriers,K_OFDMnoCP);

        x_FBMC_DFT_Antenna1 = reshape(QAM_Antenna1{i_cqi}.Bit2Symbol(CodedBits_FBMC_DFT_Antenna1),(NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC);
        x_FBMC_DFT_Antenna2 = reshape(QAM_Antenna2{i_cqi}.Bit2Symbol(CodedBits_FBMC_DFT_Antenna2),(NrSubcarriers-CP_Length_FBMC_DFT)/2,K_FBMC);

        x_FBMC_Antenna1     = reshape(PAM_Antenna1{i_cqi}.Bit2Symbol(CodedBits_Antenna1),NrSubcarriers,K_FBMC);
        x_FBMC_Antenna2     = reshape(PAM_Antenna2{i_cqi}.Bit2Symbol(CodedBits_Antenna2),NrSubcarriers,K_FBMC);


        % Generate Transmit Signal in the Time Domain
        s_OFDMnoCP_Antenna1     = OFDMnoCP.Modulation(x_OFDMnoCP_Antenna1)/sqrt(2);
        s_OFDMnoCP_Antenna2     = OFDMnoCP.Modulation(x_OFDMnoCP_Antenna2)/sqrt(2);

        s_FBMC_Antenna1         = FBMC.Modulation(x_FBMC_Antenna1)/sqrt(2);
        s_FBMC_Antenna2         = FBMC.Modulation(x_FBMC_Antenna2)/sqrt(2);

        s_DFT_OFDMnoCP_Antenna1 = OFDMnoCP.Modulation(DFTMatrix*x_OFDMnoCP_Antenna1)/sqrt(2);
        s_DFT_OFDMnoCP_Antenna2 = OFDMnoCP.Modulation(DFTMatrix*x_OFDMnoCP_Antenna2)/sqrt(2);

        s_FBMC_DFT_Antenna1     = FBMC_DFT.Modulation(C_DFTspread_TX*x_FBMC_DFT_Antenna1)/sqrt(2);
        s_FBMC_DFT_Antenna2     = FBMC_DFT.Modulation(C_DFTspread_TX*x_FBMC_DFT_Antenna2)/sqrt(2);

        s_FBMC_DFTnoCP_Antenna1  = FBMC_DFT.Modulation(C_DFTspread_TX_noCP*reshape(x_OFDMnoCP_Antenna1,NrSubcarriers/2,[]))/sqrt(2);
        s_FBMC_DFTnoCP_Antenna2  = FBMC_DFT.Modulation(C_DFTspread_TX_noCP*reshape(x_OFDMnoCP_Antenna2,NrSubcarriers/2,[]))/sqrt(2);


        % Channel
        r_OFDMnoCP_noNoise     = ChannelModel.Convolution([s_OFDMnoCP_Antenna1 s_OFDMnoCP_Antenna2]);
        r_FBMC_noNoise         = ChannelModel.Convolution([s_FBMC_Antenna1 s_FBMC_Antenna2]);

        r_DFT_OFDMnoCP_noNoise = ChannelModel.Convolution([s_DFT_OFDMnoCP_Antenna1 s_DFT_OFDMnoCP_Antenna2]);
        r_FBMC_DFT_noNoise     = ChannelModel.Convolution([s_FBMC_DFT_Antenna1 s_FBMC_DFT_Antenna2]);

        r_FBMC_DFTnoCP_noNoise     = ChannelModel.Convolution([s_FBMC_DFTnoCP_Antenna1 s_FBMC_DFTnoCP_Antenna2]);

        
        % Simulate over different noise values
        for i_SNR = 1:length(M_SNR_dB)
            SNR_dB  = M_SNR_dB(i_SNR); 
            Pn      = 10^(-SNR_dB/10);
            Pn_time = SamplingRate/(SubcarrierSpacing*NrSubcarriers)*10^(-SNR_dB/10); 

            % Add Noise
            noise_Antenna1 = sqrt(Pn_time/2)*(randn(N,1)+1j*randn(N,1));
            noise_Antenna2 = sqrt(Pn_time/2)*(randn(N,1)+1j*randn(N,1));

            r_OFDMnoCP_Antenna1     = r_OFDMnoCP_noNoise(:,1) + noise_Antenna1;
            r_OFDMnoCP_Antenna2     = r_OFDMnoCP_noNoise(:,2) + noise_Antenna2;

            r_FBMC_Antenna1         = r_FBMC_noNoise(:,1) + noise_Antenna1;
            r_FBMC_Antenna2         = r_FBMC_noNoise(:,2) + noise_Antenna2;

            r_DFT_OFDMnoCP_Antenna1 = r_DFT_OFDMnoCP_noNoise(:,1) + noise_Antenna1;
            r_DFT_OFDMnoCP_Antenna2 = r_DFT_OFDMnoCP_noNoise(:,2) + noise_Antenna2;

            r_FBMC_DFT_Antenna1     = r_FBMC_DFT_noNoise(:,1) + noise_Antenna1;
            r_FBMC_DFT_Antenna2     = r_FBMC_DFT_noNoise(:,2) + noise_Antenna2;

            r_FBMC_DFTnoCP_Antenna1     = r_FBMC_DFTnoCP_noNoise(:,1) + noise_Antenna1;
            r_FBMC_DFTnoCP_Antenna2     = r_FBMC_DFTnoCP_noNoise(:,2) + noise_Antenna2;

            % Received Symbols (Demodulation)
            y_OFDMnoCP_Antenna1     = OFDMnoCP.Demodulation(r_OFDMnoCP_Antenna1);
            y_OFDMnoCP_Antenna2     = OFDMnoCP.Demodulation(r_OFDMnoCP_Antenna2);

            y_FBMC_Antenna1         = FBMC.Demodulation(r_FBMC_Antenna1);
            y_FBMC_Antenna2         = FBMC.Demodulation(r_FBMC_Antenna2);

            y_DFT_OFDMnoCP_Antenna1 = OFDMnoCP.Demodulation(r_DFT_OFDMnoCP_Antenna1);
            y_DFT_OFDMnoCP_Antenna2 = OFDMnoCP.Demodulation(r_DFT_OFDMnoCP_Antenna2);

            y_FBMC_DFT_Antenna1     = FBMC_DFT.Demodulation(r_FBMC_DFT_Antenna1);
            y_FBMC_DFT_Antenna2     = FBMC_DFT.Demodulation(r_FBMC_DFT_Antenna2);

            y_FBMC_DFTnoCP_Antenna1     = FBMC_DFT.Demodulation(r_FBMC_DFTnoCP_Antenna1);
            y_FBMC_DFTnoCP_Antenna2     = FBMC_DFT.Demodulation(r_FBMC_DFTnoCP_Antenna2);


            % MMSE Equalizer and LLR Calculation
                % OFDM no CP
            [~,x_est_OFDMnoCP,NoiseScaling,UnbiasedScaling] = QAM_Antenna1{i_cqi}.LLR_MIMO_MMSE([y_OFDMnoCP_Antenna1(:).';y_OFDMnoCP_Antenna2(:).'],H_OFDMnoCP(:,:,:), Pn );
            LLR_OFDMnoCP_Antenna1 = QAM_Antenna1{i_cqi}.LLR_AWGN(x_est_OFDMnoCP(:,1)./UnbiasedScaling(:,1),NoiseScaling(:,1)./UnbiasedScaling(:,1).^2);
            LLR_OFDMnoCP_Antenna2 = QAM_Antenna2{i_cqi}.LLR_AWGN(x_est_OFDMnoCP(:,2)./UnbiasedScaling(:,2),NoiseScaling(:,2)./UnbiasedScaling(:,2).^2);

            [~,x_est_DFT_OFDMnoCP_BeforeDespreading,NoiseScaling_BeforeDFT,UnbiasedScaling_BeforeDFT]= QAM_Antenna1{i_cqi}.LLR_MIMO_MMSE([y_DFT_OFDMnoCP_Antenna1(:).';y_DFT_OFDMnoCP_Antenna2(:).'],H_OFDMnoCP(:,:,:), Pn );
            x_est_DFT_OFDMnoCP        = reshape(DFTMatrix'*reshape(x_est_DFT_OFDMnoCP_BeforeDespreading,OFDMnoCP.Nr.Subcarriers,OFDMnoCP.Nr.MCSymbols*2),[],2);
            NoiseScaling              = reshape(repmat(mean(reshape(NoiseScaling_BeforeDFT,OFDMnoCP.Nr.Subcarriers,OFDMnoCP.Nr.MCSymbols*2),1),OFDMnoCP.Nr.Subcarriers,1),[],2);
            UnbiasedScaling           = reshape(repmat(mean(reshape(UnbiasedScaling_BeforeDFT,OFDMnoCP.Nr.Subcarriers,OFDMnoCP.Nr.MCSymbols*2),1),OFDMnoCP.Nr.Subcarriers,1),[],2);
            LLR_OFDMnoCP_DFT_Antenna1 = QAM_Antenna1{i_cqi}.LLR_AWGN(x_est_DFT_OFDMnoCP(:,1)./UnbiasedScaling(:,1),NoiseScaling(:,1)./UnbiasedScaling(:,1).^2);
            LLR_OFDMnoCP_DFT_Antenna2 = QAM_Antenna2{i_cqi}.LLR_AWGN(x_est_DFT_OFDMnoCP(:,2)./UnbiasedScaling(:,2),NoiseScaling(:,2)./UnbiasedScaling(:,2).^2);


                % FBMC
            [~,x_est_FBMC,NoiseScaling,UnbiasedScaling] = PAM_Antenna1{i_cqi}.LLR_MIMO_MMSE([y_FBMC_Antenna1(:).';y_FBMC_Antenna2(:).'],H_FBMC(:,:,:),Pn );
            LLR_FBMC_Antenna1 = PAM_Antenna1{i_cqi}.LLR_AWGN(real(x_est_FBMC(:,1))./UnbiasedScaling(:,1),NoiseScaling(:,1)./UnbiasedScaling(:,1).^2);
            LLR_FBMC_Antenna2 = PAM_Antenna2{i_cqi}.LLR_AWGN(real(x_est_FBMC(:,2))./UnbiasedScaling(:,2),NoiseScaling(:,2)./UnbiasedScaling(:,2).^2);

            [~,x_est_DFT_FBMC_BeforeDespreading,NoiseScaling_BeforeDFT,UnbiasedScaling_BeforeDFT] = QAM_Antenna1{i_cqi}.LLR_MIMO_MMSE([y_FBMC_DFT_Antenna1(:).';y_FBMC_DFT_Antenna2(:).'],H_FBMC(:,:,:), Pn );
            x_est_DFT_FBMC        = reshape(C_DFTspread_RX'*reshape(x_est_DFT_FBMC_BeforeDespreading,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols*2),[],2)/2;
            NoiseScaling          = reshape(repmat(mean(reshape(NoiseScaling_BeforeDFT,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols*2),1),(FBMC.Nr.Subcarriers-CP_Length_FBMC_DFT)/2,1),[],2);
            UnbiasedScaling       = reshape(repmat(mean(reshape(UnbiasedScaling_BeforeDFT,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols*2),1),(FBMC.Nr.Subcarriers-CP_Length_FBMC_DFT)/2,1),[],2);
            LLR_FBMC_DFT_Antenna1 = QAM_Antenna1{i_cqi}.LLR_AWGN(x_est_DFT_FBMC(:,1)./UnbiasedScaling(:,1),NoiseScaling(:,1)./UnbiasedScaling(:,1).^2);
            LLR_FBMC_DFT_Antenna2 = QAM_Antenna2{i_cqi}.LLR_AWGN(x_est_DFT_FBMC(:,2)./UnbiasedScaling(:,2),NoiseScaling(:,2)./UnbiasedScaling(:,2).^2);

            % ML Detection
                % OFDM
            if UseSphereDecoder
                LLR_OFDMnoCP_ML = QAM_Antenna1{i_cqi}.LLR_MIMO_Sphere([y_OFDMnoCP_Antenna1(:).';y_OFDMnoCP_Antenna2(:).'],H_OFDMnoCP(:,:,:), (Pn + PI_OFDMnoCP));
            else
                LLR_OFDMnoCP_ML = QAM_Antenna1{i_cqi}.LLR_MIMO_ML([y_OFDMnoCP_Antenna1(:).';y_OFDMnoCP_Antenna2(:).'],H_OFDMnoCP(:,:,:), (Pn)*repmat(eye(2),[1 1 size(H_OFDMnoCP(:,:,:),3)]));
            end
            LLR_OFDMnoCP_ML_Antenna1 = LLR_OFDMnoCP_ML(:,1);
            LLR_OFDMnoCP_ML_Antenna2 = LLR_OFDMnoCP_ML(:,2);

                % FBMC
            y_FBMC_DFT_Despread_Antenna1 = C_DFTspread_RX'*y_FBMC_DFT_Antenna1/2;
            y_FBMC_DFT_Despread_Antenna2 = C_DFTspread_RX'*y_FBMC_DFT_Antenna2/2;
            % Rn = nan(2,2,size(H_FBMC_FrequencyMean(:,:,:),3));
            % for i_Rn = 1:size(H_FBMC_FrequencyMean(:,:,:),3)
            %     Rn(:,:,i_Rn) = diag(sum(abs(H_FBMC_FrequencyMean(:,:,i_Rn)).^2,2)*PI_DFT_Intrinsic + (Pn + PI_FBMC*2));
            % end
            Rn = (Pn + PI_DFT_Intrinsic)*repmat(eye(2),[1 1 size(H_FBMC_FrequencyMean(:,:,:),3)]);
            if UseSphereDecoder
                LLR_FBMC_DFT_ML = QAM_Antenna1{i_cqi}.LLR_MIMO_Sphere([y_FBMC_DFT_Despread_Antenna1(:).';y_FBMC_DFT_Despread_Antenna2(:).'],H_FBMC_FrequencyMean(:,:,:), (Pn + PI_FBMC*2) );
            else
                LLR_FBMC_DFT_ML = QAM_Antenna1{i_cqi}.LLR_MIMO_ML([y_FBMC_DFT_Despread_Antenna1(:).';y_FBMC_DFT_Despread_Antenna2(:).'],H_FBMC_FrequencyMean(:,:,:), Rn );
            end
            LLR_FBMC_DFT_ML_Antenna1 = LLR_FBMC_DFT_ML(:,1);
            LLR_FBMC_DFT_ML_Antenna2 = LLR_FBMC_DFT_ML(:,2);


                % FBMC no CP
            y_FBMC_DFTnoCP_Despread_Antenna1 = C_DFTspread_RX_noCP'*y_FBMC_DFTnoCP_Antenna1/2;
            y_FBMC_DFTnoCP_Despread_Antenna2 = C_DFTspread_RX_noCP'*y_FBMC_DFTnoCP_Antenna2/2;
            % Rn = nan(2,2,size(H_FBMC_FrequencyMean(:,:,:),3));
            % for i_Rn = 1:size(H_FBMC_FrequencyMean(:,:,:),3)
            %     Rn(:,:,i_Rn) = diag(sum(abs(H_FBMC_FrequencyMean(:,:,i_Rn)).^2,2)*PI_DFT_Intrinsic + (Pn + PI_FBMC*2));
            % end
            Rn = (Pn + PI_DFT_Intrinsic)*repmat(eye(2),[1 1 size(H_FBMCnoCP_FrequencyMean(:,:,:),3)]);
            if UseSphereDecoder
                LLR_FBMC_DFTnoCP_ML = QAM_Antenna1{i_cqi}.LLR_MIMO_Sphere([y_FBMC_DFTnoCP_Despread_Antenna1(:).';y_FBMC_DFTnoCP_Despread_Antenna2(:).'],H_FBMCnoCP_FrequencyMean(:,:,:), (Pn + PI_FBMC*2) );
            else
                LLR_FBMC_DFTnoCP_ML = QAM_Antenna1{i_cqi}.LLR_MIMO_ML([y_FBMC_DFTnoCP_Despread_Antenna1(:).';y_FBMC_DFTnoCP_Despread_Antenna2(:).'],H_FBMCnoCP_FrequencyMean(:,:,:), Rn );
            end
            LLR_FBMC_DFTnoCP_ML_Antenna1     = LLR_FBMC_DFTnoCP_ML(:,1);
            LLR_FBMC_DFTnoCP_ML_Antenna2     = LLR_FBMC_DFTnoCP_ML(:,2);


            % De-interleaving
            LLR_OFDMnoCP_Antenna1(BitInterleaving_Antenna1) = LLR_OFDMnoCP_Antenna1;
            LLR_OFDMnoCP_Antenna2(BitInterleaving_Antenna2) = LLR_OFDMnoCP_Antenna2;

            LLR_OFDMnoCP_DFT_Antenna1(BitInterleaving_Antenna1) = LLR_OFDMnoCP_DFT_Antenna1;
            LLR_OFDMnoCP_DFT_Antenna2(BitInterleaving_Antenna2) = LLR_OFDMnoCP_DFT_Antenna2;

            LLR_FBMC_Antenna1(BitInterleaving_Antenna1) = LLR_FBMC_Antenna1;
            LLR_FBMC_Antenna2(BitInterleaving_Antenna2) = LLR_FBMC_Antenna2;

            LLR_FBMC_DFT_Antenna1(BitInterleaving_FBMC_DFT_Antenna1) = LLR_FBMC_DFT_Antenna1;
            LLR_FBMC_DFT_Antenna2(BitInterleaving_FBMC_DFT_Antenna2) = LLR_FBMC_DFT_Antenna2;

            LLR_OFDMnoCP_ML_Antenna1(BitInterleaving_Antenna1) = LLR_OFDMnoCP_ML_Antenna1;
            LLR_OFDMnoCP_ML_Antenna2(BitInterleaving_Antenna2) = LLR_OFDMnoCP_ML_Antenna2;

            LLR_FBMC_DFT_ML_Antenna1(BitInterleaving_FBMC_DFT_Antenna1) = LLR_FBMC_DFT_ML_Antenna1;
            LLR_FBMC_DFT_ML_Antenna2(BitInterleaving_FBMC_DFT_Antenna2) = LLR_FBMC_DFT_ML_Antenna2;

            LLR_FBMC_DFTnoCP_ML_Antenna1(BitInterleaving_Antenna1) = LLR_FBMC_DFTnoCP_ML_Antenna1;
            LLR_FBMC_DFTnoCP_ML_Antenna2(BitInterleaving_Antenna2) = LLR_FBMC_DFTnoCP_ML_Antenna2;

            % Turbo-Decoding
            DecodedBits_OFDMnoCP_Antenna1 = TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_OFDMnoCP_Antenna1);
            DecodedBits_OFDMnoCP_Antenna2 = TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_OFDMnoCP_Antenna2);

            DecodedBits_OFDMnoCP_DFT_Antenna1 = TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_OFDMnoCP_DFT_Antenna1);
            DecodedBits_OFDMnoCP_DFT_Antenna2 = TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_OFDMnoCP_DFT_Antenna2);

            DecodedBits_FBMC_Antenna1 = TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_FBMC_Antenna1);
            DecodedBits_FBMC_Antenna2 = TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_FBMC_Antenna2);

            DecodedBits_FBMC_DFT_Antenna1 = FBMC_DFT_TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_FBMC_DFT_Antenna1);
            DecodedBits_FBMC_DFT_Antenna2 = FBMC_DFT_TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_FBMC_DFT_Antenna2);

            DecodedBits_OFDMnoCP_ML_Antenna1 = TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_OFDMnoCP_ML_Antenna1);
            DecodedBits_OFDMnoCP_ML_Antenna2 = TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_OFDMnoCP_ML_Antenna2);

            DecodedBits_FBMC_DFT_ML_Antenna1 = FBMC_DFT_TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_FBMC_DFT_ML_Antenna1);
            DecodedBits_FBMC_DFT_ML_Antenna2 = FBMC_DFT_TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_FBMC_DFT_ML_Antenna2);

            DecodedBits_FBMC_DFTnoCP_ML_Antenna1 = TurboCoding_Antenna1{i_cqi}.TurboDecoder(LLR_FBMC_DFTnoCP_ML_Antenna1);
            DecodedBits_FBMC_DFTnoCP_ML_Antenna2 = TurboCoding_Antenna2{i_cqi}.TurboDecoder(LLR_FBMC_DFTnoCP_ML_Antenna2);


            % Throughput
            KT_OFDMnoCP = (OFDMnoCP.Nr.MCSymbols*OFDMnoCP.PHY.TimeSpacing);
            Throughput_OFDMnoCP_OneRealization(i_SNR,i_cqi)     = (all(DecodedBits_OFDMnoCP_Antenna1==BinaryDataStream_Antenna1)*length(DecodedBits_OFDMnoCP_Antenna1)+all(DecodedBits_OFDMnoCP_Antenna2==BinaryDataStream_Antenna2)*length(DecodedBits_OFDMnoCP_Antenna2))/KT_OFDMnoCP;
            Throughput_OFDMnoCP_DFT_OneRealization(i_SNR,i_cqi) = (all(DecodedBits_OFDMnoCP_DFT_Antenna1==BinaryDataStream_Antenna1)*length(DecodedBits_OFDMnoCP_DFT_Antenna1)+all(DecodedBits_OFDMnoCP_DFT_Antenna2==BinaryDataStream_Antenna2)*length(DecodedBits_OFDMnoCP_DFT_Antenna2))/KT_OFDMnoCP;
            Throughput_OFDMnoCP_ML_OneRealization(i_SNR,i_cqi)  = (all(DecodedBits_OFDMnoCP_ML_Antenna1==BinaryDataStream_Antenna1)*length(DecodedBits_OFDMnoCP_ML_Antenna1)+all(DecodedBits_OFDMnoCP_ML_Antenna2==BinaryDataStream_Antenna2)*length(DecodedBits_OFDMnoCP_ML_Antenna2))/KT_OFDMnoCP;

            KT_FBMC = (FBMC.Nr.MCSymbols*FBMC.PHY.TimeSpacing);
            Throughput_FBMC_OneRealization(i_SNR,i_cqi)             = (all(DecodedBits_FBMC_Antenna1==BinaryDataStream_Antenna1)*length(DecodedBits_FBMC_Antenna1)+all(DecodedBits_FBMC_Antenna2==BinaryDataStream_Antenna2)*length(DecodedBits_FBMC_Antenna2))/KT_FBMC;
            Throughput_FBMC_DFT_OneRealization(i_SNR,i_cqi)         = (all(DecodedBits_FBMC_DFT_Antenna1==BinaryDataStream_FBMC_DFT_Antenna1)*length(DecodedBits_FBMC_DFT_Antenna1)+all(DecodedBits_FBMC_DFT_Antenna2==BinaryDataStream_FBMC_DFT_Antenna2)*length(DecodedBits_FBMC_DFT_Antenna2))/KT_FBMC;
            Throughput_FBMC_DFT_ML_OneRealization(i_SNR,i_cqi)      = (all(DecodedBits_FBMC_DFT_ML_Antenna1==BinaryDataStream_FBMC_DFT_Antenna1)*length(DecodedBits_FBMC_DFT_ML_Antenna1)+all(DecodedBits_FBMC_DFT_ML_Antenna2==BinaryDataStream_FBMC_DFT_Antenna2)*length(DecodedBits_FBMC_DFT_ML_Antenna2))/KT_FBMC;
            Throughput_FBMC_DFTnoCP_ML_OneRealization(i_SNR,i_cqi)  = (all(DecodedBits_FBMC_DFTnoCP_ML_Antenna1==BinaryDataStream_Antenna1)*length(DecodedBits_FBMC_DFTnoCP_ML_Antenna1)+all(DecodedBits_FBMC_DFTnoCP_ML_Antenna2==BinaryDataStream_Antenna2)*length(DecodedBits_FBMC_DFTnoCP_ML_Antenna2))/KT_FBMC;

        end

    end

    Throughput_OFDMnoCP(:,i_rep,:)     = Throughput_OFDMnoCP_OneRealization;
    Throughput_OFDMnoCP_DFT(:,i_rep,:) = Throughput_OFDMnoCP_DFT_OneRealization;
    Throughput_FBMC(:,i_rep,:)         = Throughput_FBMC_OneRealization;
    Throughput_FBMC_DFT(:,i_rep,:)     = Throughput_FBMC_DFT_OneRealization;

    Throughput_OFDMnoCP_ML(:,i_rep,:)      = Throughput_OFDMnoCP_ML_OneRealization;
    Throughput_FBMC_DFT_ML(:,i_rep,:)      = Throughput_FBMC_DFT_ML_OneRealization;
    Throughput_FBMC_DFTnoCP_ML(:,i_rep,:)  = Throughput_FBMC_DFTnoCP_ML_OneRealization;

    disp(i_rep);
    TimePassed = toc;
    disp(['Realization ' int2str(i_rep) ' took ' int2str(TimePassed/60) 'minutes, total simulation time: ' int2str(TimePassed/60 * NrRepetitions) 'minutes' ]);
end

MaxCQI_Throughput_OFDMnoCP     = max(Throughput_OFDMnoCP,[],3);
MaxCQI_Throughput_OFDMnoCP_DFT = max(Throughput_OFDMnoCP_DFT,[],3);

MaxCQI_Throughput_FBMC         = max(Throughput_FBMC,[],3);
MaxCQI_Throughput_FBMC_DFT     = max(Throughput_FBMC_DFT,[],3);

MaxCQI_Throughput_OFDMnoCP_ML      = max(Throughput_OFDMnoCP_ML,[],3);
MaxCQI_Throughput_FBMC_DFT_ML      = max(Throughput_FBMC_DFT_ML,[],3);
MaxCQI_Throughput_FBMC_DFTnoCP_ML  = max(Throughput_FBMC_DFTnoCP_ML,[],3);


%% Plot all
figure(14);
plot(M_SNR_dB,mean(MaxCQI_Throughput_OFDMnoCP_ML,2)/1e6,'red'); hold on;
plot(M_SNR_dB,mean(MaxCQI_Throughput_FBMC_DFTnoCP_ML,2)/1e6,'blue'); hold on;
plot(M_SNR_dB,mean(MaxCQI_Throughput_FBMC_DFT_ML,2)/1e6,': blue'); hold on;
plot(M_SNR_dB,mean(MaxCQI_Throughput_FBMC,2)/1e6,'magenta'); hold on;
xlabel('Signal-to-Noise Ratio');
ylabel('Throughput [Mbit/s]');
legend({'OFDM (no CP)','Pruned DFT-s FBMC (L_C_P = 0)','Pruned DFT-s FBMC (L_C_P = 2)','FBMC-OQAM (MMSE)'},'Location','NorthWest');


SaveStuff = false;
if SaveStuff 
    Name = ['.\Results\MIMO_Throughput_ML_' PowerDelayProfile '_v' int2str(Velocity_kmh) '_' int2str(NrSubcarriers) '.mat'];
    save(Name,...
        'Throughput_OFDMnoCP',...
        'Throughput_OFDMnoCP_DFT',...
        'Throughput_FBMC',...
        'Throughput_FBMC_DFT',...
        'Throughput_OFDMnoCP_ML',...
        'Throughput_FBMC_DFT_ML',...
        'Throughput_FBMC_DFTnoCP_ML',...        
        'M_SNR_dB',...
        'PowerDelayProfile',...
        'Velocity_kmh',...
        'NrSubcarriers');
end

