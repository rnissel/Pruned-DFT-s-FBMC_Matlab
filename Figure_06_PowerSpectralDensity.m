% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================     
% This script calculates the power spectal density of FBMC, OFDM and pruned
% DFT spread FBMC

% Reproduces Figure 6 of "Pruned DFT Spread FBMC: Low PAPR, Low Latency, 
% High Spectral Efficiency", R. Nissel and M. Rupp, IEEE Transactions on 
% Communications 

clear; close all;

%% Parameters
NrSubcarriers = 128;                            % Number of Subcarriers, multiple of two
CP_Length     = 0;                              % Frequency CP length, multiple of two, can often be set to zero


%% FBMC Modulation Object
FBMC = Modulation.FBMC(...          
    NrSubcarriers,...                           % Number of subcarriers
    1,...                                       % Number of FBMC symbols
    15e3,...                                    % Subcarrier spacing (Hz)
    15e3*(NrSubcarriers*7),...                  % Sampling rate (Samples/s)
    0,...                                       % Intermediate frequency first subcarrier (Hz)
    false,...                                   % Transmit real valued signal
    'HermiteCut-OQAM',...                       % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman, 'HermiteCut0.8Tukey-OQAM', 'HermiteCut1.5-OQAM','HermiteCut1.46-OQAM'
    100, ...                                    % Pseudo Overlapping factor. Chosen extremly high to get better frequency resolution!
    0, ...                                      % Initial phase shift
    true ...                                    % Polyphase implementation
    );
FBMC_O08Tukey = Modulation.FBMC(...          
    NrSubcarriers,...                           % Number of subcarriers
    1,...                                       % Number of FBMC symbols
    15e3,...                                    % Subcarrier spacing (Hz)
    15e3*(NrSubcarriers*7),...                  % Sampling rate (Samples/s)
    0,...                                       % Intermediate frequency first subcarrier (Hz)
    false,...                                   % Transmit real valued signal
    'HermiteCut0.8Tukey-OQAM',...               % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman, 'HermiteCut0.8Tukey-OQAM', 'HermiteCut1.5-OQAM','HermiteCut1.46-OQAM'
    100, ...                                    % Pseudo Overlapping factor. Chosen extremly high to get better frequency resolution!
    0, ...                                      % Initial phase shift
    true ...                                    % Polyphase implementation
    );
FBMCconventionalHermite = Modulation.FBMC(...          
    NrSubcarriers,...                           % Number of subcarriers
    1,...                                       % Number of FBMC symbols
    15e3,...                                    % Subcarrier spacing (Hz)
    15e3*(NrSubcarriers*7),...                  % Sampling rate (Samples/s)
    0,...                                       % Intermediate frequency first subcarrier (Hz)
    false,...                                   % Transmit real valued signal
    'Hermite-OQAM',...                          % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman, 'HermiteCut0.8Tukey-OQAM', 'HermiteCut1.5-OQAM','HermiteCut1.46-OQAM'
    100, ...                                    % Pseudo Overlapping factor. Chosen extremly high to get better frequency resolution!
    0, ...                                      % Initial phase shift
    true ...                                    % Polyphase implementation
    );
FBMCconventionalPHYDYAS = Modulation.FBMC(...          
    NrSubcarriers,...                           % Number of subcarriers
    1,...                                       % Number of FBMC symbols
    15e3,...                                    % Subcarrier spacing (Hz)
    15e3*(NrSubcarriers*7),...                  % Sampling rate (Samples/s)
    0,...                                       % Intermediate frequency first subcarrier (Hz)
    false,...                                   % Transmit real valued signal
    'PHYDYAS-OQAM',...                          % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman, 'HermiteCut0.8Tukey-OQAM', 'HermiteCut1.5-OQAM','HermiteCut1.46-OQAM'
    4, ...                                      % Overlapping factor. Cannot be chosen high because PHYDYAS supports only supports only O<=8 in my implementation
    0, ...                                      % Initial phase shift
    true ...                                    % Polyphase implementation
    );
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round(FBMC.PHY.SamplingRate/FBMC.PHY.SubcarrierSpacing)+1/15e3/14*FBMC.PHY.SamplingRate)*FBMC.Nr.MCSymbols)/2)/FBMC.PHY.SamplingRate;
OFDM = Modulation.OFDM(...
    NrSubcarriers,...                           % Number of  active subcarriers
    1,...                                       % Number of OFDM Symbols
    FBMC.PHY.SubcarrierSpacing,...              % Subcarrier spacing (Hz)
    FBMC.PHY.SamplingRate,...                   % Sampling rate (Samples/s)
    0,...                                       % Intermediate frequency first subcarrier (Hz)
    false,...                                   % Transmit real valued signal
    1/15e3/14, ...                              % Cyclic prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    ZeroGuardTimeLength ...                     % Zero guard length (s)
    );

%% Frequency cylic prefix (for a large number of subcarriers not necessary)
% Note that, if CP_Length==0, then T_CP and R_CP are identity matrices
T_CP                                      = zeros(NrSubcarriers,NrSubcarriers-CP_Length);
T_CP(1:CP_Length/2,end-CP_Length/2+1:end) = eye(CP_Length/2);
T_CP(CP_Length/2+1:end-CP_Length/2,:)     = eye(NrSubcarriers-CP_Length);
T_CP(end-CP_Length/2+1:end,1:CP_Length/2) = eye(CP_Length/2);

R_CP                                      = zeros(NrSubcarriers,NrSubcarriers-CP_Length);
R_CP(CP_Length/2+1:end-CP_Length/2,:)     = eye(NrSubcarriers-CP_Length);

%% Consider only one FBMC symbol to calculate the precoding matrix
FBMC.SetNrMCSymbols( 1 );
% Get transmission matrix for one FBMC symbol
Gk = FBMC.GetTXMatrix;
Gk = Gk/sqrt(Gk(:,1)'*Gk(:,1));

% DFT matrix for the spreading process
W = fft( eye( NrSubcarriers - CP_Length ) ) / sqrt( NrSubcarriers-CP_Length );

% Diagonal elements of the FBMC transmission matrix after DFT spreading and despreading
a = abs( diag( W'*R_CP'*Gk'*Gk*T_CP*W ) );
a = a + randn(size(a))*10^-10;a(1)=a(1)+10^-9; % so that sorting is unique (numerical inaccuracies)

% Sort a
a_Tilde = sort( a , 'descend' );

% Get index representing the largest values of a
alpha       = a_Tilde( (NrSubcarriers-CP_Length)/2 );
Index_Tilde = ( a >= alpha );

% Pruned DFT matrix
W_Tilde = W(:,Index_Tilde) ;

% One-tap scaling of the data symbols
b_Tilde = sqrt( 2./( a(Index_Tilde) ) );    % Compared to the paper, an addional factor of "2" so that the transmit power is the same as in FBMC-OQAM

% Final coding matrix for one FBMC symbol
C0 = W_Tilde * diag(b_Tilde);

if size(C0,2)~=(NrSubcarriers-CP_Length)/2
   error('Size should be "(Number of Subcarriers-CP Length)/2"');
end

%% Calculate the precoding matrix for FBMC_O08Tukey. Same as above.
FBMC_O08Tukey.SetNrMCSymbols( 1 );
Gk_O08Tukey = FBMC_O08Tukey.GetTXMatrix;
Gk_O08Tukey = Gk_O08Tukey/sqrt(Gk_O08Tukey(:,1)'*Gk_O08Tukey(:,1));
W           = fft( eye( NrSubcarriers - CP_Length ) ) / sqrt( NrSubcarriers-CP_Length );
a           = abs( diag( W'*R_CP'*Gk_O08Tukey'*Gk_O08Tukey*T_CP*W ) );
a           = a + randn(size(a))*10^-10;a(1)=a(1)+10^-9; % so that sorting is unique (numerical inaccuracies)
a_Tilde     = sort( a , 'descend' );
alpha       = a_Tilde( (NrSubcarriers-CP_Length)/2 );
Index_Tilde = ( a >= alpha );
W_Tilde     = W(:,Index_Tilde) ;
b_Tilde     = sqrt( 2./( a(Index_Tilde) ) );
C0_O08Tukey = W_Tilde * diag(b_Tilde);


%% Calculate the Power Spectral Densities
% OFDM
G_OFDM      = OFDM.GetTXMatrix;
G_OFDM      = G_OFDM/sqrt(G_OFDM(:,1)'*G_OFDM(:,1));
PSD_OFDM_dB = 10*log10(sum(abs(fft(G_OFDM)).^2,2));

% FBMC-OQAM Hermite Pulse
GkconventionalHermite           = FBMCconventionalHermite.GetTXMatrix;
GkconventionalHermite           = GkconventionalHermite/sqrt(GkconventionalHermite(:,1)'*GkconventionalHermite(:,1));
PSD_FBMCconventionalHermite_dB  = 10*log10(sum(abs(fft(GkconventionalHermite)).^2,2));

% FBMC-OQAM PHYDYAS Puls
GkconventionalPHYDYAS           = FBMCconventionalPHYDYAS.GetTXMatrix;
AddedZeros                      = (size(GkconventionalHermite,1)-size(GkconventionalPHYDYAS,1))/2; % We add zeros to improve the resolution in the frequency domain. For other pulses we did this implicitly with an extremly high overlapping factor
GkconventionalPHYDYAS           = [zeros(AddedZeros,NrSubcarriers);GkconventionalPHYDYAS;zeros(AddedZeros,NrSubcarriers)];
GkconventionalPHYDYAS           = GkconventionalPHYDYAS/sqrt(GkconventionalPHYDYAS(:,1)'*GkconventionalPHYDYAS(:,1));
PSD_ConventionalFBMCPHYDYAS_dB  = 10*log10(sum(abs(fft(GkconventionalPHYDYAS)).^2,2));

% Pruned DFT spread FBMC, truncated Hermite pulse (O=1.56)
PSD_pDFTsFBMC_HermiteCut_dB = 10*log10(sum(abs(fft((Gk*T_CP*C0))).^2,2));

% Pruned DFT spread FBMC, truncated Hermite pulse + Tukey window (O=0.8)
PSD_pDFTsFBMC_O08Tukey_dB = 10*log10(sum(abs(fft((Gk_O08Tukey*T_CP*C0_O08Tukey))).^2,2));

% Normalize to PSD = 0dB!
PSD_ShiftdB = mean(PSD_OFDM_dB(5:FBMC.Implementation.FrequencySpacing*(FBMC.Nr.Subcarriers-5)));


%% Plot Power Spectral Density
figure(6);
f = (0:FBMC.Nr.SamplesTotal-1)*FBMC.PHY.SamplingRate/(FBMC.Nr.SamplesTotal);
f_Normalized = f/FBMC.PHY.SubcarrierSpacing-FBMC.Nr.Subcarriers+0.5;
plot(f_Normalized, PSD_OFDM_dB                      - PSD_ShiftdB,'red');hold on;
plot(f_Normalized, PSD_ConventionalFBMCPHYDYAS_dB   - PSD_ShiftdB,'Color',[0.2 0.7 0.2]);hold on;
plot(f_Normalized, PSD_FBMCconventionalHermite_dB   - PSD_ShiftdB,'magenta');hold on;
plot(f_Normalized, PSD_pDFTsFBMC_HermiteCut_dB      - PSD_ShiftdB,'blue');hold on;
plot(f_Normalized, PSD_pDFTsFBMC_O08Tukey_dB        - PSD_ShiftdB,'Color',[0.7 0.7 0.7]);hold on;
ylim([-60 1]);
xlim([-5 20]);
plot([0,0],[-100,1],'-','Color',[1 1 1]*0.5);

xlabel('Normalized Frequency f/F');
ylabel('Power Spectral Density [dB]');
legend('CP-OFDM, SC-FDMA','FBMC (PHYDYAS)','FBMC (Hermite)','Pruned DFT-s FBMC (Hermite, O=1.56)','Pruned DFT-s FBMC (Hermite, O=0.8)');



