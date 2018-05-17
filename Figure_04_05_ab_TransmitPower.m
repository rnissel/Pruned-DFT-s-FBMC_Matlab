% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================         
% This script shows how to construct the precoding matrix C_f. Furthermore,
% the expected transmit power over time is plotted.

% Reproduces Figure 4 and 5 of "Pruned DFT Spread FBMC: Low PAPR, Low 
% Latency, High Spectral Efficiency", R. Nissel and M. Rupp, IEEE 
% Transactions on Communications

clear; close all;

%% Parameters
NrSubcarriers   = 128;                  % Number of Subcarriers, multiple of two
CP_Length       = 0;                    % Frequency CP length, multiple of two, can often be set to zero

% Only relevant for the SIR calculation
NrFBMCsymbols   = 6;                    % Number of FBMC symbols in time, 6 because we ignore the first two and last two FBMC symbols in the SIR calculation (to get rid of edge effects)


%% FBMC Modulation Object
FBMC = Modulation.FBMC(...          
    NrSubcarriers,...                   % Number of subcarriers
    NrFBMCsymbols,...                   % Number of FBMC symbols
    15e3,...                            % Subcarrier spacing (Hz)
    15e3*(NrSubcarriers*30),...         % Sampling rate (Samples/s),  30-times oversampling
    0,...                               % Intermediate frequency first subcarrier (Hz)
    false,...                           % Transmit real valued signal
    'HermiteCut-OQAM',...               % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman, 'HermiteCut0.8Tukey-OQAM', 'HermiteCut1.5-OQAM','HermiteCut1.46-OQAM'
    2, ...                              % Overlapping factor
    0, ...                              % Initial phase shift
    true ...                            % Polyphase implementation
    );


%% Frequency cylic prefix (for a large number of subcarriers not necessary!)
% Note that, if CP_Length==0, then T_CP and R_CP are identity matrices
T_CP                                        = zeros(NrSubcarriers,NrSubcarriers-CP_Length);
T_CP(1:CP_Length/2,end-CP_Length/2+1:end)   = eye(CP_Length/2);
T_CP(CP_Length/2+1:end-CP_Length/2,:)       = eye(NrSubcarriers-CP_Length);
T_CP(end-CP_Length/2+1:end,1:CP_Length/2)   = eye(CP_Length/2);

R_CP                                        = zeros(NrSubcarriers,NrSubcarriers-CP_Length);
R_CP(CP_Length/2+1:end-CP_Length/2,:)       = eye(NrSubcarriers-CP_Length);


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
b_Tilde = sqrt( 2./( a(Index_Tilde) ) );       % Compared to the paper, an additional factor of "2" so that the transmit power is the same as in FBMC-OQAM

% Final coding matrix for one FBMC symbol
Cf = W_Tilde * diag(b_Tilde);

if size(Cf,2)~=(NrSubcarriers-CP_Length)/2
   error('Size should be "(Number of Subcarriers-CP Length)/2"');
end

% Calculate the Signal-to-Interference Ratio for one FBMC symbol
SignalInterferenceMatrixk   = Cf'*R_CP'*Gk'*Gk*T_CP*Cf;
PSk                         = abs( diag( SignalInterferenceMatrixk ) ).^2;
PIk                         = sum( abs( SignalInterferenceMatrixk ).^2 ,2 ) - PSk;
SIRk_dB                     = 10*log10( sum(PSk)./sum(PIk) );


%% Consider more FBMC symbol to also account for interblock interference!
FBMC.SetNrMCSymbols( NrFBMCsymbols );

% Get transmission matrix for the whole system
G = FBMC.GetTXMatrix;
G = G/sqrt(G(:,1)'*G(:,1));

% Get Coding Matrix for the whole system
C = kron( sparse(eye(FBMC.Nr.MCSymbols)) , Cf );

% Get frequency cyclic prefix matrix for the whole system
T_CP_all = kron( sparse(eye(FBMC.Nr.MCSymbols)) , T_CP );
R_CP_all = kron( sparse(eye(FBMC.Nr.MCSymbols)) , R_CP );

% Get signal and interference power for the whole system
SignalInterferenceMatrix = C'*R_CP_all'*G'*G*T_CP_all*C;
PS = reshape(abs( diag(SignalInterferenceMatrix) ).^2,[],NrFBMCsymbols);
PI = reshape(sum( abs(SignalInterferenceMatrix).^2 , 2 ),[],NrFBMCsymbols) - PS;

% We consider only the middle time-position to get rid of edge effects
SIR_dB = 10*log10(sum(PS(:,end/2))./sum(PI(:,end/2)));


%% Display the SIR
disp(['=================================================================']);
disp(['Signal-to-Interference Ratio (one FBMC Symbol): ' num2str(SIRk_dB,4) 'dB']);
disp(['Signal-to-Interference Ratio (including ISI):   ' num2str(SIR_dB,4) 'dB']);
disp(['=================================================================']);


%% Plot the results
% Concept of finding one-tap scaling factors
figure(4); 
b_Tilde_plot = nan(size(a));
b_Tilde_plot(Index_Tilde) = b_Tilde;
plot( (1:length(a))/length(a) , a ,'red'); hold on;
plot( (1:length(a))/length(a) , b_Tilde_plot/sqrt(2) ,'blue'); % We must include a factor of two due to power normalization to FBMC-OQAM
xlabel('i-th Position / L');
legend({'a','b'});

% Calculate the transmitted power in time
% Note that: diag(A*A') = sum(A.*conj(A),2)!!!
TransmitPower_FBMC_CutHermite   = sum( Gk.*conj(Gk) , 2 );
TransmitPower_FBMC_DFTspread    = sum( (Gk*T_CP*Cf).*conj((Gk*T_CP*Cf)) , 2 );

NormalizedTime = ((0:FBMC.Implementation.FFTSize*FBMC.PrototypeFilter.OverlappingFactor-1)-(FBMC.Implementation.FFTSize*FBMC.PrototypeFilter.OverlappingFactor)/2)/FBMC.Implementation.FFTSize;
Scaling = FBMC.PHY.SamplingRate/(FBMC.PHY.SubcarrierSpacing*FBMC.Nr.Subcarriers)/2;

figure(5);
plot( NormalizedTime , TransmitPower_FBMC_CutHermite*Scaling ,'red' );hold on;
plot( NormalizedTime , TransmitPower_FBMC_DFTspread*Scaling  ,'blue');
ylim([0 1.3]);
CutPointHermite1        = NormalizedTime((sum(FBMC.PrototypeFilter.TimeDomain==0)-1)/2);
CutPointHermite2        = NormalizedTime(end-(sum(FBMC.PrototypeFilter.TimeDomain==0)-1)/2);
plot([CutPointHermite1,CutPointHermite1],[0,1.2],'--','Color',[0.8 0.8 0.8]);
plot([CutPointHermite2,CutPointHermite2],[0,1.2],'--','Color',[0.8 0.8 0.8]);
CutPointRectFilter1     = NormalizedTime(size(Gk,1)/2-floor((FBMC.Implementation.FFTSize)/2*0.8)+1);
CutPointRectFilter2     = NormalizedTime(size(Gk,1)/2+floor((FBMC.Implementation.FFTSize)/2*0.8)+1);
plot([CutPointRectFilter1,CutPointRectFilter1],[0,1.1],'--','Color',[0.8 0.8 0.8]);
plot([CutPointRectFilter2,CutPointRectFilter2],[0,1.1],'--','Color',[0.8 0.8 0.8]);
xlabel('Normalized Time, t/F');
ylabel('Transmit Power, E\{|s(t)|^2\}');
legend({'FBMC','Pruned DFT-s FBMC'});



