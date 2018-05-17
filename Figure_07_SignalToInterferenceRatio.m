% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================   
% This script calculates the signal-to-interference ratio for different
% frequency CP lengths.

% Reproduces Figure 7 of "Pruned DFT Spread FBMC: Low PAPR, Low Latency, 
% High Spectral Efficiency", R. Nissel and M. Rupp, IEEE Transactions on 
% Communications 

clear; close all; 

%% Parameters
M_CP_Length     = [0 2];                                    % Frequency cyclic prefix. Must be a multiple of two, e.g., 0, 2, 4...
M_NrSubcarriers = [4:4:100 110:40:200 200:100:300];         % Number of subcarriers, corresponds to the spreading length

% #########################################################################
% % In the paper:
% M_NrSubcarriers = [4:2:100 110:10:200 200:20:300 350:50:500 600:100:1200];  
% #########################################################################


%% For loops. The SIR calculation is the same as in "Figure_04_05_ab_TransmitPower.m".
for i_CP = 1:length(M_CP_Length)
    CP_Length = M_CP_Length(i_CP);      
    disp(['Cyclic Prefix: ' int2str(CP_Length)]);
    for i_NrSubcarriers =1:length(M_NrSubcarriers)

        NrSubcarriers = M_NrSubcarriers(i_NrSubcarriers);               % Number of Subcarriers, multiple of two
        NrFBMCsymbols = 6;                                              % Number of FBMC symbols in time, larger than 6 because we ignore the first two and last two FBMC symbols in the SIR calculation (to get rid of edge effects)

        if M_CP_Length(i_CP)==0
            PrototypeFilter = 'HermiteCut-OQAM';
        else
            PrototypeFilter = 'HermiteCut1.46-OQAM';
        end   
        %% FBMC Modulation Object        
        FBMC = Modulation.FBMC(...          
            NrSubcarriers,...               % Number of subcarriers
            NrFBMCsymbols,...               % Number of FBMC symbols
            15e3,...                        % Subcarrier spacing (Hz)
            15e3*(4400),...                 % Sampling rate (Samples/s)
            0,...                           % Intermediate frequency first subcarrier (Hz)
            false,...                       % Transmit real valued signal
            PrototypeFilter,...             % Prototype filter Hermite, PHYDYAS, InversePHYDYAS, HermiteCut, PHYDYASCut, Hann, Blackman, 'HermiteCut0.8Tukey-OQAM', 'HermiteCut1.5-OQAM','HermiteCut1.46-OQAM'
            2, ...                          % Overlapping factor
            0, ...                          % Initial phase shift
            true ...                        % Polyphase implementation
            );

        %% Frequency cylic prefix matrix (for a large number of subcarriers not necessary)
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

        % Reduced DFT matrix
        W_Tilde     = W(:,Index_Tilde) ;

        % One-tap scaling of the data symbols
        b_Tilde = sqrt( 2./( a(Index_Tilde) ) );    % an addional factor of "2" so that the transmit power is the same as in FBMC-OQAM!

        % Final coding matrix for one FBMC symbol
        C0 = W_Tilde * diag(b_Tilde);

        if size(C0,2)~=(NrSubcarriers-CP_Length)/2
           error('Size should be "(Number of Subcarriers-CP Length)/2"');
        end

        % Calculate the Signal-to-Interference Ratio for one FBMC symbol
        SignalInterferenceMatrixk   = C0'*R_CP'*Gk'*Gk*T_CP*C0;
        PSk                         = abs( diag( SignalInterferenceMatrixk ) ).^2;
        PIk                         = sum( abs( SignalInterferenceMatrixk ).^2 ,2 ) - PSk;
        SIRk_dB                     = 10*log10( sum(PSk)./sum(PIk) );


        %% Consider more FBMC symbol to also account for interblock interference!
        FBMC.SetNrMCSymbols( NrFBMCsymbols );

        % Get transmission matrix for the whole system
        G = FBMC.GetTXMatrix;
        G = G/sqrt(G(:,1)'*G(:,1));

        % Get Coding Matrix of the whole system
        C = kron( sparse(eye(FBMC.Nr.MCSymbols)) , C0 );

        % Get frequency cyclic prefix matrix for the whole system
        T_CP_all = kron( sparse(eye(FBMC.Nr.MCSymbols)) , T_CP );
        R_CP_all = kron( sparse(eye(FBMC.Nr.MCSymbols)) , R_CP );

        % Get signal and interference power for the whole system
        SignalInterferenceMatrix = C'*R_CP_all'*G'*G*T_CP_all*C;
        PS = reshape(abs( diag(SignalInterferenceMatrix) ).^2,[],NrFBMCsymbols);
        PI = reshape(sum( abs(SignalInterferenceMatrix).^2 , 2 ),[],NrFBMCsymbols) - PS;
        
        % Calculate SIR only at the middle position 
        SIR_dB      = 10*log10( sum( PS(:,end/2) ) ./ sum( PI(:,end/2) ) );
        SIR_dB_min  = 10*log10( min( PS(:,end/2) ./ PI(:,end/2) ) );
        SIR_dB_max  = 10*log10( max( PS(:,end/2) ./ PI(:,end/2) ) );
        
        
        % Save the instantaneous SIR values in a matrix
        M_SIRk_dB(i_NrSubcarriers,i_CP) = SIRk_dB;
        
        M_SIR_dB(i_NrSubcarriers,i_CP)  = SIR_dB;        
        M_SIR_dB_min(i_NrSubcarriers,i_CP)  = SIR_dB_min;
        M_SIR_dB_max(i_NrSubcarriers,i_CP)  = SIR_dB_max;   
        
        disp(['L_CP = ' int2str(CP_Length) ': ' int2str(i_NrSubcarriers/length(M_NrSubcarriers)*100) '%']);
    end
end

%% Plot results
figure(7);
plot(M_NrSubcarriers,M_SIR_dB); hold on;
plot(M_NrSubcarriers,M_SIR_dB_min,':'); hold on;
plot(M_NrSubcarriers,M_SIR_dB_max,':'); hold on;
xlabel('Number of Subcarriers, L');
ylabel('Signal-to-Interference Ratio [dB]');
ylim([0 60]);
xlim([0 max(M_NrSubcarriers)]);
legend(['L_{CP} = ' int2str(M_CP_Length(1))],['L_{CP} = ' int2str(M_CP_Length(2))],'Location','SouthEast');



