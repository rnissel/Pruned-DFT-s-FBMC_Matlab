% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================   
% This script compares the spectral efficiency of OFDM and SC-FDMA (SC-FDMA
% is equivalent to pruned DFT spread FBMC), under the assumption of
% infinitely many subcarriers (and a suboptimal receiver which ignores
% interference and noise correlation in SC-FDMA).

clear; close all;

% Parameters
M_SNR = [-10:0.1:30];

% Load the AWGN BICM capacity for {4,16,64}-QAM
Load_BICM_QAM = load('Theory\BICM_Capacity_AWGN_4_16_64_QAM');
BICM_QAM      = Load_BICM_QAM.C_max;
SNR_Ref_QAM   = Load_BICM_QAM.SNR_dB;


M_Pn = 10.^(-M_SNR/10);
for i_Pn = 1:length(M_Pn)
    Pn = M_Pn(i_Pn);

    % See (37). Note that we calculated the averagin term in Mathematica with "Integrate[ 1/(1 + Pn/(h^2))*2*h*Exp[-(h^2)], {h, 0, Infinity}, Assumptions -> Pn \[Element] Reals && Pn > 0]"
    SINR_SCFDMA(i_Pn) = 1/Pn*(exp(-Pn)/expint(Pn)-Pn);
    
    C_SCFDMA(i_Pn) = log2(1+SINR_SCFDMA(i_Pn));
    C_OFDM(i_Pn)   = exp(Pn)*expint(Pn)/log(2);
  
    BICM_SCFDMA(i_Pn) = interp1( SNR_Ref_QAM , BICM_QAM , 10*log10(SINR_SCFDMA(i_Pn)),'spline');   
    
    fun2            = @(h) interp1( SNR_Ref_QAM , BICM_QAM , 10*log10(abs(h).^2/Pn),'spline',6).*2.*h.*exp(-(h.^2));
    BICM_OFDM(i_Pn) = integral(fun2,0,inf);
      
end

%% Plot results
figure(100);
plot(M_SNR,C_OFDM',': red');hold on;
plot(M_SNR,C_SCFDMA',': blue');hold on;
plot(M_SNR,BICM_OFDM','red');hold on;
plot(M_SNR,BICM_SCFDMA','blue');hold on;
xlabel('Signal-to-Noise Ratio');
ylabel('Spectral Efficiency');
legend({'Capacity OFDM','Capacity SC-FDMA','BICM OFDM','BICM SC-FDMA'});
title('For infinitely many subcarriers, L=>inf.');






