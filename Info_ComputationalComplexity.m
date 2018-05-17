% =========================================================================   
% (c) 2018 Ronald Nissel, https://www.linkedin.com/in/ronaldnissel/
% =========================================================================   
% This script compares the computational complexity of pruned DFT spread
% FBMC, SC-FDMA and "Low PAPR FBMC".


clear; close all;

% Parameters
N_FFT =  128:16:2048;                   % FFT size for OFDM and FBMC
NoverL = [1 2 4 8 16];                  % Factor: "N_FFT/L"                      
O = 0.8;                                % Overlapping factor

% Start calculation
for i_NoverL = 1:length(NoverL)
    
    L = N_FFT/NoverL(i_NoverL);

    CFBMC = 2*(L/2+L.*log2(L/2)+N_FFT.*log2(N_FFT)+O*N_FFT);
    COFDM = L.*log2(L)+N_FFT.*log2(N_FFT);
    
    ClowPAPRFBMC =  L.*log2(L) + 2*N_FFT.*log2(N_FFT)+ 4*4*N_FFT;  % O=4; See the paper "Low PAPR FBMC"

    ComplexityIncrease(:,i_NoverL) = CFBMC./COFDM;    
    ComplexityIncreaseLowPAPRFBMC(:,i_NoverL) = ClowPAPRFBMC./CFBMC;
   
end

% Plot the computational complexity relative to the reference
figure(1);
plot(N_FFT,ComplexityIncrease);
xlabel('FFT Size for FBMC and OFDM'); 
ylabel('Relative Complexity');
title('Pruned DFT-s FBMC vs SC-FDMA (Reference)');
legend([repmat('N_{FFT}/L = ',length(NoverL),1) int2str(NoverL')]);
xlim([min(N_FFT) max(N_FFT)]);

figure(2);
plot(N_FFT,ComplexityIncreaseLowPAPRFBMC);
xlabel('FFT Size for FBMC and OFDM'); 
ylabel('Relative Complexity');
title('Low PAPR FBMC vs Pruned DFT-s FBMC (Reference)')
legend([repmat('N_{FFT}/L = ',length(NoverL),1) int2str(NoverL')]);
xlim([min(N_FFT) max(N_FFT)]);





