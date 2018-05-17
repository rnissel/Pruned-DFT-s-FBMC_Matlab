% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.tc.tuwien.ac.at    

% This script calculates the ergodic Bit-Interleaved Coded Modulation(BICM)
% capacity for an AWGN channel. It requires the class 
% "../Modulation/SignalConstellation".

clear; close all;

addpath('..');

%% Parameters
PlotResult                      = true; 
SaveResult                      = true;

SNR_dB                          = -30:0.25:50;
MonteCarloRepetitionNoise       = 50000; % 1000
MonteCarloRepetitionChannel     = 1; %20000

Method                          = 'PAM';
M_ModulationOrder               = [2 4 8];
% Method = 'QAM';
% M_ModulationOrder = [4 16 64];


%% Start Calculation
for i_Modulation = 1:length(M_ModulationOrder)
    SignalConstellationTemp = Modulation.SignalConstellation(M_ModulationOrder(i_Modulation),Method);
    M_SignalConstellationSymbolMapping{i_Modulation}= SignalConstellationTemp.SymbolMapping;
    M_SignalConstallationBitMapping{i_Modulation} = SignalConstellationTemp.BitMapping;
end

if strcmp(Method,'PAM')
    MethodIsPAM = true;
else
    MethodIsPAM = false;
end

cluster = parcluster('local');
NrWorkers = cluster.NumWorkers; 
M_C_all = nan(length(SNR_dB),size(M_SignalConstallationBitMapping,2),MonteCarloRepetitionChannel);
% parfor i_rep = 1:MonteCarloRepetitionChannel
for i_rep = 1:MonteCarloRepetitionChannel   
    channel = 1; % AWGN channel!
    Ref_noise = 1/sqrt(2)*randn(MonteCarloRepetitionNoise,2)*[1;1j];

    C_all_temp = nan(length(SNR_dB),size(M_SignalConstallationBitMapping,2));
    for i_ModulationOrder = 1:size(M_SignalConstallationBitMapping,2)

        SignalConstellationSymbolMapping = M_SignalConstellationSymbolMapping{i_ModulationOrder};
        SignalConstallationBitMapping = M_SignalConstallationBitMapping{i_ModulationOrder};

        M = size(SignalConstellationSymbolMapping,1);
        all_symbol_indices = repmat((1:M)',[MonteCarloRepetitionNoise 1 ]);
        DataSymbols = SignalConstellationSymbolMapping(all_symbol_indices);

        noise = reshape(repmat(Ref_noise.',size(SignalConstellationSymbolMapping)),[],1);

        temp = channel*bsxfun(@minus,DataSymbols,SignalConstellationSymbolMapping.');
        for i_SNR = 1:length(SNR_dB)
            Pn = 10^(-SNR_dB(i_SNR)/10); 
            
            %Taking the real part in PAM detection reduces the noise power by 2: 
            if MethodIsPAM
                Pn = Pn*2;
            end
            
            pmf = exp(-abs(bsxfun(@plus,temp,noise*sqrt(Pn))).^2/Pn);
            I_b_m = nan(log2(M),2,size(SignalConstellationSymbolMapping,1)/2*MonteCarloRepetitionNoise);
            for i_b_index = 1:log2(M)
                for bit_value = 0:1  
                    nominator = sum(pmf(SignalConstallationBitMapping(all_symbol_indices,i_b_index)==bit_value,:),2);
                    denominator = sum(pmf(SignalConstallationBitMapping(all_symbol_indices,i_b_index)==bit_value,SignalConstallationBitMapping(:,i_b_index)==bit_value),2);
                    I_b_m(i_b_index,bit_value+1,:) = log2(nominator./denominator);
                end
            end
            C_all_temp(i_SNR,i_ModulationOrder) = mean(reshape(mean(reshape(log2(M) - sum(mean(I_b_m,2),1),size(SignalConstellationSymbolMapping,1)/2,[]),1),MonteCarloRepetitionNoise,[]),1);
            
            disp(['Modulation ' int2str(i_ModulationOrder) ' of ' int2str(size(M_SignalConstallationBitMapping,2))  ', SNR ' int2str(i_SNR) ' of '  int2str(length(SNR_dB))]);
        end
    end

    M_C_all(:,:,i_rep) = C_all_temp;

end


M_C_max = squeeze(max(M_C_all,[],2));
C_max = mean(M_C_max,2);
C = mean(M_C_all,3);

if PlotResult
    figure();
    plot(SNR_dB,C_max,'black');
    hold on;
    plot(SNR_dB,C);
    LegendStringArray{1} = 'Max. over modulation order';
    for i_ModOrd=1:length(M_ModulationOrder)
        LegendStringArray{i_ModOrd+1} = [int2str(M_ModulationOrder(i_ModOrd)) ' - ' Method];
    end
    legend(LegendStringArray,'Location','NorthWest');
    ylabel('BICM Capacity');
    xlabel('SNR (dB)');
end

if SaveResult
    StringTemp ='';
    for i_ModOrd=1:length(M_ModulationOrder)
        StringTemp = [StringTemp int2str(M_ModulationOrder(i_ModOrd)) '_'];
    end
    save(['BICM_Capacity_AWGN_' StringTemp Method '.mat'], ...
         'C_max',...
         'C',...
         'SNR_dB',...
         'MonteCarloRepetitionNoise',...
         'MonteCarloRepetitionChannel',...
         'Method',...
         'M_ModulationOrder');
end



