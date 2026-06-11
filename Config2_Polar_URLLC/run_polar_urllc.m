function run_polar_urllc()
% run_polar_urllc - Моделирование полярного кода для URLLC с адаптивным остановом
%   K = 512 бит, N = 1024, скорость 1/2, SCL-8, канал AWGN + QPSK

clear; close all; clc;

%% Параметры моделирования
K = 512; N = 1024; R = K/N; L = 8;
EbNo_dB = 0:0.25:4;
numFramesMax = 100000;
targetFER = 200;

M = 4; bps = log2(M);

resultsDir = './Results';
figuresDir = './Figures';
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end

results = zeros(length(EbNo_dB),4);
fprintf('\n========== Запуск симуляции ==========\n');
fprintf('K=%d, N=%d, R=%.2f, L=%d (UCI, nmax=10)\n',K,N,R,L);

parfor idx = 1:length(EbNo_dB)
    EbNo = EbNo_dB(idx);
    fprintf('\n--- Eb/N0 = %.2f dB ---\n',EbNo);
    
    totalBitErrors = 0; totalBits = 0;
    totalFrameErrors = 0; totalFrames = 0;
    
    while (totalFrameErrors < targetFER) && (totalFrames < numFramesMax)
        data = randi([0 1], K, 1);
        encData = nrPolarEncode(data, N, 10, false);
        modData = pskmod(encData, M, pi/4, 'gray', 'InputType','bit');
        SNR = EbNo + 10*log10(R) + 10*log10(bps);
        rxSig = awgn(modData, SNR, 'measured');
        rxLLR = pskdemod(rxSig, M, pi/4, 'gray', 'OutputType','llr');
        rxBits = nrPolarDecode(rxLLR, K, N, L, 10, false, 24);
        
        bitErrors = sum(rxBits ~= data);
        frameError = bitErrors > 0;
        
        totalBitErrors = totalBitErrors + bitErrors;
        totalBits = totalBits + K;
        totalFrameErrors = totalFrameErrors + frameError;
        totalFrames = totalFrames + 1;
        
        if totalFrameErrors >= targetFER, break; end
        if totalFrames >= 15000 && totalFrameErrors == 0, break; end
        if totalFrames >= 30000 && totalFrameErrors < 50, break; end
        if totalFrames >= 60000 && totalFrameErrors < 100, break; end
    end
    
    BER = totalBitErrors / totalBits;
    FER = totalFrameErrors / totalFrames;
    results(idx,:) = [EbNo, BER, FER, totalFrames];
    fprintf('BER=%.2e, FER=%.2e, frames=%d\n', BER, FER, totalFrames);
end

csvwrite(fullfile(resultsDir, 'polar_urllc_results.csv'), results);

figure('Position',[100,100,800,600]);
semilogy(results(:,1), results(:,2), 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
semilogy(results(:,1), results(:,3), 'r--s', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER / FER');
title('Полярный код, URLLC (K=512, R=1/2, L=8, UCI, nmax=10)');
legend('BER','FER','Location','southwest');
saveas(gcf, fullfile(figuresDir, 'polar_urllc_ber_fer.png'));

fprintf('\nРезультаты сохранены в %s\n', fullfile(resultsDir, 'polar_urllc_results.csv'));
fprintf('График сохранён в %s\n', fullfile(figuresDir, 'polar_urllc_ber_fer.png'));
end