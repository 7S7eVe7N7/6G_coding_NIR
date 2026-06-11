function run_polar_mmtc()
% mMTC: K=256, N=512, скорость 0.5, SCL-8
clear; close all; clc;

K = 256; N = 512; R = K/N; L = 8;
EbNo_dB = 0:0.25:4;
numFramesMax = 100000;
targetFER = 200;
M = 4; bps = log2(M);

resultsDir = './Results';
figuresDir = './Figures';
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end

results = zeros(length(EbNo_dB),4);
fprintf('\n========== mMTC (Polar) ==========\n');

for idx = 1:length(EbNo_dB)
    EbNo = EbNo_dB(idx);
    totalBitErrors = 0; totalBits = 0;
    totalFrameErrors = 0; totalFrames = 0;
    fprintf('\n--- Eb/N0 = %.2f dB ---\n', EbNo);
    
    while (totalFrameErrors < targetFER) && (totalFrames < numFramesMax)
        data = randi([0 1], K, 1);
        encData = nrPolarEncode(data, N, 10, false);
        modData = pskmod(encData, M, pi/4, 'gray', 'InputType','bit');
        SNR = EbNo + 10*log10(R) + 10*log10(bps);
        rxSig = awgn(modData, SNR, 'measured');
        rxLLR = pskdemod(rxSig, M, pi/4, 'gray', 'OutputType','llr');
        rxBits = nrPolarDecode(rxLLR, K, N, L, 10, false, 24);
        
        errs = sum(rxBits ~= data);
        totalBitErrors = totalBitErrors + errs;
        totalBits = totalBits + K;
        totalFrameErrors = totalFrameErrors + (errs > 0);
        totalFrames = totalFrames + 1;
        
        if totalFrameErrors >= targetFER, break; end
        if totalFrames >= 20000 && totalFrameErrors == 0, break; end
        if totalFrames >= 50000 && totalFrameErrors < 100, break; end
    end
    
    BER = totalBitErrors / totalBits;
    FER = totalFrameErrors / totalFrames;
    results(idx,:) = [EbNo, BER, FER, totalFrames];
    fprintf('BER=%.2e, FER=%.2e, frames=%d\n', BER, FER, totalFrames);
end

csvwrite(fullfile(resultsDir, 'polar_mmtc.csv'), results);
figure;
semilogy(results(:,1), results(:,2), 'b-o'); hold on;
semilogy(results(:,1), results(:,3), 'r--s'); grid;
xlabel('Eb/N0 (dB)'); ylabel('BER/FER');
legend('BER','FER');
title('mMTC: Polar (K=256,N=512,SCL-8)');
saveas(gcf, fullfile(figuresDir, 'polar_mmtc.png'));
fprintf('\nГотово.\n');
end