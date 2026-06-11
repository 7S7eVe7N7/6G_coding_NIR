function run_cascade_bch_polar_urllc()
% Каскад BCH(511,493) + полярный (512,1024) для URLLC
clear; close all; clc;

resultsDir = 'Results';
figuresDir = 'Figures';
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end

n_bch = 511; k_bch = 493;
bchEnc = comm.BCHEncoder('BitInput', true, 'CodewordLength', n_bch, 'MessageLength', k_bch);
bchDec = comm.BCHDecoder('BitInput', true, 'CodewordLength', n_bch, 'MessageLength', k_bch);

K_polar = 512; N_polar = 1024; L = 8;
nmax = 10; iil = false; CRClen = 24;
R_polar = K_polar / N_polar;
R_total = k_bch / N_polar;  

EbNo_dB = 0:0.5:4;
maxFrames = 20000;
targetFER = 200;
M = 4; bps = log2(M);

results = zeros(length(EbNo_dB),4);

fprintf('\n========== Каскад BCH+Polar (URLLC) ==========\n');
fprintf('BCH(%d,%d) + Polar(%d,%d), L=%d\n', n_bch, k_bch, K_polar, N_polar, L);
fprintf('Общая скорость = %.3f\n', R_total);

for idx = 1:length(EbNo_dB)
    EbNo = EbNo_dB(idx);
    totalBitErrors = 0; totalBits = 0;
    totalFrameErrors = 0; totalFrames = 0;
    fprintf('\n--- Eb/N0 = %.2f dB ---\n', EbNo);
    
    while (totalFrameErrors < targetFER) && (totalFrames < maxFrames)
        data_bch_in = randi([0 1], k_bch, 1);
        enc_bch = bchEnc(data_bch_in);   % 511 бит
        
        enc_polar_in = [enc_bch; 0];
        
        enc_polar = nrPolarEncode(enc_polar_in, N_polar, nmax, iil);
        
        tx = pskmod(enc_polar, M, pi/4, 'gray', 'InputType','bit');
        
        SNR = EbNo + 10*log10(R_total) + 10*log10(bps);
        rx = awgn(tx, SNR, 'measured');
        
        rxLLR = pskdemod(rx, M, pi/4, 'gray', 'OutputType','llr');
        
        dec_polar = nrPolarDecode(rxLLR, K_polar, N_polar, L, nmax, iil, CRClen);
        
        dec_bch_in = dec_polar(1:n_bch);
        
        dec_bch = bchDec(dec_bch_in);
        
        errs = sum(dec_bch ~= data_bch_in);
        totalBitErrors = totalBitErrors + errs;
        totalBits = totalBits + k_bch;
        totalFrameErrors = totalFrameErrors + (errs > 0);
        totalFrames = totalFrames + 1;
        
        if totalFrameErrors >= targetFER, break; end
        if totalFrames >= 10000 && totalFrameErrors == 0, break; end
    end
    
    BER = totalBitErrors / totalBits;
    FER = totalFrameErrors / totalFrames;
    results(idx,:) = [EbNo, BER, FER, totalFrames];
    fprintf('BER = %.2e, FER = %.2e, frames = %d\n', BER, FER, totalFrames);
end

csvwrite(fullfile(resultsDir, 'cascade_bch_polar.csv'), results);
figure;
semilogy(results(:,1), results(:,2), 'b-o'); hold on;
semilogy(results(:,1), results(:,3), 'r--s'); grid;
xlabel('Eb/N0 (dB)'); ylabel('BER / FER');
legend('BER', 'FER');
title('Каскад BCH(511,493)+Полярный (512,1024)');
saveas(gcf, fullfile(figuresDir, 'cascade_bch_polar.png'));
fprintf('\nРезультаты сохранены в %s и %s\n', resultsDir, figuresDir);
end