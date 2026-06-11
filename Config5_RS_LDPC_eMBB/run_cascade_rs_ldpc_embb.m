function run_cascade_rs_ldpc_embb()
% Каскад RS(255,239) + LDPC (eMBB) для MATLAB R2025b

clear; close all; clc;

resultsDir = 'Results';
figuresDir = 'Figures';
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end


n_rs = 255; k_rs = 239;
m_rs = 8;                     
prim_poly = 285;              
genpoly = rsgenpoly(n_rs, k_rs, prim_poly);  


blockSize = 81;
P = [16 17 22 24 9 3 14 -1 4 2 7 -1 26 -1 2 -1 21 -1 1 0 -1 -1 -1 -1;
     25 12 12 3 3 26 6 21 -1 15 22 -1 15 -1 4 -1 -1 16 -1 0 0 -1 -1 -1;
     25 18 26 16 22 23 9 -1 0 -1 4 -1 4 -1 8 23 11 -1 -1 -1 0 0 -1 -1;
     9 7 0 1 17 -1 -1 7 3 -1 3 23 -1 16 -1 -1 21 -1 0 -1 -1 0 0 -1;
     24 5 26 7 1 -1 -1 15 24 15 -1 8 -1 13 -1 13 -1 11 -1 -1 -1 -1 0 0;
     2 2 19 14 24 1 15 19 -1 21 -1 2 -1 24 -1 3 -1 2 1 -1 -1 -1 -1 0];
H_num = ldpcQuasiCyclicMatrix(blockSize, P);
H = logical(sparse(H_num));
[Nrows, Ncols] = size(H);
K_ldpc = Ncols - Nrows;   
N_ldpc = Ncols;           
cfgEnc = ldpcEncoderConfig(H);
cfgDec = ldpcDecoderConfig(H);


R_total = (k_rs / n_rs) * (K_ldpc / N_ldpc);   
fprintf('Каскад RS(%d,%d) + LDPC(%d,%d)\n', n_rs, k_rs, K_ldpc, N_ldpc);
fprintf('Общая скорость = %.3f\n', R_total);


EbNo_dB = 0:0.5:4;
maxFrames = 20000;
targetFER = 200;
maxIter = 50;
M = 4; bps = log2(M);     

results = zeros(length(EbNo_dB), 4);


for idx = 1:length(EbNo_dB)
    EbNo = EbNo_dB(idx);
    totalBitErrors = 0; totalBits = 0;
    totalFrameErrors = 0; totalFrames = 0;
    fprintf('\n--- Eb/N0 = %.2f dB ---\n', EbNo);
    
    SNR = EbNo + 10*log10(R_total) + 10*log10(bps);
    snr_lin = 10^(SNR/10);
    noiseVar = 1/(2*snr_lin);
    
    while (totalFrameErrors < targetFER) && (totalFrames < maxFrames)
        data_sym = randi([0, 255], 1, k_rs);
        data_gf = gf(data_sym, m_rs, prim_poly);
        code_gf = rsenc(data_gf, n_rs, k_rs, genpoly);
        code_sym = double(code_gf.x);            
        
        code_bits = reshape(de2bi(code_sym, m_rs, 'left-msb')', [], 1);
        
        ldpc_data = code_bits(1:K_ldpc);        
        
        enc_ldpc = ldpcEncode(ldpc_data, cfgEnc);   
        
        tx = pskmod(enc_ldpc, M, pi/4, 'gray', 'InputType','bit');
        rx = awgn(tx, SNR, 'measured');
        rxLLR = pskdemod(rx, M, pi/4, 'gray', 'OutputType','llr', 'NoiseVariance', noiseVar);
        
        dec_ldpc = ldpcDecode(rxLLR, cfgDec, maxIter);   
        
       
        dec_bits_full = [dec_ldpc; zeros(n_rs*m_rs - K_ldpc, 1)];
        
        dec_sym_col = bi2de(reshape(dec_bits_full, m_rs, [])', 'left-msb');
        dec_sym = dec_sym_col';      
        
        rx_gf = gf(dec_sym, m_rs, prim_poly);
        dec_gf = rsdec(rx_gf, n_rs, k_rs, genpoly);
        if isempty(dec_gf)
            dec_sym_out = zeros(1, k_rs);
        else
            dec_sym_out = double(dec_gf.x);
        end
        dec_bits = reshape(de2bi(dec_sym_out, m_rs, 'left-msb')', [], 1);
        
        orig_bits = reshape(de2bi(data_sym, m_rs, 'left-msb')', [], 1);
        errs = sum(dec_bits ~= orig_bits);
        
        totalBitErrors = totalBitErrors + errs;
        totalBits = totalBits + length(orig_bits);
        totalFrameErrors = totalFrameErrors + (errs > 0);
        totalFrames = totalFrames + 1;
        
        if totalFrameErrors >= targetFER, break; end
        if totalFrames >= 10000 && totalFrameErrors == 0, break; end
    end
    
    BER = totalBitErrors / totalBits;
    FER = totalFrameErrors / totalFrames;
    results(idx, :) = [EbNo, BER, FER, totalFrames];
    fprintf('BER = %.2e, FER = %.2e, frames = %d\n', BER, FER, totalFrames);
end


csvwrite(fullfile(resultsDir, 'cascade_rs_ldpc.csv'), results);

figure;
semilogy(results(:,1), results(:,2), 'b-o', 'LineWidth', 1.5); hold on;
semilogy(results(:,1), results(:,3), 'r--s', 'LineWidth', 1.5);
grid on;
xlabel('Eb/N0 (дБ)');
ylabel('BER / FER');
legend('BER', 'FER', 'Location', 'southwest');
title('Каскад RS(255,239) + LDPC (eMBB)');
saveas(gcf, fullfile(figuresDir, 'cascade_rs_ldpc.png'));

fprintf('\nРезультаты сохранены в %s и %s\n', resultsDir, figuresDir);
end