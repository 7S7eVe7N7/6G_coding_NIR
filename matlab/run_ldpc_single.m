function run_ldpc_single()
% LDPC с параметрами K_target=1912, N_target=2552
% Используется реальный код (K_real=1458, N_real=1944) с дополнением нулями.

clear; close all; clc;

% Создание рабочей проверочной матрицы H (размер 486x1944)
blockSize = 81;
P = [16 17 22 24 9 3 14 -1 4 2 7 -1 26 -1 2 -1 21 -1 1 0 -1 -1 -1 -1;
     25 12 12 3 3 26 6 21 -1 15 22 -1 15 -1 4 -1 -1 16 -1 0 0 -1 -1 -1;
     25 18 26 16 22 23 9 -1 0 -1 4 -1 4 -1 8 23 11 -1 -1 -1 0 0 -1 -1;
     9 7 0 1 17 -1 -1 7 3 -1 3 23 -1 16 -1 -1 21 -1 0 -1 -1 0 0 -1;
     24 5 26 7 1 -1 -1 15 24 15 -1 8 -1 13 -1 13 -1 11 -1 -1 -1 -1 0 0;
     2 2 19 14 24 1 15 19 -1 21 -1 2 -1 24 -1 3 -1 2 1 -1 -1 -1 -1 0];
H_num = ldpcQuasiCyclicMatrix(blockSize, P);
H = logical(sparse(H_num));   % разреженная логическая матрица
[Nrows, Ncols] = size(H);     % 486 x 1944
K_real = Ncols - Nrows;       % 1458
N_real = Ncols;               % 1944
fprintf('Реальный LDPC код: K_real = %d, N_real = %d, скорость = %.3f\n', K_real, N_real, K_real/N_real);

% Конфигурации для реального кода
cfgEnc = ldpcEncoderConfig(H);
cfgDec = ldpcDecoderConfig(H);

% Целевые параметры
K_target = 1912;
N_target = 2552;
pad_info = K_target - K_real;   % 454
pad_code = N_target - N_real;   % 608

% Параметры симуляции
EbNo_dB = 0:0.5:4;
numFramesMax = 20000;
targetFER = 200;
maxIter = 50;
M = 4; bps = log2(M);
R_total = K_target / N_target;   % ~0.749

% Создание папок для результатов 
if ~exist('../Results','dir'), mkdir('../Results'); end
if ~exist('../Figures','dir'), mkdir('../Figures'); end

results = zeros(length(EbNo_dB),4);

for idx = 1:length(EbNo_dB)
    EbNo = EbNo_dB(idx);
    totalBitErrors = 0; totalBits = 0;
    totalFrameErrors = 0; totalFrames = 0;
    fprintf('\n--- Eb/N0 = %.2f dB ---\n', EbNo);
    
    % Расчёт дисперсии шума для QPSK 
    EsN0_lin = 10^((EbNo + 10*log10(R_total) + 10*log10(bps))/10);
    sigma2 = 1/(2*EsN0_lin);
    
    while (totalFrameErrors < targetFER) && (totalFrames < numFramesMax)
        % Генерация K_target бит (первые K_real случайные, остальные нули)
        data_real = randi([0 1], K_real, 1);
        data_target = [data_real; zeros(K_target - K_real, 1)];
        
        % LDPC кодирование реальных бит (длина N_real)
        encData_real = ldpcEncode(data_real, cfgEnc);
        
        % Дополнение кодового слова нулями до N_target
        encData_target = [encData_real; zeros(N_target - N_real, 1)];
        
        % QPSK модуляция
        tx = pskmod(encData_target, M, pi/4, 'gray', 'InputType','bit');
        
        % AWGN канал
        noise = sqrt(sigma2)*(randn(size(tx)) + 1j*randn(size(tx)));
        rx = tx + noise;
        
        % Демодуляция LLR (мягкие решения)
        rxLLR_target = pskdemod(rx, M, pi/4, 'gray', 'OutputType','llr', 'NoiseVariance', sigma2);
        
        % Обрезка до реальной длины N_real
        rxLLR_real = rxLLR_target(1:N_real);
        
        % LDPC декодирование реальных бит
        decData_real = ldpcDecode(rxLLR_real, cfgDec, maxIter);
        
        % Восстановление целевой длины (добавление нулей)
        decData_target = [decData_real; zeros(K_target - K_real, 1)];
        
        % Подсчёт ошибок (по всем K_target битам)
        errs = sum(decData_target ~= data_target);
        totalBitErrors = totalBitErrors + errs;
        totalBits = totalBits + K_target;
        totalFrameErrors = totalFrameErrors + (errs > 0);
        totalFrames = totalFrames + 1;
        
        if mod(totalFrames, 500) == 0
            fprintf('  frames: %d, текущая FER: %.4f\n', totalFrames, totalFrameErrors/totalFrames);
        end
    end
    
    BER = totalBitErrors / totalBits;
    FER = totalFrameErrors / totalFrames;
    results(idx,:) = [EbNo, BER, FER, totalFrames];
    fprintf('  >> BER = %.2e, FER = %.2e\n', BER, FER);
end

% Сохранение и построение графиков 
csvwrite('Results/ldpc_embb.csv', results);
figure;
semilogy(results(:,1), results(:,2), 'b-o', 'LineWidth', 1.5); hold on;
semilogy(results(:,1), results(:,3), 'r--s', 'LineWidth', 1.5);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER / FER');
legend('BER', 'FER', 'Location', 'southwest');
title(sprintf('LDPC (eMBB) K=%d, N=%d', K_target, N_target));
saveas(gcf, 'Figures/ldpc_embb.png');

fprintf('\n Результаты для K=%d, N=%d сохранены.\n', K_target, N_target);
end