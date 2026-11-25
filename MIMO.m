clc;
clear;

N = 1e5;        
snr_db = -30:5:20;         
numTx = 32;                
numRx = 64;                
ber = zeros(size(snr_db));
for i = 1:length(snr_db)
    data_bits = randi([0 1], numTx, N);
    tx_symbols = (2*data_bits - 1); 
    snr_linear = 10^(snr_db(i)/10);
    noise_var =  numRx /(snr_linear); 
    err = 0;
    for n = 1:N
        H = (randn(numRx, numTx) + 1j*randn(numRx, numTx))/sqrt(2); 
        noise = sqrt(noise_var/2) * (randn(numRx,1) + 1j*randn(numRx,1)); 
        y = H * tx_symbols(:, n) + noise;     
        x_hat = pinv(H) * y;                   
        detected = real(x_hat) > 0;              
        err = err + sum(detected ~= data_bits(:, n));
    end
    ber(i) = err / (numTx * N);
    fprintf('SNR = %d dB, BER = %e\n', snr_db(i), ber(i));
end
figure;
semilogy(snr_db, ber, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlim([1 20]);
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('BER vs SNR for 8x4 MIMO with Zero-Forcing Detection', ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontAngle', 'italic');
axis([min(snr_db) max(snr_db) 1e-6 1]);  
set(gca, 'FontSize', 12, 'LineWidth', 1.2);