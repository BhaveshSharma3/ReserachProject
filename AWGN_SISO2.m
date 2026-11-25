clc;
clear;
N = 1e5;
SNR_dB = -30:5:40;
numTx = 2;
numRx = 2;

% Initialize BER arrays
BER_SISO = zeros(size(SNR_dB));
BER_SISO_theoretical = zeros(size(SNR_dB));
BER_MIMO_AWGN = zeros(size(SNR_dB));
BER_MIMO_AWGN_theoretical = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    snr_linear = 10^(SNR_dB(i)/10);
    noise_power = 1 / snr_linear;
    
    %% SISO Rayleigh Fading Channel
    data = randi([0 1], 1, N);
    tx = 2*data - 1;  % BPSK modulation
    
    % Rayleigh fading channel
    h = (randn(1, N) + 1j * randn(1, N)) / sqrt(2);
    rx_channel = h .* tx;
    
    % Add noise
    noise = (randn(1, N) + 1j * randn(1, N)) * sqrt(noise_power/2);
    rx1 = rx_channel + noise;
    
    % Equalization
    rx_eq = rx1 ./ h;
    
    % Detection
    detected_bits = real(rx_eq) > 0;
    BER_SISO(i) = sum(data ~= detected_bits) / N;
    
    % Theoretical BER for SISO Rayleigh
    BER_SISO_theoretical(i) = 0.5*(1-sqrt(snr_linear/(1 + snr_linear)));
    
    %% 2×2 MIMO AWGN System
    data_bits_mimo = randi([0 1], numTx, N);
    tx_symbols = (2*data_bits_mimo - 1);  % BPSK modulation
    
    err_awgn = 0;
    
    for n = 1:N
        x = tx_symbols(:, n);
        
        % Pure AWGN channel (2×2)
        H_awgn = eye(numRx, numTx);
        noise_awgn = sqrt(noise_power/2) * (randn(numRx,1) + 1j*randn(numRx,1));
        y_awgn = H_awgn * x + noise_awgn;
        
        % Zero-Forcing detection for AWGN
        x_hat_awgn = pinv(H_awgn) * y_awgn;
        detected_awgn = real(x_hat_awgn) > 0;
        err_awgn = err_awgn + sum(detected_awgn ~= data_bits_mimo(:, n));
    end
    
    BER_MIMO_AWGN(i) = err_awgn / (numTx * N);
    
    % Theoretical BER for 2×2 MIMO AWGN with ZF
    % For 2×2 MIMO AWGN with ZF detection, each stream experiences same SNR as SISO AWGN
    BER_MIMO_AWGN_theoretical(i) = qfunc(sqrt(2*snr_linear));
    
    % Display progress
    fprintf('SNR = %d dB: SISO_Rayleigh = %e, MIMO_AWGN = %e\n', ...
        SNR_dB(i), BER_SISO(i), BER_MIMO_AWGN(i));
end

% Handle zero values for logarithmic plot
min_ber = 1e-8;
BER_SISO(BER_SISO == 0) = min_ber;
BER_SISO_theoretical(BER_SISO_theoretical == 0) = min_ber;
BER_MIMO_AWGN(BER_MIMO_AWGN == 0) = min_ber;
BER_MIMO_AWGN_theoretical(BER_MIMO_AWGN_theoretical == 0) = min_ber;

% Plotting
figure('Position', [100, 100, 800, 600]);
semilogy(SNR_dB, BER_SISO, 'm-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'SISO Rayleigh (Empirical)');
hold on;
semilogy(SNR_dB, BER_SISO_theoretical, 'y--', 'LineWidth', 2, ...
    'DisplayName', 'SISO Rayleigh (Theoretical)');
semilogy(SNR_dB, BER_MIMO_AWGN, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', '2×2 MIMO AWGN (Empirical)');
semilogy(SNR_dB, BER_MIMO_AWGN_theoretical, 'b--', 'LineWidth', 2, ...
    'DisplayName', '2×2 MIMO AWGN (Theoretical)');
hold off;

grid on;
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('BER Performance: SISO Rayleigh vs 2×2 MIMO AWGN', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
axis([min(SNR_dB) max(SNR_dB) 1e-4 1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);


