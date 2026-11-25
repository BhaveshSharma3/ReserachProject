N = 1e5;
SNR_dB = -30:5:40;
BER_empirical = zeros(size(SNR_dB));
BER_theoretical = zeros(size(SNR_dB));
BER_SISO = zeros(size(SNR_dB));
BER_SISO_theoretical = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    data = randi([0 1], 1, N);
    tx = 2*data - 1;
    snr_linear = 10^(SNR_dB(i)/10);
    noise_power = 1 / snr_linear;
    noise = sqrt(noise_power/2) * (randn(size(tx)) + 1j*randn(size(tx)));
    rx = tx + noise;
    data_hat = rx > 0;
    BER_empirical(i) = sum(data ~= data_hat)/N;
    
    % Theoretical BER for AWGN
    BER_theoretical(i) = qfunc(sqrt(2 *snr_linear));
    
    % SISO channel
    h = (randn(1, N) + 1j * randn(1, N)) / sqrt(2);
    rx_channel = h .* tx;
    noise = (randn(1, N) + 1j * randn(1, N)) * sqrt(noise_power/2);
    rx1 = rx_channel + noise;
    rx_eq = rx1 ./ h;
    detected_bits = real(rx_eq) > 0;
    BER_SISO(i) = sum(data ~= detected_bits) / N;

    %Theoretical BER for SISO
    BER_SISO_theoretical(i) = 0.5*(1-sqrt(snr_linear/(1 + snr_linear))); 


end
BER_empirical(BER_empirical == 0) = 3.87211e-06;
figure;
hold on;
semilogy(SNR_dB, BER_empirical, 'o-', 'LineWidth', 2, 'DisplayName', 'Empirical AWGN');
semilogy(SNR_dB, BER_theoretical, 's--', 'LineWidth', 2, 'DisplayName', 'Theoretical AWGN');
semilogy(SNR_dB, BER_SISO, 'd-', 'LineWidth', 1.5, 'DisplayName', 'Empirical SISO');
semilogy(SNR_dB, BER_SISO_theoretical, 'x--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical SISO');
hold off;
grid on;
set(gca, 'YScale', 'log');
ylim([1e-4 1]);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for BPSK');
legend('show');