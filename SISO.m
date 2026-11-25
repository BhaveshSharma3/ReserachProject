N = 1e7;
snr_db = -30:5:15;
ber = zeros(size(snr_db));
for i = 1:length(snr_db)
data_bits = randi([0 1], 1, N);
tx = 2*data_bits - 1;
snr_linear = 10^(snr_db(i)/10);
h = (randn(1, N) + 1j * randn(1, N)) / sqrt(2);
rx_channel = h .* tx;
noise = (randn(1, N) + 1j * randn(1, N)) / sqrt(2 * snr_linear);
rx = rx_channel + noise;
rx_eq = rx ./ h;
detected_bits = real(rx_eq) > 0;
ber(i) = sum(data_bits ~= detected_bits) / N;
end
ber(ber == 0) = 1e-10;
semilogy(snr_db, ber, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('BER vs SNR for BPSK over Rayleigh Fading Channel', 'FontSize', 1);
axis([min(snr_db) max(snr_db) 1e-3 1]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);