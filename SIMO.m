%% Combined AWGN + SISO + SIMO BER Simulation (BPSK)
clc;
clear;

%% Simulation parameters
N = 1e5;                     % Number of bits/symbols per trial
snr_db = -30:5:40;           % SNR range in dB (combined range)
snr_lin = 10.^(snr_db/10);

%% Preallocate arrays for AWGN & SISO
BER_empirical_awgn = zeros(size(snr_db));
BER_theoretical_awgn = zeros(size(snr_db));
BER_empirical_siso = zeros(size(snr_db));
BER_theoretical_siso = zeros(size(snr_db));

%% AWGN and SISO simulations (vectorized per SNR)
for ii = 1:length(snr_db)
    % ----- AWGN (no fading) -----
    data = randi([0 1], 1, N);
    tx = 2*data - 1;                       % BPSK: -1 / +1
    snr_linear = snr_lin(ii);
    noise_power = 1 / snr_linear;         % assuming signal power = 1
    noise = sqrt(noise_power/2) * (randn(size(tx)) + 1j*randn(size(tx)));
    rx = tx + noise;
    data_hat = real(rx) > 0;
    BER_empirical_awgn(ii) = sum(data ~= data_hat) / N;
    
    % theoretical AWGN BER for BPSK (coherent)
    BER_theoretical_awgn(ii) = qfunc(sqrt(2 * snr_linear));
    
    % ----- SISO (flat Rayleigh fading) -----
    h = (randn(1, N) + 1j * randn(1, N)) / sqrt(2);
    noise_siso = sqrt(noise_power/2) * (randn(1, N) + 1j * randn(1, N));
    rx_siso = h .* tx + noise_siso;
    rx_eq = rx_siso ./ h;                  % simple equalization (zero-forcing)
    detected_siso = real(rx_eq) > 0;
    BER_empirical_siso(ii) = sum(data ~= detected_siso) / N;
    
    % theoretical SISO (Rayleigh single-branch) closed-form BE R
    BER_theoretical_siso(ii) = 0.5 * (1 - sqrt(snr_linear ./ (1 + snr_linear)));
end

% avoid exact zeros for plotting (replace zeros with numerical floor)
min_nonzero = 1/(N);  % smallest resolvable BER with N bits
BER_empirical_awgn(BER_empirical_awgn == 0) = min_nonzero;
BER_empirical_siso(BER_empirical_siso == 0) = min_nonzero;

%% SIMO configurations
simo_configs = [
    1, 1;
    1, 2;   % 1x2
    1, 4;   % 1x4
    1, 8;   % 1x8
    1, 16;  % 1x16
];

numConfigs = size(simo_configs,1);
BER_simo = zeros(numConfigs, length(snr_db));

% Loop through SIMO configurations
for cfg = 1:numConfigs
    numTx = simo_configs(cfg,1);
    numRx = simo_configs(cfg,2);
    fprintf('Simulating %dx%d SIMO ...\n', numTx, numRx);
    
    for ii = 1:length(snr_db)
        snr_linear = snr_lin(ii);
        % The original script used noise_var = numRx / snr_linear.
        % We'll follow that normalization (interpreting received energy scaling),
        % but ensure consistency: noise variance per complex component:
        noise_var = numRx / snr_linear;
        
        % Generate tx bits/symbols for numTx x N (here numTx = 1)
        data_bits = randi([0 1], numTx, N);
        tx_symbols = 2 * data_bits - 1;   % values -1 / +1
        
        % Vectorized SIMO Monte Carlo:
        % Generate channel H: size numRx x N (each column is channel for one symbol)
        H = (randn(numRx, N) + 1j*randn(numRx, N)) / sqrt(2);
        
        % Generate noise matrix: numRx x N
        noise = sqrt(noise_var/2) * (randn(numRx, N) + 1j*randn(numRx, N));
        
        % Received signal: y = H .* tx + noise  (tx is 1 x N, broadcast across rows)
        y = H .* repmat(tx_symbols, numRx, 1) + noise;
        
        % MRC / matched filter combining:
        % x_hat (1 x N) = sum(conj(H) .* y, over rows) ./ sum(abs(H).^2, over rows)
        numerator = sum(conj(H).*y, 1);            % 1 x N
        denom = sum(abs(H).^2, 1);                 % 1 x N
        x_hat = numerator ./ denom;                % 1 x N (complex)
        
        % Decisions (real part)
        detected = real(x_hat) > 0;                % logical 1 x N
        
        % Count bit errors (for single tx, compare with data_bits(1,:))
        bit_errors = sum(detected ~= data_bits(1,:));
        BER_simo(cfg, ii) = bit_errors / N;
        
        % optional progress print (comment out if verbose)
        % fprintf('Config %dx%d: SNR=%2d dB, BER=%e\n', numTx, numRx, snr_db(ii), BER_simo(cfg,ii));
    end
    
    % replace zeros if any
    BER_simo(cfg, BER_simo(cfg,:) == 0) = min_nonzero;
end

%% Plotting - combine everything on one figure
figure('Position', [100, 100, 1000, 650]);
hold on;
grid on;
box on;

% AWGN empirical & theoretical
semilogy(snr_db, BER_empirical_awgn, 'o-', 'LineWidth', 1.8, 'DisplayName', 'Empirical AWGN');
semilogy(snr_db, BER_theoretical_awgn, 's--', 'LineWidth', 1.6, 'DisplayName', 'Theoretical AWGN');

% SISO empirical & theoretical
semilogy(snr_db, BER_empirical_siso, 'd-', 'LineWidth', 1.6, 'DisplayName', 'Empirical SISO');
semilogy(snr_db, BER_theoretical_siso, 'x--', 'LineWidth', 1.6, 'DisplayName', 'Theoretical SISO');

% SIMO curves with distinct markers/colors
colors = {'b','r','g','m','c','k','y'};
markers = {'o','s','d','^','v','<'};
for cfg = 1:numConfigs
    numTx = simo_configs(cfg,1);
    numRx = simo_configs(cfg,2);
    semilogy(snr_db, BER_simo(cfg,:), ...
        'Color', colors{mod(cfg-1,length(colors))+1}, ...
        'Marker', markers{mod(cfg-1,length(markers))+1}, ...
        'LineWidth', 1.8, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', colors{mod(cfg-1,length(colors))+1}, ...
        'DisplayName', sprintf('%d×%d SIMO (empirical)', numTx, numRx));
end

% Axis & labels
set(gca, 'YScale', 'log');
xlabel('SNR (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER) [Log Scale]', 'FontSize', 13, 'FontWeight', 'bold');
title('BER vs SNR: AWGN, SISO, and SIMO (BPSK)', 'FontSize', 15, 'FontWeight', 'bold');

xlim([min(snr_db) max(snr_db)]);
ylim([1e-4 1]);   % expanded lower limit for better SIMO visibility
yticks([1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1]);
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'});

legend('Location', 'southwest', 'FontSize', 10);
set(gca, 'FontSize', 12, 'LineWidth', 1.1);
grid on;
hold off;

%% Final note
fprintf('\nSimulation complete. Plotted AWGN (empirical/theoretical), SISO (empirical/theoretical), and SIMO empirical for %d configurations.\n', numConfigs);
