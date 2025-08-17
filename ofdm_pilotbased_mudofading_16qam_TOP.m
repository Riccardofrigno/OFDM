clc; clear; close all;

%% OFDM with Pilot-Based Channel Estimation in Multipath+Doppler Fading (16-QAM, no LDPC)
% Parameters
M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 10;     % make this larger than the channel delay spread channel in samples
padType = 'ZP';  % this example requires ZP for ISI mitigation
bps      = 4;              % bits per symbol (16-QAM)
modOrder = 2^bps;
EbN0_vec = 2:0.5:10;        % Eb/N0 range (dB)
Num      = 100;            % number of OFDM frames to simulate

%% OFDM bitstream length
totalBitsPerFrame = M * bps * N;

%% Channel: multipath + Doppler
chanParams.pathDelays       = [0, 5, 8];
chanParams.pathGains        = [1, 0.7, 0.5];
chanParams.pathDopplerFreqs = [0, -100, 150];

%% Precompute pilot OFDM signal
pilotGrid = exp(1i*pi/4) * ones(M, N);
txPilots  = ofdmmod(pilotGrid, M, padLen);

BER_vec = zeros(size(EbN0_vec));

for idx = 1:length(EbN0_vec)
    EbN0dB = EbN0_vec(idx);
    fprintf('Eb/N0 = %.1f dB\n', EbN0dB);
    SNRdB  = EbN0dB + 10*log10(bps);
    n0     = 1/(10^(SNRdB/10));

    %% Channel estimation via pilots
    rxP = dopplerChannel(txPilots, fsamp, chanParams);
    rxP = awgn(rxP, SNRdB, 'measured');
    rxPilots = ofdmdemod(rxP(1:length(txPilots)), M, padLen);
    H_est = rxPilots ./ pilotGrid;

    totalErr  = 0;
    totalBits = 0;

    for fr = 1:Num
        %% Generate random bits, map to QAM and OFDM modulate
        bits = randi([0 1], totalBitsPerFrame, 1);
        symIdx = bi2de(reshape(bits, bps, []).','left-msb');
        qamSyms = qammod(symIdx, modOrder, 'gray', 'UnitAveragePower', true);
        dataGrid = reshape(qamSyms, M, N);
        txData   = ofdmmod(dataGrid, M, padLen);

        %% Transmit through channel + AWGN
        rxD = dopplerChannel(txData, fsamp, chanParams);
        rxD = awgn(rxD, SNRdB, 'measured');
        rxGrid = ofdmdemod(rxD(1:length(txData)), M, padLen);

        %% Equalize
        eqd = (conj(H_est) .* rxGrid) ./ (abs(H_est).^2 + n0);

        %% Demap to bits
        rxBits = de2bi(qamdemod(eqd(:), modOrder, 'gray', 'OutputType','integer', ...
                     'UnitAveragePower', true, 'NoiseVariance', n0), bps, 'left-msb');
        rxBits = rxBits.';
        rxBits = rxBits(:);

        %% BER count
        totalErr  = totalErr + sum(rxBits ~= bits);
        totalBits = totalBits + numel(bits);
    end

    BER_vec(idx) = totalErr / totalBits;
    fprintf('BER = %.2e\n\n', BER_vec(idx));
end

%% Plot BER
figure;
semilogy(EbN0_vec, BER_vec, '-o','LineWidth',1.5);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('OFDM + 16-QAM + Pilot-FDE over Multipath/Doppler');
grid on;
