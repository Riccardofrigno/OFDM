clc; clear; %close all;

%% OFDM + LDPC DVB-S2 + Pilot-Based Channel Estimation in Multipath+Doppler Fading (16-QAM)
% Parameters
M        = 256;            % number of subcarriers
N        = 64;             % number of OFDM symbols per frame
padLen   = 16;             % cyclic prefix length (samples)
bps      = 4;              % bits per symbol (16-QAM)
modOrder = 2^bps;
EbN0_vec = 2:0.5:10;        % Eb/N0 range (dB)
Num      = 10;             % number of LDPC-coded frames

%% LDPC DVB-S2 configuration
rate      = 3/4;
dvbs2H    = dvbs2ldpc(rate);
encCfg    = ldpcEncoderConfig(dvbs2H);
decCfg    = ldpcDecoderConfig(dvbs2H);
K         = encCfg.NumInformationBits;        % info bits per codeword
N_ldpc    = encCfg.NumParityCheckBits + K;    % total codeword length

%% OFDM grid padding
ofdm_bits = M * bps * N;
padOFDM   = ofdm_bits - N_ldpc;

%% Channel: multipath + Doppler
df       = 15e3;
fsamp    = M * df;
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
        %% 1) LDPC encoding
        txBits = randi([0 1], K, 1);             % generate information bits
        coded  = ldpcEncode(txBits, encCfg);     % encode to N_ldpc bits
        codedP = [coded; zeros(padOFDM,1)];      % pad to fill OFDM grid

        %% 2) QAM mapping and OFDM modulation
        symIdx  = reshape(codedP, bps, []).';
        qamSyms = qammod(bi2de(symIdx,'left-msb'), modOrder, 'gray', 'UnitAveragePower', true);
        dataGrid= reshape(qamSyms, M, N);
        txData  = ofdmmod(dataGrid, M, padLen);

        %% 3) Transmit through channel + AWGN
        rxD = dopplerChannel(txData, fsamp, chanParams);
        rxD = awgn(rxD, SNRdB, 'measured');
        rxGrid = ofdmdemod(rxD(1:length(txData)), M, padLen);

        %% 4) Equalization
        eqd = (conj(H_est) .* rxGrid) ./ (abs(H_est).^2 + n0);

        %% 5) Soft demapping to LLR
        llr = qamdemod(eqd(:), modOrder, 'gray', ...
               'OutputType','approxllr','UnitAveragePower',true, ...
               'NoiseVariance', n0);
        llr = llr(1:N_ldpc);

        %% 6) LDPC decoding
        decBits = ldpcDecode(llr, decCfg, 50);
        rxBits  = decBits(1:K);

        %% 7) BER count
        totalErr  = totalErr + sum(rxBits ~= txBits);
        totalBits = totalBits + K;
    end

    BER_vec(idx) = totalErr / totalBits;
    fprintf('BER = %.2e\n\n', BER_vec(idx));
end

%% Plot BER
figure;
semilogy(EbN0_vec, BER_vec, '-o','LineWidth',1.5);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('OFDM + LDPC(3/4) + 16-QAM over Multipath + Doppler');
grid on;
