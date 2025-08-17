clc; clear; %close all;

%% Parametri OFDM
bps = 4; % bit per simbolo (16-QAM)
M_mod = 2^bps; % ordine modulazione
nFFT = 256; % numero di sottoportanti OFDM
ofdm_bits = nFFT * bps; % 1024 bit per simbolo OFDM
EbN0_vec = 2:0.2:8; % Eb/N0 in dB 
BER_vec = zeros(size(EbN0_vec));
Num = 10; % numero di pacchetti OFDM totali

%% Parametri LDPC DVB-S2 rate 3/4
rate = 3/4; % rate LDPC DVB-S2 (numerico)
dvbs2H = dvbs2ldpc(rate); % matrice di parità DVB-S2
cfgEnc = ldpcEncoderConfig(dvbs2H);
cfgDec = ldpcDecoderConfig(dvbs2H);
K = cfgEnc.NumInformationBits; % 48600 bit informazione

% Calcolo numero di blocchi OFDM necessari
num_pck = floor(K / ofdm_bits); % 47 blocchi OFDM

%% Ciclo BER
for idx = 1:length(EbN0_vec)
    EbN0_dB = EbN0_vec(idx);
    fprintf('Eb/N0 = %d dB\n', EbN0_dB);
    SNR_dB  = EbN0_dB + 10*log10(bps);

    total_err = 0;
    total_bits = 0;
    bufferBits = [];

    for i = 1:Num
        % Genera pacchetto casuale di bits
        tx_bits = randi([0 1], K, 1);

            % Codifica LDPC e padding per OFDM
            coded = ldpcEncode(tx_bits, cfgEnc);
            padO  = ceil(numel(coded)/ofdm_bits) * ofdm_bits - numel(coded);
            coded_p = [coded; zeros(padO,1)];

            % Mapping 16-QAM Gray
            symIdx = bi2de(reshape(coded_p, bps, []).', 'left-msb');
            qamSyms = qammod(symIdx, M_mod, 'gray', 'UnitAveragePower', true);

            % OFDM IFFT
            txGrid = reshape(qamSyms, nFFT, []);
            txTime = ifft(txGrid, nFFT) * sqrt(nFFT);
            txSig  = txTime(:);

            % Canale AWGN
            rxSig = awgn(txSig, SNR_dB, 'measured');

            % Ricezione OFDM: FFT e demapping
            rxTime = reshape(rxSig, nFFT, []);
            rxGrid = fft(rxTime, nFFT) / sqrt(nFFT);

            % soft‐demapping diretto in LLR continui
            noiseVar = 1/(10^(SNR_dB/10));  
            rx_llr   = qamdemod(rxGrid(:), M_mod, ...
                         'gray','OutputType','approxllr', ...
                         'UnitAveragePower',true, ...
                         'NoiseVariance',noiseVar);
            
            % prendo i primi BlockLength LLR e decodifico
            rx_llr = rx_llr(1:cfgDec.BlockLength);
            dec    = ldpcDecode(rx_llr, cfgDec, 50);

            % Estrai dati e calcola errori
            rx_data = dec(1 : numel(tx_bits));
            [err, ~] = biterr(tx_bits, rx_data);
            total_err  = total_err + err;
            total_bits = total_bits + numel(tx_bits);
        
    end

    ber = total_err / total_bits;
    BER_vec(idx) = ber;
    ber = 0;
    fprintf('BER = %.8f (%d errori su %d bit)\n\n', BER_vec(idx), total_err, total_bits);

end

%% 36 ogni 1440 simboli dopo qam piloti per stiam canale

%% Plot BER
figure;
semilogy(EbN0_vec, BER_vec, '-o');
xlabel('E_b/N_0 (dB)');
ylabel('BER');
title('FER con LDPC DVB-S2 rate 3/4, OFDM 256 subcarriers, 16QAM');
grid on;
hold on;

% Visualizza in tabella tutti i valori di Eb/N0 e BER
T = table(EbN0_vec.', BER_vec.', 'VariableNames', {'EbN0_dB', 'BER'});
disp(T);