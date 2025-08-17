clc; clear; close all;

%% Parametri OFDM
bps       = 4;               
M         = 2^bps;           
nFFT      = 256;             
ofdm_bits = nFFT * bps;      
EbN0_vec  = 1:7;            
BER_vec   = zeros(size(EbN0_vec));

%% Parametri LDPC DVB-S2 rate 3/4
rate      = 3/4;
dvbs2H    = dvbs2ldpc(rate);
cfgEnc    = ldpcEncoderConfig(dvbs2H);
cfgDec    = ldpcDecoderConfig(dvbs2H);
K         = cfgEnc.NumInformationBits;
num_pck   = floor(K / ofdm_bits);

%% Caricamento immagine e vettorializzazione (una sola volta)
img = imread('lena.png');
if size(img,3)==3, img = rgb2gray(img); end
[H,W]      = size(img);
imgVec     = img(:);
byteMatrix = de2bi(imgVec,8,'left-msb');
txBitsImage= reshape(byteMatrix.', [], 1);

totalBits  = numel(txBitsImage);
numFrames  = ceil(totalBits/(num_pck*ofdm_bits));
padImg     = numFrames*num_pck*ofdm_bits - totalBits;
txBitsImage= [txBitsImage; zeros(padImg,1)];

% **Calcolo dinamico del numero di frame OFDM**
Num = numel(txBitsImage) / ofdm_bits;  % ora Num è intero

%% Pre-allocazione per la ricostruzione immagine
rxBitsImg  = zeros(totalBits,1);

frameCount = 0;

%% Ciclo su Eb/N0
for idx = 1:length(EbN0_vec)
    EbN0_dB = EbN0_vec(idx);
    fprintf('Eb/N0 = %d dB\n', EbN0_dB);
    SNR_dB  = EbN0_dB + 10*log10(bps*rate);

    total_err  = 0;
    total_bits = 0;
    bufferBits = [];
    
    % Reset contatore frame per ogni Eb/N0
    frameCount = 0;

    for i = 1:Num
        startBit = (i-1)*ofdm_bits + 1;
        if startBit + ofdm_bits - 1 > numel(txBitsImage)
            break;
        end
        tx_bits  = txBitsImage(startBit : startBit + ofdm_bits - 1);
        bufferBits = [bufferBits; tx_bits];

        if numel(bufferBits) >= num_pck*ofdm_bits
            % Preleva blocco dati
            data_in = bufferBits(1:num_pck*ofdm_bits);
            bufferBits(1:num_pck*ofdm_bits) = [];
            %% INIZIO MODULAZIONE
            % Padding fino a K bit e conversione a double
            dataBits = double([data_in; zeros(K - numel(data_in),1)]);
                     
            % Codifica LDPC + zero-pad per quadrare OFDM
            coded = ldpcEncode(dataBits, cfgEnc);
            padO  = ceil(numel(coded)/ofdm_bits)*ofdm_bits - numel(coded);
            coded_p = [coded; zeros(padO,1)];

            % 16-QAM + OFDM IFFT
            symIdx = bi2de(reshape(coded_p, bps, []).','left-msb');
            qamSyms= qammod(symIdx, M, 'gray','UnitAveragePower',true);
            txGrid = reshape(qamSyms,nFFT,[]);
            txTime = ifft(txGrid,nFFT)*sqrt(nFFT);
            txSig  = txTime(:);

            % Canale AWGN
            rxSig = awgn(txSig, SNR_dB, 'measured');

            % FFT per tornare al dominio frequenza
            rxMat   = reshape(rxSig, nFFT, []);
            rxGrid  = fft(rxMat,nFFT)/sqrt(nFFT);
            
            % Soft‐demapping con varianza del rumore
            noiseVar = 10^(-SNR_dB/10);
            rx_llr = qamdemod(rxGrid(:), M, ...
                       'gray', ...
                       'OutputType','llr', ...
                       'UnitAveragePower',true, ...
                       'NoiseVariance',noiseVar);
            
            % Soft‐decoding LDPC
            rx_llr = rx_llr(1:cfgDec.BlockLength);  % 64800 LLR
            rx_decoded = ldpcDecode(rx_llr, cfgDec, 50);
            
            % Estrai dati e rimuovi padding
            rx_frame = rx_decoded(1:num_pck*ofdm_bits);


            % Accumula errori
            [err, ~] = biterr(data_in, rx_frame);
            total_err  = total_err + err;
            total_bits = total_bits + numel(data_in);

            % Accumulo bit nel vettore immagine ricevuta
            frameCount = frameCount + 1;
            idxStart   = (frameCount-1)*numel(data_in) + 1;
            idxEnd     = idxStart + numel(data_in) - 1;
            rxBitsImg(idxStart:idxEnd) = rx_frame;
        end
    end

    % Calcola BER e visualizza
    BER_vec(idx) = total_err/total_bits;
    fprintf('BER = %.8f (%d errori su %d bit)\n\n', BER_vec(idx), total_err, total_bits);

    % Ricostruzione e visualizzazione immagine
    imgBits = rxBitsImg(1:totalBits);
    byteM   = reshape(imgBits,8,[]).';
    rxImgVec= bi2de(byteM,'left-msb');
    rxImg   = reshape(uint8(rxImgVec), H, W);
    figure; imshow(rxImg);
    title(sprintf('Received image at Eb/N0 = %d dB', EbN0_dB));
end

%% Plot BER
figure;
semilogy(EbN0_vec,BER_vec,'-o');
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('OFDM+LDPC DVB-S2 Rate 3/4, 16QAM');
grid on; ylim([1e-6,1]);
