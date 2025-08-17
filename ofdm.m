clc; 
clear; 
%close all;

%% Parametri OFDM
bps = 4;              % Bit per simbolo -> 16-QAM
M = 2^bps;            % Ordine modulazione
nFFT = 256;            % Numero sottoportanti OFDM
%cp_len = 16;          % Lunghezza del Cyclic Prefix
%EbN0_dB = 12;         % Eb/N0 in dB
EbN0_dB_min = 3;
EbN0_dB_max = 15;
EbN0_vec = EbN0_dB_min:1:EbN0_dB_max;
BER_vec = zeros(1,length(EbN0_vec));
N = 5000;              % Numero di pacchetti inviati

%% CICLIO BER
for j=1:length(EbN0_vec)
    EbN0_dB = EbN0_vec(j);
    fprintf('Eb/N0 = %f \n', EbN0_dB);

    %% Calcolo fattori e relazione tra SNR e Eb/N0
    bps_term = 10*log10(bps);
    %cp_term  = 10*log10(nFFT / (nFFT + cp_len));
    SNR_dB = EbN0_dB + bps_term ;
    
    % Inizializzazione contatori errori
    total_err = 0;
    total_bits = 0;
    
    %% CICLO PER INVIO DI N PACCHETTI
    for i=1:N
        %% Trasmissione
        % 1. Generazione di bit casuali
        num_bits = nFFT * bps;              % Numero totale di bit: 1024
        tx_bits = randi([0 1], num_bits, 1);  % Vettore di bit casuali
        % 2. Raggruppamento dei bit in simboli QAM
        bit_matrix = reshape(tx_bits, bps, []).';   % Ogni riga è un simbolo
        tx_symbols = bi2de(bit_matrix, 'left-msb'); % Conversione binario → decimale (simboli QAM)
        % 3. Modulazione QAM
        tx_grid = qammod(tx_symbols, M, 'UnitAveragePower', true);
        % 4. IFFT -> dominio tempo
        tx_ifft = ifft(tx_grid, nFFT)*sqrt(nFFT);                             
        % 5. Aggiunta del Cyclic Prefix
        %tx_out_cp = [tx_ifft(end - cp_len + 1:end); tx_ifft];
    
        %% Canale AWGN
        rx_in = awgn(tx_ifft, SNR_dB, 'measured');
    
        %% Ricezione
        % Rimozione del Cyclic Prefix
        %rx_in = rx_in_cp(cp_len+1:end);
        % FFT -> ritorno al dominio della frequenza
        rx_grid = fft(rx_in, nFFT)/sqrt(nFFT);
        % Demodulazione QAM (ottenimento dei simboli)
        rx_symbols = qamdemod(rx_grid, M, 'UnitAveragePower', true);
        % Da simboli a bit
        rx_bits_matrix = de2bi(rx_symbols, bps, 'left-msb');
        rx_bits = reshape(rx_bits_matrix.', [], 1);   % Vettore colonna
       
        %% Accumulo errori
        [err, ~] = biterr(tx_bits, rx_bits);
        total_err = total_err + err;
        total_bits = total_bits + num_bits;
    end
    
    %% Calcolo finale della BER
    ber = total_err / total_bits;
    fprintf("SNR: %f\n", SNR_dB);
    fprintf("BER = %.5f (%d errori su %d bit totali)\n", ber, total_err, total_bits);
    BER_vec(j) = ber ;
end

%% PLOT UTILI
figure (2); % Plot BER
semilogy(EbN0_vec, BER_vec);
xlabel(' E_{b} / N_{0} (dB) ');
ylabel('BER');
grid on;
hold on;

figure; % Plot costellazione a tx
scatter(real(tx_grid), imag(tx_grid), 'filled');
xlabel('Re'); ylabel('Im');
title('Diagramma di costellazione QAM inviato');
grid on; axis equal;

figure; % Plot costellazione a rx
scatter(real(rx_grid), imag(rx_grid), 'filled');
xlabel('Re'); ylabel('Im');
title('Diagramma di costellazione QAM ricevuto');
grid on; axis equal;
