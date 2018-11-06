clc; close all; clear global; clearvars;

load('Input.mat','Npx','in','M');
in_qam = qammod(in,16,'UnitAveragePower',true).';		% QAM modulation
SNR_db_vect = 0:10;					% SNR vector (in dB)
sigma_a = var(in_qam);

%% OFDM TRANSMITTER
a_pad = in_qam;
a_matrix = reshape(a_pad, M, []);	% S/P converter
A_matrix = ifft(a_matrix,M);
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix];	% Adding CP
r = sqrt(M)*reshape(A_matrix, [], 1);		% P/S converter
num_realizations = 50;
Ser = zeros(length(SNR_db_vect),num_realizations);

%% OFDM CHANNEL AND RECEIVER
for i=1:length(SNR_db_vect)
	snr_db = SNR_db_vect(i);
	snr_lin = 10^(snr_db/10);
	for k=1:num_realizations
		h = channel(22);
		filt = filter(h,1,r);
		E_tot = h * conj(h).';
		sigma_w = sigma_a/M * E_tot / snr_lin;
		noise_wgn = wgn(length(filt),1,10*log10(sigma_w),'complex');
		filt = filt + noise_wgn;
		b_matrix = reshape(filt, M+Npx, []);		% S/P converter
		%b_matrix = reshape(filt(1:end-mod(length(filt),M+Npx)), M+Npx, []);	% S/P converter
		rn = b_matrix(Npx+1:end,:);
		x_k = fft(rn,M);
		G = fft(h,M).';
		K_i = 1./G;
		y_matrix = x_k.*K_i;
		rec = reshape(y_matrix,[],1)/sqrt(M);		% P/S converter
		dec_bits = qamdemod(rec,16,'UnitAveragePower',true).';		% QAM demodulation
		Ser(i,k) = length(find(in(1:length(dec_bits))~=dec_bits))/length(dec_bits);
		% Ber(i,k) = Ser(i,k)/log2(16);			% Bit error rate is SER / log2(cardinality)
	end
end
Ser_mean = mean(Ser,2);		% Compute mean by row

%% NICE PLOTTING
figure();
semilogy(SNR_db_vect,Ser_mean,'Color','r');
title('SER versus SNR \Lambda');
grid on;
legend('SER for OFDM');
xlabel('SNR \Lambda'); ylabel('P_{s}');
xlim([SNR_db_vect(1) SNR_db_vect(end)]); ylim([10^-7 10^-0.3010]);