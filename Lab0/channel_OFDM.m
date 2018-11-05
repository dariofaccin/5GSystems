clc; close all; clear global; clearvars;

load('Input.mat');
in_qam = qammod(in,16,'UnitAveragePower',true).';		% QAM modulation
SNR_db_vect = 10;					% SNR vector (in dB)

% a_pad = [in_qam; ones(M - mod(length(in_qam), M), 1) * (1+1i)];
a_pad = in_qam;
a_matrix = reshape(a_pad, M, []);	% S/P converter

A_matrix = ifft(a_matrix,M);
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix];	% Adding CP

r = reshape(A_matrix, [], 1);		% P/S converter

in_filt = filter(h,1,r);	% Signal after passing the channel

sigma_a = var(in_qam);				% Random number here, do computations
h = h ./ sqrt((h * conj(h).'));
E_tot = h * conj(h).';

Ser = zeros(length(SNR_db_vect),1);
Ber = zeros(length(SNR_db_vect),1);

for i=1:length(SNR_db_vect)
	snr_db = SNR_db_vect(i);
	snr_lin = 10^(snr_db/10);
	sigma_w = sigma_a/M * E_tot / snr_lin;
	samples = in_filt(1:end);
	noise_wgn = wgn(length(samples),1,10*log10(sigma_w),'complex');
	samples = samples + noise_wgn;
	b_matrix = reshape(samples(1:end-mod(length(samples),M+Npx)), M+Npx, []);	% S/P converter
	rn = b_matrix(Npx+1:end,:);
	x_k = fft(rn,M);
	G = fft(h,M).';
	K_i = 1./G;
	y_matrix = x_k.*K_i;
	rec = reshape(y_matrix,[],1);		% P/S converter
	dec_bits = qamdemod(rec,16,'UnitAveragePower',true).';		% QAM demodulation
	Ser(i) = length(find(in(1:length(dec_bits))~=dec_bits))/length(dec_bits);
	Ber(i) = Ser(i)/log2(16);			% Bit error rate is SER / log2(cardinality)
end

% figure();
% semilogy(SNR_db_vect,Ser);
% title('P_{symbol} versus SNR \Lambda');
% grid on;
% xlabel('SNR \Lambda'); ylabel('P_{bit}');
% xlim([SNR_db_vect(1) SNR_db_vect(10)]);% ylim([10^-2 10^-0.3010]);