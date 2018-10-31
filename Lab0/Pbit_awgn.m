set(0,'DefaultTextInterpreter','latex');

clc; close all; clear global; clearvars;

N = 10000000;
in = rand(1,N);

in(in>0.5) = 1;
in(in<=0.5) = -1;

Pbit = zeros(1,15);

for i=1:15
	snr_lin = 10^(i/10);
	sigma_w = 2 / snr_lin;
	wgn = sqrt(sigma_w)*randn(1,N);
	in_after_noise = in+wgn;
	dec = zeros(1,N);
	dec(in_after_noise>0) = 1;
	dec(in_after_noise<=0) = -1;
	Pbit(i) = length(find(in(1:length(dec))~=dec)) / length(dec);
end

figure()
semilogy(1:15,Pbit);
title('Bit error probability');
grid on;
xlim([1 15]);
xlabel('SNR $\Lambda$'); ylabel('$P_{bit}$');