clc; close all; clear global; clearvars;

ch_length = 21;
h = zeros(1,ch_length);
alpha = (1-exp(-1/3))/(1-exp(-22/3));

for k=0:ch_length
	sigma = alpha * exp(-k/3);
	tap = sqrt(sigma/2)*randn(1,1) + 1i*sqrt(sigma/2)*randn(1,1);
	h(k+1) = tap;
end

Npx = length(h);
M = 2^(ceil(log2(5*Npx)));			% n. of subcarriers
in = randi(15,1,(M+Npx)*2^8);		% Input generation

save('Input.mat','M','in','Npx','h');