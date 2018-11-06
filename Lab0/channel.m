function h = channel(ch_length)
% Generates channel with exponentially decaying power profile
h = zeros(1,ch_length);
alpha = (1-exp(-1/3))/(1-exp(-22/3));

for i=1:ch_length
	sigma = alpha * exp(-(i-1)/3);
	h(i) = sqrt(sigma/2)*randn(1,1) + 1i*sqrt(sigma/2)*randn(1,1);
end
h = h ./ sqrt((h * conj(h).'));
end