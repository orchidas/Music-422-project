
[input,fs]  = audioread('inputs/test.wav');

inputL = fft(input(:,1));

inputR = fft(input(:,2));


coherence = abs(inputL.*inputR).^2 ./ (abs(inputL).^2 .* abs(inputR).^2);

plot(coherence)
