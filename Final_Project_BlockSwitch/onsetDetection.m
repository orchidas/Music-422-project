

% Method-1 Transient Detection 

[x,fs] = audioread('inputs/german.wav');

blockSize = 1024
numOfBlocks = length(x(:,1)/2);

sig = x(:,1);

level = 0.2;
blocks = 0;

number_onsets = 0
while blocks <= numOfBlocks

start_pos = 1 + blocks*blockSize
end_pos = start_pos + blockSize
if(end_pos >= length(sig))
    
    break
end
signal = sig(start_pos : end_pos);

threshold_energy = 0.14
X = fft(signal);
N_half = floor(length(signal)/2);

X_half = X(1:N_half);

weights = zeros(1,N_half)';

weights(floor(length(weights))/2 : end) = 1.0;

E = sum(weights.*abs(X_half))/(N_half);



% plot(signal(:,1))

if( E > threshold_energy)
    
    disp('Onset Found')
    number_onsets = number_onsets + 1;
    yL = get(gca,'YLim');
    pos = (start_pos + end_pos)/2
    line([pos pos],yL,'Color','r');
    hold on
    
    
end

blocks = blocks + 1;


end

disp("Number of Onsets Detected : " + number_onsets)
plot(x(:,1))
hold off
