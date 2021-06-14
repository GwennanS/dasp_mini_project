%% Set up sounds
[s, fs] = audioread('data/clean_speech.wav');
[n1, fsn1] = audioread('data/babble_noise.wav');
[n2, fsn2] = audioread('data/aritificial_nonstat_noise.wav');
[n3, fsn3] = audioread('data/Speech_shaped_noise.wav');

n1(numel(s)) = 0;
n2(numel(s)) = 0;
n3(numel(s)) = 0;
n1 = n1(1 : numel(s))*0.05;
n2 = n2(1 : numel(s))*0.35;
n3 = n3(1 : numel(s))*0.35;
n = n2+n3+n1;
y = s+n;
%y = y(:,1);

%% Set frames
N = numel(y);
frame_size = 0.01 * fs; 
hann = hanning(frame_size*2, 'periodic');
l = floor(N/frame_size) -2;
FRAMES = zeros(l, frame_size*2);
for i = 1:l
    frame = y(i*frame_size-0.5*frame_size+1 : (i+1)*frame_size+0.5*frame_size);
    FRAMES(i, :) = hann .* frame;
end

%% Find FFT, mag, phase and estimate PSD for each frame
L = frame_size;

FRAMESFFT = zeros(l, frame_size*2);
FRAMESMAG = zeros(l, frame_size*2);
FRAMESPHA = zeros(l, frame_size*2);
FRAMESPSD = zeros(l, frame_size*2);
for i = 1:l
    framefft = fft(FRAMES(i, :));
    FRAMESFFT(i, :) = framefft;
    FRAMESMAG(i, :) = abs(framefft);
    FRAMESPHA(i, :) = angle(framefft);
    FRAMESPSD(i, :) = (FRAMESMAG(i, :)).^2;
end

%% Calculate bartlett for each frame
M = 3;
FRAMESbartlett = zeros(l, frame_size*2);
for i = 1:l
    sumFRAMES = FRAMESPSD(i, :);
    for ii = i-M+1:i
        if ii < 1
            ii = 1;
        end
        sumFRAMES = sumFRAMES + FRAMESPSD(ii, :);
    end
    FRAMESbartlett(i, :) = 1/M * sumFRAMES;
end

expSmoother = FRAMESPSD;
alph = 0.98;
for i = 2:l
    expSmoother(i, :) = alph * expSmoother(i-1, :) + (1-alph) * expSmoother(i, :);
end
%% Find noise in a way???
FRAMESNOISE = zeros(l, frame_size*2);
%FRAMESNOISE = movmin((expSmoother), 10, 2);
%filtered = FRAMESPSD - FRAMESNOISE;

sumNOISE = zeros(1, frame_size*2);
for i = 1:10 
    sumNOISE = sumNOISE + FRAMESPSD(i, :);
end


vadNOISE = sumNOISE * 1/10;
for i = 1:l
    FRAMESNOISE(i, :) = vadNOISE;
end


%% power spectral subtraction
pss = max((1-(FRAMESNOISE  ./ FRAMESPSD)), 0.2) .^0.5;

%% Calculate SNR
priori_speech = zeros(1, frame_size*2); 
weiner = zeros(l, frame_size*2);
alp = 0.96;
for i = 1:l
    first_part = (priori_speech.^2) ./ (FRAMESNOISE(i, :));
    second_part = max((FRAMESPSD(i, :)) ./ (FRAMESNOISE(i, :)) - 1, 0);
    SNR = (alp * first_part) + (1-alp) * second_part;
    weiner(i, :) =  SNR ./ (SNR + 1) ;
    priori_speech = weiner(i, :) .* FRAMESMAG(i, :);
end

%% Estimate Clean speech

gain = weiner;
CLEAN = zeros(size(y));
NOISE = zeros(size(y));
for i = 1:l
    noise = real(ifft(sqrt(FRAMESNOISE(i, :)) .*exp(1j*FRAMESPHA(i, :))));
    new_mag = gain(i, :) .* FRAMESMAG(i, :) ;
    clean = new_mag .* exp(1j*FRAMESPHA(i, :));
    clean = real(ifft(clean));
    CLEAN(i*frame_size+1:(i+2)*frame_size) = CLEAN(i*frame_size+1:(i+2)*frame_size) + clean.';
    NOISE(i*frame_size+1:(i+2)*frame_size) = NOISE(i*frame_size+1:(i+2)*frame_size) + noise.';
end

%% Show it
sound(CLEAN(1:100000), fs)

figure;
plot(y)
figure;
plot(CLEAN)
%figure;
%plot(noise)

%% Evaluation
%figure;
%plot(CLEAN.^2-s.^2)

snr1 = 10*log10(mean(s.^2)/mean(n.^2));
snr2 = 10*log10(mean(CLEAN.^2)/mean((CLEAN-s).^2));
snr2 = 10*log10(mean(CLEAN.^2)/mean((n).^2));

frame_size_snr = 0.04 * fs;
num_snr = floor(N/frame_size_snr) -1;
sum_snr = 0;
for i = 0:num_snr
    first_part_snr = sum(abs(y(i*frame_size_snr + 1 : (i+1)*frame_size_snr )), 'all')^2;
    second_part_snr = sum(abs(y(i*frame_size_snr +1 : (i+1)*frame_size_snr)-CLEAN(i*frame_size_snr +1 : (i+1)*frame_size_snr)))^2;
    sum_snr = sum_snr + 10*log10(first_part_snr/second_part_snr);
end
sum_snr = 1/num_snr * sum_snr;

stoiy = stoi(s, y, fs)
stoic = stoi(s, CLEAN, fs)
stoin = stoi(s, NOISE, fs)



