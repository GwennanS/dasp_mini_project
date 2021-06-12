%% Set up sounds
[s, fs] = audioread('data/clean_speech.wav');
[n1, fsn1] = audioread('data/babble_noise.wav');
[n2, fsn2] = audioread('data/aritificial_nonstat_noise.wav');
[n3, fsn3] = audioread('data/Speech_shaped_noise.wav');

n1(numel(s)) = 0;
n2(numel(s)) = 0;
n3(numel(s)) = 0;
n1 = n1(1 : numel(s))*0.05;
n2 = n2(1 : numel(s))*0.05;
n3 = n3(1 : numel(s))*0.05;
n = n2+n3;
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
    FRAMESPSD(i, :) = (1/L) * (FRAMESMAG(i, :)).^2;
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

%% Find noise in a way???
FRAMESNOISE = movmin((FRAMESbartlett), 20, 2);
%filtered = FRAMESPSD - FRAMESNOISE;

%% Calculate SNR
priori_speech = zeros(1, frame_size*2); 
gain = zeros(size(FRAMES));
alp = 0.96;
for i = 1:l
    first_part = (priori_speech.^2) ./ (FRAMESNOISE(i, :));
    second_part = max((FRAMESPSD(i, :)) ./ (FRAMESNOISE(i, :)) - 1, 0);
    SNR = (alp * first_part) + (1-alp) * second_part;
    gain(i, :) =  SNR ./ (SNR + 1) ;
    priori_speech = gain(i, :) .* FRAMESMAG(i, :);
end

%% Estimate Clean speech
CLEAN = zeros(size(y));
for i = 1:l
    new_mag = gain(i, :) .* FRAMESMAG(i, :) ;
    clean = new_mag .* exp(1j*FRAMESPHA(i, :));
    clean = real(ifft(clean));
    CLEAN(i*frame_size+1:(i+2)*frame_size) = CLEAN(i*frame_size+1:(i+2)*frame_size) + clean.';
end

%sound(CLEAN, fs)
