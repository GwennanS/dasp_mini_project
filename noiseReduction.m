function res = noiseReduction(l, parameters)

estimateNoise = noisePSDestimator(l, parameters)
estimateSpeech = speechPSDestimator(l, parameters)

res = gainfunction(l, estimateNoise, estimateSpeech, parameters)

res = l;
end