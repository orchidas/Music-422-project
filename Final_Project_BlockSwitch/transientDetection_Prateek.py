import numpy as np 

def transient_detection(data):

	N_half = len(data)/2

	fftData = np.fft.fft(data)[ : N_half]

	threshold_energy = 0.14

	weights = np.ones(N_half)

	weights[ : len(weights)/2] = 0.0

	# print weights

	E = np.sum(weights*abs(fftData))/ (N_half)

	return E > threshold_energy



# x = np.ones(20)

# print transient_detection(x)