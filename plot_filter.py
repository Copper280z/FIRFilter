from matplotlib import pyplot as plt
import numpy as np

fname = "fir_filter.txt"

coeffs = np.loadtxt(fname)
sample_rate = coeffs[0]
coeffs = coeffs[1:]

fig, axs = plt.subplots(3)
axs[0].set_title("Filter kernel")
axs[0].plot(coeffs)

resp = np.fft.rfft(coeffs)
freqs = np.fft.rfftfreq(len(coeffs),1/sample_rate)

axs[1].set_title("Filter Frequency response")
axs[1].semilogx(freqs, 20*np.log10(np.abs(resp)))

axs[2].set_title("Filter Phase response")
axs[2].semilogx(freqs, np.angle(resp))
plt.tight_layout()
plt.show()