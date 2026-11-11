import numpy as np
r10 = np.load("r10a.npy"); r9 = np.load("line_counts.npz")["r9a"]
Ns = np.arange(1, 1 << 15, 2)
for k in range(6, 15):
    idx = np.arange((1 << k) // 2, (1 << (k + 1)) // 2)
    print("[2^%d,2^%d): share r_9 %.5f  share r_10 %.5f  mean r_9 %.2f mean r_10 %.2f" % (k, k + 1, np.mean(r9[idx] / (Ns[idx] / 6)), np.mean(r10[idx] / (Ns[idx] / 6)), r9[idx].mean(), r10[idx].mean()))
