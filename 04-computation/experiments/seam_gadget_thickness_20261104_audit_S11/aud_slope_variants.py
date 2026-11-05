"""aud_slope_variants.py -- audit S11: the after-the-fact statistic "D_8 averaged over the realised slopes" (claimed 0.8825).
Inputs (regenerate with: AUD_DUMP=band_states8.txt ./aud_bfs plane 5000 aud_plane_seeds.txt aud_mats1..8.txt, run in DIR):
DIR/band_states8.txt (order<=8 plane BFS states with 801<=H<=1799) and DIR/aud_plane_maps_dump.txt (the 21,023 maps).
Usage: python3 aud_slope_variants.py DIR"""
import sys, os
os.chdir(sys.argv[1])
import numpy as np
S = np.loadtxt('band_states8.txt', dtype=np.int64); A = np.loadtxt('aud_plane_maps_dump.txt', dtype=np.int64)
a = A[:, 0].astype(float); b = A[:, 1].astype(float)
def D(t):
    t = np.atleast_1d(t).astype(float); out = np.empty(len(t))
    for i in range(0, len(t), 2000):
        tt = t[i:i+2000][:, None]; out[i:i+2000] = (1.0 / (a[None, :] + b[None, :] * tt) ** 2).sum(axis=1)
    return out
H = S[:, 0]; R = S[:, 1]; t = R / H
Dt = D(t)
print('(a) mean D_8 over all 54994 states with 801<=H<=1799:', round(Dt.mean(), 4))
u = np.unique(t); print('(b) mean over distinct slopes (%d):' % len(u), round(D(u).mean(), 4))
m = t <= 2.75; print('(c) states with slope <= 2.75 (%d):' % m.sum(), round(Dt[m].mean(), 4))
hist = {2: 36235, 3: 12363, 4: 5221, 5: 1150, 6: 25}
print('(d) histogram bin centres:', round(sum(c * D(k + 0.5)[0] for k, c in hist.items()) / sum(hist.values()), 4))
print('(e) D at median slope %.4f: %.4f; at mean slope %.4f: %.4f' % (np.median(t), D(np.median(t))[0], t.mean(), D(t.mean())[0]))
for w in (0.01, 0.05, 0.1):
    bc = (np.floor(t / w) + 0.5) * w; print('(f) slopes binned at %.2f, centre: %.4f' % (w, D(bc).mean()))
tr = np.round(t, 2); print('(g) slopes rounded to 2 decimals:', round(D(tr).mean(), 4))
