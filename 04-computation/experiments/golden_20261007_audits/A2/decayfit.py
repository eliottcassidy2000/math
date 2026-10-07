"""Fit per-octave decay of log-Fourier modes (k=16..28) and compare with the exact root of
   e^{2 pi i theta + w log2(3/2)} + e^{-w} = 2 (linear Haar transfer psi(u)=psi(u+log2(3/2))/2+psi(u-1)/2)
   and with the Gaussian approximation 2 pi^2 * 27.97 * theta^2."""
import re, math, cmath
L3 = math.log2(3); al = L3 - 1
def root(th):
    w = 0j; steps = 4000
    for i in range(1, steps + 1):
        t = th * i / steps
        w += 2j * math.pi * (th / steps) / (1 - al)
        for _ in range(30):
            F = cmath.exp(2j*math.pi*t + al*w) + cmath.exp(-w) - 2
            dF = al*cmath.exp(2j*math.pi*t + al*w) - cmath.exp(-w)
            d = F/dF; w -= d
            if abs(d) < 1e-15: break
    return w
def parse(files):
    amp = {}
    for fn in files:
        txt = open(fn).read()
        for blk in re.finditer(r"k=(\d+) mode=\d.*?\n  list:(.*?)\n", txt):
            k = int(blk.group(1))
            for f, a in re.findall(r" (\d+):([0-9.]+)", blk.group(2)):
                amp.setdefault(int(f), {})[k] = float(a)
    return amp
for name, files in [("-1 basin", ["spec_neg24.out", "spec_neg28.out"]), ("entry via 85", ["spec_ent24.out", "spec_ent28.out"])]:
    amp = parse(files)
    print(f"== {name}: per-octave decay rate, least squares on k=16..28 vs exact root vs Gaussian 552.7 th^2")
    for f in sorted(amp):
        ks = sorted(amp[f]); ys = [math.log(max(amp[f][k], 1e-9)) for k in ks]
        if min(amp[f][k] for k in ks) < 0.0015: continue
        n = len(ks); mk = sum(ks)/n; my = sum(ys)/n
        slope = sum((k-mk)*(y-my) for k, y in zip(ks, ys)) / sum((k-mk)**2 for k in ks)
        th = f*L3 - round(f*L3); w = root(th)
        print(f"  f={f:4d} theta={th:+.5f}  observed {-slope:.4f}  exact-root {-w.real:.4f} (freq shift {w.imag/2/math.pi:+.4f})  Gaussian {552.7*th*th:.4f}   amp k16->k28: {amp[f][16]:.4f}->{amp[f][28]:.4f}")
