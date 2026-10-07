import cmath, math
L3 = math.log2(3); al = L3 - 1
# root w of e^{2 pi i th + al w} + e^{-w} = 2, continued from th=0 (w=0)
def F(w, th): return cmath.exp(2j*math.pi*th + al*w) + cmath.exp(-w) - 2
def dF(w, th): return al*cmath.exp(2j*math.pi*th + al*w) - cmath.exp(-w)
w = 0j; out = {}
ths = [i*1e-5 for i in range(0, 8001)]
for th in ths:
    w = w + 2j*math.pi*1e-5/(1-al)  # predictor
    for _ in range(50):
        d = F(w, th)/dF(w, th); w -= d
        if abs(d) < 1e-15: break
    out[round(th, 5)] = w
for th in [0.003, 0.006, 0.0105, 0.0135, 0.0165, 0.0196, 0.0226, 0.0361, 0.0391, 0.0556, 0.0587, 0.075]:
    w = out[round(th,5)]
    print(f"theta={th:.4f}  mu={w.real:+.4f}  nu_shift={w.imag/(2*math.pi):+.4f}  effective C=-mu/theta^2={-w.real/th**2:7.1f}  gauss 552.7*th^2={552.7*th*th:.4f}")
# scan for any other roots with -0.6<mu<0 in nu-window [f-0.5, f+0.5] for f=12,17,24
for f in [12, 17, 24, 36, 53]:
    th = f*L3 - round(f*L3)
    roots = set()
    for mu0 in [x*0.05 for x in range(-20, 2)]:
        for nu0 in [x*0.05 for x in range(-10, 11)]:
            w = complex(mu0, 2*math.pi*nu0)
            ok = False
            for _ in range(100):
                d = F(w, th)/dF(w, th); w -= d
                if abs(d) < 1e-13: ok = True; break
            if ok and abs(w.imag/(2*math.pi)) <= 0.5 and -1.5 < w.real < 0.5:
                roots.add((round(w.real, 4), round(w.imag/(2*math.pi), 4)))
    print(f, sorted(roots))
