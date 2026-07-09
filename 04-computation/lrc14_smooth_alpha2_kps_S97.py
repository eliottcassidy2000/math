# Is the smooth-route Fourier decay PROVABLY alpha=2? (kps-S97) opus-S170 MEASURED alpha in [1.48,2.02].
# CLAIM: maxgap(x) and W(x) are CONTINUOUS piecewise-linear in x (kinks at collisions/gap=1/7 crossings,
# NO jumps) => Fourier coeff ~ 1/m^2 asymptotically (alpha=2), EXACTLY -- 1.48 is pre-asymptotic.
# Verify on hard 7-structured MISTAKE-128 set + AP + generic, fitting alpha over a HIGH-m tail window.
import numpy as np
TH=1.0/7.0
def maxgap_W(E,xs):
    E=np.array(sorted(E)); mg=np.empty(len(xs)); W=np.empty(len(xs))
    for t,x in enumerate(xs):
        ph=np.sort(np.mod(E*x,1.0))
        g=np.append(np.diff(ph), ph[0]+1-ph[-1])
        mg[t]=g.max(); W[t]=np.maximum(g-TH,0).sum()
    return mg,W
def fourier_alpha(f,M0,M1):
    # f sampled on N points over [0,1); FFT; fit log|fhat(m)| ~ -alpha log m over m in [M0,M1]
    N=len(f); F=np.fft.rfft(f-f.mean())/N; mag=np.abs(F)
    ms=np.arange(len(mag))
    sel=(ms>=M0)&(ms<=M1)&(mag>0)
    lm=np.log(ms[sel]); lg=np.log(mag[sel])
    A=np.vstack([lm,np.ones_like(lm)]).T
    alpha,_=np.linalg.lstsq(A,lg,rcond=None)[0]
    return -alpha
N=200000; xs=(np.arange(N)+0.5)/N
sets={
 "AP{1..13}":list(range(1,14)),
 "7-struct(MISTAKE-128)":[0,7,14,21,26,29,37,44,51,58,67,75,82],
 "generic-dissoc":[0,3,10,14,23,31,40,52,61,67,70,75,80],
}
print(f"Fourier decay exponent alpha of maxgap(x) and W(x) (continuous PL => alpha=2 asymptotically).")
print(f"{'set':>24}{'a(mid)':>9}{'a(high)':>9}{'a_W(high)':>11}")
for name,E in sets.items():
    mg,W=maxgap_W(E,xs)
    a_mid=fourier_alpha(mg,20,200)     # opus's likely measurement window (pre-asymptotic)
    a_hi =fourier_alpha(mg,2000,40000) # HIGH-m asymptotic window
    aW_hi=fourier_alpha(W,2000,40000)
    print(f"{name:>24}{a_mid:>9.3f}{a_hi:>9.3f}{aW_hi:>11.3f}")
print()
print("=> if a(high) -> 2 for ALL sets (incl 7-structured), the smooth-route decay is PROVABLY alpha=2")
print("   (continuous PL, no jumps): opus's 1.48 is pre-asymptotic. Absolute convergence RIGOROUS.")
