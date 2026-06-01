#!/usr/bin/env python3
"""
lrc_gabor_frame_bounds_s521.py   claudebox-2026-06-01-S521

The Gabor / time-frequency angle on LRC (reflection:
07-reflections/lrc-gabor-time-frequency-angle-s521.md).

Runners = frequencies v_i; safe band chi=1_[1/n,1-1/n] is the window; runner i sees
chi dilated by v_i. N(t)=sum_i 1_{B_i}(t) is the diagonal of the time-frequency
frame operator. Gabor FRAME BOUNDS = (A,B)=(min_t N, max_t N); LRC <=> A=0 (a hole
in the overcomplete system, density 2(n-1)/n>1). Bridge: B = congestion = clique
number of the danger graph. Doubling pair v_j=2v_i = time-frequency commensurability
(ambiguity peak, overlap ratio n/4).
"""
from fractions import Fraction as F
def fr(x): return x % 1
def dist(x):
    x = x % 1; return min(x, 1 - x)
def cells(sp, n):
    W = set([F(0)])
    for v in sp:
        for k in range(v):
            W.add(fr(F(k*n-1, n*v))); W.add(fr(F(k*n+1, n*v)))
    W = sorted(w for w in W if 0 <= w < 1); W2 = W + [F(1)]
    return W + [(a+b)/2 for a, b in zip(W, W2[1:])]
def N(sp, t, n): return sum(1 for v in sp if dist(F(v)*t) < F(1, n))
def congestion(sp, t, n):
    pts = sorted([F(0)] + [fr(F(v)*t) for v in sp]); k = len(pts); best = 0
    for i in range(k):
        c = 1
        for d in range(1, k):
            if (pts[(i+d) % k]-pts[i]) % 1 < F(1, n): c += 1
            else: break
        best = max(best, c)
    return best
def overlap(a, b, n):
    W = set([F(0), F(1)])
    for v in (a, b):
        for k in range(v):
            W.add(fr(F(k*n-1, n*v))); W.add(fr(F(k*n+1, n*v)))
    W = sorted(w for w in W if 0 <= w <= 1); tot = F(0)
    for x, y in zip(W, W[1:]):
        m = (x+y)/2
        if dist(F(a)*m) < F(1, n) and dist(F(b)*m) < F(1, n): tot += y-x
    return tot

def main():
    print("Gabor frame bounds (A=min N, B=max N) of the runner-window time-frequency system:")
    print(f"{'speeds':22} {'A=min N':8} {'B=max N':8} {'congestion(omega)':18} {'density 2(n-1)/n':16} {'LRC(A=0)':9}")
    for sp in [(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7,11)]:
        n = len(sp)+1; ts = cells(list(sp), n)
        Ns = [N(list(sp), t, n) for t in ts]; A = min(Ns); B = max(Ns)
        omega = max(congestion(list(sp), t, n) for t in ts)
        print(f"{str(sp):22} {A:8} {B:8} {omega:18} {float(2*(n-1)/n):16.3f} {str(A==0):9}")
    print("\n  LRC <=> lower frame bound A=0 (overcomplete density ~2 yet a time-frequency HOLE).")
    print("  Doubling pair = commensurable resonance (ambiguity peak):")
    for n in (5, 7, 9):
        print(f"   n={n}: overlap(1,2)*n^2/4 = {float(overlap(1,2,n)*n*n/4):.3f}  (=n/4)")

if __name__ == "__main__":
    main()
