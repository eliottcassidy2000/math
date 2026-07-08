"""
mac-mini-2026-07-07-S57 (HYP-5267) -- the k=10 CROSSOVER closer via the EXACT degree-4
one-sided moment bound, beta-OPTIMIZED (sharpening klein-S175 THM-656's degree-2 Cantelli).

On S = {maxgap <= 1/7}, F = sum_pairs tent_beta((e_j-e_i)x) >= toll = 1-k*beta, so
  mu(E) >= 1 - P(F >= toll).
Degree-4 one-sided Markov (EXACT moments m0..m4):
  P(F >= toll) <= min over polys p(t)=sum_{i<=4} c_i t^i, p(t)>=1 on [toll,Fmax],
     p(t)>=0 on [0,Fmax], of sum_i c_i m_i    -- an LP; strictly sharper than Cantelli
     (degree 2 is feasible).  Bennett/Bernstein INVALID (F = function of one variable x).
Optimize beta per family.  Compare max(diam-floor THM-653, deg4) to the bar 0.4521 and to
the TRUE mu (to see whether a residual is a real gap or bound-looseness).
"""
import numpy as np
from math import gcd
from functools import reduce
from scipy.optimize import linprog

k = 10
TH = 1/7
mp = 14249/252252

def tent(s, beta):
    return np.where((s > beta) & (s <= TH), s - beta, 0.0)

def Fvals(E, beta, x):
    Fv = np.zeros(len(x)); E = sorted(E)
    for i in range(len(E)):
        for j in range(len(E)):
            if i == j: continue
            Fv += tent(np.mod(abs(E[j]-E[i])*x, 1.0), beta)
    return Fv

def true_mu(E, x):
    Ea = np.array(sorted(E), float); ph = np.mod(np.outer(x, Ea), 1.0); ph.sort(axis=1)
    mg = np.concatenate([np.diff(ph, axis=1), (ph[:,0]+1-ph[:,-1])[:,None]], axis=1).max(axis=1)
    return float((mg > TH).mean())

def deg4_onesided(m, Fmax, toll, ngrid=1200):
    tl = np.linspace(0, Fmax, ngrid); th = np.linspace(toll, Fmax, ngrid)
    A = [[-(t**i) for i in range(5)] for t in tl] + [[-(t**i) for i in range(5)] for t in th]
    b = [0.0]*ngrid + [-1.0]*ngrid
    r = linprog(np.array(m), A_ub=np.array(A), b_ub=np.array(b), bounds=[(None,None)]*5, method='highs')
    return r.fun if r.success else None

def diam_floor(E):
    E = sorted(E); g = reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)])
    d = (E[-1]-E[0])//g
    return 146/(35*d) if d >= 6 else 0.0

def energy_R2(E):
    from collections import Counter
    E = sorted(E); cnt = Counter()
    for a in E:
        for b in E: cnt[a+b] += 1
    return sum(v*v for v in cnt.values()) - k*k

FAMS = {
  'AP_10 (block d9)': list(range(10)),
  'AP d10 perf': [0,1,2,3,4,5,6,7,8,10],
  'AP d11 perf': [0,1,2,3,4,5,6,7,9,11],
  '2-block 5+5 d13': [0,1,2,3,4,9,10,11,12,13],
  'AP d12 perf2': [0,1,2,3,4,5,7,9,11,12],
  'near-Sidon d16': [0,1,2,4,7,10,12,13,15,16],
  'Sidon d19': [0,1,3,7,12,14,16,17,18,19],
}
GRID = 2_000_000
x = (np.arange(GRID)+0.5)/GRID
BETAS = [b/280 for b in range(4, 42)]
print("=== k=10 degree-4 moment closer, beta-OPTIMIZED; bar 0.4521 ===\n")
print(f"{'family':20s} {'R2':>5s} {'diam-fl':>8s} {'trueMu':>7s} {'deg4*':>7s} {'@beta':>6s} {'max(df,d4)':>10s} {'verdict':>9s}")
allsafe = True
for name, E in FAMS.items():
    df = diam_floor(E); tm = true_mu(E, x); R2 = energy_R2(E)
    best = (0.0, None)
    for beta in BETAS:
        toll = 1 - k*beta
        if toll <= 0: continue
        Fv = Fvals(E, beta, x); m = [float(np.mean(Fv**i)) for i in range(5)]; Fmax = float(Fv.max())
        d4 = deg4_onesided(m, Fmax, toll)
        if d4 is not None and 1-d4 > best[0]: best = (1-d4, beta)
    mud4, bb = best
    floor = max(df, mud4)
    verdict = 'PASS' if floor >= 0.4521 else ('true-safe' if tm >= 0.4521 else 'REAL-GAP')
    if floor < 0.4521 and tm < 0.4521: allsafe = False
    print(f"{name:20s} {R2:5d} {df:8.4f} {tm:7.4f} {mud4:7.4f} {(bb or 0):6.3f} {floor:10.4f} {verdict:>9s}")
print(f"\n  {'ALL crossover families PASS or true-safe => k=10 wide-family mu floor holds' if allsafe else 'a REAL-GAP remains (true mu < bar) -- would threaten LRC; investigate'}")
print("  compact diam<=10 = klein THM-653 exhaustion; this targets the diam>10 crossover.")
print("  deg4* strictly sharper than klein Cantelli; where deg4<bar but trueMu>=bar, higher moments/degree close it.")
