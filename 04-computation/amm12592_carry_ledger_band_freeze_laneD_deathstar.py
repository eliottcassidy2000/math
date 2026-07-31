#!/usr/bin/env python3
"""Lane D referee: carry-chain ledger for gamma < 1 (AMM 12592 / HYP-9061).

Frame (THM-2966 spine normal form, cited not re-proved): an extractor with
deadline T(m) = m + 1 + d_m, d_m = floor(gamma*m) + D0, exists iff doubled
deficits delta_{m,side,k} (integers, |delta| <= binom(d_m,k), parity
delta == binom(d_m,k) mod 2) can be chosen so that

    D_M(p) = (1/2) * sum_{m<=M} sum_cells delta * p^z q^o  ->  0

pointwise on (0,1), where the 0-side cell (m,k) has monomial
z = m + d_m - k, o = k + 1 and the 1-side cell has z = k + 1,
o = m + d_m - k.  All cells of row m lie on the anti-diagonal
z + o = A_m := m + d_m + 1 = T(m).

This script:
  A. PROVED geometry checks (exact, all m <= 400, several gamma/D0):
     cone coverage o in [1, d_m+1] U [m, A_m - 1]; band = [d_m+2, m-1]
     nonempty iff m >= d_m + 3; cross-side same-monomial pairing exists
     iff m <= d_m + 1 ("pairing horizon"); the corner monomial
     (z,o) = (m, d_m+1) is attained by NO other cell of any row once
     m >= d_m + 2 (corner isolation), and binom(d_m,d_m) = 1 is odd, so
     every row emits an isolated forced half-token at its corner.
  B. PROVED forced-parity closed form, verified against the simulation at
     every step: the homogenized doubled ledger's parity vector is
     choice-independent and equals
        beta_M(x) = (1+x)^{A_M} + (1 + x^{M+1}) (1+x)^{d_M}   mod 2.
     (Derivation: 2 D_M = 2 S_M - 1 + p^{M+1} + q^{M+1}, and 2 S_M has
     even integer coefficients; homogenize to degree A_M.)
  C. Exact greedy carry-chain simulation (Task 2): gamma in
     {1/2, 3/5, 3/4, 9/10, 19/20, 99/100, 1}, exact integer/Fraction
     arithmetic, M up to 150 rows (>= 30 shells), with referee asserts:
     box + parity at every cell, and independent recomputation of D_M(p)
     from the raw choice list vs the descended ledger (must agree exactly).
  D. Local sign-search upgrade on one failing run, with held-out biases.
  E. Realization cross-check: for gamma = 1 (control), rebuild the actual
     stopping rule from the chosen w's, enumerate all stopped words, and
     verify the polynomial identity  P_H(p) = 1/2 - (p^{M+1}+q^{M+1})/2
     + D_M(p) exactly, plus the pathwise deadline.
  F. Certificate (27) re-verification (exact) and slope reverse-engineering
     numerics (reported, not asserted).
"""

from fractions import Fraction as Fr
from math import comb, log
import sys

def require(c, msg):
    if not c:
        raise RuntimeError(msg)

# ---------------------------------------------------------------- A. geometry
def geometry_audit():
    for gn, gd, D0 in [(1,2,0),(1,2,6),(3,5,0),(3,5,6),(3,4,0),(9,10,0),
                       (19,20,0),(2457,6592,0),(2457,6592,10)]:
        d = lambda m: (gn*m)//gd + D0
        A = lambda m: m + d(m) + 1
        # anti-diagonals strictly increase -> cells of different rows never share
        # a monomial; within row m the 0-side corner (m, d+1) coincides with a
        # 1-side cell iff d_m + 1 >= m.
        for m in range(1, 401):
            require(A(m+1) > A(m), "A not strictly increasing")
            zs0 = set(range(1, d(m)+2))                 # 0-side o-positions
            zs1 = set(range(m, A(m)))                   # 1-side o-positions
            require(zs0 == {k+1 for k in range(d(m)+1)}, "cone0 drift")
            require(zs1 == {m + d(m) - k for k in range(d(m)+1)}, "cone1 drift")
            band = set(range(d(m)+2, m))
            require((len(band) > 0) == (m >= d(m)+3), "band birth drift")
            corner_shared = (d(m)+1) in zs1
            require(corner_shared == (m <= d(m)+1), "pairing horizon drift")
            require(comb(d(m), d(m)) % 2 == 1, "corner parity drift")
    print("A_geometry=PROVED_checked_m<=400_9_depth_laws"
          " (corner_isolated_iff_m>=d_m+2; band_birth_iff_m>=d_m+3)")

# ------------------------------------------------------- B/C. ledger machine
class Ledger:
    def __init__(self, gn, gd, D0):
        self.gn, self.gd, self.D0 = gn, gd, D0
        self.b = {}        # doubled homogenized coeffs on current anti-diag
        self.A = None
        self.rows = []     # (m, A_m, [(o, cap, par, delta)])
    def d(self, m):  return (self.gn*m)//self.gd + self.D0
    def Ad(self, m): return m + self.d(m) + 1
    def cells(self, m):
        dm, A = self.d(m), self.Ad(m)
        cel = {}
        for k in range(dm+1):
            c = comb(dm, k)
            for o in (k+1, A-1-k):
                cap, par = cel.get(o, (0, 0))
                cel[o] = (cap + c, (par + c) % 2)
        return cel
    def descend(self, t):
        if t == 0 or not self.b:
            return
        row = [comb(t, i) for i in range(t+1)]
        nb = {}
        for o, c in self.b.items():
            for i, r in enumerate(row):
                nb[o+i] = nb.get(o+i, 0) + c*r
        self.b = {o: c for o, c in nb.items() if c}
    def step(self, m, policy):
        Anew = self.Ad(m)
        if self.A is not None:
            self.descend(Anew - self.A)
        self.A = Anew
        rec = []
        for o, (cap, par) in sorted(self.cells(m).items()):
            delta = policy(self, o, cap, par)
            require(abs(delta) <= cap and (delta - par) % 2 == 0,
                    "box/parity violation")
            nv = self.b.get(o, 0) + delta
            if nv: self.b[o] = nv
            else:  self.b.pop(o, None)
            rec.append((o, cap, par, delta))
        self.rows.append((m, Anew, rec))
    # exact value of D_M at rational p, from the descended ledger
    def val(self, p):
        q = 1 - p
        return sum(c * p**(self.A - o) * q**o for o, c in self.b.items()) / 2
    # independent recomputation from the raw choices (referee identity)
    def val_raw(self, p):
        q = 1 - p
        tot = Fr(0)
        for m, Am, rec in self.rows:
            for o, cap, par, delta in rec:
                if delta:
                    tot += Fr(delta) * p**(Am - o) * q**o
        return tot / 2
    def parity_formula(self):
        # beta_M = (1+x)^{A} + (1+x^{M+1})(1+x)^{d_M} mod 2, M = last row
        M = self.rows[-1][0]
        dM, A = self.d(M), self.A
        par = {}
        for o in range(A+1):
            v = 1 if (o & A) == o else 0
            if o <= dM and (o & dM) == o:
                v ^= 1
            j = o - (M+1)
            if 0 <= j <= dM and (j & dM) == j:
                v ^= 1
            if v: par[o] = 1
        return par

def greedy(led, o, cap, par):
    t = -led.b.get(o, 0)
    t = max(-cap, min(cap, t))
    if (t - par) % 2:
        b0 = led.b.get(o, 0)
        cand = [x for x in (t-1, t+1) if -cap <= x <= cap]
        t = min(cand, key=lambda x: (abs(b0+x), abs(x)))
    return t

def lazy(led, o, cap, par):   # hostile control: never absorbs backlog
    return par

def run(gn, gd, D0, M, policy, plist, checkpoints, tag, band_stats=True):
    led = Ledger(gn, gd, D0)
    print(f"--- run {tag}: gamma={gn}/{gd} D0={D0} M={M} "
          f"policy={policy.__name__}")
    hist = {}
    for m in range(1, M+1):
        led.step(m, policy)
        # forced-parity closed form must hold at every step (claim B)
        got = {o: c % 2 for o, c in led.b.items() if c % 2}
        require(got == led.parity_formula(), f"parity closed form fails m={m}")
        if m % 10 == 0:
            require(led.val(Fr(1,3)) == led.val_raw(Fr(1,3)),
                    "raw-vs-descended value drift")
        if m in checkpoints:
            row = []
            for p in plist:
                v = led.val(p)
                row.append(v)
            hist[m] = row
            dm = led.d(m)
            band = {o: c for o, c in led.b.items() if dm+2 <= o <= m-1}
            bmax = max((abs(c) for c in band.values()), default=0)
            mval = " ".join(f"{float(v):+.3e}" for v in row)
            print(f"  m={m:4d} A={led.A:4d} band[{dm+2},{m-1}] "
                  f"maxband={bmax} |D_M(p)|: {mval}")
    return led, hist

# ------------------------------------------------- D. local sign-flip search
def local_search(led, fit_ps, hold_ps, m_hi, cap_hi):
    # candidate moves: change delta at small-capacity cells by +-2 (stays in
    # parity class), respecting the box; objective = sum of squares of D(p)
    # over fit_ps (floats); exact re-eval at the end.
    moves = []
    for ri, (m, Am, rec) in enumerate(led.rows):
        if m > m_hi: continue
        for ci, (o, cap, par, delta) in enumerate(rec):
            if cap <= cap_hi:
                moves.append((ri, ci, Am, o))
    def massf(p, Am, o):
        return float(p)**(Am-o) * float(1-p)**o / 2
    resid = {p: float(led.val_raw(p)) for p in list(fit_ps) + list(hold_ps)}
    deltas = {}
    improved, rounds = True, 0
    while improved and rounds < 400:
        improved = False; rounds += 1
        best = None
        for mv in moves:
            ri, ci, Am, o = mv
            o_, cap, par, delta = led.rows[ri][2][ci]
            cur = delta + deltas.get(mv, 0)
            for e in (-2, 2):
                if abs(cur + e) > cap: continue
                obj = sum((resid[p] + e*massf(p, Am, o))**2 for p in fit_ps)
                if best is None or obj < best[0]:
                    best = (obj, mv, e)
        base = sum(resid[p]**2 for p in fit_ps)
        if best and best[0] < base * (1 - 1e-12):
            obj, mv, e = best
            deltas[mv] = deltas.get(mv, 0) + e
            ri, ci, Am, o = mv
            for p in resid:
                resid[p] += e * massf(p, Am, o)
            improved = True
    # apply exactly and re-report
    for (ri, ci, Am, o), e in deltas.items():
        o_, cap, par, delta = led.rows[ri][2][ci]
        led.rows[ri][2][ci] = (o_, cap, par, delta + e)
        require(abs(delta+e) <= cap and (delta+e-par) % 2 == 0, "ls box")
    return deltas, rounds

# ------------------------------------------------- E. realization control
def realization_control(M=4):
    # gamma = 1, D0 = 0: d_m = m, T(m) = 2m+1. Build the scheme from greedy
    # choices, enumerate every stopped word, and check the exact identity.
    led = Ledger(1, 1, 0)
    for m in range(1, M+1):
        led.step(m, greedy)
    # decompose merged position-deltas into per-side cell values
    # 0-side cell (m,k): o = k+1, cap comb(d,k); 1-side (m,k'): o = A-1-k'.
    schemes = {}   # (m, side) -> list over k of w_{m,k}
    for m, Am, rec in led.rows:
        dm = m
        caps0 = {k+1: comb(dm, k) for k in range(dm+1)}
        caps1 = {Am-1-k: comb(dm, k) for k in range(dm+1)}
        w0 = [0]*(dm+1); w1 = [0]*(dm+1)
        for o, cap, par, delta in rec:
            c1 = caps0.get(o, 0); c2 = caps1.get(o, 0)
            require(c1 + c2 == cap, "cap decomposition drift")
            found = False
            for d1 in range(-c1, c1+1):
                if (d1 - c1) % 2: continue
                d2 = delta - d1
                if abs(d2) <= c2 and (d2 - c2) % 2 == 0:
                    found = True; break
            require(found, "delta split failed")
            if c1: w0[o-1] = (c1 + d1)//2
            if c2: w1[Am-1-o] = (c2 + d2)//2
        schemes[m] = (w0, w1)
    # enumerate stopped words; polynomial arithmetic over Fr, poly in p
    def pmul(a, b):
        out = [Fr(0)]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i+j] += x*y
        return out
    def padd(a, b):
        n = max(len(a), len(b))
        return [ (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(n) ]
    Pp, Qq = [Fr(0), Fr(1)], [Fr(1), Fr(-1)]
    def mono(z, o):
        out = [Fr(1)]
        for _ in range(z): out = pmul(out, Pp)
        for _ in range(o): out = pmul(out, Qq)
        return out
    from itertools import product as iproduct
    PH = [Fr(0)]
    for m in range(1, M+1):
        dm = m
        w0, w1 = schemes[m]
        for side in (0, 1):
            seen = {}
            for ext in iproduct((0,1), repeat=dm):
                wgt = (sum(1 for b in ext if b == 1) if side == 0
                       else sum(1 for b in ext if b == 0))
                r = seen.get(wgt, 0); seen[wgt] = r+1
                wv = (w0 if side == 0 else w1)[wgt]
                if r < wv:
                    word = ([side]*m + [1-side] + list(ext))
                    require(len(word) == 2*m+1, "deadline drift")
                    z = sum(1 for b in word if b == 0)
                    PH = padd(PH, mono(z, len(word)-z))
    # target: 1/2 - (p^{M+1}+q^{M+1})/2 + D_M
    tgt = padd([Fr(1,2)], [Fr(-1,2)*c for c in padd(mono(M+1,0), mono(0,M+1))])
    DM = [Fr(0)]
    for m, Am, rec in led.rows:
        for o, cap, par, delta in rec:
            if delta:
                DM = padd(DM, [Fr(delta,2)*c for c in mono(Am-o, o)])
    tgt = padd(tgt, DM)
    diff = padd(PH, [-c for c in tgt])
    require(all(c == 0 for c in diff), "realization identity failed")
    print(f"E_realization_control=gamma1_M{M}_all_stopped_words_enumerated;"
          " P_H == 1/2 - (p^{M+1}+q^{M+1})/2 + D_M as exact polynomial;"
          " pathwise stop length == T(m)=2m+1")

# ------------------------------------------------- F. certificate (27)
def certificate_recheck():
    def upper(t): return 2*(t + t**3/3 + t**5/(5*(1-t**2)))
    def lower(t): return 2*(t + t**3/3 + t**5/5)
    tA, tB = Fr(389,2181), Fr(5872957,11821757)
    F = Fr(391926968594914200867482400554891567498742649630277,
           82738859282193417287303438726081463937219800938169600)
    require(Fr(2457,6592)*lower(tB) - upper(tA) == Fr(1,25) + F and F > 0,
            "certificate identity drift")
    pA, pB = (1+tA)/2, (1+tB)/2
    print("F_certificate27=EXACT_recheck_ok: (2457/6592)L(tB)-U(tA)=1/25+F,F>0")
    # slope reverse-engineering numerics (reported, not asserted)
    g = Fr(2457, 6592)
    C = 1 + g
    print(f"  C_candidate=9049/6592={float(C):.6f} gamma={float(g):.6f}")
    print(f"  band_edges_at_gamma: theta=({float(g/(1+g)):.5f},"
          f"{float(1/(1+g)):.5f}) binding p=({float(1/(1+g)):.5f},"
          f"{float(g/(1+g)):.5f})")
    print(f"  p_A={float(pA):.5f} q_A={float(1-pA):.5f} "
          f"p_B={float(pB):.5f} q_B={float(1-pB):.5f}")
    la, lb = log(1285/896), log(8847357/2974400)
    print(f"  rho_A=log(pA/qA)={la:.6f} rho_B={lb:.6f} "
          f"(2457/6592)*rho_B-rho_A={g*Fr(1)*0 + (2457/6592)*lb-la:.6f}"
          f" vs 1/25=0.04")
    print(f"  slope_test: -log pA/log qA={-log(1285/2181)/log(896/2181):.6f}"
          f" vs 2457/(6592-2457)={2457/4135:.6f} (near-miss, flag numerology)")
    print(f"  1/(1+gamma)=6592/9049={6592/9049:.5f} vs p_B={float(pB):.5f}"
          f" (2.7% off); q_B/p_B={2974400/8847357:.5f} -> gamma'=0.336,"
          f" C'=1.336 (alt reading)")

# ---------------------------------------------------------------- main
def main():
    sys.setrecursionlimit(10000)
    geometry_audit()
    certificate_recheck()
    realization_control()

    plist = [Fr(1,2), Fr(1,3), Fr(2,3), Fr(1,4), Fr(3,4)]
    cps = [10, 20, 30, 40, 60, 80]
    print("columns: D_M(p) at p = 1/2, 1/3, 2/3, 1/4, 3/4")

    # control: gamma = 1 (no band ever) and hostile lazy policy
    run(1, 1, 0, 60, greedy, plist, [10, 30, 60], "control-gamma1")
    run(1, 2, 0, 40, lazy, plist, [10, 20, 40], "hostile-lazy-g0.5")

    # main sweep
    results = {}
    for gn, gd, D0, M in [(1,2,0,80), (1,2,6,80), (3,5,0,80), (3,5,6,80),
                          (3,4,0,80), (3,4,6,80), (9,10,0,120),
                          (19,20,0,150), (99,100,0,60)]:
        cp = sorted(set([10,20,30,40,60,80,100,120,150]) & set(range(1, M+1)))
        led, hist = run(gn, gd, D0, M, greedy, plist, cp,
                        f"g{gn}/{gd}-D{D0}")
        results[(gn,gd,D0)] = (led, hist)
        ms = sorted(hist)
        if len(ms) >= 2:
            a, b = hist[ms[-2]][0], hist[ms[-1]][0]
            frozen = (a != 0 and b != 0 and abs(float(a-b)) <
                      1e-3*max(abs(float(a)), abs(float(b))))
            print(f"  verdict[g={gn}/{gd},D0={D0}]: "
                  + ("FROZEN residual (no amortization): D stuck at "
                     f"{float(b):.3e} at p=1/2" if frozen else
                     ("decaying" if abs(float(b)) < abs(float(a))
                      else "growing/oscillating")))

    # frozen-residual law: log|D(1/2)| vs band-birth scale m* = (D0+2)/(1-g)
    print("--- frozen-residual law: log2|D_M(1/2)| vs (1+gamma)*m*")
    for (gn, gd, D0), (led, hist) in results.items():
        g = gn/gd
        if g >= 1: continue
        mstar = (D0 + 2) / (1 - g)
        last = float(hist[max(hist)][0])
        if last != 0:
            print(f"  g={g:.2f} D0={D0}: m*={mstar:6.1f} "
                  f"log2|D|={log(abs(last),2):8.2f} "
                  f"-(1+g)m*={-(1+g)*mstar:8.2f}")
        else:
            print(f"  g={g:.2f} D0={D0}: m*={mstar:6.1f} D==0 at final")

    # local search on gamma=3/5 D0=0
    led = results[(3,5,0)][0]
    fit = [Fr(1,2), Fr(1,3), Fr(2,3), Fr(2,5), Fr(3,4)]
    hold = [Fr(5,12), Fr(7,16), Fr(3,5), Fr(5,8), Fr(1,4), Fr(9,16)]
    before_f = [led.val_raw(p) for p in fit]
    before_h = [led.val_raw(p) for p in hold]
    deltas, rounds = local_search(led, fit, hold, m_hi=36, cap_hi=9)
    after_f = [led.val_raw(p) for p in fit]
    after_h = [led.val_raw(p) for p in hold]
    print(f"--- local sign-search g=3/5 D0=0: {len(deltas)} flips, "
          f"{rounds} rounds")
    print("  fit  before:", " ".join(f"{float(v):+.2e}" for v in before_f))
    print("  fit  after :", " ".join(f"{float(v):+.2e}" for v in after_f))
    print("  hold before:", " ".join(f"{float(v):+.2e}" for v in before_h))
    print("  hold after :", " ".join(f"{float(v):+.2e}" for v in after_h))
    print("status=LANE_D_LEDGER_COMPLETE")

if __name__ == "__main__":
    main()
