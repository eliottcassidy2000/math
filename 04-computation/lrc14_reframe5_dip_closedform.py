#!/usr/bin/env python3
"""
REFRAME 5 (mac-mini-2026-06-17): self-similar recursion with EXPLICIT dip control.

Level-n covering family for LRC(N), N=14.  The "core" of a covering-N set is
covering-(N-1), whose core is covering-(N-2), ... down to covering-7 (PROVEN base).

At level n we look at the family:
    A_n = {1,...,n} \ {q}          (an "easy core": uncovered at modulus q)
    S   = A_n  u  { w }            where w = lcm(q, n) * k  parks a runner
                                   that re-covers modulus q (perfect-middle runner).
M_n = M(A_n) is the easy-core gap; dip_n(k) = M(A_n) - M(S) is the resonance dip.

REFRAME 5 questions:
  (1) per-level dip closed form generalizing M = 7m/(84m+5)  (the N=14, q=12 family);
  (2) THE KEY non-tautological question: is dip_n bounded by an explicit function NOT
      requiring M(S) >= 1/N as input? (e.g. dip <= c/lcm, or dip <= flank/D, flank<=n);
  (3) do the per-level dips telescope below M(base) - 1/14 INDEPENDENTLY?

Exact M tool exactly as supplied in the prompt.  stdlib only.
"""
from fractions import Fraction as F
from math import gcd


# ---------------------------------------------------------------- exact M tool
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C


def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at


def lcm(a, b):
    return a * b // gcd(a, b)


# ----------------------------------------------- binding-pair diagnostics at tau*
def binding(S, tstar):
    """Return (gap_value, list of binding runners, list of (kind,a,b,D) crossings)."""
    S = sorted(set(S)); val = g(S, tstar)
    mr = [v for v in S if nrm(v * tstar) == val]
    pairs = []
    for i in range(len(mr)):
        for j in range(i + 1, len(mr)):
            a, b = mr[i], mr[j]
            if nrm((a + b) * tstar) == 0:
                pairs.append(('sum', a, b, a + b))
            if b != a and nrm((b - a) * tstar) == 0:
                pairs.append(('diff', a, b, b - a))
    return val, mr, pairs


N = 14
print("=" * 78)
print("REFRAME 5 - per-level dip closed forms + NON-CIRCULARITY of the dip bound")
print("=" * 78)

# =====================================================================
# PART 1: canonical N=14, q=12 family  M = 7m/(84m+5)  re-derived exactly.
# =====================================================================
print("\n" + "-" * 78)
print("(1) Canonical N=14, q=12 family  A={1..11,13}, w=84k, re-derived exactly")
print("-" * 78)
A = list(range(1, 12)) + [13]            # {1..11,13}; uncovered at 12
MA, tA = M(A)
print(f"  A={A}")
print(f"  M(A) = {MA}  (= 1/12 = {F(1,12)})   tau_A = {tA}")
for k in range(1, 8):
    w = 84 * k
    S = A + [w]
    MS, tS = M(S)
    dip = MA - MS
    val, mr, pairs = binding(S, tS)
    pred = F(7 * k, 84 * k + 5)
    Ds = [p[3] for p in pairs]
    print(f"  k={k} w={w:4d}: M(S)={str(MS):8s} (=7k/(84k+5)? {MS==pred})  "
          f"dip={str(dip):10s}  bind={mr}  D(crossings)={Ds}")

# =====================================================================
# PART 2: GENERAL family A_n={1..n}\{q}, park w = lcm(q,n)*k.
# =====================================================================
print("\n" + "-" * 78)
print("(2) GENERAL family  A_n={1..n}\\{q},  w = lcm(q,n)*k   (q-resonant runner)")
print("-" * 78)

records = []
for n in range(7, 15):
    for q in range(2, n + 1):
        A = [v for v in range(1, n + 1) if v != q]
        MA, tA = M(A)
        L = lcm(q, n)
        rowdips = []
        for k in range(1, 7):
            w = L * k
            S = A + [w]
            MS, tS = M(S)
            dip = MA - MS
            val, mr, pairs = binding(S, tS)
            flank = None; Dval = None
            for (kind, a, b, D) in pairs:
                if w in (a, b):
                    flank = a if b == w else b; Dval = D; break
            rowdips.append((k, w, MS, dip, flank, Dval, mr))
        maxdip = max(r[3] for r in rowdips)
        records.append((n, q, A, MA, L, rowdips, maxdip))

print(f"\n  {'n':>2} {'q':>2} {'M(A_n)':>9} {'lcm':>5} "
      f"{'M(S)k=1':>10} {'dip k=1':>12} {'flank':>5} {'D':>6} {'maxdip k<=6':>12}")
for (n, q, A, MA, L, rowdips, maxdip) in records:
    k, w, MS, dip, flank, Dval, mr = rowdips[0]
    print(f"  {n:>2} {q:>2} {str(MA):>9} {L:>5} {str(MS):>10} "
          f"{str(dip):>12} {str(flank):>5} {str(Dval):>6} {str(maxdip):>12}")

# =====================================================================
# PART 3: dip CLOSED FORM as rational function of k.
# =====================================================================
print("\n" + "-" * 78)
print("(3) CLOSED FORM:  M(S) = p*k/(k+r)  fit per (n,q); dip(k) closed form")
print("-" * 78)


def fit_rational(ks, Ms):
    """M(k) = p k/(k+r); solve from first two distinct, verify all."""
    k1, k2 = ks[0], ks[1]
    M1, M2 = Ms[0], Ms[1]
    denom = M1 * k2 - M2 * k1
    if denom == 0:
        return None
    r = (M2 - M1) * k1 * k2 / denom
    if k1 + r == 0:
        return None
    p = M1 * (k1 + r) / k1
    for k, m in zip(ks, Ms):
        if k + r == 0 or p * k / (k + r) != m:
            return None
    return p, r


closedforms = {}
for (n, q, A, MA, L, rowdips, maxdip) in records:
    ks = [r[0] for r in rowdips]
    Ms = [r[2] for r in rowdips]
    fit = fit_rational(ks, Ms)
    if fit:
        p, r = fit
        diplim = MA - p
        closedforms[(n, q)] = (p, r, diplim)
        print(f"  n={n:>2} q={q:>2}: M(S)= {p}*k/(k+{r})  "
              f"k->inf M={p}  dip_limit={diplim}  M(A)={MA}")
    else:
        print(f"  n={n:>2} q={q:>2}: NOT p k/(k+r) form. Ms={Ms}")

# =====================================================================
# PART 4: NON-CIRCULARITY tests.
# =====================================================================
print("\n" + "-" * 78)
print("(4) NON-CIRCULARITY: bounds NOT using M(S)>=1/N")
print("-" * 78)

print("\n  (B1)  dip <= 1/lcm(q,n) ?   (purely arithmetic in q,n)")
viol_b1 = 0; tight_b1 = None
for (n, q, A, MA, L, rowdips, maxdip) in records:
    for (k, w, MS, dip, flank, Dval, mr) in rowdips:
        if dip > F(1, L): viol_b1 += 1
        ratio = dip / F(1, L)
        if tight_b1 is None or ratio > tight_b1[0]:
            tight_b1 = (ratio, n, q, k, dip, F(1, L))
print(f"     violations of dip<=1/lcm: {viol_b1}")
print(f"     tightest dip/(1/lcm): {float(tight_b1[0]):.4f} at "
      f"n={tight_b1[1]} q={tight_b1[2]} k={tight_b1[3]} dip={tight_b1[4]} 1/lcm={tight_b1[5]}")

print("\n  (B4)  flank <= n  for binding pair (flank,w)?")
viol_flank = 0; maxflank = None
for (n, q, A, MA, L, rowdips, maxdip) in records:
    for (k, w, MS, dip, flank, Dval, mr) in rowdips:
        if flank is None: continue
        if flank > n: viol_flank += 1
        if maxflank is None or flank > maxflank[0]:
            maxflank = (flank, n, q, k)
print(f"     violations of flank<=n: {viol_flank}")
if maxflank:
    print(f"     max flank: {maxflank[0]} at n={maxflank[1]} q={maxflank[2]} k={maxflank[3]}")

print("\n  (CRUX) binding D=flank+w, M(S)=j/D.  D<=N*j  <=>  M(S)>=1/N (tautology).")
print("         INDEPENDENT route would be: j (crossing index) bounded below by")
print("         (flank,w) arithmetic so that j/D >= 1/N follows WITHOUT assuming it.")
print("\n  Detail (N=14 level: A={1..13}\\{q}, re-cover modulus 14, w=lcm(q,14)):")
for q in range(2, 14):
    A = [v for v in range(1, 14) if v != q]
    MA, _ = M(A)
    L = lcm(q, 14)
    w = L
    S = A + [w]
    MS, tS = M(S)
    dip = MA - MS
    val, mr, pairs = binding(S, tS)
    flank = None; Dval = None
    for (kind, a, b, D) in pairs:
        if w in (a, b):
            flank = a if b == w else b; Dval = D; break
    geN = MS >= F(1, 14)
    Dok = (Dval is not None and Dval <= 14 * MS.numerator)
    print(f"   q={q:>2}: w=lcm({q},14)={L:>3}, M(A)={str(MA):>6}, M(S)={str(MS):>7}, "
          f"flank={flank}, D={Dval}, j={MS.numerator}, D<=14j?{Dok}, M>=1/14:{geN}")

print("\nDONE.")
