#!/usr/bin/env python3
"""(c) Paper 1: Niemeier root systems, Coxeter numbers, equal-h pairs, genus-1 theta, Table 2.1,
plus consistency tests of the genus-3 / genus-4 claims on root-type Fourier coefficients."""
import sys, time, itertools
sys.path.insert(0, '.')
import numpy as np
from fractions import Fraction
from rootsys import roots, gram, cartan, count_tuples, direct_sum

NIEMEIER = {
 'Leech': [], 'A1^24': [('A',1,24)], 'A2^12': [('A',2,12)], 'A3^8': [('A',3,8)], 'A4^6': [('A',4,6)],
 'A5^4D4': [('A',5,4),('D',4,1)], 'D4^6': [('D',4,6)], 'A6^4': [('A',6,4)], 'A7^2D5^2': [('A',7,2),('D',5,2)],
 'A8^3': [('A',8,3)], 'A9^2D6': [('A',9,2),('D',6,1)], 'D6^4': [('D',6,4)], 'A11D7E6': [('A',11,1),('D',7,1),('E',6,1)],
 'E6^4': [('E',6,4)], 'A12^2': [('A',12,2)], 'D8^3': [('D',8,3)], 'A15D9': [('A',15,1),('D',9,1)],
 'A17E7': [('A',17,1),('E',7,1)], 'D10E7^2': [('D',10,1),('E',7,2)], 'D12^2': [('D',12,2)], 'A24': [('A',24,1)],
 'D16E8': [('D',16,1),('E',8,1)], 'E8^3': [('E',8,3)], 'D24': [('D',24,1)]}
def cox(t, k): return {'A': k+1, 'D': 2*k-2}.get(t, {6: 12, 7: 18, 8: 30}.get(k))
def nroots(t, k): return {'A': k*(k+1), 'D': 2*k*(k-1)}.get(t, {6: 72, 7: 126, 8: 240}.get(k))

print("number of Niemeier root systems:", len(NIEMEIER))
H = {}
for name, comps in NIEMEIER.items():
    if not comps:
        H[name] = 0; continue
    hs = {cox(t, k) for t, k, m in comps}
    assert len(hs) == 1, name
    rank = sum(k*m for t, k, m in comps)
    nr = sum(nroots(t, k)*m for t, k, m in comps)
    h = hs.pop()
    assert rank == 24 and nr == 24*h, (name, rank, nr, h)
    H[name] = h
vals = sorted(set(H.values()))
print("distinct Coxeter numbers:", len(vals), vals)
pairs = {h: [n for n in H if H[n] == h] for h in vals if list(H.values()).count(h) > 1}
print("equal-h classes:", pairs)
print("all have rank 24 and #roots = 24h: OK")

# genus-1: Theta = E12 + (24h - 65520/691) Delta = E4^3 + (24h-720) Delta
assert Fraction(65520 + 432000, 691) == 720
print("E12 = E4^3 - (432000/691) Delta  =>  E12 + (24h - 65520/691)Delta = E4^3 + (24h - 720)Delta: OK")
# Leech: E4^3 - 720 Delta, q^2 coefficient
N = 6
def sigma(n, k): return sum(d**k for d in range(1, n+1) if n % d == 0)
E4 = [1] + [240*sigma(n, 3) for n in range(1, N)]
def mul(a, b): return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(N)]
P = [1] + [0]*(N-1)
for n in range(1, N):
    for _ in range(24): P = [P[i] - (P[i-n] if i >= n else 0) for i in range(N)]
Delta = [0] + P[:N-1]
E43 = mul(mul(E4, E4), E4)
leech = [E43[i] - 720*Delta[i] for i in range(N)]
print("Leech theta E4^3 - 720 Delta:", leech[:5], "(196560 minimal vectors:", leech[2] == 196560, ")")

# build root systems of the 10 lattices in the 5 pairs
def build(name):
    parts = []
    reps = []   # (index of a root in component, #roots in component)
    for t, k, m in NIEMEIER[name]:
        R, sc = roots(t, k)
        for _ in range(m):
            parts.append((R, sc))
    RR, L, lab = direct_sum(parts)
    IP = (RR @ RR.T) // L
    assert (np.diag(IP) == 2).all()
    # orbit representatives: one per isomorphism type of component
    seen = {}
    ci = 0
    for t, k, m in NIEMEIER[name]:
        first = int(np.where(lab == ci)[0][0])
        seen[(t, k)] = (first, nroots(t, k) * m)
        ci += m
    return IP, list(seen.values())

def count_L(IPreps, G):
    IP, reps = IPreps
    return sum(mult * count_tuples(IP, G, first_candidates=[idx]) for idx, mult in reps)

PAIRS = [(6, 'D4^6', 'A5^4D4'), (10, 'D6^4', 'A9^2D6'), (12, 'E6^4', 'A11D7E6'), (18, 'D10E7^2', 'A17E7'), (30, 'E8^3', 'D16E8')]
t0 = time.time()
LAT = {n: build(n) for _, a, b in PAIRS for n in (a, b)}
print(f"built root systems ({time.time()-t0:.1f}s)")

STAR = np.ones((5, 5), dtype=np.int64) + np.eye(5, dtype=np.int64)
printed = {6: ('A5', -5760), 10: ('D5', 69120), 12: ('D5', 230400), 18: ('D5', 2419200), 30: ('D5', 41287680)}
printed_star = {6: -5760, 10: -466560, 12: -1497600, 18: -15120000}
print("\nTable 2.1 / Theorem 2.8 / Theorem 4.11 (genus-5 root-frame coefficients a(T, F_h^(5)) = #frames(L1) - #frames(L2)):")
for h, a, b in PAIRS:
    X, val = printed[h]
    G = cartan(X[0], 5)
    d = count_L(LAT[a], G) - count_L(LAT[b], G)
    s = count_L(LAT[a], STAR) - count_L(LAT[b], STAR)
    ps = printed_star.get(h)
    print(f"  h={h}: a(T_{X}) = {d:,} (printed {val:,}) {'OK' if d == val else 'MISMATCH'};  star-frame a(T_cyc) = {s:,}" +
          (f" (printed {ps:,}) {'OK' if s == ps else 'MISMATCH'}" if ps is not None else ""))

# canonical classes of root-type Gram matrices (diag 2, offdiag in {-1,0,1}, positive definite)
def classes(m):
    reps = {}
    idx = [(i, j) for i in range(m) for j in range(i+1, m)]
    perms = list(itertools.permutations(range(m)))
    signs = list(itertools.product((1, -1), repeat=m))
    for vals in itertools.product((-1, 0, 1), repeat=len(idx)):
        G = 2*np.eye(m, dtype=np.int64)
        for (i, j), v in zip(idx, vals):
            G[i, j] = G[j, i] = v
        if np.linalg.eigvalsh(G.astype(float)).min() <= 1e-9:
            continue
        best = None
        for p in perms:
            Gp = G[np.ix_(p, p)]
            for s in signs:
                S = np.array(s)
                key = tuple((Gp * np.outer(S, S)).flatten())
                if best is None or key < best: best = key
        reps.setdefault(best, G)
    return list(reps.values())

C3 = classes(3); C4 = classes(4)
print(f"\nroot-type positive-definite Gram classes (up to permutation/sign): rank 3: {len(C3)}, rank 4: {len(C4)}")
print("Genus-3 test (Erokhin/Freitag claim: equal h => equal degree-3 theta): differences on all rank-3 root-type T:")
g3ok = True
for h, a, b in PAIRS:
    diffs = [count_L(LAT[a], G) - count_L(LAT[b], G) for G in C3]
    g3ok &= all(x == 0 for x in diffs)
    print(f"  h={h}: {diffs}")
print("  all zero:", g3ok)

print("\nGenus-4 test: F_h^(4) proportional to F_30^(4) = E4*(Theta_{E8^2}-Theta_{D16+}) on all rank-4 root-type T?")
D30 = [count_L(LAT['E8^3'], G) - count_L(LAT['D16E8'], G) for G in C4]
for h, a, b in PAIRS:
    D = [count_L(LAT[a], G) - count_L(LAT[b], G) for G in C4]
    ratios = set()
    ok = True
    for x, y in zip(D, D30):
        if y == 0:
            ok &= (x == 0)
        else:
            ratios.add(Fraction(x, y))
    print(f"  h={h}: ratios F_h/F_30 over {len(C4)} classes = {sorted(ratios)}; zeros consistent: {ok}")
GD4 = cartan('D', 4)
for h, a, b in PAIRS:
    dm = (count_L(LAT[a], GD4) - count_L(LAT[b], GD4)) // 1152
    print(f"  h={h}: Delta m(D4) = {dm}, Delta m(D4)/5 = {Fraction(dm,5)}, c_h = Delta m(D4)/4480 = {Fraction(dm, 4480)}")
print(f"(total {time.time()-t0:.1f}s)")
