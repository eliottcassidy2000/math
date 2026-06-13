# ocf_neg_candidates_kpo2_verify.py
# ADVERSARIAL VERIFIER, thread C (HYP-2380), claim C7 (refutation table of 31
# pre-registered closed-form candidates for I(Omega,-2) etc.). FRESH code.
#
# Candidates re-stated here BEFORE testing (same 31 forms the worker lists):
#   For r0 in {3/2, 1, 5/2, -3/2, -1}:                              (20 forms)
#     I(-2) = +W(r0), -W(r0), +2^(n-1) W(r0), -2^(n-1) W(r0)
#   I(-2) = det(J-2A);  |I(-2)| = det(J-2A)                          (2 forms)
#   I(-2) = +det(I-A), -det(I-A), +det(I+A), -det(I+A),
#           +per(I-A), -per(I-A)                                     (6 forms)
#   I(-1) = det(I-A);  I(1) = det(I+A);  I(1) = per(I+A)             (3 forms)
#
# W(r) = sum over all n! permutations P of prod_{i=0}^{n-2} (r + s_i),
#        s_i = +1/2 if P(i)->P(i+1) else -1/2 (THM-059 convention).
# Computed EXACTLY via integers: 2(r+s_i) = 2r+-1, so W = S / 2^(n-1) with S integer.
# det/per by Leibniz with exact integers. Scan: FULL labeled n=3,4,5.
# Sanity anchors: W(1/2)=H (THM-059) and W(-1/2)=(-1)^(n-1) H (THM-061) asserted.

import itertools
from fractions import Fraction
from collections import Counter

def pairs_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, pairs):
    adj = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(pairs):
        if (mask >> b) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def odd_cycle_supports(n, adj):
    res = []
    for k in range(3, n + 1, 2):
        for sub in itertools.combinations(range(n), k):
            v0 = sub[0]
            for perm in itertools.permutations(sub[1:]):
                cyc = (v0,) + perm
                if all(adj[cyc[i]][cyc[(i + 1) % k]] for i in range(k)):
                    m = 0
                    for v in sub:
                        m |= 1 << v
                    res.append(m)
    return res

def alphas(cycles):
    a1 = len(cycles)
    a2 = sum(1 for i in range(a1) for j in range(i + 1, a1)
             if cycles[i] & cycles[j] == 0)
    return a1, a2

def W_scaled(n, adj, two_r):
    """returns S = sum_P prod (2r + (1 if fwd else -1)); W(r) = S / 2^(n-1)."""
    S = 0
    for p in itertools.permutations(range(n)):
        prod = 1
        for i in range(n - 1):
            prod *= two_r + (1 if adj[p[i]][p[i + 1]] else -1)
        S += prod
    return S

def det_int(M):
    n = len(M)
    tot = 0
    for p in itertools.permutations(range(n)):
        seen = [False] * n
        sgn = 1
        for i in range(n):
            if not seen[i]:
                j = i; clen = 0
                while not seen[j]:
                    seen[j] = True; j = p[j]; clen += 1
                if clen % 2 == 0:
                    sgn = -sgn
        prod = 1
        for i in range(n):
            prod *= M[i][p[i]]
            if prod == 0:
                break
        tot += sgn * prod
    return tot

def per_int(M):
    n = len(M)
    tot = 0
    for p in itertools.permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i][p[i]]
            if prod == 0:
                break
        tot += prod
    return tot

R0S = [Fraction(3, 2), Fraction(1), Fraction(5, 2), Fraction(-3, 2), Fraction(-1)]

candidates = []  # (name, function(n, data) -> bool holds)
for r0 in R0S:
    for sgn in (1, -1):
        candidates.append((f"I(-2) == {'+' if sgn>0 else '-'}W({r0})",
                           ('W', r0, sgn, False)))
        candidates.append((f"I(-2) == {'+' if sgn>0 else '-'}2^(n-1)*W({r0})",
                           ('W', r0, sgn, True)))
candidates += [
    ("I(-2) == det(J-2A)", ('M', 'Jm2A', 1, 'Im2')),
    ("|I(-2)| == det(J-2A)", ('M', 'Jm2A', 1, 'absIm2')),
    ("I(-2) == +det(I-A)", ('M', 'ImA_det', 1, 'Im2')),
    ("I(-2) == -det(I-A)", ('M', 'ImA_det', -1, 'Im2')),
    ("I(-2) == +det(I+A)", ('M', 'IpA_det', 1, 'Im2')),
    ("I(-2) == -det(I+A)", ('M', 'IpA_det', -1, 'Im2')),
    ("I(-2) == +per(I-A)", ('M', 'ImA_per', 1, 'Im2')),
    ("I(-2) == -per(I-A)", ('M', 'ImA_per', -1, 'Im2')),
    ("I(-1) == det(I-A)", ('M', 'ImA_det', 1, 'Im1')),
    ("I(1) == det(I+A)", ('M', 'IpA_det', 1, 'I1')),
    ("I(1) == per(I+A)", ('M', 'IpA_per', 1, 'I1')),
]
print(f"number of pre-registered candidates: {len(candidates)} (expect 31)")

first_cex = {}   # name -> (n, mask, lhs, rhs)
alive = {name for name, _ in candidates}

for n in (3, 4, 5):
    pairs = pairs_list(n)
    for mask in range(1 << len(pairs)):
        if not alive and n > 3:
            break
        adj = adj_from_mask(n, mask, pairs)
        cyc = odd_cycle_supports(n, adj)
        a1, a2 = alphas(cyc)
        Im2 = 1 - 2 * a1 + 4 * a2
        Im1 = 1 - a1 + a2
        I1 = 1 + a1 + a2
        H = 1 + 2 * a1 + 4 * a2  # OCF (verified in census script)
        # sanity anchors
        S_half = W_scaled(n, adj, 1)        # 2r=1 -> W(1/2)
        assert S_half == (2 ** (n - 1)) * H // (2 ** (n - 1)) * (2 ** (n - 1)) // (2 ** (n - 1)) or True
        assert Fraction(S_half, 2 ** (n - 1)) == H, ("THM-059 anchor", n, mask)
        S_mhalf = W_scaled(n, adj, -1)
        assert Fraction(S_mhalf, 2 ** (n - 1)) == (-1) ** (n - 1) * H, ("THM-061 anchor", n, mask)
        # W values
        Wv = {}
        for r0 in R0S:
            S = 0
            two_r = 2 * r0  # Fraction, integer-valued for our r0 list
            assert two_r.denominator == 1
            S = W_scaled(n, adj, int(two_r))
            Wv[r0] = Fraction(S, 2 ** (n - 1))
        # matrices
        J2A = [[1 - 2 * adj[i][j] for j in range(n)] for i in range(n)]
        ImA = [[(1 if i == j else 0) - adj[i][j] for j in range(n)] for i in range(n)]
        IpA = [[(1 if i == j else 0) + adj[i][j] for j in range(n)] for i in range(n)]
        Mv = {'Jm2A': det_int(J2A), 'ImA_det': det_int(ImA),
              'IpA_det': det_int(IpA), 'ImA_per': per_int(ImA),
              'IpA_per': per_int(IpA)}
        for name, spec in candidates:
            if name not in alive:
                continue
            if spec[0] == 'W':
                _, r0, sgn, scaled = spec
                rhs = sgn * Wv[r0] * (2 ** (n - 1) if scaled else 1)
                lhs = Fraction(Im2)
            else:
                _, key, sgn, which = spec
                rhs = sgn * Mv[key]
                lhs = {'Im2': Im2, 'absIm2': abs(Im2), 'Im1': Im1, 'I1': I1}[which]
            if lhs != rhs:
                alive.discard(name)
                first_cex[name] = (n, mask, lhs, rhs)

print("\nREFUTATION TABLE (my first counterexample, my bit convention bit=>i->j):")
for name, _ in candidates:
    if name in first_cex:
        n, mask, lhs, rhs = first_cex[name]
        print(f"  REFUTED  {name:36s} at n={n} mask={mask}: lhs={lhs} rhs={rhs}")
    else:
        print(f"  SURVIVED {name} on full labeled n=3..5  <-- WOULD CONTRADICT WORKER")

print(f"\nrefuted: {len(first_cex)}/{len(candidates)}; survivors: {len(alive)}")

# Specific spot values for the n=3 transitive tournament (mask=7: 0->1,0->2,1->2)
n = 3
pairs = pairs_list(3)
adj = adj_from_mask(3, 0b111, pairs)
print("\nn=3 transitive (arcs 0->1,0->2,1->2): alpha=(0,0), I(-2)=1")
for r0 in R0S:
    S = W_scaled(3, adj, int(2 * r0))
    print(f"  W({r0}) = {Fraction(S,4)}   2^(n-1) W = {Fraction(S,1)}")
J2A = [[1 - 2 * adj[i][j] for j in range(3)] for i in range(3)]
ImA = [[(1 if i == j else 0) - adj[i][j] for j in range(3)] for i in range(3)]
IpA = [[(1 if i == j else 0) + adj[i][j] for j in range(3)] for i in range(3)]
print(f"  det(J-2A)={det_int(J2A)}  det(I-A)={det_int(ImA)}  det(I+A)={det_int(IpA)}"
      f"  per(I-A)={per_int(ImA)}  per(I+A)={per_int(IpA)}")
print("done")
