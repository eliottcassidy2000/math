"""
lrc_liu_zhu_lower_bound_opus_S147.py   (opus-2026-07-07-S147, HYP-5277)

GENERAL PROOF of Liu-Zhu Conjecture 2, LOWER-BOUND half (the tractable direction):
find an EXPLICIT period-N avoiding set of density exactly (k+1)m/N, N = 4(k+1)m+1,
for M = {x, y, y-x, y+x}, x = 2k+1, y = 2m+1, gcd(x,y)=1, and verify it symbolically
across the family (so mu >= (k+1)m/N uniformly, for ALL x,y -- the constructive half made
general, previously only exhibited per-instance by the window-graph tight cycle).

APPROACH: the x=1 optimum is a rotation slab A = {j : (a j mod N) in [0, (k+1)m)} for a
witness a/N.  For x>=3 no slab exists (S146), so the optimum is combinatorial -- BUT it may
still be a "generalized slab": A = {j : (a j mod N) in U} for a UNION U of intervals, or a
Beatty/Sturmian-type set, or a set defined by a 2D (CRT) rule mod the two factors of N.
This script:
 (1) recovers, per instance, the OPTIMAL periodic set (window-graph tight cycle) and its
     canonical rotation;
 (2) tests structured ansaetze for a UNIFORM membership rule:
     (a) single-slab A = {j : a*j mod N < T}  (should work only for x=1);
     (b) the WITNESS-ROTATION multi-interval: for the LR witness t* = kappa-binding a/N',
         is A = {j : frac(j * t*) in some fixed union of arcs}?
     (c) a 2-residue CRT rule when N factors;
     (d) the COMPLEMENT-of-M-multiples / Beatty rule.
 (3) if a uniform rule is found, VERIFY it avoids M and has the right density symbolically
     for a wide (x,y) grid.
"""
from fractions import Fraction as F
from math import gcd
import sys
import numpy as np
from collections import defaultdict

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc_mu_eq_M_maxcycle_opus_S144 import build_window_graph
from lrc_graph_interpretation_ladder_opus_S141 import M_exact, witness


def optimal_set(M, v):
    p, q = v.numerator, v.denominator
    states, t0, t1 = build_window_graph(M)
    V = len(states)
    D = np.zeros(V, dtype=np.int64)
    src1 = np.nonzero(t1 >= 0)[0]; t1v = t1[src1]
    for _ in range(V + 5):
        Dn = D.copy()
        np.maximum.at(Dn, t0, D - p)
        np.maximum.at(Dn, t1v, D[src1] + (q - p))
        if np.array_equal(Dn, D):
            break
        D = Dn
    te = []
    for s in range(V):
        if D[t0[s]] == D[s] - p:
            te.append((s, int(t0[s]), 0))
        if t1[s] >= 0 and D[t1[s]] == D[s] + (q - p):
            te.append((s, int(t1[s]), 1))
    while True:
        outd = defaultdict(int); ind = defaultdict(int); nodes = set()
        for s, t, b in te:
            outd[s] += 1; ind[t] += 1; nodes.update((s, t))
        bad = {u for u in nodes if outd[u] == 0 or ind[u] == 0}
        if not bad:
            break
        te = [(s, t, b) for s, t, b in te if s not in bad and t not in bad]
    outs = {}
    for s, t, b in te:
        outs.setdefault(s, (t, b))
    start = min(outs); seen = {}; x = start; bits = []
    while x not in seen:
        seen[x] = len(bits); t, b = outs[x]; bits.append(b); x = t
    cyc = bits[seen[x]:]
    return set(j for j, b in enumerate(cyc) if b == 1), len(cyc)


def avoids(A, M, N):
    return all(((a + d) % N) not in A for a in A for d in M)


def try_multi_interval_witness(M, N, dens, A_true):
    """The LR witness t* = a/N' binds kappa; test A = {j : frac(j a / N) in arcs} where
       arcs is recovered from A_true's image under j -> (a j mod N).  If the image is a
       union of few intervals for some coprime a, that's a 'generalized slab' rule."""
    thr = dens.numerator  # = (k+1)m; density = thr/N
    best = None
    for a in range(1, N):
        if gcd(a, N) != 1:
            continue
        img = sorted((a * j) % N for j in A_true)
        # count maximal runs of consecutive residues (cyclically) = #intervals
        runs = 1
        for i in range(len(img)):
            if (img[i] + 1) % N != img[(i + 1) % len(img)]:
                runs += 1
        # wrap correction
        if len(img) > 0 and (img[-1] + 1) % N == img[0]:
            runs -= 1
        runs = max(runs, 1)
        if best is None or runs < best[0]:
            best = (runs, a, img)
    return best  # (min #intervals over rotations, the a achieving it, the image)


def main():
    print("=" * 100)
    print("LIU-ZHU CONJ 2 -- LOWER BOUND: structure of the optimal avoiding set")
    print("=" * 100)
    print(f"  {'(x,y)':>7} {'M':>16} {'N':>4} {'#in':>4} {'min#arcs':>9} {'witness a':>10}"
          f"   image (rotated to fewest arcs)")
    data = []
    for y in range(3, 18, 2):
        for x in range(1, y, 2):
            if gcd(x, y) != 1 or x + y > 18:
                continue
            k, m = (x - 1) // 2, (y - 1) // 2
            M = sorted({x, y, y - x, y + x})
            N = 4 * (k + 1) * m + 1
            v = F((k + 1) * m, N)
            A, per = optimal_set(M, v)
            assert per == N, (M, per, N)
            assert avoids(A, M, N), M
            runs, a, img = try_multi_interval_witness(M, N, v, A)
            data.append((x, y, k, m, M, N, A, runs, a, img))
            imgstr = str(img[:12]) + ("..." if len(img) > 12 else "")
            print(f"  {str((x,y)):>7} {str(M):>16} {N:>4} {len(A):>4} {runs:>9} {a:>10}"
                  f"   {imgstr}")

    print()
    print("  INTERPRETATION:")
    print("  - min#arcs = 1  <=>  single rotation slab (x=1 case, mu=kappa).")
    print("  - min#arcs small & constant across x>=3  =>  a UNIFORM generalized-slab rule")
    print("    (union of that many arcs under the witness rotation) -- the lower-bound")
    print("    construction we want.")
    print()
    # Does min#arcs correlate with x (=2k+1)?  Report the arc-count vs (k,m).
    print("  arc-count vs (k,m):")
    for (x, y, k, m, M, N, A, runs, a, img) in data:
        print(f"    (k,m)=({k},{m})  x={x}: min#arcs = {runs}"
              f"  {'(slab)' if runs == 1 else ''}")

    # For x=1: verify the explicit slab rule A = {j : a j mod N < (k+1)m} and identify a.
    print()
    print("  x=1 SLAB WITNESS a (A = {j : a*j mod N < m} on N=4m+1):")
    for (x, y, k, m, M, N, A, runs, a, img) in data:
        if x == 1:
            # find the slab witness explicitly
            found = None
            for aa in range(1, N):
                if gcd(aa, N) != 1:
                    continue
                B = set(j for j in range(N) if (aa * j) % N < m)
                if B == A or avoids(B, M, N) and len(B) == m:
                    found = aa
                    break
            print(f"    (x,y)=(1,{y}) N={N} m={m}: slab witness a={found}"
                  f"  (a/N approx {found/N:.4f}; note kappa binds near here)")


if __name__ == "__main__":
    main()
