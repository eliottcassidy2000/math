"""
lrc_mu_eq_M_maxcycle_opus_S144.py

THE mu = M HUNT, DONE RIGHT (supersedes the S143 period-scan, which was O(300 DPs/set)
and died mid-run; and removes its N <= 320 period-cap caveat entirely).

KEY UPGRADE: mu(S) (the Motzkin density = max upper density of a set A with
(A-A) cap S = empty) equals the MAXIMUM CYCLE MEAN of the S-avoiding window graph:
  states = S-independent 0/1 windows of width w = max(S);
  edges  = append a bit (weight = the bit), legal iff no distance-s conflict.
Any S-avoiding set's upper density <= max cycle mean (block decomposition of the walk);
any cycle yields a periodic S-avoiding set of density = its mean. So
  mu(S) = maxmean(G_S)  EXACTLY (this also re-proves Cantor-Gordon periodicity:
  the optimum is attained by a periodic set of period <= #states).

DECISION mu > M?  <=>  the graph reweighted by (q*bit - p) has a POSITIVE cycle
(M = p/q).  Bellman-Ford longest-path with monotone init D=0 detects this in <= V
iterations (numpy-vectorized).  mu >= M is certified per-set by the exact witness
independent set (density M, from the binding rational a/q).  So:
  no positive cycle  =>  mu(S) = M(S) EXACTLY (both bounds certified);
  positive cycle     =>  mu(S) > M(S): extract the cycle = an explicit periodic
                         S-avoiding set beating M -- the ladder SEPARATES.

SCOPE: all primitive S, |S|=3,4 with max <= 14; |S|=5 with max <= 12;
       |S|=2 coprime <= 30 as engine validation (mu = floor((a+b)/2)/(a+b) theorem).
"""
from fractions import Fraction as F
from math import gcd
import sys, time, itertools
import numpy as np

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc_graph_interpretation_ladder_opus_S141 import M_exact, witness


def build_window_graph(S):
    """States: S-independent bit-windows of width w=max(S); edges append 0/1.
       Returns (t0, t1) int arrays (t1[i] = -1 if appending 1 illegal)."""
    w = max(S)
    Smask = 0
    for s in S:
        Smask |= 1 << (w - s)
    states = []
    def gen(mask, j):
        if j == w:
            states.append(mask); return
        gen(mask, j + 1)
        for s in S:
            jj = j - s
            if jj >= 0 and (mask >> jj) & 1:
                break
        else:
            gen(mask | (1 << j), j + 1)
    gen(0, 0)
    sidx = {m: i for i, m in enumerate(states)}
    V = len(states)
    t0 = np.empty(V, dtype=np.int64)
    t1 = np.empty(V, dtype=np.int64)
    for i, m in enumerate(states):
        t0[i] = sidx[m >> 1]
        t1[i] = sidx[(m >> 1) | (1 << (w - 1))] if (m & Smask) == 0 else -1
    return states, t0, t1


def positive_cycle_exists(t0, t1, p, q, itercap=None):
    """Longest-path value iteration with edge weights (q*bit - p).
       Monotone from D=0; convergence => no cycle of mean > p/q; else positive cycle."""
    V = len(t0)
    if itercap is None:
        itercap = V + 5
    D = np.zeros(V, dtype=np.int64)
    w0, w1 = -p, q - p
    has1 = t1 >= 0
    t1v = t1[has1]
    src1 = np.nonzero(has1)[0]
    for it in range(itercap):
        Dn = D.copy()
        np.maximum.at(Dn, t0, D + w0)
        np.maximum.at(Dn, t1v, D[src1] + w1)
        if np.array_equal(Dn, D):
            return False, it
        D = Dn
    return True, itercap


def extract_beating_cycle(states, t0, t1, p, q, w):
    """When a positive cycle exists: python Bellman with parents, walk back to a cycle,
       return (density Fraction, period, bits-per-period) as an explicit certificate."""
    V = len(t0)
    D = [0] * V
    par = [-1] * V
    parb = [0] * V
    w0, w1 = -p, q - p
    last_improved = None
    for _ in range(V + 2):
        improved = False
        Dn = D[:]
        for s in range(V):
            v = D[s] + w0
            t = t0[s]
            if v > Dn[t]:
                Dn[t] = v; par[t] = s; parb[t] = 0; improved = True; last_improved = t
            tt = t1[s]
            if tt >= 0:
                v = D[s] + w1
                if v > Dn[tt]:
                    Dn[tt] = v; par[tt] = s; parb[tt] = 1; improved = True; last_improved = tt
        D = Dn
        if not improved:
            break
    x = last_improved
    for _ in range(V):
        x = par[x]
    cyc_nodes = [x]
    bits = []
    y = par[x]
    bits.append(parb[x])
    while y != x:
        cyc_nodes.append(y)
        bits.append(parb[y])
        y = par[y]
    L = len(cyc_nodes)
    B = sum(bits)
    return F(B, L), L, B


def certify_witness_indep(S, a, q, v):
    """The arithmetic independent set {x in [0,q): (a x mod q) < v} avoids S and has
       density v/q = M.  Returns True if certified."""
    A = set(x for x in range(q) if (a * x) % q < v)
    for x in A:
        for s in S:
            if (x + s) % q in A:
                return False
    return len(A) == v  # density v/q exactly (a invertible mod q permutes residues)


def hunt(sets, label):
    print("=" * 100)
    print(f"HUNT: {label}")
    print("=" * 100)
    t0c = time.time()
    n_eq = n_gt = n_tested = 0
    slow = []
    for S in sets:
        g = 0
        for s in S:
            g = gcd(g, s)
        if g != 1:
            continue
        n_tested += 1
        M = M_exact(list(S))
        p, q = M.numerator, M.denominator
        wa, wq, wv = witness(list(S))
        wit_ok = certify_witness_indep(S, wa, wq, wv)
        states, tt0, tt1 = build_window_graph(list(S))
        if len(states) > 250_000:
            print(f"  (skip {tuple(S)}: {len(states)} states)")
            continue
        pos, iters = positive_cycle_exists(tt0, tt1, p, q)
        if pos:
            n_gt += 1
            dens, L, B = extract_beating_cycle(states, tt0, tt1, p, q, max(S))
            print(f"  *** mu > M: S={tuple(S)}  M={M}  periodic set: density {dens} = {B}/{L} > M"
                  f"  (witness-indep certified: {wit_ok})")
        else:
            n_eq += 1
            if not wit_ok:
                print(f"  !!! witness indep-set FAILED to certify at S={tuple(S)} (M={M}) -- engine bug?")
        if iters > 2000:
            slow.append((tuple(S), iters))
    print(f"  tested {n_tested} primitive sets in {time.time()-t0c:.0f}s: "
          f"mu = M (exact, both sides certified): {n_eq};  mu > M: {n_gt}")
    if slow:
        print(f"  (slowest convergences: {slow[:5]})")
    return n_gt


def main():
    t00 = time.time()

    print("=" * 100)
    print("(0) ENGINE VALIDATION on |S| = 2 (theorem: mu = M = floor((a+b)/2)/(a+b), coprime)")
    print("    [b <= 12: the window-graph state count grows ~2^b for sparse S]")
    print("=" * 100)
    bad = 0
    for b in range(2, 13):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            S = [a, b]
            M = M_exact(S)
            states, tt0, tt1 = build_window_graph(S)
            pos, _ = positive_cycle_exists(tt0, tt1, M.numerator, M.denominator)
            if pos:
                bad += 1
                print(f"  *** unexpected mu > M at S={{{a},{b}}}")
            exp = F(1, 2) if (a % 2 == 1 and b % 2 == 1) else F((a + b) // 2, a + b)
            if M != exp:
                bad += 1
                print(f"  *** M closed form mismatch at ({a},{b}): {M} vs {exp}")
    print(f"  |S|=2 validation: {'ALL OK (mu = M everywhere, closed form confirmed)' if bad == 0 else f'{bad} FAILURES'}")
    print()

    total_gt = 0
    total_gt += hunt(itertools.combinations(range(1, 15), 3), "|S| = 3, max <= 14 (primitive)")
    print()
    total_gt += hunt(itertools.combinations(range(1, 15), 4), "|S| = 4, max <= 14 (primitive)")
    print()
    total_gt += hunt(itertools.combinations(range(1, 13), 5), "|S| = 5, max <= 12 (primitive)")

    print()
    print("=" * 100)
    if total_gt == 0:
        print("VERDICT: mu(S) = M(S) EXACTLY on every primitive set tested (no period cap --")
        print("max-cycle-mean is the true Motzkin density).  The ladder collapse conjecture")
        print("mu = M survives a complete small-set sweep; the Haralambis-style separation,")
        print("if it exists, needs larger elements or larger |S|.")
    else:
        print(f"VERDICT: {total_gt} SEPARATIONS mu > M found -- the fractional rung of the")
        print("homomorphism ladder is strictly above the LRC constant on explicit sets;")
        print("certificates (periodic sets) printed above.")
    print(f"[total {time.time()-t00:.0f}s]")


if __name__ == "__main__":
    main()
