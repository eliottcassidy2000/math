"""
lrc_liu_zhu_pattern_structure_opus_S146.py   (opus-2026-07-07-S146, HYP-5217 part 2)

PROOF-MINING Liu-Zhu Conjecture 2 (confirmed exactly on 27 instances in part 1) toward a
GENERAL proof, via the window-graph optimal-pattern structure.

M = {x, y, y-x, y+x}, x = 2k+1, y = 2m+1, gcd(x,y) = 1;  mu(M) = (k+1)m / (4(k+1)m + 1).
Period of the optimal Motzkin set = N := 4(k+1)m + 1 = the denominator; elements per
period = (k+1)m.  Note gmM signature: x,y ODD; y-x, y+x EVEN -- the composite/small-gcd
structure the fleet's density-floor side keeps hitting (mac-mini-S56).

QUESTIONS:
 (1) NORMALIZED PATTERN: rotate each optimal set to canonical form (lexicographically
     minimal 0/1 necklace); is there a closed-form membership rule A = {a : phi(a) < mu}
     for an explicit witness rotation t = c/N?  (The witness-independent-set: A(t) =
     {a in Z_N : (a*t mod N) in [0, mu*N)} for the binding a/N.)
 (2) THE WITNESS RATIONAL: the max-cycle-mean optimum should be attained at t = a/N for a
     specific a; recover it (M(M) side gives the LR witness; is the Motzkin witness the
     SAME rational?  If mu = kappa the sandwich forces the same t).  Report kappa(M) too:
     Liu-Zhu Thm 5.7 says mu = kappa iff (x=1,y odd) or (distinct parity, none == 0 mod4);
     for x,y BOTH ODD the family has mu > kappa GENERICALLY -- so the Motzkin optimum is
     NOT a rotation slab (chi_c < 1/kappa territory) => the optimal set is genuinely
     combinatorial.  Confirm mu > kappa on the family and quantify the gap.
 (3) THE LOWER-BOUND FORMULA: guess an explicit periodic set of density (k+1)m/N and
     verify it symbolically across the family (the constructive half, made uniform).
 (4) POTENTIAL RANGE C(M) of the window graph (the Haralambis-lemma upper-bound constant):
     is C uniform / a simple function of (k,m)?  If S(N-window) <= (N + C)/... with C
     predictable, the open UPPER half generalizes.
"""
from fractions import Fraction as F
from math import gcd
import sys
import numpy as np
from collections import defaultdict

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc_mu_eq_M_maxcycle_opus_S144 import build_window_graph
from lrc_graph_interpretation_ladder_opus_S141 import M_exact


def optimal_set_and_potential(M, v):
    """Return (offsets of one optimal period-N set, N, potential range C)."""
    p, q = v.numerator, v.denominator
    states, t0, t1 = build_window_graph(M)
    V = len(states)
    D = np.zeros(V, dtype=np.int64)
    src1 = np.nonzero(t1 >= 0)[0]
    t1v = t1[src1]
    for _ in range(V + 5):
        Dn = D.copy()
        np.maximum.at(Dn, t0, D - p)
        np.maximum.at(Dn, t1v, D[src1] + (q - p))
        if np.array_equal(Dn, D):
            break
        D = Dn
    C = int(D.max() - D.min())
    # tight cycle
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
    start = min(outs)
    seen = {}; x = start; bits = []
    while x not in seen:
        seen[x] = len(bits); t, b = outs[x]; bits.append(b); x = t
    cyc = bits[seen[x]:]
    offs = [j for j, b in enumerate(cyc) if b == 1]
    return offs, len(cyc), C


def kappa_and_witness(M):
    """kappa(M) = M_exact, plus the binding (a, q) witness."""
    from lrc_graph_interpretation_ladder_opus_S141 import witness
    kap = M_exact(M)
    w = witness(M)
    return kap, w


def canonical_necklace(offs, N):
    """Lexicographically minimal rotation of the 0/1 indicator (as a tuple of gaps)."""
    ind = [0] * N
    for o in offs:
        ind[o % N] = 1
    best = None
    for r in range(N):
        rot = tuple(ind[(r + i) % N] for i in range(N))
        if best is None or rot < best:
            best = rot
    pos = [i for i, b in enumerate(best) if b]
    gaps = [(pos[(i + 1) % len(pos)] - pos[i]) % N for i in range(len(pos))]
    return pos, gaps


def try_witness_slab(M, N, dens):
    """Is there a rotation t = a/N s.t. A = {j : (a*j mod N) < dens*N} avoids M?
       (That would make the Motzkin optimum a slab = mu = kappa, chi_c = 1/kappa.)"""
    thr = int(dens * N)  # dens = (k+1)m/N so dens*N = numerator exactly
    for a in range(1, N):
        if gcd(a, N) != 1:
            continue
        A = set(j for j in range(N) if (a * j) % N < thr)
        if len(A) != thr:
            continue
        ok = all(((j + d) % N) not in A for j in A for d in M)
        if ok:
            return a
    return None


def main():
    print("=" * 104)
    print("LIU-ZHU CONJ 2 PATTERN STRUCTURE (proof-mining): M={x,y,y-x,y+x}, both odd")
    print("=" * 104)
    print(f"  {'(x,y)':>7} {'M':>16} {'mu=v':>8} {'kappa':>8} {'mu>kap?':>7} {'N':>5} "
          f"{'#elts':>5} {'C':>3} {'slab?':>6}  canonical gaps")
    rows = []
    for y in range(3, 21, 2):
        for x in range(1, y, 2):
            if gcd(x, y) != 1 or x + y > 20:
                continue
            k, m = (x - 1) // 2, (y - 1) // 2
            M = sorted({x, y, y - x, y + x})
            v = F((k + 1) * m, 4 * (k + 1) * m + 1)
            offs, N, C = optimal_set_and_potential(M, v)
            kap, w = kappa_and_witness(M)
            slab = try_witness_slab(M, N, v)
            pos, gaps = canonical_necklace(offs, N)
            gt = "yes" if v > kap else "EQ"
            rows.append((x, y, k, m, M, v, kap, N, len(offs), C, slab, gaps))
            gapstr = str(gaps) if len(gaps) <= 10 else str(gaps[:10]) + "..."
            print(f"  {str((x,y)):>7} {str(M):>16} {str(v):>8} {str(kap):>8} {gt:>7} "
                  f"{N:>5} {len(offs):>5} {C:>3} {str(slab):>6}  {gapstr}")

    print()
    print("STRUCTURAL OBSERVATIONS:")
    # C vs (k,m)
    print("  potential range C as a function of (k,m):")
    for (x, y, k, m, M, v, kap, N, ne, C, slab, gaps) in rows:
        print(f"    (k,m)=({k},{m}): C={C}, y+x={y+x} (largest diff), 2*(y+x)-1={2*(y+x)-1}"
              f"  C vs y+x: {'C = y+x' if C == y + x else ('C=2(y+x)-?' )}")
    # gap alphabet
    print("  gap alphabets (distinct gaps in the optimal necklace):")
    for (x, y, k, m, M, v, kap, N, ne, C, slab, gaps) in rows:
        alpha = sorted(set(gaps))
        print(f"    (x,y)=({x},{y}) M={M}: gaps use {alpha}"
              f"  {'(two-gap = Steinhaus/CF!)' if len(alpha) <= 2 else ''}")
    # mu > kappa quantification
    print("  mu - kappa gap (both-odd family): always positive?")
    allgt = all(v > kap for (_,_,_,_,_,v,kap,_,_,_,_,_) in rows if not (rows and False))
    xeq1 = [r for r in rows if r[0] == 1]
    xge3 = [r for r in rows if r[0] >= 3]
    print(f"    x=1 rows: mu vs kappa: {['EQ' if r[5]==r[6] else 'GT' for r in xeq1]}")
    print(f"    x>=3 rows: mu vs kappa: {['EQ' if r[5]==r[6] else 'GT' for r in xge3]}")
    print(f"    (Liu-Zhu Thm 5.7: mu=kappa iff x=1,y odd OR distinct-parity-none-mult-4;")
    print(f"     both-odd x>=3 => mu > kappa => optimum is NOT a rotation slab => genuinely")
    print(f"     combinatorial Motzkin set; slab column above should be None for x>=3.)")


if __name__ == "__main__":
    main()
